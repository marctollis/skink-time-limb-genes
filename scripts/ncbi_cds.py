#!/usr/bin/env python3
"""Datasets gene-list downloads and longest RefSeq CDS per species."""
import csv
import hashlib
import json
import re
import subprocess
import sys
import zipfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from Bio.Seq import Seq
from pipeline import fasta

SCOPE = 'datasets_refseq_lepidosauria'


def gene_list(path):
    genes = []
    for line in Path(path).read_text().splitlines():
        gene = line.split('#', 1)[0].strip().upper()
        if not gene:
            continue
        if not re.fullmatch(r'[A-Z0-9][A-Z0-9_.-]*', gene):
            raise ValueError(f'Invalid gene symbol in list: {gene!r}')
        if gene in genes:
            raise ValueError(f'Duplicate gene in list: {gene}')
        genes.append(gene)
    if not genes:
        raise ValueError('Gene list is empty')
    return genes


def package_files(directory):
    data = directory / 'downloads/package/ncbi_dataset/data'
    return data / 'cds.fna', data / 'data_report.jsonl'


def download(directory, gene):
    directory = Path(directory)
    dest = directory / 'downloads'
    if dest.exists():
        raise ValueError('downloads already exists; move the old gene directory aside before restarting')
    dest.mkdir(parents=True)
    archive = dest / f'{gene}_lepidosauria.zip'
    command = ['datasets', 'download', 'gene', 'symbol', gene, '--taxon', 'Lepidosauria',
               '--include', 'cds,product-report', '--filename', str(archive)]
    manifest = dict(gene=gene, scope=SCOPE, complete=False, command=command,
                    retrieved_utc=datetime.now(timezone.utc).isoformat())
    manifest_path = directory / 'download_manifest.json'
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    with (dest / 'download.log').open('w') as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode:
        raise ValueError(f'Datasets failed; see {dest}/download.log')
    with zipfile.ZipFile(archive) as package:
        target = (dest / 'package').resolve()
        for member in package.infolist():
            resolved = (target / member.filename).resolve()
            if target != resolved and target not in resolved.parents:
                raise ValueError('Unsafe path in downloaded ZIP')
        if package.testzip():
            raise ValueError('Downloaded ZIP failed CRC validation')
        package.extractall(target)
    cds, report = package_files(directory)
    if not cds.is_file() or not report.is_file():
        raise ValueError('Package lacks cds.fna or data_report.jsonl; inspect download.log')
    count = sum(1 for _ in fasta(cds))
    if not count:
        raise ValueError('Package contains no CDS sequences; inspect download.log and metadata')
    # Reading every JSON line detects truncated/malformed metadata before completion.
    reports = [json.loads(line) for line in report.read_text().splitlines() if line.strip()]
    if not reports:
        raise ValueError('Empty gene report')
    manifest.update(complete=True, downloaded_cds_count=count, gene_report_count=len(reports),
                    cds_sha256=hashlib.sha256(cds.read_bytes()).hexdigest(),
                    report_sha256=hashlib.sha256(report.read_bytes()).hexdigest())
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    print(f'{gene}: downloaded {count} CDS records', flush=True)


def completed(directory, gene):
    path = directory / 'download_manifest.json'
    if not path.exists():
        return False
    manifest = json.loads(path.read_text())
    cds, report = package_files(directory)
    return (manifest.get('complete') and manifest.get('scope') == SCOPE
            and manifest.get('gene') == gene and cds.is_file() and report.is_file()
            and hashlib.sha256(cds.read_bytes()).hexdigest() == manifest.get('cds_sha256')
            and hashlib.sha256(report.read_bytes()).hexdigest() == manifest.get('report_sha256'))


def download_list(root, path):
    root = Path(root)
    genes = gene_list(path)
    (root / 'genes').mkdir(parents=True, exist_ok=True)
    statuses, failed = [], False
    for gene in genes:
        directory = root / 'genes' / gene
        directory.mkdir(exist_ok=True)
        print(f'\nDownloading {gene}: RefSeq Lepidosauria CDS', flush=True)
        try:
            if completed(directory, gene):
                status, note = 'already_complete', 'Verified existing matching download'
                print(f'{gene}: already downloaded; skipping', flush=True)
            else:
                download(directory, gene)
                status, note = 'downloaded', ''
        except Exception as error:
            status, note, failed = 'failed', str(error), True
            print(f'{gene}: {note}', file=sys.stderr, flush=True)
        statuses.append(dict(gene=gene, status=status, notes=note))
        with (root / 'genes/download_batch_status.tsv').open('w') as handle:
            writer = csv.DictWriter(handle, fieldnames=['gene', 'status', 'notes'], delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(statuses)
    if failed:
        raise ValueError('Some genes failed; inspect genes/download_batch_status.tsv. Successful genes were retained.')


FIELDS = ['species', 'organism', 'label', 'taxid', 'source', 'transcript_accession',
          'gene_symbol', 'gene_id', 'CDS_length', 'eligible', 'selected', 'notes', 'original_header']


def table(path, rows):
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows({key: row.get(key, '') for key in FIELDS} for row in rows)


def extract(directory, gene):
    directory = Path(directory)
    output = directory / f'{gene}_longest_cds.fasta'
    if output.exists() or (directory / 'candidate_cds.tsv').exists():
        raise ValueError('Extraction output already exists; move old extraction outputs aside before restarting')
    if not completed(directory, gene):
        raise ValueError('Missing, incomplete, modified or old-scope package; run step 01 in a fresh gene directory')
    cds, report = package_files(directory)
    reports = [json.loads(line) for line in report.read_text().splitlines() if line.strip()]
    metadata = {str(row['geneId']): row for row in reports}
    candidates, groups = [], defaultdict(list)
    for header, sequence in fasta(cds):
        attributes = dict(re.findall(r'\[([^=\]]+)=([^\]]*)\]', header))
        geneid = attributes.get('GeneID', '')
        info = metadata.get(geneid, {})
        organism = attributes.get('organism', info.get('taxname', ''))
        words = organism.split()
        species = ' '.join(words[:2])
        accession = header.split()[0].split(':')[0]
        symbol = info.get('symbol', '')
        problems, warnings = [], []
        if not re.fullmatch(r'(?:NM|XM)_\d+\.\d+', accession):
            problems.append('not_RefSeq_coding_transcript')
        if not info:
            problems.append('missing_gene_metadata')
        if info.get('taxname') and info['taxname'] != organism:
            problems.append('organism_metadata_mismatch')
        if not re.fullmatch(r'[A-Z][a-z]+ [a-z]+', species) or (len(words) > 1 and words[1] in {'sp', 'cf', 'aff', 'hybrid'}):
            problems.append('unresolved_species')
        synonyms = {str(x).upper() for x in info.get('synonyms', [])}
        if symbol.upper() != gene:
            if gene in synonyms:
                warnings.append('matched_synonym_review')
            else:
                problems.append('gene_symbol_mismatch')
        if info.get('type') and info['type'] != 'PROTEIN_CODING':
            problems.append('not_protein_coding_gene')
        if not re.fullmatch(r'[ACGTRYSWKMBDHVN]+', sequence) or not set(sequence) & set('ACGT'):
            problems.append('invalid_or_uninformative_sequence')
        elif len(sequence) % 3:
            problems.append('length_not_multiple_of_3')
        else:
            protein = str(Seq(sequence).translate())
            if '*' in protein[:-1]:
                problems.append('internal_stop')
            if not protein.endswith('*'):
                warnings.append('no_terminal_stop_review')
            if sequence[:3] != 'ATG':
                warnings.append('non_ATG_start_review')
            if re.search(r'[^ACGT]', sequence):
                warnings.append('ambiguous_bases_review')
        row = dict(species=species, organism=organism, label=species.replace(' ', '_'),
                   taxid=str(info.get('taxId', '')), source='RefSeq',
                   transcript_accession=accession, gene_symbol=symbol, gene_id=geneid,
                   CDS_length=len(sequence), eligible='no' if problems else 'yes',
                   selected='no', notes=';'.join(problems + warnings), original_header=header,
                   sequence=sequence)
        candidates.append(row)
        if not problems:
            groups[species].append(row)
    selections = []
    for species, rows in sorted(groups.items()):
        if len({r['gene_id'] for r in rows}) > 1:
            for row in rows:
                row['eligible'] = 'no'
                row['notes'] = ';'.join(filter(None, [row['notes'], 'ambiguous_gene_ids_review']))
            continue
        winner = min(rows, key=lambda r: (-r['CDS_length'], r['transcript_accession'], r['original_header']))
        winner['selected'] = 'yes'
        selections.append(winner)
    table(directory / 'candidate_cds.tsv', candidates)
    table(directory / 'selection.tsv', selections)
    (directory / 'extraction_summary.json').write_text(json.dumps(dict(
        downloaded_cds_count=len(candidates), candidate_species_count=len({r['species'] for r in candidates}),
        species_selected=len(selections)), indent=2) + '\n')
    if not selections:
        raise ValueError('No eligible CDS; inspect candidate_cds.tsv')
    with output.open('w') as handle:
        for row in selections:
            handle.write(f'>{row["label"]}\n{row["sequence"]}\n')
    print(f'Selected {len(selections)} species from {len(candidates)} downloaded CDS records')


if __name__ == '__main__':
    try:
        mode, directory, gene = sys.argv[1:]
        {'download': download, 'download-list': download_list, 'extract': extract}[mode](directory, gene)
    except Exception as error:
        sys.exit(f'Error: {error}')
