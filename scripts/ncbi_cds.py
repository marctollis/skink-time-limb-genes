#!/usr/bin/env python3
"""Datasets gene-list downloads and longest RefSeq CDS per species."""
import csv
import hashlib
import json
import os
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
import re
import subprocess
import sys
import zipfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from Bio.Seq import Seq
from pipeline import fasta

SCOPE = 'datasets_refseq_squamata_gene_ids_v4'


def normalize_report(value):
    """Older CLI reports use snake_case; newer documentation uses camelCase."""
    if isinstance(value, list):
        return [normalize_report(item) for item in value]
    if isinstance(value, dict):
        result = {}
        for key, item in value.items():
            parts = key.split('_')
            normalized = parts[0] + ''.join(part[:1].upper() + part[1:] for part in parts[1:])
            result[normalized] = normalize_report(item)
        return result
    return value


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
    ids, query = discover_gene_ids(gene)
    ids_path = dest / 'gene_ids.txt'
    ids_path.write_text(''.join(identifier + '\n' for identifier in ids))
    command = ['datasets', 'download', 'gene', 'gene-id', '--inputfile', str(ids_path),
               '--include', 'cds,product-report', '--filename', str(archive)]
    manifest = dict(gene=gene, scope=SCOPE, complete=False, command=command,
                    discovery_query=query, gene_id_count=len(ids),
                    include_sphenodon=include_outgroup(),
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
    reports = [normalize_report(json.loads(line)) for line in report.read_text().splitlines() if line.strip()]
    if not reports:
        raise ValueError('Empty gene report')
    manifest.update(complete=True, downloaded_cds_count=count, gene_report_count=len(reports),
                    cds_sha256=hashlib.sha256(cds.read_bytes()).hexdigest(),
                    report_sha256=hashlib.sha256(report.read_bytes()).hexdigest())
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    print(f'{gene}: downloaded {count} CDS records', flush=True)


def search_gene(params):
    params = dict(params, db='gene', tool='skink-time-limb-genes')
    if os.environ.get('NCBI_EMAIL'):
        params['email'] = os.environ['NCBI_EMAIL']
    if os.environ.get('NCBI_API_KEY'):
        params['api_key'] = os.environ['NCBI_API_KEY']
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?' + urllib.parse.urlencode(params)
    for attempt in range(3):
        time.sleep(0.4)
        try:
            request = urllib.request.Request(url)
            with urllib.request.urlopen(request, timeout=30) as handle:
                root = ET.fromstring(handle.read())
            errors = root.findall('.//ERROR') + root.findall('.//ErrorList/*')
            if errors:
                raise ValueError('NCBI Gene search error: ' + '; '.join(e.text or '' for e in errors))
            return root
        except (OSError, ET.ParseError):
            if attempt == 2:
                raise ValueError('NCBI Gene discovery failed after three attempts') from None
            time.sleep(2 * (attempt + 1))


def include_outgroup():
    value = os.environ.get('SKINK_INCLUDE_OUTGROUP', '0')
    if value not in ('0', '1'):
        raise ValueError('SKINK_INCLUDE_OUTGROUP must be 0 or 1')
    return value == '1'


def discover_gene_ids(gene):
    # Search genes directly across the clade; no species/taxonomy enumeration.
    # RefSeq CDS identity is enforced from package accessions during step 02.
    taxa = '(txid8509[Organism:exp] OR Sphenodon[Organism])' if include_outgroup() else 'txid8509[Organism:exp]'
    query = f'"{gene}"[Gene Name] AND {taxa}'
    first = search_gene(dict(term=query, retmax=1000, usehistory='y'))
    count = int(first.findtext('Count', default='0'))
    if not count:
        raise ValueError(f'No matching NCBI Gene records for {gene} in Squamata/outgroup scope')
    ids = [node.text for node in first.findall('./IdList/Id')]
    for start in range(1000, count, 1000):
        key, history = first.findtext('QueryKey'), first.findtext('WebEnv')
        if not key or not history:
            raise ValueError('Missing NCBI history needed to retrieve all Gene IDs')
        page = search_gene(dict(term='#' + key, WebEnv=history, retstart=start, retmax=1000))
        ids.extend(node.text for node in page.findall('./IdList/Id'))
    if len(ids) != count or len(set(ids)) != count or not all(re.fullmatch(r'\d+', x or '') for x in ids):
        raise ValueError('Incomplete or invalid Gene ID discovery; download stopped')
    print(f'{gene}: found {len(ids)} matching Squamata/outgroup Gene IDs', flush=True)
    return ids, query


def completed(directory, gene):
    path = directory / 'download_manifest.json'
    if not path.exists():
        return False
    manifest = json.loads(path.read_text())
    cds, report = package_files(directory)
    return (manifest.get('complete') and manifest.get('scope') == SCOPE
            and manifest.get('gene') == gene and manifest.get('include_sphenodon') == include_outgroup() and cds.is_file() and report.is_file()
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
        print(f'\nDownloading {gene}: RefSeq squamate CDS', flush=True)
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
    reports = [normalize_report(json.loads(line)) for line in report.read_text().splitlines() if line.strip()]
    metadata = {str(row['geneId']): row for row in reports}
    candidates, groups = [], defaultdict(list)
    seen = set()
    for header, sequence in fasta(cds):
        if (header, sequence) in seen:
            continue
        seen.add((header, sequence))
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
