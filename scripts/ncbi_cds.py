#!/usr/bin/env python3
"""Retrieve annotated RefSeq CDS across Lepidosauria for a gene list."""
import csv
import json
import os
import re
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from Bio import Entrez, SeqIO
from Bio.SeqFeature import ExactPosition


def read_records(path):
    with path.open() as handle:
        yield from SeqIO.parse(handle, 'genbank')


def names_for(gene):
    # Batch genes must not accidentally share aliases belonging to another gene.
    names = [gene]
    names = sorted({x.strip().upper() for x in names if x.strip()})
    if not all(re.fullmatch(r'[A-Z0-9][A-Z0-9_.-]*', x) for x in names):
        raise ValueError('Gene aliases must be comma-separated symbols, not search expressions')
    return names


def query_for(names):
    # Organism expansion includes descendant species; RefSeq includes predictions.
    terms = ' OR '.join(f'"{name}"[Gene Name]' for name in names)
    return f'({terms}) AND Lepidosauria[Organism] AND srcdb_refseq[PROP]'


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


def download_list(root, path):
    root = Path(root)
    genes = gene_list(path)  # Validate the entire list before any network calls.
    if not os.environ.get('NCBI_EMAIL', '').strip():
        raise ValueError('Set NCBI_EMAIL to your contact email')
    (root / 'genes').mkdir(parents=True, exist_ok=True)
    statuses = []
    failed = False
    for gene in genes:
        directory = root / 'genes' / gene
        directory.mkdir(exist_ok=True)
        print(f'\nDownloading RefSeq Lepidosauria CDS records for {gene}', flush=True)
        try:
            manifest_path = directory / 'download_manifest.json'
            existing = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
            if (existing.get('complete') and existing.get('scope') == 'refseq_lepidosauria'
                    and existing.get('gene') == gene
                    and all((directory / 'downloads' / b['file']).is_file() for b in existing['batches'])):
                status, note = 'already_complete', 'Existing matching download retained'
                print(f'{gene}: already downloaded; skipping', flush=True)
            else:
                download(directory, gene)
                status, note = 'downloaded', ''
        except Exception as error:
            status, note, failed = 'failed', str(error), True
            print(f'{gene}: {note}', file=sys.stderr, flush=True)
        statuses.append({'gene': gene, 'status': status, 'notes': note})
        with (root / 'genes/download_batch_status.tsv').open('w') as handle:
            writer = csv.DictWriter(handle, fieldnames=['gene', 'status', 'notes'], delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(statuses)
    if failed:
        raise ValueError('Some genes failed; inspect genes/download_batch_status.tsv. Successful genes were retained.')


def download(directory, gene):
    directory = Path(directory)
    if (directory / 'downloads').exists():
        raise ValueError('downloads already exists; move the old gene directory aside before restarting')
    email = os.environ.get('NCBI_EMAIL', '').strip()
    if not re.fullmatch(r'[^\s@]+@[^\s@]+\.[^\s@]+', email):
        raise ValueError('Set NCBI_EMAIL to your contact email')
    Entrez.email, Entrez.tool = email, 'skink-time-limb-genes'
    Entrez.api_key = os.environ.get('NCBI_API_KEY') or None
    Entrez.max_tries, Entrez.sleep_between_tries = 3, 5
    names = names_for(gene)
    query = query_for(names)
    with Entrez.esearch(db='nuccore', term=query, usehistory='y', retmax=0) as handle:
        result = Entrez.read(handle)
    if result.get('ErrorList'):
        raise ValueError(f'NCBI search error: {result["ErrorList"]}')
    count = int(result['Count'])
    if not count:
        raise ValueError('No matching RefSeq Lepidosauria records; check the gene symbol')
    dest = directory / 'downloads'
    dest.mkdir(parents=True)
    manifest = {'gene': gene, 'gene_names': names, 'scope': 'refseq_lepidosauria', 'query': query,
                'query_translation': result.get('QueryTranslation', ''),
                'database': 'nuccore', 'search_record_count': count,
                'retrieved_record_count': 0, 'batches': [], 'complete': False,
                'retrieved_utc': datetime.now(timezone.utc).isoformat()}
    manifest_path = directory / 'download_manifest.json'
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    print(f'Found {count} records; retrieving all in batches of 20.', flush=True)
    for start in range(0, count, 20):
        expected = min(20, count - start)
        target = dest / f'batch_{start:08d}.gb'
        temporary = target.with_suffix('.part')
        with Entrez.efetch(db='nuccore', query_key=result['QueryKey'],
                          WebEnv=result['WebEnv'], rettype='gbwithparts',
                          retmode='text', retstart=start, retmax=expected) as handle:
            # Stream, since annotated chromosomes can be large.
            with temporary.open('w') as out:
                while True:
                    chunk = handle.read(1024 * 1024)
                    if not chunk:
                        break
                    out.write(chunk)
        actual = sum(1 for _ in read_records(temporary))
        if actual != expected:
            raise ValueError(f'Incomplete batch at {start}: expected {expected}, got {actual}; extraction is blocked')
        temporary.rename(target)
        manifest['batches'].append({'file': target.name, 'record_count': actual})
        manifest['retrieved_record_count'] += actual
        manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
        print(f'Retrieved {manifest["retrieved_record_count"]}/{count}', flush=True)
    manifest['complete'] = True
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')


FIELDS = ['species', 'organism', 'label', 'taxid', 'source', 'nucleotide_accession',
          'protein_accession', 'gene_annotation', 'gene_ids', 'location', 'CDS_length',
          'eligible', 'selected', 'notes']


def table(path, rows):
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows({key: row.get(key, '') for key in FIELDS} for row in rows)


def feature_names(feature):
    names = list(feature.qualifiers.get('gene', []))
    for entry in feature.qualifiers.get('gene_synonym', []):
        names.extend(re.split(r'[;,]', entry))
    return {name.strip().upper() for name in names if name.strip()}


def extract(directory, gene):
    directory = Path(directory)
    output = directory / f'{gene}_longest_cds.fasta'
    if output.exists() or (directory / 'candidate_cds.tsv').exists():
        raise ValueError('Extraction output already exists; move the old gene directory aside before restarting')
    manifest = json.loads((directory / 'download_manifest.json').read_text())
    if not manifest.get('complete') or manifest.get('gene') != gene:
        raise ValueError('Download is incomplete or belongs to another gene; rerun step 01 in a fresh gene directory')
    if manifest.get('scope') != 'refseq_lepidosauria':
        raise ValueError('Old retrieval scope: move the gene directory aside and run the new RefSeq downloader')
    names = set(manifest['gene_names'])
    candidates, groups, seen, unmatched = [], defaultdict(list), set(), 0
    total = 0
    for batch in manifest['batches']:
        path = directory / 'downloads' / batch['file']
        actual = 0
        for record in read_records(path):
            actual += 1
            organism = record.annotations.get('organism', '')
            words = organism.split()
            species = ' '.join(words[:2])
            valid_species = len(words) >= 2 and bool(re.fullmatch(r'[A-Z][a-z]+ [a-z]+', species))
            valid_species = valid_species and words[1] not in {'sp', 'cf', 'aff', 'hybrid', 'unidentified'}
            taxonomy = record.annotations.get('taxonomy', [])
            valid_taxon = 'Lepidosauria' in taxonomy or 'Squamata' in taxonomy or (words and words[0] == 'Sphenodon')
            label = species.replace(' ', '_')
            sources = [f for f in record.features if f.type == 'source']
            taxids = sorted({x.split(':', 1)[1] for f in sources
                             for x in f.qualifiers.get('db_xref', []) if x.startswith('taxon:')})
            matched = False
            for feature in record.features:
                if feature.type != 'CDS' or not feature_names(feature) & names:
                    continue
                matched = True
                location = str(feature.location)
                identity = (record.id, location)
                if identity in seen:
                    continue
                seen.add(identity)
                refs = feature.qualifiers.get('db_xref', [])
                geneids = sorted({x.split(':', 1)[1] for x in refs if x.startswith('GeneID:')})
                row = dict(species=species, organism=organism, label=label,
                           taxid=','.join(taxids), source='RefSeq' if re.match(r'^[A-Z]{2}_', record.id) else 'GenBank/INSDC',
                           nucleotide_accession=record.id,
                           protein_accession=','.join(feature.qualifiers.get('protein_id', [])),
                           gene_annotation=','.join(feature.qualifiers.get('gene', [])),
                           gene_ids=','.join(geneids), location=location,
                           eligible='no', selected='no', CDS_length='', notes='', sequence='')
                problems, warnings = [], []
                if row['source'] != 'RefSeq':
                    problems.append('not_RefSeq')
                if not valid_species:
                    problems.append('unresolved_species')
                if not valid_taxon:
                    problems.append('outside_target_taxonomy')
                if 'pseudo' in feature.qualifiers or 'pseudogene' in feature.qualifiers:
                    problems.append('pseudogene')
                if feature.location is None:
                    problems.append('missing_location')
                else:
                    if any(part.ref for part in feature.location.parts):
                        problems.append('remote_location_needs_manual_extraction')
                    if any(not isinstance(end, ExactPosition) for part in feature.location.parts
                           for end in (part.start, part.end)) or 'partial' in feature.qualifiers:
                        problems.append('partial_CDS')
                if feature.qualifiers.get('codon_start', ['1']) != ['1']:
                    problems.append('partial_reading_frame')
                if 'transl_except' in feature.qualifiers or 'exception' in feature.qualifiers:
                    problems.append('translation_exception_review')
                if not problems:
                    try:
                        seq = feature.extract(record.seq)
                        sequence = str(seq).upper()
                        row['CDS_length'] = len(sequence)
                        row['sequence'] = sequence
                        if not re.fullmatch(r'[ACGTRYSWKMBDHVN]+', sequence) or not set(sequence) & set('ACGT'):
                            problems.append('invalid_or_uninformative_sequence')
                        elif len(sequence) % 3:
                            problems.append('length_not_multiple_of_3')
                        else:
                            code = int(feature.qualifiers.get('transl_table', ['1'])[0])
                            if code != 1:
                                problems.append('nonstandard_genetic_code')
                            protein = str(seq.translate(table=code))
                            if '*' in protein[:-1]:
                                problems.append('internal_stop')
                            if not protein.endswith('*'):
                                warnings.append('no_terminal_stop_review')
                            if sequence[:3] != 'ATG':
                                warnings.append('non_ATG_start_review')
                            if re.search(r'[^ACGT]', sequence):
                                warnings.append('ambiguous_bases_review')
                    except (ValueError, TypeError) as error:
                        problems.append(f'extraction_error:{error}')
                row['notes'] = ';'.join(problems + warnings)
                row['eligible'] = 'no' if problems else 'yes'
                candidates.append(row)
                if not problems:
                    groups[species].append(row)
            if not matched:
                unmatched += 1
        if actual != batch['record_count']:
            raise ValueError(f'{path.name}: record count differs from the download manifest')
        total += actual
    if total != manifest['retrieved_record_count']:
        raise ValueError('Total retrieved records differ from manifest')
    selections = []
    for species, rows in sorted(groups.items()):
        # Longest first; prefer RefSeq only when lengths tie, then accession/location.
        winner = min(rows, key=lambda r: (-int(r['CDS_length']), r['source'] != 'RefSeq',
                                         r['nucleotide_accession'], r['location']))
        winner['selected'] = 'yes'
        geneids = {x for row in rows for x in row['gene_ids'].split(',') if x}
        if len(geneids) > 1:
            winner['notes'] = ';'.join(filter(None, [winner['notes'], 'multiple_gene_ids_review']))
        selections.append(winner)
    table(directory / 'candidate_cds.tsv', candidates)
    table(directory / 'selection.tsv', selections)
    (directory / 'extraction_summary.json').write_text(json.dumps({
        'retrieved_records': total, 'records_without_matching_CDS': unmatched,
        'matching_CDS_candidates': len(candidates), 'species_selected': len(selections),
        'selected_RefSeq': sum(r['source'] == 'RefSeq' for r in selections),
        'selected_GenBank_INSDC': sum(r['source'] != 'RefSeq' for r in selections),
    }, indent=2) + '\n')
    if not selections:
        raise ValueError('No eligible CDS; inspect candidate_cds.tsv and extraction_summary.json')
    with output.open('w') as handle:
        for row in selections:
            handle.write(f'>{row["label"]}\n{row["sequence"]}\n')
    print(f'Selected {len(selections)} species from {len(candidates)} matching CDS candidates')


if __name__ == '__main__':
    try:
        mode, directory, gene = sys.argv[1:]
        {'download': download, 'download-list': download_list, 'extract': extract}[mode](directory, gene)
    except Exception as error:
        sys.exit(f'Error: {error}')
