#!/usr/bin/env python3
"""Small standard-library helpers for the interactive training scripts."""
import csv
import json
import re
import sys
from pathlib import Path


def fasta(path):
    header, parts = None, []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(parts).upper()
                header, parts = line[1:], []
            else:
                if header is None:
                    raise ValueError(f'{path}: sequence before first FASTA header')
                parts.append(line)
    if header is not None:
        yield header, ''.join(parts).upper()


def species(path):
    rows, labels = [], set()
    for line in Path(path).read_text().splitlines():
        name = ' '.join(line.split('#', 1)[0].split())
        if not name:
            continue
        if not re.fullmatch(r'[A-Za-z]+ [A-Za-z]+', name):
            raise ValueError(f'Expected a scientific binomial, got: {name!r}')
        label = name.replace(' ', '_')
        if label in labels:
            raise ValueError(f'Duplicate species: {name}')
        labels.add(label)
        rows.append((name, label))
    if not rows:
        raise ValueError('No species in species.txt')
    for name, label in rows:
        print(f'{name}\t{label}')


def extract(directory, gene):
    directory = Path(directory)
    selected, report = [], []
    with (directory / 'download_status.tsv').open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    for row in rows:
        name, label = row['species'], row['label']
        candidates, reason = [], row['status']
        if row['status'] == 'downloaded':
            package = directory / 'downloads' / label / 'package'
            gene_ids = set()
            for path in package.rglob('data_report.jsonl'):
                for line in path.read_text().splitlines():
                    if line.strip():
                        record = json.loads(line)
                        if record.get('geneId') is not None:
                            gene_ids.add(str(record['geneId']))
            # Symbol queries also match synonyms. Do not mix different genes.
            if len(gene_ids) > 1:
                reason = 'ambiguous_gene_ids:resolve_package_manually'
            else:
                for path in sorted(package.rglob('cds.fna')):
                    for header, seq in fasta(path):
                        accession = re.search(r'(?:NM|XM|NP|XP)_\d+(?:\.\d+)?', header)
                        symbol = re.search(r'\[gene=([^\]]+)\]', header)
                        if not accession or (symbol and symbol[1].upper() != gene):
                            continue
                        if not seq or not re.fullmatch(r'[ACGTRYSWKMBDHVN]+', seq):
                            continue
                        candidates.append((header, seq, accession[0]))
                reason = 'no_matching_refseq_cds' if not candidates else 'selected'
        if candidates:
            # Longest across curated and predicted RefSeqs; header breaks ties.
            header, seq, accession = sorted(candidates, key=lambda x: (-len(x[1]), x[0]))[0]
            selected.append((label, seq))
            report.append((gene, name, label, len(candidates), accession, len(seq),
                           'yes', 'length_not_multiple_of_3' if len(seq) % 3 else '', header))
        else:
            report.append((gene, name, label, 0, '', '', 'no', reason, ''))
    with (directory / 'selection.tsv').open('w') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['gene', 'species', 'label', 'candidate_count', 'accession',
                         'CDS_length', 'included', 'notes', 'original_header'])
        writer.writerows(report)
    if not selected:
        raise ValueError('No usable RefSeq CDS; inspect selection.tsv and download logs')
    with (directory / f'{gene}_longest_refseq.fasta').open('w') as handle:
        for label, seq in selected:
            handle.write(f'>{label}\n{seq}\n')
    print(f'Selected {len(selected)} of {len(rows)} requested species')


def check(path, minimum, aligned=False):
    records = list(fasta(path))
    if len(records) < minimum:
        raise ValueError(f'Need at least {minimum} sequences; found {len(records)}')
    names = [h.split()[0] for h, _ in records]
    if len(names) != len(set(names)):
        raise ValueError('Duplicate sequence labels')
    alphabet = r'[ACGTRYSWKMBDHVN-]+' if aligned else r'[ACGTRYSWKMBDHVN]+'
    for header, seq in records:
        if not re.fullmatch(alphabet, seq) or not set(seq) & set('ACGT'):
            raise ValueError(f'Invalid or uninformative sequence: {header}')
    lengths = {len(seq) for _, seq in records}
    if aligned and (len(lengths) != 1 or next(iter(lengths)) % 3):
        raise ValueError('Expected equal alignment lengths divisible by three')
    print(f'Validated {len(records)} sequences; lengths: {sorted(lengths)}')


if __name__ == '__main__':
    try:
        mode, *args = sys.argv[1:]
        if mode == 'species':
            species(*args)
        elif mode == 'extract':
            extract(*args)
        elif mode == 'check':
            check(args[0], int(args[1]), len(args) == 3 and args[2] == 'aligned')
        else:
            raise ValueError(f'Unknown mode: {mode}')
    except (ValueError, OSError, KeyError) as error:
        sys.exit(f'Error: {error}')
