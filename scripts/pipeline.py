#!/usr/bin/env python3
"""Small standard-library helpers for the interactive training scripts."""
import re
import sys


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
        if mode == 'check':
            check(args[0], int(args[1]), len(args) == 3 and args[2] == 'aligned')
        else:
            raise ValueError(f'Unknown mode: {mode}')
    except (ValueError, OSError, KeyError) as error:
        sys.exit(f'Error: {error}')
