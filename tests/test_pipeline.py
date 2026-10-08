"""Offline Datasets package and workflow checks; no live downloads."""
import csv
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
import zipfile
from contextlib import redirect_stdout, redirect_stderr
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import ncbi_cds as ncbi
import pipeline


def content(records=None):
    records = records or [('NM_1.1', 'Anolis carolinensis', '1', 'SHH', 'ATGAAATAA')]
    fasta, reports = [], {}
    for accession, species, geneid, symbol, seq in records:
        fasta.append(f'>{accession}:1-{len(seq)} {symbol} [organism={species}] [GeneID={geneid}] [region=cds]\n{seq}\n')
        reports[geneid] = dict(geneId=geneid, symbol=symbol, taxId='123', taxname=species, type='PROTEIN_CODING')
    return ''.join(fasta), ''.join(json.dumps(r) + '\n' for r in reports.values())


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def download_fixture(self, records=None, directory=None):
        directory = directory or self.root
        text, report = content(records)
        def fake_run(command, **kwargs):
            archive = Path(command[command.index('--filename') + 1])
            with zipfile.ZipFile(archive, 'w') as out:
                out.writestr('ncbi_dataset/data/cds.fna', text)
                out.writestr('ncbi_dataset/data/data_report.jsonl', report)
            return subprocess.CompletedProcess(command, 0)
        with patch.object(ncbi.subprocess, 'run', side_effect=fake_run) as run, redirect_stdout(io.StringIO()):
            ncbi.download(directory, 'SHH')
        self.assertEqual(run.call_args.args[0][0:5], ['datasets', 'download', 'gene', 'symbol', 'SHH'])
        self.assertEqual(run.call_args.args[0][6], 'Lepidosauria')
        return directory

    def test_package_download_and_longest(self):
        self.download_fixture([('NM_1.1', 'Anolis carolinensis', '1', 'SHH', 'ATGAAATAA'),
                               ('XM_2.1', 'Anolis carolinensis', '1', 'SHH', 'ATGAAAAAATAA'),
                               ('NM_3.1', 'Sphenodon punctatus', '2', 'SHH', 'ATGAAATAA')])
        with redirect_stdout(io.StringIO()):
            ncbi.extract(self.root, 'SHH')
        with (self.root / 'selection.tsv').open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]['transcript_accession'], 'XM_2.1')
        self.assertEqual(rows[1]['species'], 'Sphenodon punctatus')

    def test_exclusions_and_ambiguous_gene_ids(self):
        self.download_fixture([('AB1.1', 'Anolis carolinensis', '1', 'SHH', 'ATGAAATAA'),
                               ('NM_2.1', 'Pogona vitticeps', '2', 'SHH', 'ATGTAATAA'),
                               ('NM_3.1', 'Python bivittatus', '3', 'SHH', 'ATGAAATAA'),
                               ('NM_4.1', 'Python bivittatus', '4', 'SHH', 'ATGAAATAA'),
                               ('NM_5.1', 'Sphenodon punctatus', '5', 'SHH', 'ATGAAATAA')])
        with redirect_stdout(io.StringIO()):
            ncbi.extract(self.root, 'SHH')
        with (self.root / 'candidate_cds.tsv').open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual([r['eligible'] for r in rows], ['no', 'no', 'no', 'no', 'yes'])
        self.assertIn('ambiguous_gene_ids', rows[2]['notes'])

    def test_modified_or_old_package_blocked(self):
        self.download_fixture()
        cds, _ = ncbi.package_files(self.root)
        cds.write_text(cds.read_text() + '\n')
        self.assertFalse(ncbi.completed(self.root, 'SHH'))
        with self.assertRaisesRegex(ValueError, 'modified or old-scope'):
            ncbi.extract(self.root, 'SHH')

    def test_zero_sequences_not_complete(self):
        def fake_run(command, **kwargs):
            with zipfile.ZipFile(command[command.index('--filename') + 1], 'w') as out:
                out.writestr('ncbi_dataset/data/cds.fna', '')
                out.writestr('ncbi_dataset/data/data_report.jsonl', '{}\n')
            return subprocess.CompletedProcess(command, 0)
        with patch.object(ncbi.subprocess, 'run', side_effect=fake_run):
            with self.assertRaisesRegex(ValueError, 'no CDS'):
                ncbi.download(self.root, 'SHH')
        self.assertFalse(json.loads((self.root / 'download_manifest.json').read_text())['complete'])

    def test_list_validation_and_continue_after_failure(self):
        path = self.root / 'genes.txt'
        path.write_text('shh\nFGF8 # comment\n')
        self.assertEqual(ncbi.gene_list(path), ['SHH', 'FGF8'])
        with patch.object(ncbi, 'download', side_effect=[ValueError('no results'), None]) as run, \
             redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
            with self.assertRaisesRegex(ValueError, 'Some genes failed'):
                ncbi.download_list(self.root, path)
        self.assertEqual(run.call_count, 2)
        self.assertIn('FGF8\tdownloaded', (self.root / 'genes/download_batch_status.tsv').read_text())
        for text in ['SHH\nshh\n', '../escape\n', '# empty\n']:
            path.write_text(text)
            with self.assertRaises(ValueError):
                ncbi.gene_list(path)

    def test_completed_skip(self):
        path = self.root / 'genes.txt'
        path.write_text('SHH\n')
        self.download_fixture(directory=self.root / 'genes/SHH')
        with patch.object(ncbi, 'download') as run, redirect_stdout(io.StringIO()):
            ncbi.download_list(self.root, path)
        run.assert_not_called()
        self.assertIn('already_complete', (self.root / 'genes/download_batch_status.tsv').read_text())

    def test_alignment_validation(self):
        path = self.root / 'input.fasta'
        for text in ['>a\nATG\n>a\nATG\n', '>a\nAT!\n>b\nATG\n', '>a\nATG\n>b\nATGATG\n']:
            path.write_text(text)
            with self.assertRaises(ValueError):
                pipeline.check(path, 2, True)

    def test_shell_stages(self):
        repo = self.root / 'repo with spaces'
        shutil.copytree(ROOT / 'scripts', repo / 'scripts')
        self.download_fixture([(f'NM_{i}.1', species, str(i), 'SHH', 'ATGAAATAA') for i, species in enumerate(
            ['Anolis carolinensis', 'Pogona vitticeps', 'Python bivittatus', 'Sphenodon punctatus'], 1)], repo / 'genes/SHH')
        bindir = self.root / 'bin'
        bindir.mkdir()
        (bindir / 'python').symlink_to(sys.executable)
        mock = '''#!/usr/bin/env python3
import pathlib, sys
args = sys.argv[1:]
def value(key): return args[args.index(key) + 1]
tool = pathlib.Path(sys.argv[0]).name
if tool == 'macse':
    source = value('-seq') if value('-prog') == 'alignSequences' else value('-align')
    pathlib.Path(value('-out_NT')).write_text(pathlib.Path(source).read_text())
    if '-out_AA' in args: pathlib.Path(value('-out_AA')).write_text('>mock\\nMK*\\n')
elif tool == 'iqtree2':
    assert value('-T') == '2'
    pathlib.Path(value('-pre')+'.treefile').write_text('(a,b);\\n')
'''
        for tool in ('macse', 'iqtree2', 'seqkit'):
            path = bindir / tool
            path.write_text(mock)
            path.chmod(0o755)
        env = dict(os.environ, PATH=str(bindir) + os.pathsep + os.environ['PATH'], SLURM_CPUS_PER_TASK='2')
        for stage in ['02_extract_longest_refseq.sh', '03_macse.sh', '04_build_tree.sh']:
            command = ['bash', str(repo / 'scripts' / stage)]
            result = subprocess.run(command, input='shh\n', text=True, capture_output=True, cwd=self.root, env=env)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            again = subprocess.run(command, input='SHH\n', text=True, capture_output=True, cwd=self.root, env=env)
            self.assertNotEqual(again.returncode, 0)
            self.assertIn('already exists', again.stderr)


if __name__ == '__main__':
    unittest.main()
