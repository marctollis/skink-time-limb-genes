"""Offline checks; mock tools test orchestration, not scientific algorithms."""
import csv
import importlib.util
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('pipeline', ROOT / 'scripts/pipeline.py')
pipeline = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pipeline)


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_longest_tie_exclusions_and_ambiguity(self):
        rows = [('One species', 'One_species', 'downloaded'),
                ('Two species', 'Two_species', 'downloaded'),
                ('Three species', 'Three_species', 'download_failed')]
        with (self.root / 'download_status.tsv').open('w') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['species', 'label', 'status'])
            writer.writerows(rows)
        for _, label, _ in rows[:2]:
            package = self.root / 'downloads' / label / 'package/ncbi_dataset/data'
            package.mkdir(parents=True)
            (package / 'cds.fna').write_text(
                '>NM_2.1 [gene=SHH]\nATG\n'
                '>XM_3.1 [gene=SHH]\nATGATG\n'
                '>NM_1.1 [gene=SHH]\nATGATG\n'
                '>NM_4.1 [gene=OTHER]\nATGATGATG\n'
                '>ABC123\nATGATGATG\n'
                '>NM_5.1 [gene=SHH]\nINVALID\n')
            ids = [1] if label == 'One_species' else [1, 2]
            (package / 'data_report.jsonl').write_text(
                ''.join(json.dumps({'geneId': x}) + '\n' for x in ids))
        with redirect_stdout(io.StringIO()):
            pipeline.extract(self.root, 'SHH')
        self.assertEqual((self.root / 'SHH_longest_refseq.fasta').read_text(),
                         '>One_species\nATGATG\n')
        with (self.root / 'selection.tsv').open() as handle:
            report = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual(report[0]['accession'], 'NM_1.1')
        self.assertEqual(report[0]['candidate_count'], '3')
        self.assertIn('ambiguous_gene_ids', report[1]['notes'])
        self.assertEqual(report[2]['notes'], 'download_failed')

    def test_alignment_validation(self):
        path = self.root / 'input.fasta'
        for contents in ['>a\nATG\n>a\nATG\n', '>a\nAT!\n>b\nATG\n',
                         '>a\nATG\n>b\nATGATG\n', '>a\nAT\n>b\nAT\n',
                         '>a\nNNN\n>b\nATG\n']:
            path.write_text(contents)
            with self.assertRaises(ValueError):
                pipeline.check(path, 2, True)
        path.write_text('>a\nATGNNN\n>b\nATG---\n')
        with redirect_stdout(io.StringIO()):
            pipeline.check(path, 2, True)

    def test_species_input(self):
        path = self.root / 'species.txt'
        path.write_text('# comment\nAnolis carolinensis # note\nSphenodon punctatus\n')
        output = io.StringIO()
        with redirect_stdout(output):
            pipeline.species(path)
        self.assertIn('Anolis carolinensis\tAnolis_carolinensis', output.getvalue())
        path.write_text('Anolis carolinensis\nAnolis carolinensis\n')
        with self.assertRaises(ValueError):
            pipeline.species(path)

    def test_four_stages_and_overwrite_guards(self):
        repo = self.root / 'repo with spaces'
        shutil.copytree(ROOT / 'scripts', repo / 'scripts')
        (repo / 'species.txt').write_text(
            'Anolis carolinensis\nSphenodon punctatus\nPogona vitticeps\nPython bivittatus\n')
        bindir = self.root / 'bin'
        bindir.mkdir()
        (bindir / 'python').symlink_to(sys.executable)
        mock = '''#!/usr/bin/env python3
import json, pathlib, sys, zipfile
tool = pathlib.Path(sys.argv[0]).name
args = sys.argv[1:]
def value(key): return args[args.index(key) + 1]
if tool == 'datasets':
    assert args[:3] == ['download', 'gene', 'symbol']
    assert value('--include') == 'cds'
    with zipfile.ZipFile(value('--filename'), 'w') as z:
        z.writestr('ncbi_dataset/data/cds.fna', '>NM_1.1 [gene=SHH]\\nATGATGATG\\n')
        z.writestr('ncbi_dataset/data/data_report.jsonl', json.dumps({'geneId': 1})+'\\n')
elif tool == 'macse':
    source = value('-seq') if value('-prog') == 'alignSequences' else value('-align')
    pathlib.Path(value('-out_NT')).write_text(pathlib.Path(source).read_text())
    if '-out_AA' in args: pathlib.Path(value('-out_AA')).write_text('>mock\\nMMM\\n')
elif tool == 'iqtree2':
    assert value('-T') == '2' and value('-st') == 'DNA'
    assert value('-B') == '1000' and value('-seed') == '12345'
    pathlib.Path(value('-pre')+'.treefile').write_text('(Anolis_carolinensis,Python_bivittatus);\\n')
elif tool == 'seqkit':
    assert args[0] == 'stats'
'''
        for tool in ('datasets', 'macse', 'iqtree2', 'seqkit'):
            path = bindir / tool
            path.write_text(mock)
            path.chmod(0o755)
        env = dict(os.environ, PATH=str(bindir) + os.pathsep + os.environ['PATH'],
                   SLURM_CPUS_PER_TASK='2')
        stages = ['01_download_cds.sh', '02_extract_longest_refseq.sh',
                  '03_macse.sh', '04_build_tree.sh']
        for stage in stages:
            command = ['bash', str(repo / 'scripts' / stage)]
            result = subprocess.run(command, input='shh\n', text=True,
                                    capture_output=True, cwd=self.root, env=env)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            again = subprocess.run(command, input='SHH\n', text=True,
                                   capture_output=True, cwd=self.root, env=env)
            self.assertNotEqual(again.returncode, 0)
            self.assertIn('already exists', again.stderr)
        self.assertTrue((repo / 'genes/SHH/iqtree/SHH.treefile').exists())
        result = subprocess.run(['bash', str(repo / 'scripts' / stages[0])],
                                input='../escape\n', text=True, capture_output=True, env=env)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse((repo / 'escape').exists())


if __name__ == '__main__':
    unittest.main()
