"""Offline checks: synthetic annotated records and mock external tools."""
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
from unittest.mock import patch

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation, BeforePosition
from Bio.SeqRecord import SeqRecord

ROOT = Path(__file__).resolve().parents[1]

def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'scripts' / f'{name}.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

ncbi = load('ncbi_cds')
pipeline = load('pipeline')


def record(accession, organism='Anolis carolinensis', seq='ATGAAATAA', gene='SHH', location=None, **qualifiers):
    rec = SeqRecord(Seq(seq), id=accession, name=accession.split('.')[0], description='synthetic test')
    rec.annotations = dict(molecule_type='DNA', organism=organism,
                           taxonomy=['Eukaryota', 'Sauropsida', 'Squamata'])
    rec.features = [SeqFeature(SimpleLocation(0, len(seq)), type='source',
                               qualifiers={'organism': [organism], 'db_xref': ['taxon:123']}),
                    SeqFeature(location if location is not None else SimpleLocation(0, len(seq)), type='CDS',
                               qualifiers=dict(gene=[gene], protein_id=['TEST123.1'], **qualifiers))]
    return rec


def genbank(records):
    out = io.StringIO()
    SeqIO.write(records, out, 'genbank')
    return out.getvalue()


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def fixture(self, records, directory=None):
        directory = directory or self.root
        (directory / 'downloads').mkdir(parents=True)
        (directory / 'downloads/batch_00000000.gb').write_text(genbank(records))
        (directory / 'download_manifest.json').write_text(json.dumps(dict(
            gene='SHH', gene_names=['SHH'], scope='refseq_lepidosauria', complete=True, retrieved_record_count=len(records),
            batches=[dict(file='batch_00000000.gb', record_count=len(records))])))

    def test_refseq_only_longest_and_unrelated_genes(self):
        self.fixture([record('NM_1.1'), record('NM_000001.1', seq='ATGAAAAAATAA'),
                      record('NM_000002.1', organism='Pogona vitticeps'),
                      record('NM_000003.1', gene='IHH', seq='ATGAAAAAATAA'),
                      record('AB000004.1', seq='ATGAAAAAAAAATAA')])
        with redirect_stdout(io.StringIO()):
            ncbi.extract(self.root, 'SHH')
        with (self.root / 'selection.tsv').open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]['nucleotide_accession'], 'NM_000001.1')
        self.assertEqual(rows[0]['source'], 'RefSeq')
        self.assertEqual(json.loads((self.root / 'extraction_summary.json').read_text())['records_without_matching_CDS'], 1)
        with self.assertRaisesRegex(ValueError, 'already exists'):
            ncbi.extract(self.root, 'SHH')

    def test_exon_join_and_reverse_strand(self):
        joined = record('NM_1.1', seq='ATGCCCCCCAAATAA',
                        location=CompoundLocation([SimpleLocation(0, 3), SimpleLocation(9, 15)]))
        reverse = record('NM_2.1', organism='Pogona vitticeps', seq='TTATTTCAT',
                         location=SimpleLocation(0, 9, strand=-1))
        self.fixture([joined, reverse])
        with redirect_stdout(io.StringIO()):
            ncbi.extract(self.root, 'SHH')
        sequences = list(pipeline.fasta(self.root / 'SHH_longest_cds.fasta'))
        self.assertEqual([s for _, s in sequences], ['ATGAAATAA', 'ATGAAATAA'])

    def test_exclusions_and_refseq_tie(self):
        partial = record('NM_1.1', location=SimpleLocation(BeforePosition(0), 9))
        pseudo = record('NM_2.1', pseudo=[''])
        stop = record('NM_3.1', seq='ATGTAATAA')
        unrelated_taxon = record('NM_4.1', organism='Mus musculus')
        unrelated_taxon.annotations['taxonomy'] = ['Mammalia']
        self.fixture([partial, pseudo, stop, unrelated_taxon, record('NM_5.1'), record('NM_6.1')])
        with redirect_stdout(io.StringIO()):
            ncbi.extract(self.root, 'SHH')
        with (self.root / 'candidate_cds.tsv').open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual([r['eligible'] for r in rows], ['no', 'no', 'no', 'no', 'yes', 'yes'])
        self.assertEqual([r['nucleotide_accession'] for r in rows if r['selected'] == 'yes'], ['NM_5.1'])
        self.assertIn('partial_CDS', rows[0]['notes'])

    def test_aliases_and_geneid_review(self):
        a = record('NM_1.1', gene='ALIAS', db_xref=['GeneID:1'])
        b = record('NM_2.1', db_xref=['GeneID:2'])
        self.fixture([a, b])
        manifest = json.loads((self.root / 'download_manifest.json').read_text())
        manifest['gene_names'].append('ALIAS')
        (self.root / 'download_manifest.json').write_text(json.dumps(manifest))
        with redirect_stdout(io.StringIO()):
            ncbi.extract(self.root, 'SHH')
        self.assertIn('multiple_gene_ids_review', (self.root / 'selection.tsv').read_text())

    def test_download_pagination(self):
        records = [record(f'NM_{i:06d}.1') for i in range(21)]
        result = dict(Count='21', QueryKey='1', WebEnv='mock', QueryTranslation='mock')
        with patch.dict(os.environ, {'NCBI_EMAIL': 'test@example.org', 'NCBI_GENE_ALIASES': ''}), \
             patch.object(ncbi.Entrez, 'esearch', return_value=io.StringIO('mock')) as search, \
             patch.object(ncbi.Entrez, 'read', return_value=result), \
             patch.object(ncbi.Entrez, 'efetch', side_effect=[io.StringIO(genbank(records[:20])), io.StringIO(genbank(records[20:]))]) as fetch, \
             redirect_stdout(io.StringIO()):
            ncbi.download(self.root, 'SHH')
        manifest = json.loads((self.root / 'download_manifest.json').read_text())
        self.assertTrue(manifest['complete'])
        self.assertEqual(manifest['retrieved_record_count'], 21)
        self.assertEqual([c.kwargs['retstart'] for c in fetch.call_args_list], [0, 20])
        self.assertIn('Lepidosauria[Organism]', search.call_args.kwargs['term'])
        self.assertIn('srcdb_refseq[PROP]', search.call_args.kwargs['term'])

    def test_incomplete_download_blocks_extraction(self):
        result = dict(Count='2', QueryKey='1', WebEnv='mock')
        with patch.dict(os.environ, {'NCBI_EMAIL': 'test@example.org', 'NCBI_GENE_ALIASES': ''}), \
             patch.object(ncbi.Entrez, 'esearch', return_value=io.StringIO('mock')), \
             patch.object(ncbi.Entrez, 'read', return_value=result), \
             patch.object(ncbi.Entrez, 'efetch', return_value=io.StringIO(genbank([record('NM_1.1')]))), \
             redirect_stdout(io.StringIO()):
            with self.assertRaisesRegex(ValueError, 'Incomplete batch'):
                ncbi.download(self.root, 'SHH')
        with self.assertRaisesRegex(ValueError, 'incomplete'):
            ncbi.extract(self.root, 'SHH')

    def test_gene_list_and_batch_failures(self):
        path = self.root / 'genes.txt'
        path.write_text('# genes\nshh\nFGF8 # comment\n')
        self.assertEqual(ncbi.gene_list(path), ['SHH', 'FGF8'])
        with patch.dict(os.environ, {'NCBI_EMAIL': 'test@example.org'}), \
             patch.object(ncbi, 'download', side_effect=[ValueError('No matching records'), None]) as fetch, \
             redirect_stdout(io.StringIO()):
            with self.assertRaisesRegex(ValueError, 'Some genes failed'):
                ncbi.download_list(self.root, path)
        self.assertEqual(fetch.call_count, 2)
        self.assertIn('FGF8\tdownloaded', (self.root / 'genes/download_batch_status.tsv').read_text())
        path.write_text('SHH\nshh\n')
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            ncbi.gene_list(path)
        path.write_text('../SHH\n')
        with self.assertRaisesRegex(ValueError, 'Invalid'):
            ncbi.gene_list(path)

    def test_completed_batch_skip_and_old_scope(self):
        path = self.root / 'genes.txt'
        path.write_text('SHH\n')
        directory = self.root / 'genes/SHH'
        self.fixture([record('NM_1.1')], directory)
        with patch.dict(os.environ, {'NCBI_EMAIL': 'test@example.org'}), \
             patch.object(ncbi, 'download') as fetch, redirect_stdout(io.StringIO()):
            ncbi.download_list(self.root, path)
        fetch.assert_not_called()
        self.assertIn('already_complete', (self.root / 'genes/download_batch_status.tsv').read_text())
        manifest = json.loads((directory / 'download_manifest.json').read_text())
        manifest['scope'] = 'old_genbank_scope'
        (directory / 'download_manifest.json').write_text(json.dumps(manifest))
        with self.assertRaisesRegex(ValueError, 'Old retrieval scope'):
            ncbi.extract(directory, 'SHH')

    def test_alignment_validation(self):
        path = self.root / 'input.fasta'
        for contents in ['>a\nATG\n>a\nATG\n', '>a\nAT!\n>b\nATG\n',
                         '>a\nATG\n>b\nATGATG\n', '>a\nAT\n>b\nAT\n', '>a\nNNN\n>b\nATG\n']:
            path.write_text(contents)
            with self.assertRaises(ValueError):
                pipeline.check(path, 2, True)

    def test_shell_stages_and_overwrite_guards(self):
        repo = self.root / 'repo with spaces'
        shutil.copytree(ROOT / 'scripts', repo / 'scripts')
        self.fixture([record(f'NM_{i}.1', organism=species) for i, species in enumerate(
            ['Anolis carolinensis', 'Pogona vitticeps', 'Python bivittatus', 'Sphenodon punctatus'])], repo / 'genes/SHH')
        bindir = self.root / 'bin'
        bindir.mkdir()
        (bindir / 'python').symlink_to(sys.executable)
        mock = '''#!/usr/bin/env python3
import pathlib, sys
tool = pathlib.Path(sys.argv[0]).name
args = sys.argv[1:]
def value(key): return args[args.index(key) + 1]
if tool == 'macse':
    source = value('-seq') if value('-prog') == 'alignSequences' else value('-align')
    pathlib.Path(value('-out_NT')).write_text(pathlib.Path(source).read_text())
    if '-out_AA' in args: pathlib.Path(value('-out_AA')).write_text('>mock\\nMK*\\n')
elif tool == 'iqtree2':
    assert value('-T') == '2' and value('-st') == 'DNA'
    pathlib.Path(value('-pre')+'.treefile').write_text('(Anolis_carolinensis,Python_bivittatus);\\n')
elif tool == 'seqkit':
    assert args[0] == 'stats'
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
        self.assertTrue((repo / 'genes/SHH/iqtree/SHH.treefile').exists())
        bad = repo / 'bad_genes.txt'
        bad.write_text('../escape\n')
        env['NCBI_EMAIL'] = 'test@example.org'
        result = subprocess.run(['bash', str(repo / 'scripts/01_download_cds.sh'), str(bad)],
                                text=True, capture_output=True, env=env)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse((repo / 'escape').exists())


if __name__ == '__main__':
    unittest.main()
