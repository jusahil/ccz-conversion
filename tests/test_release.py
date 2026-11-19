"""Fast mathematical property, data-schema and failure-mode tests."""
import csv, hashlib, itertools, random, subprocess, sys, tempfile, unittest
from pathlib import Path
from src.ccz_research import aligned_yield, has_accessible_cubic
from verification.check_manifest import verify
from verification.branch_records import parse
from verification.reproduce import check_pdf
ROOT=Path(__file__).resolve().parents[1]
class ResearchTests(unittest.TestCase):
    def test_cubic_witness(self):
        self.assertFalse(has_accessible_cubic(3,[7],[1,2]));self.assertTrue(has_accessible_cubic(3,[7],[1,2,4]));self.assertFalse(has_accessible_cubic(3,[7,7],[1,2,4]))
    def test_cubic_random_restrictions(self):
        rng=random.Random(26)
        for n in range(1,6):
            for _ in range(12):
                terms=[rng.randrange(1<<n) for _ in range(6)];terms=[m for m in terms if m.bit_count()<=3];g=[rng.randrange(1<<n) for _ in range(3)]
                span={0}
                for v in g: span|={x^v for x in span}
                def f(x):return sum(x&m==m for m in terms)&1
                expected=any(sum(f(x) for x in (0,a,b,c,a^b,a^c,b^c,a^b^c))&1 for a,b,c in itertools.combinations(span,3))
                self.assertEqual(has_accessible_cubic(n,terms,g),expected)
    def test_input_guards(self):
        for args in [(True,[7],[]),(3,[8],[]),(4,[15],[]),(3,[7],[-1])]:
            with self.assertRaises(ValueError):has_accessible_cubic(*args)
        for args in [(0,0),(1,4),(True,1),(1,True)]:
            with self.assertRaises(ValueError):aligned_yield(*args)
        self.assertTrue(has_accessible_cubic(100,[7],[1,2,4]))
    def test_yield(self):
        with (ROOT/'data/batch_yield.csv').open() as f:
            for row in csv.DictReader(f):n,e,y=map(int,row.values());self.assertEqual(aligned_yield(n,e),y)
    def test_snapshot(self):self.assertGreater(verify(ROOT),0)
    def test_manifest_negative(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'a').write_text('x');(p/'SHA256SUMS').write_text(hashlib.sha256(b'x').hexdigest()+'  a\n');self.assertEqual(verify(p),1)
            (p/'unlisted').write_text('x')
            with self.assertRaises(ValueError):verify(p)
            (p/'unlisted').unlink();(p/'a').write_text('y')
            with self.assertRaises(ValueError):verify(p)
            (p/'a').unlink()
            with self.assertRaises(ValueError):verify(p)
    def test_duplicate_or_invalid_branch(self):
        for s in ['BR 0 1 -1 6 128\nBR 0 1 -1 6 128\n','nullity 6 1\nnullity 6 2\n','BR 0 0 1 6 128\n','BR 1 0 0 6 128\n','unknown\n']:
            with self.assertRaises(ValueError):parse(s)
    def test_truncated_pdf(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'broken.pdf';p.write_bytes(b'%PDF-1.4\n')
            with self.assertRaises((subprocess.CalledProcessError,ValueError)):check_pdf(p)
    def test_optimized_runner_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'run';r=subprocess.run([sys.executable,'-O',str(ROOT/'verification/reproduce.py'),'--out',str(p)],capture_output=True,text=True)
            self.assertNotEqual(r.returncode,0);self.assertFalse(p.exists());self.assertIn('assertions',r.stderr)
if __name__=='__main__':unittest.main()
