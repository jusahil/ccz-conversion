"""Strict release integrity verification (not cryptographic authenticity)."""
import hashlib, re, sys
from pathlib import Path, PurePosixPath
ROOT=Path(__file__).resolve().parents[1]

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def ignored(rel):
    p=PurePosixPath(rel)
    return p.parts[0] in {'.git','.venv','reproduction'} or '__pycache__' in p.parts or p.suffix=='.pyc'
def verify(root=ROOT):
    expected={}
    for line in (root/'SHA256SUMS').read_text().splitlines():
        m=re.fullmatch(r'([0-9a-f]{64})  (.+)',line)
        if not m: raise ValueError('malformed manifest row')
        digest,rel=m.groups(); p=PurePosixPath(rel)
        if p.is_absolute() or '..' in p.parts or str(p)!=rel or rel in expected or rel=='SHA256SUMS' or ignored(rel): raise ValueError('invalid or duplicate manifest path: '+rel)
        expected[rel]=digest
    actual={p.relative_to(root).as_posix() for p in root.rglob('*') if (p.is_file() or p.is_symlink()) and not ignored(p.relative_to(root).as_posix()) and p.relative_to(root).as_posix()!='SHA256SUMS'}
    if actual!=set(expected): raise ValueError(f'missing={sorted(set(expected)-actual)} unexpected={sorted(actual-set(expected))}')
    for rel,digest in expected.items():
        f=root/rel
        if f.is_symlink() or sha(f)!=digest:raise ValueError('checksum mismatch: '+rel)
    return len(expected)
if __name__=='__main__':
    try: print(f'{verify()} release files verified')
    except (OSError,ValueError) as e: print(e,file=sys.stderr);sys.exit(1)
