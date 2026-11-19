"""Regenerate the release manifest after an intentional, reviewed change."""
from pathlib import Path
import argparse,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from verification.check_manifest import ignored,sha

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--confirm-reviewed-change',action='store_true');a=p.parse_args()
    if not a.confirm_reviewed_change:p.error('review all changed files before updating the manifest')
    files=[]
    for f in sorted(ROOT.rglob('*')):
        rel=f.relative_to(ROOT).as_posix()
        if ignored(rel) or rel=='SHA256SUMS':continue
        if f.is_symlink():raise ValueError('unexpected symlink: '+rel)
        if f.is_file():files.append(f'{sha(f)}  {rel}\n')
    (ROOT/'SHA256SUMS').write_text(''.join(files));print(len(files),'files inventoried')
if __name__=='__main__':main()
