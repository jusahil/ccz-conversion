"""Run finite scientific certificates in a new directory; fail closed."""
from pathlib import Path
import argparse, concurrent.futures, csv, hashlib, json, os, platform, shutil, subprocess, sys, time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from verification.check_manifest import verify,sha
from verification.branch_records import reconcile,require
from scripts.export_data import export
from scripts.plot_data import plot

def read_csv(p):
    with p.open(newline='') as f:return list(csv.reader(f))
def check_pdf(p):
    from PIL import Image,ImageStat
    info=subprocess.run(['pdfinfo',str(p)],capture_output=True,text=True,check=True,timeout=30).stdout
    import re
    require(re.search(r'^Pages:\s+1\s*$',info,re.M),'not a one-page PDF')
    raster=p.with_suffix('.png');subprocess.run(['pdftoppm','-f','1','-l','1','-r','120','-singlefile','-png',str(p),str(raster.with_suffix(''))],check=True,capture_output=True,timeout=30)
    with Image.open(raster) as im: require(min(im.size)>100 and ImageStat.Stat(im.convert('L')).stddev[0]>2,'blank/invalid rendered plot')
    return raster

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--profile',choices=['smoke','core','full'],default='full');p.add_argument('--out',type=Path,required=True);p.add_argument('--jobs',type=int,default=3);p.add_argument('--timeout',type=int,default=600);a=p.parse_args()
    if not __debug__:p.error('assertions must be enabled')
    if a.jobs<1 or a.timeout<1:p.error('positive jobs and timeout required')
    out=a.out.resolve()
    if out.exists():p.error('choose a new output directory')
    if out.is_relative_to(ROOT) and not out.is_relative_to(ROOT/'reproduction'):p.error('output must be external or inside reproduction/')
    for cmd in ('g++','pdfinfo','pdftoppm'):
        if not shutil.which(cmd):p.error(cmd+' required')
    verify(ROOT)
    selected=[r for r in json.loads((ROOT/'verification/programs.json').read_text()) if a.profile=='full' or (a.profile=='core' and r['profile']=='core') or (a.profile=='smoke' and r['id'] in {'C08','C09','S7-F'})]
    out.mkdir(parents=True);results=out/'results';data=out/'data';figures=out/'figures'
    for d in (results,data,figures):d.mkdir()
    env=os.environ.copy();env.pop('PYTHONOPTIMIZE',None);env.update(PYTHONHASHSEED='0',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONPATH=os.pathsep.join([str(ROOT/'src'),str(ROOT)]),PYTHONUNBUFFERED='1',MPLBACKEND='Agg',PYTHONDONTWRITEBYTECODE='1',MPLCONFIGDIR=str(out/'mpl-cache'),LC_ALL='C.UTF-8')
    start=time.monotonic();checks=[]
    def run(r):
        id=r['id'];src=ROOT/r['path'];record={'id':id,'source':r['path']};t=time.monotonic()
        try:
            require(sha(src)==r['source_sha256'],'source hash mismatch');cmd=[sys.executable,str(src)]
            if src.suffix=='.cpp':
                binary=results/(id+'.bin');build=['g++','-O3','-std=c++17','-UNDEBUG',str(src),'-o',str(binary)]
                with (results/(id+'.build.log')).open('wb') as log:subprocess.run(build,check=True,stdout=log,stderr=subprocess.STDOUT,timeout=a.timeout,env=env)
                cmd=[str(binary)]
            with (results/(id+'.out')).open('wb') as o,(results/(id+'.err')).open('wb') as e:
                subprocess.run(cmd+r['args'],cwd=out,check=True,stdout=o,stderr=e,env=env,timeout=a.timeout)
            raw=(results/(id+'.out')).read_bytes();record['stdout_sha256']=hashlib.sha256(raw).hexdigest()
            if r['reference_mode'] in ('exact','aggregate'):
                if r['reference_mode']=='aggregate':raw=b''.join(l for l in raw.splitlines(keepends=True) if not l.startswith((b'BR ',b'CERT ')))
                require(hashlib.sha256(raw).hexdigest()==r['expected_stdout_sha256'],'reference stdout differs')
            record['status']='PASS'
        except Exception as e:record.update(status='FAIL',error=str(e))
        record['seconds']=round(time.monotonic()-t,3);print(id,record['status'],record['seconds'],flush=True);(results/(id+'.json')).write_text(json.dumps(record,indent=2)+'\n');return record
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as ex:runs=list(ex.map(run,selected))
    good=all(r['status']=='PASS' for r in runs);branch=None
    if good:
        try:
            branch=reconcile(results,data);checks.append({'check':'signed branch equality','status':'PASS'})
            for name in ('signed_pauli_branches.csv','fano_histograms.csv'):
                require(read_csv(ROOT/'data'/name)==read_csv(data/name),'signed certificate snapshot mismatch: '+name)
                checks.append({'check':name+' snapshot','status':'PASS'})
            if a.profile!='smoke':
                export(results,data,branch['maximum_bound'])
                for reference in sorted((ROOT/'data').glob('*.csv')):
                    if reference.name in ('signed_pauli_branches.csv','fano_histograms.csv'):continue
                    require(read_csv(reference)==read_csv(data/reference.name),'table mismatch '+reference.name);checks.append({'check':reference.name,'status':'PASS'})
                with (results/'table_semantics.log').open('wb') as log:subprocess.run([sys.executable,str(ROOT/'checks/table_semantics.py'),str(data)],check=True,stdout=log,stderr=subprocess.STDOUT,env=env,timeout=a.timeout)
                checks.append({'check':'independent table semantics','status':'PASS'});plot(data,figures)
                for f in sorted(figures.glob('*.pdf')):
                    check_pdf(f);checks.append({'check':f.name+' parse and render','status':'PASS'})
                require(len(list(figures.glob('*.pdf')))==2,'missing figure')
        except Exception as e:good=False;checks.append({'check':'downstream','status':'FAIL','error':str(e)});print('FAIL downstream:',e,flush=True)
    summary={'status':'PASS' if good else 'FAIL','profile':a.profile,'python':sys.version,'compiler':subprocess.check_output(['g++','--version'],text=True).splitlines()[0],'platform':platform.platform(),'seconds':round(time.monotonic()-start,3),'runs':runs,'checks':checks,'branch_certificate':branch}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary['status'],summary['seconds'],flush=True);return int(not good)
if __name__=='__main__':sys.exit(main())
