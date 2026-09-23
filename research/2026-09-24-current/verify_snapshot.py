"""Verify full bytes and optionally restore included evidence; never execute it."""
from pathlib import Path,PurePosixPath
import argparse,hashlib,json,zipfile
HERE=Path(__file__).resolve().parent
def sha_stream(f):
    h=hashlib.sha256();size=0
    for b in iter(lambda:f.read(1048576),b''):h.update(b);size+=len(b)
    return h.hexdigest(),size
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def verify():
    manifest=read(HERE/'MANIFEST.json');cat=read(HERE/'CATALOG.json');roster=read(HERE/'SOURCE_ROSTER.json')
    for n,r in manifest['files'].items():
        p=(HERE/n).resolve();assert p.is_relative_to(HERE.resolve())
        with p.open('rb') as f:h,size=sha_stream(f)
        assert h==r['sha256'] and size==r['bytes'],n
    verified={};bundle_count=0
    for row in cat['bundles']:
        p=HERE/row['path']
        with p.open('rb') as f:h,size=sha_stream(f)
        assert (h,size)==(row['sha256'],row['bytes'])
        with zipfile.ZipFile(p) as z:
            assert len(z.namelist())==len(set(z.namelist()))
            for name in z.namelist():
                with z.open(name) as f:bh,bs=sha_stream(f)
                assert name=='blobs/'+bh and bh not in verified
                br=cat['blobs'][bh];assert br['bytes']==bs and br['bundle']==row['path'] and br['member']==name
                verified[bh]=bs
        bundle_count+=1
    assert set(verified)==set(cat['blobs'])
    for name,r in cat['files'].items():assert verified[r['sha256']]==r['bytes'],name
    excluded={r['path']:r for r in cat['excluded_administrative_files']}
    archives={r['path']:r for r in cat['original_archives']}
    dispositions={**cat['files'],**excluded,**archives}
    assert len(dispositions)==len(cat['files'])+len(excluded)+len(archives)
    for n,r in roster['files'].items():
        dest=dispositions['workspace/'+n];assert dest['sha256']==r['sha256'] and dest['bytes']==r['bytes'],n
    for n,r in archives.items():
        for member in r['members']:
            dest=dispositions[n+'/members/'+member['path']]
            assert (dest['sha256'],dest['bytes'])==(member['sha256'],member['bytes'])
    for n,r in cat['latest_files'].items():
        with (HERE/n).open('rb') as f:h,size=sha_stream(f)
        assert h==r['sha256'] and size==r['bytes']
    for r in cat['stages']:
        if r['report']:
            with (HERE/r['report']).open('rb') as f:h,size=sha_stream(f)
            assert h==r['report_sha256']
    return dict(passed=True,manifest_files=len(manifest['files']),stages=len(cat['stages']),
        pinned_source_paths=len(roster['files']),published_logical_files=len(cat['files']),
        unique_full_blobs=len(verified),bundles=bundle_count,archives_accounted=len(archives),
        excluded_administrative_paths=len(excluded),truncated_research_files=0),cat
def restore(cat,directory,prefix):
    root=Path(directory).resolve();root.mkdir(parents=True,exist_ok=False);count=0
    for name,r in cat['files'].items():
        if not name.startswith(prefix):continue
        rel=PurePosixPath(name);assert not rel.is_absolute() and '..' not in rel.parts
        dest=(root/str(rel)).resolve();assert dest.is_relative_to(root)
        dest.parent.mkdir(parents=True,exist_ok=True);br=cat['blobs'][r['sha256']]
        with zipfile.ZipFile(HERE/br['bundle']) as z, z.open(br['member']) as src, dest.open('xb') as out:
            for block in iter(lambda:src.read(1048576),b''):out.write(block)
        with dest.open('rb') as f:h,size=sha_stream(f)
        assert h==r['sha256'] and size==r['bytes'];count+=1
    assert count>0,'No matching files'
    return count
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--restore');ap.add_argument('--prefix',default='workspace/')
    args=ap.parse_args();result,cat=verify()
    if args.restore:result['restored_files']=restore(cat,args.restore,args.prefix)
    print(json.dumps(result,sort_keys=True,indent=2))
