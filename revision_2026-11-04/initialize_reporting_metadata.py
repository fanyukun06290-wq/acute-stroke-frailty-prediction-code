from pathlib import Path
import shutil,hashlib,json
ROOT=Path(__file__).resolve().parent
def main():
 source=ROOT/'reference_metadata';target=ROOT/'phenotype_v2_results'
 files=list(source.glob('*.csv'))
 assert len(files)==5
 assert not any((target/p.name).exists() for p in files),'Existing dictionaries preserved'
 target.mkdir(exist_ok=True)
 for p in files:shutil.copy2(p,target/p.name)
 print(json.dumps({'copied_static_dictionaries':len(files),'raw_data_read':False}))
if __name__=='__main__':main()
