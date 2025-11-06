import json
import traceback
from pathlib import Path
from rgt import RGT

repo_root = Path(__file__).resolve().parents[1]
settings_path = repo_root / 'tests' / 'data' / 'sample_settings.json'
fastq_path = repo_root / 'tests' / 'data' / 'sample.fastq.txt'

with open(settings_path) as fh:
    settings = json.load(fh)

print('settings loaded, fastq:', fastq_path)

rgt_inst = RGT(settings, str(repo_root/'tmp_input'), str(repo_root/'tmp_output'))
try:
    res = rgt_inst.rgt(str(fastq_path))
    print('result:', res)
except Exception as e:
    print('Exception:', repr(e))
    traceback.print_exc()
