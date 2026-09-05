"""File interface for the exact exhaustive C++ slice engine (run in WSL)."""
import hashlib
import json
from pathlib import Path
import subprocess

HERE=Path(__file__).resolve().parent


def orbit(b,word,p,directory,identifier):
    directory=Path(directory);directory.mkdir(exist_ok=True)
    source=directory/(identifier+'.txt');target=directory/(identifier+'.json')
    source.write_text(' '.join(map(str,[len(word),len(b)]+[x for row in b for x in row]+list(word)+list(p)))+'\n')
    subprocess.run([str(HERE/'bin/slice_orbit'),str(source),str(target)],check=True)
    result=json.loads(target.read_text());assert result['complete'] and result['states']==result['processed']
    result['engine_source_sha256']=hashlib.sha256((HERE/'slice_orbit.cpp').read_bytes()).hexdigest()
    result['engine_binary_sha256']=hashlib.sha256((HERE/'bin/slice_orbit').read_bytes()).hexdigest()
    result['input_sha256']=hashlib.sha256(source.read_bytes()).hexdigest()
    target.write_text(json.dumps(result,indent=2)+'\n')
    return result
