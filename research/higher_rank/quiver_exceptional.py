"""Reuse the exact exceptional-class search for unresolved small quivers."""
import json
from pathlib import Path
import signal
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'catalogue'))
import quiver_more

def main(rank):
    directory=HERE/f'rank{rank}';quiver_more.HERE=directory
    records=json.loads((directory/'base-records.json').read_text())
    quiver_more.dynkin_path=lambda *a,**kw:None
    signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError()))
    for record in records:
        q=json.loads((directory/'quiver-data.json').read_text())[record['id']]
        if q['status'].startswith('certified'):continue
        quiver_more.records=[record]
        try:
            signal.alarm(25);quiver_more.main()
        except TimeoutError:print(record['id'],'exceptional search remains unresolved',flush=True)
        finally:signal.alarm(0)
if __name__=='__main__':main(int(sys.argv[1]))
