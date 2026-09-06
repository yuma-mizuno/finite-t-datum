"""Replay the complete upper/partner partition and canonical union from the final archives."""
import hashlib,json,math,sys,zipfile
from pathlib import Path
HERE=Path(__file__).resolve().parent
def main(rank):
    directory=HERE/f'rank{rank}';constants=json.loads((directory/'constant_candidates.json').read_text())
    assert constants['enumeration_complete'];count=math.factorial(rank+1);union=set();seen=set();blocks=0
    with zipfile.ZipFile(directory/'constant_tasks.zip') as tasks,zipfile.ZipFile(directory/'constant_parts.zip') as parts:
        for name in tasks.namelist():
            task=json.loads(tasks.read(name));a=task['upper_index'];assert a not in seen and task['completed'] and task['rank']==rank;seen.add(a)
            keys={tuple(k) for k in task['keys']};intervals=task.get('partner_blocks')
            if intervals:
                expected=a;from_parts=set()
                for lo,hi in intervals:
                    assert lo==expected and lo<hi<=count;expected=hi
                    item=json.loads(parts.read(f'upper-{a}-part-{lo}-{hi}.json'))
                    assert item['completed'] and item['upper_index']==a and item['rank']==rank
                    assert item['partner_interval']==[lo,hi]
                    assert item['source_sha256'] in constants['task_source_sha256']
                    from_parts.update(tuple(k) for k in item['keys']);blocks+=1
                assert expected==count and from_parts==keys
            union.update(keys)
        assert seen==set(range(count))
    expected={tuple(c['symmetrizer'])+tuple(x for sign in ('N_plus_1','N_minus_1') for row in c[sign] for x in row) for c in constants['candidates']}
    assert union==expected and len(union)==constants['count']
    result={'rank':rank,'complete_upper_tasks':len(seen),'replayed_partner_blocks':blocks,'canonical_union':len(union),'all_partner_intervals_cover_their_full_range':True,'archived_task_union_equals_final_candidates':True,
            'archives':{name:hashlib.sha256((directory/name).read_bytes()).hexdigest() for name in ('constant_tasks.zip','constant_parts.zip')}}
    (directory/'constant-archive-audit.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8',newline='\n');print(result,flush=True)
if __name__=='__main__':main(int(sys.argv[1]))
