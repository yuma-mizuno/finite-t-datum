"""Standalone scientific SVG plots; the HTML reader also has interactive plots."""
from pathlib import Path
import json
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

HERE=Path(__file__).resolve().parent
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'svg.fonttype':'none','svg.hashsalt':'finite-t-data-v2'})

def main():
 (HERE/'plots').mkdir(exist_ok=True)
 spectra=json.loads((HERE/'spectral-data.json').read_text())
 for id,item in spectra.items():
  e=item['exponents'];N=e['root_order'];fig,ax=plt.subplots(figsize=(6,6))
  ax.add_patch(Circle((0,0),1,facecolor='#edf4ee',edgecolor='#8aaba0',linewidth=1))
  ax.axhline(0,color='#b6c8bd',linewidth=.6);ax.axvline(0,color='#b6c8bd',linewidth=.6)
  for m in range(N):ax.plot(math.cos(2*math.pi*m/N),math.sin(2*math.pi*m/N),'o',markersize=3,color='#b6c8bd')
  for term in e['multiplicities']:
   m,a=term['m'],term['multiplicity'];x,y=math.cos(2*math.pi*m/N),math.sin(2*math.pi*m/N)
   ax.scatter([x],[y],s=85 if a==1 else 230,c='#146452',zorder=3)
   if a>1:ax.text(x,y,str(a),ha='center',va='center',color='white',fontsize=10,zorder=4)
   ax.text(x*1.12,y*1.12,str(m),ha='center',va='center',fontsize=10)
  ax.set(xlim=(-1.28,1.28),ylim=(-1.28,1.28),aspect='equal',xlabel='Real part',ylabel='Imaginary part',xticks=[-1,0,1],yticks=[-1,0,1])
  ax.set_title(f'{id}: exponents modulo {N}\n'+r'$\zeta^m=\exp(2\pi i m/N)$'+f'   (total multiplicity {e["degree"]})',pad=15)
  for spine in ax.spines.values():spine.set_visible(False)
  fig.text(.5,.015,'Outer labels: m. Numbers inside points: multiplicities > 1.',ha='center',fontsize=9,color='#586c66')
  fig.tight_layout(rect=(0,.03,1,1));fig.savefig(HERE/'plots'/f'{id}.svg',metadata={'Date':None,'Creator':'Finite T-data catalogue; Matplotlib'});plt.close(fig)
 print('Generated 61 standalone SVG exponent plots.')

if __name__=='__main__':main()
