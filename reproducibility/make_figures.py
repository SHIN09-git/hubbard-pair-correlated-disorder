"""Publication-style figures made only from saved evidence and analytic formulas."""
from pathlib import Path
import json,time,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm
from src.physics import born_rate

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'docs/figures';OUT.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,
 'figure.dpi':140,'savefig.dpi':180,'svg.fonttype':'none','axes.labelsize':11})
COL=['#2166ac','#b35806','#5e3c99','#1b7837']
def save(fig,name):
 fig.savefig(OUT/(name+'.png'),bbox_inches='tight');fig.savefig(OUT/(name+'.svg'),bbox_inches='tight');plt.close(fig)
def latest(pattern):return sorted(ROOT.glob(pattern))[-1]
def main():
 start=time.perf_counter();cpu=time.process_time()
 data_path=latest('results/analysis/*/summary.json');data=json.loads(data_path.read_text())
 diag_path=latest('results/sample_diagnostic/*/summary.json');diag=json.loads(diag_path.read_text())
 theory_path=ROOT/'results/theory_audit/closed_kernel_initial.json';theory=json.loads(theory_path.read_text())
 def group(U,W,f,norm='variance',n=256):return next(r for r in data['groups'] if (r['U'],r['W'],r['f'],r['normalization'],r['n'])==(U,W,f,norm,n))
 # Analytic parameter classification and normalization controls.
 fig,ax=plt.subplots(1,3,figsize=(14.3,4.3),layout='constrained')
 kk=np.linspace(.0001,.9999,500)*np.pi;uu=np.linspace(.0001,5.5,450);K,U=np.meshgrid(kk,uu)
 E=np.sqrt(U**2+16*np.cos(K/2)**2);Q=2*np.arccos(np.clip(E/4,0,1));z=np.where(E>=4,0,np.where(K+Q>=np.pi,1,2))
 cmap=ListedColormap(['#e1e1e1','#efbd95','#8dbadd'])
 ax[0].pcolormesh(kk/np.pi,uu,z,cmap=cmap,vmin=-.5,vmax=2.5,rasterized=True)
 ax[0].plot(kk/np.pi,4*np.sin(kk/2),'k-',lw=1)
 kz=kk[kk>np.pi/2];ax[0].plot(kz/np.pi,4*np.sqrt(-np.cos(kz)),'k--',lw=1)
 ax[0].text(.32,4.8,'Energy closed',ha='center',fontsize=10)
 ax[0].text(.32,.72,'Finite cutoff exists',ha='center',fontsize=9)
 ax[0].text(.84,1.2,'No nonempty\nhigh-pass closure',ha='center',fontsize=9)
 for u in [2,3.5,5]:ax[0].plot(.75,u,'ko',ms=4)
 ax[0].set(xlabel=r'$|K|/\pi$',ylabel=r'$U/J$',title='(a) Leading-order dissociation support',xlim=(0,1),ylim=(0,5.5))
 fs=np.linspace(0,.995,260)
 for U,c in [(2,COL[0]),(3.5,COL[1])]:
  rates=np.array([born_rate(.75*np.pi,U,f*np.pi)[0] for f in fs])
  ax[1].plot(fs,rates,color=c,label=f'U/J = {U}')
  ax[2].plot(fs,rates*(1-fs),color=c,label=f'U/J = {U}')
 for a,title in [(ax[1],'(b) Fixed site variance'),(ax[2],'(c) Fixed pass-band height')]:
  a.axvline(.5,color='gray',ls=':',lw=1,label=r'$q_B/\pi=0.5$')
  a.axvline(.9416570958,color=COL[1],ls='--',lw=.9)
  a.set(xlabel=r'$q_0/\pi$',ylabel=r'$c=\Gamma/W^2$',title=title,xlim=(0,1),ylim=(0,None));a.legend(frameon=False,fontsize=9)
 save(fig,'01_support_and_normalization')
 # Frozen endpoints, with finite-W stress check kept visible.
 fig,ax=plt.subplots(1,2,figsize=(10.5,4.3),layout='constrained')
 for i,(f,norm,label) in enumerate([(0.,'variance','No gap'),(.8,'variance','Gap 0.8, fixed variance'),(.8,'passband','Gap 0.8, fixed height')]):
  rr=[group(2,w,f,norm) for w in [.02,.04,.08]];x=np.array([r['W'] for r in rr])
  y=np.array([r['main']['mean'] for r in rr]);err=np.array([[r['main']['mean']-r['main']['ci95'][0],r['main']['ci95'][1]-r['main']['mean']] for r in rr]).T
  ax[0].errorbar(x,y,yerr=err,fmt='o-',color=COL[i],capsize=3,label=label)
  ax[0].axhline(rr[0]['coefficient'],color=COL[i],ls=':',lw=1)
 ax[0].axvspan(.065,.09,color='#efefef',zorder=0)
 ax[0].text(.0775,.13,'Stress check',ha='center',fontsize=8)
 ax[0].set(xlabel='Input W/J',ylabel=r'Loss slope / $W^2$',title='(a) U/J=2, N=256, frozen t=12-24',xlim=(.012,.088),ylim=(0,2.4));ax[0].legend(frameon=True,facecolor='white',edgecolor='white',framealpha=1,fontsize=8)
 for i,f in enumerate([0.,.8,.9,.97]):
  r=group(3.5,.02,f);ts=np.array(r['times']);mean=np.array(r['mean_loss'])/.02**2;sem=np.array(r['sem_loss'])/.02**2
  pred=[r for r in diag['rows'] if (r['n'],r['U'],r['W'],r['f'],r['normalization'])==(256,3.5,.02,f,'variance')]
  yy=np.mean([r['predicted_loss'] for r in pred],axis=0)/.02**2
  ax[1].plot(ts,mean,color=COL[i],label=f'gap {f}');ax[1].fill_between(ts,mean-sem,mean+sem,color=COL[i],alpha=.13)
  ax[1].plot(ts,yy,color=COL[i],ls=':',lw=1.8)
 ax[1].set(xlabel=r'Time $tJ/\hbar$',ylabel=r'Bound-band loss / $W^2$',title='(b) U/J=3.5: solid dynamics, dotted 2nd order',xlim=(0,32),ylim=(0,None));ax[1].legend(frameon=False,fontsize=8)
 save(fig,'02_dynamics_and_virtual_mixing')
 # Endpoint powers, no finite-W scaling claim.
 fig,ax=plt.subplots(1,2,figsize=(10.5,4.2),layout='constrained')
 for i,r in enumerate(theory['edges']):
  d=np.array([x['delta'] for x in r['rows']]);v=np.array([x['c'] for x in r['rows']]);p=np.array([x['passband'] for x in r['rows']])
  label_v={'wrap_open':'U/J=2: nonzero limit','critical':r'$U=U_*:$ exponent 1/2','internal':'U/J=3.5: exponent 3/2'}[r['mode']]
  label_h={'wrap_open':'U/J=2: exponent 1','critical':r'$U=U_*:$ exponent 3/2','internal':'U/J=3.5: exponent 3/2'}[r['mode']]
  ax[0].loglog(d,v,'o-',color=COL[i],label=label_v,ms=4)
  ax[1].loglog(d,p,'o-',color=COL[i],label=label_h,ms=4)
 for a,title in zip(ax,['(a) Fixed variance','(b) Fixed pass-band height']):
  a.set(xlabel=r'$\Delta=q_c-q_0$ (radians)',ylabel=r'$c$',title=title);a.legend(frameon=False,fontsize=8);a.grid(alpha=.12,which='major')
 save(fig,'03_endpoint_laws')
 # Leakage is the only robustness extension.
 leak_run=latest('results/leakage/*/manifest.json').parent;meta=json.loads((leak_run/'manifest.json').read_text());assert meta['status']=='PASS'
 rr=[json.loads(p.read_text()) for p in leak_run.glob('epsilon*.json')]
 fig,ax=plt.subplots(1,2,figsize=(10.5,4.2),layout='constrained');eps=np.logspace(-5,-1,250)
 ax[0].semilogx(eps,eps/(.03+.97*eps),color=COL[0]);ax[0].axhline(.1,color='gray',ls=':')
 tol=.1*.03/(1-.1*.97);ax[0].axvline(tol,color=COL[1],ls='--');ax[0].text(tol*1.15,.35,f'10% target\nfloor <= {tol:.4f}',fontsize=9)
 ax[0].set(xlabel=r'Gap / pass-band power ratio $\epsilon$',ylabel=r'$c_\epsilon/c_{white}$',title='(a) Narrow pass-band amplifies leakage',ylim=(0,.9))
 for i,e in enumerate([0.,.003,.01]):
  rows=[r for r in rr if r['epsilon']==e];ts=rows[0]['times'];y=np.array([r['loss'] for r in rows])/.02**2;pred=np.array([r['predicted_loss'] for r in rows])/.02**2
  ax[1].plot(ts,y.mean(axis=0),color=COL[i],label=f'floor = {e}');ax[1].plot(ts,pred.mean(axis=0),color=COL[i],ls=':',lw=1.8)
 ax[1].set(xlabel=r'Time $tJ/\hbar$',ylabel=r'Bound-band loss / $W^2$',title='(b) Finite-time check: U=3.5, gap=0.97',xlim=(0,32),ylim=(0,None));ax[1].legend(frameon=False,fontsize=8)
 save(fig,'04_leakage_robustness')
 sources=[data_path,diag_path,theory_path,leak_run/'manifest.json',Path(__file__)]
 (OUT/'manifest.json').write_text(json.dumps(dict(inputs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},outputs={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.glob('*') if p.suffix in ['.png','.svg']},wall_seconds=time.perf_counter()-start,process_cpu_seconds=time.process_time()-cpu),indent=2),encoding='utf-8')
 print('Saved four PNG and SVG figures',OUT)
if __name__=='__main__':main()
