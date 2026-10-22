"""Export the numerical arrays used by make_figures.py, without new trajectories."""
from pathlib import Path
import sys, json, csv, hashlib
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.physics import born_rate

def latest(pattern):
    return sorted(ROOT.glob(pattern))[-1]

def main():
    out = ROOT / 'data/figure_data'
    out.mkdir(parents=True, exist_ok=True)
    inputs = [latest('results/analysis/*/summary.json'),
              latest('results/sample_diagnostic/*/summary.json'),
              ROOT/'results/theory_audit/closed_kernel_initial.json',
              latest('results/leakage/*/manifest.json')]
    analysis, diagnostic, theory, _ = [json.loads(p.read_text('utf-8')) for p in inputs]
    def save(name, rows):
        rows = list(rows)
        with (out/name).open('w', encoding='utf-8', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    kk = np.linspace(.0001,.9999,500)*np.pi
    uu = np.linspace(.0001,5.5,450)
    K,U = np.meshgrid(kk,uu)
    E=np.sqrt(U**2+16*np.cos(K/2)**2)
    Q=2*np.arccos(np.clip(E/4,0,1))
    z=np.where(E>=4,0,np.where(K+Q>=np.pi,1,2))
    np.savez_compressed(out/'01_phase_grid.npz',K_over_pi=kk/np.pi,U_over_J=uu,phase=z)
    save('01_boundaries.csv', (dict(K_over_pi=float(k/np.pi),U_energy=float(4*np.sin(k/2)),
             U_star=float(4*np.sqrt(-np.cos(k))) if k>np.pi/2 else '') for k in kk))
    save('01_rates.csv', (dict(U=U,f=float(f),c_variance=born_rate(.75*np.pi,U,f*np.pi)[0],
             c_passband=born_rate(.75*np.pi,U,f*np.pi)[0]*(1-f))
             for U in (2.,3.5) for f in np.linspace(0,.995,260)))
    groups=analysis['groups']
    def group(U,W,f,norm='variance',n=256):
        return next(r for r in groups if (r['U'],r['W'],r['f'],r['normalization'],r['n'])==(U,W,f,norm,n))
    rows=[]
    for f,norm in [(0.,'variance'),(.8,'variance'),(.8,'passband')]:
        for W in [.02,.04,.08]:
            r=group(2,W,f,norm)
            rows.append(dict(f=f,normalization=norm,W=W,n=256,U=2,
               slope=r['main']['mean'],ci95_lo=r['main']['ci95'][0],ci95_hi=r['main']['ci95'][1],
               coefficient=r['coefficient'],seeds=len(r['seeds'])))
    save('02_slopes.csv',rows)
    rows=[]
    for f in [0.,.8,.9,.97]:
        r=group(3.5,.02,f)
        pred=[v['predicted_loss'] for v in diagnostic['rows']
              if (v['n'],v['U'],v['W'],v['f'],v['normalization'])==(256,3.5,.02,f,'variance')]
        for t,mean,sem,pr in zip(r['times'],r['mean_loss'],r['sem_loss'],np.mean(pred,axis=0)):
            rows.append(dict(f=f,time=t,mean_loss_over_W2=mean/.02**2,
                       sem_loss_over_W2=sem/.02**2,prediction_over_W2=pr/.02**2))
    save('02_curves.csv',rows)
    save('03_endpoints.csv',(dict(mode=e['mode'],U=e['U'],qc_over_pi=e['qc_over_pi'],**r)
          for e in theory['edges'] for r in e['rows']))
    save('04_leakage_formula.csv',(dict(epsilon=e,coefficient_ratio=e/(.03+.97*e))
          for e in np.logspace(-5,-1,250)))
    leak=[json.loads(p.read_text()) for p in inputs[3].parent.glob('epsilon*.json')]
    rows=[]
    for e in [0.,.003,.01]:
        rs=[r for r in leak if r['epsilon']==e]
        mean=np.mean([r['loss'] for r in rs],axis=0)/.02**2
        pred=np.mean([r['predicted_loss'] for r in rs],axis=0)/.02**2
        for t,y,pr in zip(rs[0]['times'],mean,pred):
            rows.append(dict(epsilon=e,time=t,mean_loss_over_W2=y,prediction_over_W2=pr))
    save('04_leakage_curves.csv',rows)
    save('all_group_statistics.csv',(dict(n=r['n'],U=r['U'],W=r['W'],f=r['f'],
        normalization=r['normalization'],seeds=len(r['seeds']),window_start=r['main_window'][0],
        window_end=r['main_window'][1],coefficient=r['coefficient'],slope=r['main']['mean'],
        ci95_lo=r['main']['ci95'][0],ci95_hi=r['main']['ci95'][1],loss_at_24=r['mean_loss_at_24']['mean']) for r in groups))
    save('paired_comparisons.csv',(dict(n=r['n'],W=r['W'],normalization=r['normalization'],ratio=r['ratio'],
          ratio_ci95_lo=r['bootstrap_ratio_ci95'][0],ratio_ci95_hi=r['bootstrap_ratio_ci95'][1],
          difference=r['difference']['mean'],difference_ci95_lo=r['difference']['ci95'][0],
          difference_ci95_hi=r['difference']['ci95'][1],theory_ratio=r['theory_ratio']) for r in analysis['comparisons']))
    names={
       '01_phase_grid.npz':'图1a：K_over_pi[500]、U_over_J[450]、phase[450,500]；0能量关闭，1非空高通不能关闭，2可有限缺口关闭。',
       '01_boundaries.csv':'图1a：能量边界与U*边界；空U_star表示不适用。',
       '01_rates.csv':'图1b/c：K=.75π，J=1；f=q₀/π；两种归一化的系数。',
       '02_slopes.csv':'历史报告图3a：误差棒为95% Student区间。',
       '02_curves.csv':'历史报告图3b：色带为均值标准误；点线为逐样本二阶预测平均。',
       '03_endpoints.csv':'本稿图2／历史图2：delta以弧度计，c为固定方差，passband为固定高度。',
       '04_leakage_formula.csv':'历史图4a：ε为功率比，固定总方差。',
       '04_leakage_curves.csv':'历史图4b：实线均值与点线二阶预测。',
       'all_group_statistics.csv':'全部确认组别；window为实际拟合窗，coefficient不是有限W寿命。',
       'paired_comparisons.csv':'配对差、比值和描述性重采样区间；不重复计数种子。'}
    (out/'INDEX.md').write_text('# 图源数据与统计表\n\nJ=ħ=1；时间单位ħ/J，能量单位J。CSV采用UTF-8。原始逐样本轨迹仍保存在results目录；这些表只是可读导出。图2/4编号中的02/04对应历史报告，最小研究稿只使用01和03。\n\n'+
       '\n'.join(f'- [{n}]({n})：{desc}' for n,desc in names.items())+
       '\n\n生成入口：[export_figure_data.py](../../tools/export_figure_data.py)。重建说明：[REPRODUCTION.md](../../REPRODUCTION.md)。\n',encoding='utf-8')
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    manifest=dict(inputs={p.relative_to(ROOT).as_posix():sha(p) for p in inputs},
        source={p.relative_to(ROOT).as_posix():sha(p) for p in [Path(__file__),ROOT/'make_figures.py',ROOT/'src/physics.py']},
        outputs={n:sha(out/n) for n in names})
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print('Exported',len(names),'figure and table data files to',out)
if __name__=='__main__':main()
