"""Post-hoc analysis of saved trajectories only; no dynamics or theory recomputation."""
from pathlib import Path
from collections import defaultdict
from datetime import datetime, timezone
import argparse, csv, hashlib, json, platform, sys, time
import numpy as np


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def key(r):
    return tuple(r[k] for k in ('n','U','W','f','normalization'))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--input-root',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();root=args.input_root.resolve();out=args.output.resolve()
    out.mkdir(parents=True,exist_ok=False)
    start=time.perf_counter();cpu=time.process_time();inputs={}
    def read(p):
        inputs[p.relative_to(root).as_posix()]=sha(p)
        return json.loads(p.read_text('utf-8'))
    diagdir=root/'results/sample_diagnostic/20261013T072213Z'
    diagnostic=read(diagdir/'summary.json')
    frozen=read(root/'results/analysis/20261013T071959Z/summary.json')
    oldgroups={key(r):r for r in frozen['groups']}
    dr={key(r)+(r['seed'],):r for r in diagnostic['rows']}
    assert len(dr)==len(diagnostic['rows'])
    groups=defaultdict(list);seen=set();max_curve_error=0.;max_slope_error=0.;max_frozen_error=0.
    for folder in sorted((root/'results/confirmation').iterdir()):
        if not folder.is_dir():continue
        status=read(folder/'manifest.json')
        assert status['status']=='PASS'
        for p in sorted(folder.glob('N*.json')):
            r=read(p);k=key(r);identity=k+(r['seed'],)
            assert identity not in seen;seen.add(identity)
            assert identity in dr,('Missing existing prediction',identity)
            d=dr[identity];ts=np.asarray(r['times']);pop=np.asarray(r['bound_populations'])
            assert np.all(pop>0),('Nonpositive population',identity)
            modefile=diagdir/f"modes_N{r['n']}_U{r['U']}.npz"
            inputs[modefile.relative_to(root).as_posix()]=sha(modefile)
            with np.load(modefile) as z:assert np.array_equal(z['times'],ts)
            L2=np.asarray(d['predicted_loss']);actual=1-pop
            assert L2.shape==actual.shape
            err=float(np.max(np.abs(actual-np.asarray(d['actual_loss']))));max_curve_error=max(max_curve_error,err)
            assert err<1e-12
            groups[k].append((r,d,ts,pop,L2,p.relative_to(root).as_posix()))
    assert seen==set(dr),('Unmatched predictions',len(set(dr)-seen))
    assert set(groups)==set(oldgroups)
    B=10000;rng=np.random.default_rng(2026101303);group_results=[];seed_results=[]
    for gi,(k,records) in enumerate(sorted(groups.items())):
        records.sort(key=lambda item:item[0]['seed']);n,U,W,f,norm=k;nseed=len(records)
        ts=records[0][2];assert all(np.array_equal(r[2],ts) for r in records)
        pop=np.stack([r[3] for r in records]);L2=np.stack([r[4] for r in records])
        assert np.min(L2)>-1e-12
        logpop=-np.log(pop)
        curve_exp_loss=-np.expm1(-L2)
        indices=rng.integers(0,nseed,size=(B,nseed))
        counts=np.zeros((B,nseed),dtype=np.int16)
        np.add.at(counts,(np.repeat(np.arange(B),nseed),indices.ravel()),1)
        weights=counts.astype(float)/nseed
        windows=[(12,24),(8,20),(16,28)] if U==2 else [(16,28),(12,24),(20,32)]
        for wi,window in enumerate(windows):
            mask=(ts>=window[0])&(ts<=window[1]);t=ts[mask]
            centered=t-t.mean();a=centered/(centered@centered)/W**2
            act=(1-pop[:,mask])@a;pred=L2[:,mask]@a
            exp_curve=curve_exp_loss[:,mask]@a
            logs=logpop[:,mask]@a
            const_valid=bool(np.all(pred>0))
            exp_const=None
            if const_valid:
                const_population=np.exp(-W**2*pred[:,None]*t[None,:])
                exp_const=(1-const_population)@a
                assert np.allclose((-np.log(const_population))@a,pred,rtol=1e-11,atol=1e-12)
            else:const_population=None
            bpred=weights@pred;bact=weights@act;bcv=weights@exp_curve;blogs=weights@logs
            log_mean=float((-np.log(pop[:,mask].mean(axis=0)))@a)
            blog_mean=(-np.log(weights@pop[:,mask]))@a
            metrics={}
            def metric(name,value,boot):
                metrics[name]={'estimate':float(value),'bootstrap_ci95':np.quantile(boot,[.025,.975]).tolist()}
            def ratio(name,x,y,bx,by):
                ci=np.quantile(by,[.025,.975])
                unstable=bool(ci[0]<=0<=ci[1] or abs(y)<1e-12)
                metrics[name]={'estimate':float(x/y) if abs(y)>1e-12 else None,
                    'bootstrap_ci95':None if unstable else np.quantile(bx/by,[.025,.975]).tolist(),
                    'denominator_interval_includes_zero':unstable}
            metric('actual_loss_slope',act.mean(),bact)
            metric('predicted_second_order_slope',pred.mean(),bpred)
            metric('curve_exponential_loss_slope',exp_curve.mean(),bcv)
            metric('slope_log_mean_population',log_mean,blog_mean)
            metric('mean_slope_log_population',logs.mean(),blogs)
            metric('difference_two_log_estimators',logs.mean()-log_mean,blogs-blog_mean)
            ratio('actual_over_pred',act.mean(),pred.mean(),bact,bpred)
            ratio('actual_over_curve_exponential',act.mean(),exp_curve.mean(),bact,bcv)
            ratio('log_mean_over_pred',log_mean,pred.mean(),blog_mean,bpred)
            ratio('mean_log_over_pred',logs.mean(),pred.mean(),blogs,bpred)
            coefficient=float(oldgroups[k]['coefficient'])
            if abs(coefficient)>1e-12:
                metric('actual_over_infinite_coefficient',act.mean()/coefficient,bact/coefficient)
                metric('pred_over_infinite_coefficient',pred.mean()/coefficient,bpred/coefficient)
            if const_valid:
                bexp=weights@exp_const
                metric('constant_exponential_loss_slope',exp_const.mean(),bexp)
                ratio('actual_over_constant_exponential',act.mean(),exp_const.mean(),bact,bexp)
            else:
                metrics['constant_exponential_loss_slope']={'estimate':None,'bootstrap_ci95':None,'reason':'At least one predicted slope is nonpositive; no seed was dropped.'}
            if wi==0:
                stored_pred=np.array([r[1]['predicted_slope'] for r in records])
                stored_act=np.array([r[1]['actual_slope'] for r in records])
                err=float(max(np.max(abs(pred-stored_pred)),np.max(abs(act-stored_act))))
                max_slope_error=max(max_slope_error,err);assert err<2e-12
                ferr=abs(float(act.mean())-oldgroups[k]['main']['mean'])
                max_frozen_error=max(max_frozen_error,ferr);assert ferr<2e-12
            actual_window=1-pop[:,mask]
            curve_model=curve_exp_loss[:,mask]
            denom=np.linalg.norm(actual_window,axis=1)
            relcurve=np.linalg.norm(actual_window-curve_model,axis=1)/denom
            deltapop=np.diff(pop[:,mask],axis=1)
            log_linear=logpop[:,mask].mean(axis=1)[:,None]+(W**2*logs)[:,None]*centered[None,:]
            logrms=np.sqrt(np.mean((logpop[:,mask]-log_linear)**2,axis=1))
            gr={'n':n,'U':U,'W':W,'f':f,'normalization':norm,'window':list(window),
                'role':'original_main_window' if wi==0 else 'original_sensitivity_window',
                'n_seeds':nseed,'seeds':[r[0]['seed'] for r in records],
                'time_points':t.tolist(),'coefficient':coefficient,'metrics':metrics,
                'constant_exponential_all_slopes_positive':const_valid,
                'predicted_nonpositive_slopes':int(np.sum(pred<=0)),
                'actual_negative_slopes':int(np.sum(act<0)),
                'seeds_with_population_increase_in_window':int(np.any(deltapop>1e-10,axis=1).sum()),
                'curve_exponential_relative_window_error_mean':float(relcurve.mean()),
                'curve_exponential_relative_window_error_median':float(np.median(relcurve)),
                'log_linear_fit_rms_mean':float(logrms.mean()),
                'actual_mean_loss_at_window_end':float(actual_window[:,-1].mean())}
            if const_valid:
                consterr=np.linalg.norm(actual_window-(1-const_population),axis=1)/denom
                gr['constant_exponential_relative_window_error_mean']=float(consterr.mean())
                gr['constant_exponential_relative_window_error_median']=float(np.median(consterr))
            group_results.append(gr)
            for i,r in enumerate(records):
                seed_results.append(dict(n=n,U=U,W=W,f=f,normalization=norm,seed=r[0]['seed'],
                    window_start=window[0],window_end=window[1],role=gr['role'],
                    actual_slope=float(act[i]),predicted_slope=float(pred[i]),
                    curve_exponential_slope=float(exp_curve[i]),
                    constant_exponential_slope=float(exp_const[i]) if const_valid else None,
                    log_population_slope=float(logs[i]),
                    curve_exponential_relative_window_error=float(relcurve[i]),
                    trajectory=r[5]))
    result={'status':'PASS','analysis_type':'POST_HOC_EXISTING_DATA_ONLY',
        'bootstrap':{'repetitions':B,'seed':2026101303,'unit':'whole matched seed trajectory','interval':'percentile 95%','time_points_resampled':False},
        'alignment':{'trajectory_count':len(seen),'prediction_count':len(dr),'groups':len(groups),
            'missing':[],'duplicate':[],'max_actual_curve_difference':max_curve_error,
            'max_saved_slope_difference':max_slope_error,'max_frozen_group_mean_difference':max_frozen_error},
        'groups':group_results,
        'interpretation':'Exponentials are unproved diagnostic ansatzes. Neither log estimator is asserted to be an asymptotic lifetime. Frozen primary results remain unchanged.'}
    (out/'summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    with (out/'per_seed_estimates.csv').open('w',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(seed_results[0]));writer.writeheader();writer.writerows(seed_results)
    flat=[]
    for gr in group_results:
        for metric_name,m in gr['metrics'].items():
            ci=m.get('bootstrap_ci95')
            flat.append({**{k:gr[k] for k in ['n','U','W','f','normalization','n_seeds','role']},
                'window':f"{gr['window'][0]}-{gr['window'][1]}",'metric':metric_name,
                'estimate':m['estimate'],'ci_low':None if ci is None else ci[0],'ci_high':None if ci is None else ci[1]})
    with (out/'group_estimates.csv').open('w',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(flat[0]));writer.writeheader();writer.writerows(flat)
    manifest={'status':'PASS','created_utc':datetime.now(timezone.utc).isoformat(),
        'python':sys.version,'numpy':np.__version__,'platform':platform.platform(),
        'source_sha256':sha(Path(__file__)),'input_sha256':inputs,
        'output_sha256':{p.name:sha(p) for p in out.iterdir() if p.is_file()},
        'wall_seconds':time.perf_counter()-start,'process_cpu_seconds':time.process_time()-cpu,
        'exit_code':0,'no_new_dynamics':True,'no_regenerated_predictions':True}
    (out/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':'PASS','alignment':result['alignment'],'group_windows':len(group_results),
        'wall_seconds':manifest['wall_seconds'],'process_cpu_seconds':manifest['process_cpu_seconds']},ensure_ascii=False))
    for gr in group_results:
        if gr['n']==256 and gr['U']==2 and gr['W']==.04 and gr['window']==[12,24]:
            print(json.dumps(gr,ensure_ascii=False))


if __name__=='__main__':main()
