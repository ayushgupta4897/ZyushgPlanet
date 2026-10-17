import os; os.environ['OMP_NUM_THREADS']='1'
import numpy as np, multiprocessing as mp, time, warnings; warnings.filterwarnings('ignore')
def work(use_threads):
    from transitleastsquares import transitleastsquares
    rng=np.random.default_rng(0); t=np.arange(0,27,1/48); f=1+rng.normal(0,1e-3,len(t))
    f[(t%2.5)<0.08]-=0.01
    try:
        r=transitleastsquares(t,f).power(period_min=1,period_max=5,use_threads=use_threads,show_progress_bar=False,verbose=False)
        return (use_threads,'ok',round(float(r.period),4),round(float(r.SDE),1))
    except Exception as e:
        return (use_threads,'ERR',repr(e)[:200])
if __name__=='__main__':
    t0=time.time()
    with mp.Pool(2) as pool:
        for res in pool.imap_unordered(work,[1,2]): print(res, round(time.time()-t0,1),'s', flush=True)
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(2) as ex:
        for res in ex.map(work,[2]): print('PPE',res, round(time.time()-t0,1),'s')
