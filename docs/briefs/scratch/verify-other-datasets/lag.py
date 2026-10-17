import pandas as pd, re, numpy as np
toi=pd.read_csv('dl/toi.csv'); ctoi=pd.read_csv('dl/ctoi.csv')
print('TOI rows',len(toi),'CTOI rows',len(ctoi))
print('TOI T>13.5',(toi['TESS Mag']>13.5).sum(),'>15',(toi['TESS Mag']>15).sum(),'>16',(toi['TESS Mag']>16).sum())
print(toi[toi['TESS Mag']>13.5]['TFOPWG Disposition'].value_counts(dropna=False).to_dict())
print('CTOI T>13.5',(ctoi['TESS Mag']>13.5).sum(),'>15',(ctoi['TESS Mag']>15).sum(),'>16',(ctoi['TESS Mag']>16).sum())
print('CTOI median T',ctoi['TESS Mag'].median(),'75pct',ctoi['TESS Mag'].quantile(.75))
print('CTOI promoted',ctoi['Promoted to TOI'].notna().sum(), ctoi['Promoted to TOI'].astype(str).str.strip().replace('nan','').ne('').sum())
print('faint CTOI users',ctoi[ctoi['TESS Mag']>13.5]['User'].value_counts().head(8).to_dict())
print('top CTOI users',ctoi['User'].value_counts().head(10).to_dict())
# sector end dates (HEASARC, fetched 2026-10-17)
ends={86:'2025-01-14',87:'2025-02-11',88:'2025-03-12',89:'2025-04-09',90:'2025-05-07',91:'2025-06-01',92:'2025-06-29',93:'2025-07-25',94:'2025-08-20',95:'2025-08-20',96:'2025-09-15',97:'2025-11-09',98:'2026-01-05',99:'2026-02-02',100:'2026-03-01',101:'2026-03-27',102:'2026-04-21',103:'2026-05-17',104:'2026-06-13',105:'2026-07-11',106:'2026-08-09',107:'2026-09-07'}
# note: 94 end date placeholder
toi['alert']=pd.to_datetime(toi['Date TOI Alerted (UTC)'],errors='coerce')
src=toi['Source'].astype(str)
print('latest alert',toi['alert'].max())
print(toi.sort_values('alert',ascending=False)[['TOI','Source','Sectors','TESS Mag','alert']].head(15).to_string())
# single-sector spoc
m=src.str.extract(r'^spoc-s(\d{2,3})-b0\d$')[0]
toi['spoc_sec']=pd.to_numeric(m)
g=toi.dropna(subset=['spoc_sec']).groupby('spoc_sec')['alert'].min()
for s in range(90,110):
    if s in g.index:
        e=pd.Timestamp(ends.get(s,'NaT')) if s in ends else pd.NaT
        print('SPOC single sector',s,'first alert',g[s].date(),'lag d',(g[s]-e).days if s in ends else 'NA', 'n',(toi['spoc_sec']==s).sum())
m2=src.str.extract(r'^qlp-s(\d{2,3})-faintsearch')[0]
toi['qlpf']=pd.to_numeric(m2)
g2=toi.dropna(subset=['qlpf']).groupby('qlpf')['alert'].min()
for s in range(85,110):
    if s in g2.index:
        print('QLP faint',s,'first alert',g2[s].date(),'lag',(g2[s]-pd.Timestamp(ends[s])).days if s in ends else 'NA','n',(toi['qlpf']==s).sum())
m3=src.str.extract(r'^qlp-s(\d{2,3})-ffi')[0]
toi['qlpffi']=pd.to_numeric(m3)
g3=toi.dropna(subset=['qlpffi']).groupby('qlpffi')['alert'].min()
for s in range(85,110):
    if s in g3.index: print('QLP ffi',s,'first alert',g3[s].date(),'n',(toi['qlpffi']==s).sum())
print('Source values 2026:',toi[toi['alert']>='2026-01-01']['Source'].value_counts().head(25).to_dict())
print('2026 TOIs by month',toi[toi['alert']>='2026-01-01']['alert'].dt.month.value_counts().sort_index().to_dict())
print('median T 2026',toi[toi['alert']>='2026-01-01']['TESS Mag'].median())
# multi-sector spoc faint
ms=toi[src.str.contains('spoc-s14-s86')]
print('spoc-s14-s86 n',len(ms),'Tmax',ms['TESS Mag'].max(),'alert range',ms['alert'].min(),ms['alert'].max())
# any TOI with sectors including 103-107?
def hassec(x,lo=103,hi=107):
    try: return any(lo<=int(t)<=hi for t in str(x).split(','))
    except: return False
print('TOIs whose Sectors field includes 103-107:',toi['Sectors'].apply(hassec).sum(), toi[toi['Sectors'].apply(hassec)][['TOI','Source','Sectors','alert']].head(10).to_string())
