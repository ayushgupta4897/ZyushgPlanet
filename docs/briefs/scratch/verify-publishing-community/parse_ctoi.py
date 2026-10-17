import csv,sys,collections
rows=list(csv.DictReader(open(sys.argv[1],encoding='utf-8',errors='replace')))
print('CTOI rows',len(rows))
prom=[r for r in rows if r['Promoted to TOI'].strip()]
print('promoted',len(prom), '%.1f%%'%(100*len(prom)/len(rows)))
print('Discovery sources',collections.Counter(r['Discovery Data Source'] for r in rows).most_common(5))
print('TFOPWG disp',collections.Counter(r['TFOPWG Disposition'] for r in rows).most_common(8))
withpaper=[r for r in rows if r['Paper'].strip()]
print('with paper',len(withpaper))
bydate=collections.Counter(r['CTOI lastmod'][:10] for r in withpaper)
print('paper rows by lastmod date',sorted(bydate.items()))
pre=[r for r in withpaper if r['CTOI lastmod'][:10]<'2026-08-19']
print('paper rows pre 2026-08-19',len(pre))
for r in withpaper:
    if r['CTOI lastmod'][:10]>='2026-09-01':
        print({k:r[k] for k in ['TIC ID','CTOI','Discovery Data Source','Candidate Name','User','Group','Tag','Paper','Notes','CTOI lastmod','Period (days)','Transit Epoch (BJD)','Depth ppm','Duration (hrs)']})
users=collections.Counter(r['User'] for r in rows)
print('top users',users.most_common(12))
for u in ['eisner','olmschenk','montalto','mason','kaplan']:
    rs=[r for r in rows if r['User'].lower()==u]
    print(u,len(rs),'promoted',sum(1 for r in rs if r['Promoted to TOI'].strip()))
for tic in ['4206066','370612437','252721481','120686465','22760643','392664683']:
    print(tic,[ (r['CTOI'],r['User'],r['CTOI lastmod'],r['Promoted to TOI'],r['Paper']) for r in rows if r['TIC ID']==tic])
print('lastmod dates after 2026-08-01', sorted(collections.Counter(r['CTOI lastmod'][:10] for r in rows if r['CTOI lastmod'][:10]>='2026-03-01').items()))
# max lastmod
print('max lastmod',max(r['CTOI lastmod'] for r in rows))
