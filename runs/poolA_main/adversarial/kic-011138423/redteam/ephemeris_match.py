"""Ephemeris / proximity match of the kic-011138423 signal against data/known/ catalogues (red team, own code)."""
import numpy as np, pandas as pd, json
P = 4.894451; t0 = 125.33663  # BKJD
ra0, dec0 = 292.3076266, 48.7511453
K = '/home/user/ZyushgPlanet/data/known/'
out = {'period_matches': [], 'within_5_arcmin': {}}

def sep_arcmin(ra, dec):
    return np.hypot((ra - ra0) * np.cos(np.radians(dec0)), dec - dec0) * 60

koi = pd.read_csv(K + 'koi_cumulative.csv'); tce = pd.read_csv(K + 'dr25_tce.csv')
eb = pd.read_csv(K + 'J_AJ_151_68_catalog.csv'); ps = pd.read_csv(K + 'ps_keplerfield.csv')

def check(name, df, pcol, tcol, racol, deccol, idcols):
    p = df[pcol].values.astype(float)
    for num, den in [(1, 1), (1, 2), (2, 1), (1, 3), (3, 1), (2, 3), (3, 2), (1, 4), (4, 1)]:
        m = np.abs(p / (P * num / den) - 1) < 2e-3
        for i in np.where(m)[0]:
            r = df.iloc[i]
            tt = float(r[tcol]) if tcol and np.isfinite(float(r[tcol])) else np.nan
            if name == 'EB':  # BJD0 is BJD-2400000 in Kirk+2016 -> BKJD = BJD - 2454833
                tt = tt + 2400000 - 2454833 if np.isfinite(tt) else np.nan
            ph = ((tt - t0) / P) % 1 if np.isfinite(tt) else np.nan
            out['period_matches'].append({'cat': name, 'harmonic': f'{num}/{den}', 'id': ' '.join(str(r[c]) for c in idcols),
                                          'P': float(p[i]), 'dP_over_P': float(p[i] / (P * num / den) - 1), 't0_bkjd': tt,
                                          'phase_in_signal_ephemeris': float(ph) if np.isfinite(ph) else None,
                                          'sep_arcmin': float(sep_arcmin(float(r[racol]), float(r[deccol])))})
check('KOI', koi, 'koi_period', 'koi_time0bk', 'ra', 'dec', ['kepid', 'kepoi_name', 'koi_disposition'])
check('TCE', tce, 'tce_period', 'tce_time0bk', 'ra', 'dec', ['kepid', 'tce_plnt_num'])
check('EB', eb, 'Per', 'BJD0', '_RA', '_DE', ['KIC'])
for name, df, rc, dc in [('KOI', koi, 'ra', 'dec'), ('TCE', tce, 'ra', 'dec'), ('EB', eb, '_RA', '_DE'), ('PS', ps, 'ra', 'dec')]:
    s = sep_arcmin(df[rc].values.astype(float), df[dc].values.astype(float))
    out['within_5_arcmin'][name] = int(np.sum(s < 5))
pm = pd.DataFrame(out['period_matches'])
print('objects within 5 arcmin:', out['within_5_arcmin'])
print('period/harmonic matches (|dP/P|<2e-3):', len(pm))
if len(pm):
    pm['coherent_over_300_epochs'] = np.abs(pm['dP_over_P']) < 1e-4
    pm['phase_within_0.02_of_0_or_0.5'] = pm['phase_in_signal_ephemeris'].apply(lambda x: x is not None and (min(abs(x), abs(1 - x)) < 0.02 or abs(x - 0.5) < 0.02))
    print(pm.sort_values(['cat', 'harmonic', 'sep_arcmin']).to_string())
    print('\nmatches that are BOTH period-coherent (|dP/P|<1e-4) AND phase-aligned (within 0.02 of 0 or 0.5):',
          int((pm['coherent_over_300_epochs'] & pm['phase_within_0.02_of_0_or_0.5']).sum()))
    print('nearest period match of any kind: %.1f arcmin' % pm['sep_arcmin'].min())
json.dump(out, open('/home/user/ZyushgPlanet/runs/poolA_main/adversarial/kic-011138423/redteam/ephemeris_match.json', 'w'), indent=1, default=float)
