#!/usr/bin/env python3
"""
sky_scan.py — pre-analisi automatica di una giornata SKY

Uso:  python3 sky_scan.py ~/Desktop/SKY_2026-09-23
Legge: TUTTE_<data>.csv (+ NOAA_*.json se presenti)
Scrive nella cartella: SCAN_<data>.txt  (report da incollare all'AI)
Aggiorna nella cartella madre: SKY_REGISTRO.csv (una riga per giorno)

Regola: il report elenca OSSERVAZIONI e SEGNALAZIONI, mai verdetti.
"""
import sys, json, glob, math, pathlib
import numpy as np
import pandas as pd

EU      = ['WIC', 'LON', 'THY', 'CLF', 'BFO', 'NGK', 'BEL', 'IZN', 'DUR']
AURORAL = ['ABK', 'TRO']
GLOBAL  = ['KAK', 'FRD', 'OTT', 'BRD']
DUR_ROT_DEG = 4.27
# stazioni che pubblicano in un sistema ruotato (verificato con IGRF 2026): angolo da applicare
ROT = {'DUR': 4.27, 'FRD': -10.49}
SYNC_MIN_ST = 5          # stazioni EU per un evento "sincrono"
ROC_THR     = 8.0        # nT/30min (come sky.py)
WIN22 = ('21:45', '21:55', '22:05', '22:30')   # baseline / finestra fascia 22 UT
L1_DEFAULT_KM = 1.5e6
DERIVED_F = set()
NOAA_COVER = None
DRIFT_EST = set()
ACTIVE_SH = {}
CULM = {}
REG_LAGS = {}
SYNC_LOG = []          # eventi sincroni della giornata (per registro eventi e analisi aggiuntive)

# coordinate (lat, lon) — SD = punto di osservazione
COORD = {'SD': (46.155, 13.005), 'WIC': (47.93, 15.87), 'LON': (45.41, 16.66), 'THY': (46.90, 17.89),
         'CLF': (48.02, 2.26), 'BFO': (48.33, 8.33), 'NGK': (52.07, 12.68), 'BEL': (51.84, 20.79),
         'IZN': (40.50, 29.73), 'DUR': (41.65, 14.47), 'ABK': (68.36, 18.82), 'TRO': (69.66, 18.94),
         'KAK': (36.23, 140.19), 'FRD': (38.20, -77.37), 'OTT': (45.40, -75.55), 'BRD': (49.87, -99.97)}

# radianti (AR, Dec in gradi) per nome come compare nel RIEPILOGO
RADIANTS = {'Perseidi': (48, 58), 'Eps Perseidi': (48, 40), 'Alpha Aurigidi': (91, 39),
            'Kappa Cignidi': (286, 59), 'Delta Aquaridi': (340, -16), 'Alpha Capric.': (307, -10),
            'Sextantidi': (152, 0), 'Draconidi': (262, 54), 'Orionidi': (95, 16),
            'Tauridi Sud': (32, 9), 'Tauridi Nord': (58, 22), 'Leonidi': (152, 22),
            'Geminidi': (112, 33), 'Quadrantidi': (230, 49), 'Liridi': (271, 34), 'Eta Aquaridi': (338, -1)}
SIDEREAL_SHIFT_MIN = -3.93    # min/giorno in UT per un fenomeno legato alle stelle
LUNAR_SHIFT_MIN    = 50.5     # min/giorno in UT per un fenomeno legato alla Luna

# ─────────────────────────────────────────────────────────────────────────────
def load_all(folder):
    f = sorted(glob.glob(str(folder / 'TUTTE_*.csv')))
    if not f:
        sys.exit('TUTTE_*.csv non trovato nella cartella')
    df = pd.read_csv(f[0], comment='#')
    df['t'] = pd.to_datetime(df.time_UT)
    df = df.set_index('t').drop(columns='time_UT')
    day = df.index[0].strftime('%Y-%m-%d')
    D = {}
    for col in df.columns:
        st, comp = col.split('_')
        if st == 'DURg':
            continue
        D.setdefault(st, pd.DataFrame(index=df.index))[comp] = df[col]
    for st, ang in ROT.items():                 # riporto nel sistema geografico
        if st in D:
            th = math.radians(ang); x, y = D[st].X.copy(), D[st].Y.copy()
            D[st]['X'] = x * math.cos(th) - y * math.sin(th); D[st]['Y'] = x * math.sin(th) + y * math.cos(th)
    D = {k: v for k, v in D.items() if v[['X', 'Y', 'Z']].notna().sum().sum() > 0}
    return day, D

def runs(mask):
    """lista (start, end, lunghezza) dei tratti True consecutivi"""
    out, s, prev, n = [], None, None, 0
    for t, v in mask.items():
        if v:
            if s is None: s, n = t, 0
            n += 1
        elif s is not None:
            out.append((s, prev, n)); s = None
        prev = t
    if s is not None: out.append((s, prev, n))
    return out

def hm(t): return t.strftime('%H:%M')

# ─────────────────────────────────────────────────────────────────────────────
def sec_integrity(D, L, reg):
    L.append('\n══ 1. INTEGRITÀ DEI DATI ══')
    for s, d in D.items():
        n = int(d.X.notna().sum()); nan = len(d) - n
        msg = [f'{s}: {n}/{len(d)} minuti']
        gaps = [r for r in runs(d.X.isna()) if r[2] >= 3]
        if gaps: msg.append('buchi ' + ' '.join(f'{hm(a)}–{hm(b)}' for a, b, _ in gaps[:6]))
        # outlier grossolani
        for c in 'XYZF':
            if c in d and d[c].notna().any():
                bad = (d[c] - d[c].median()).abs() > (3000 if s in AURORAL else 500)
                if bad.any(): msg.append(f'⚠ {c} fuori scala in {int(bad.sum())} min ({hm(d[c][bad].index[0])})')
        # risoluzione e tratti piatti
        frac = ((d.Y.dropna() * 100) % 10).round().astype(int)
        res = 0.1 if frac.nunique() <= 2 else 0.01
        for c in 'XYZ':
            fl = [r for r in runs(d[c].diff() == 0) if r[2] + 1 >= (20 if res == 0.1 else 8)]
            if fl:
                a, b, k = max(fl, key=lambda r: r[2])
                msg.append(f'{c} piatto {k+1} min fino alle {hm(b)}')
        # F misurato o derivato
        if 'F' in d and d.F.notna().sum() > 100:
            dF = d.F - np.sqrt(d.X**2 + d.Y**2 + d.Z**2)
            m, sd = dF.mean(), dF.std()
            tag = ' ← F probabilmente CALCOLATO (scarto troppo perfetto)' if sd < 0.06 and abs(m) < 0.05 else ''
            msg.append(f'F−|B| {m:+.2f}±{sd:.2f}{tag}')
            reg[f'{s}_dF'] = round(m, 2)
            if tag: DERIVED_F.add(s)
        msg.append(f'risoluzione ~{res} nT')
        L.append('  ' + ' | '.join(msg))

def sec_jumps(D, L, reg):
    L.append('\n══ 2. SALTI SU SINGOLA STAZIONE (1 min, con controllo su F) ══')
    L.append('  (F segue = variazione reale vicino al sensore; F fermo = strumentale; F derivato = test non valido)')
    for s, d in D.items():
        if s not in EU: continue     # aurorali ed extra-EU: salti propri attesi, non 'locali'
        hits = []
        for c in 'XYZ':
            d1 = d[c].diff()
            mad = np.nanmedian(np.abs(d1 - np.nanmedian(d1))) * 1.4826
            thr = max(6 * mad, 1.5)
            for t, v in d1[d1.abs() > thr].items():
                others = [o for o in D if o != s and o in EU and abs(D[o][c].diff().get(t, 0) or 0) > thr * 0.5]
                if len(others) >= 2:
                    continue                      # è un evento di rete, non locale
                fv = d.F.diff().get(t, np.nan) if 'F' in d else np.nan
                ftag = '' if np.isnan(fv) else (f' F{fv:+.1f}' + ('✓' if abs(fv) > 0.5 * abs(v) else '✗'))
                if s in DERIVED_F: ftag = ' (F derivato)'
                hits.append(f'{hm(t)} {c}{v:+.1f}{ftag}')
        if hits:
            L.append(f'  {s}: ' + ' '.join(hits[:15]) + (' …' if len(hits) > 15 else ''))
        reg[f'{s}_jumps'] = len(hits)

def sync_events(D, stations, comp, kind):
    idx = D[stations[0]].index
    if kind == '1min':
        M = pd.concat([(D[s][comp].diff().abs() > D[s][comp].diff().abs().quantile(.99)).astype(int)
                       for s in stations], axis=1).fillna(0)
        n = M.sum(axis=1)
        return [(t, int(v)) for t, v in n[n >= SYNC_MIN_ST].items()]
    out = []
    R = {s: D[s][comp].diff(30).reindex(idx) for s in stations}
    for sg, f in (('+', lambda v: v > ROC_THR), ('−', lambda v: v < -ROC_THR)):
        m = sum(f(R[s]).astype(int) for s in stations) >= SYNC_MIN_ST
        for a, b, k in runs(m):
            out.append((a, b, sg))
    return sorted(out)

def sec_sync(D, L, reg, arrivals):
    eu = [s for s in EU if s in D]
    ext = [s for s in AURORAL + GLOBAL if s in D]
    L.append(f'\n══ 3. EVENTI SINCRONI (≥{SYNC_MIN_ST} stazioni EU: {", ".join(eu)}) ══')
    tot = 0
    for c in 'XYZ':
        ev = sync_events(D, eu, c, '1min')
        # raggruppa minuti consecutivi
        groups = []
        for t, n in ev:
            if groups and (t - groups[-1][1]).total_seconds() <= 120:
                groups[-1][1] = t; groups[-1][2] = max(groups[-1][2], n)
            else:
                groups.append([t, t, n])
        tot += len(groups)
        if groups:
            desc = []
            for a, b, n in groups:
                big = max(eu, key=lambda s: abs(D[s][c].diff().loc[a:b].abs().max() or 0))
                val = D[big][c].diff().loc[a:b]
                v = val.loc[val.abs().idxmax()]
                x = ''
                for e in ext:                      # risposta fuori Europa
                    dv = D[e][c].diff().loc[a - pd.Timedelta(minutes=3):b + pd.Timedelta(minutes=3)]
                    if dv.notna().any(): x += f' {e}{dv.loc[dv.abs().idxmax()]:+.1f}'
                arr = ''
                for ta, lab in arrivals:           # arrivo di discontinuità dal vento solare
                    if abs((ta - a).total_seconds()) <= 900: arr = f' ⟵ L1 {lab} (arrivo {hm(ta)})'
                covered = NOAA_COVER is not None and NOAA_COVER[0] <= a <= NOAA_COVER[1]
                if not covered: arr = arr or ' (vento solare: fuori copertura)'
                SYNC_LOG.append({'t': a, 'end': b, 'comp': c, 'n': n, 'st': big, 'val': round(float(v), 2), 'orphan': (int(not arr) if covered else -1)})
                desc.append(f'{hm(a)}{"–"+hm(b) if b != a else ""}({n}) max {big}{v:+.1f}' + (f' |fuori EU:{x}' if x else '') + arr)
            L.append(f'  {c} al minuto: ' + ' ; '.join(desc))
    for c in 'XYZ':
        ev = sync_events(D, eu, c, 'roc30')
        if ev:
            L.append(f'  {c} su 30 min (|Δ|>{ROC_THR}): ' + ' '.join(f'{hm(a)}–{hm(b)}{sg}' for a, b, sg in ev))
    reg['sync_1min'] = tot

def sec_diurnal(D, L, reg):
    eu = [s for s in EU if s in D]
    L.append('\n══ 4. CURVA DIURNA ══')
    amp = {s: D[s].Y.max() - D[s].Y.min() for s in eu}
    L.append('  Ampiezza Y: ' + ' '.join(f'{s}{amp[s]:.0f}' for s in eu))
    ymin = {s: D[s].Y.idxmin() for s in eu}
    ymax = {s: D[s].Y.idxmax() for s in eu}
    L.append('  Min Y: ' + ' '.join(f'{s}{hm(ymin[s])}' for s in sorted(eu, key=lambda s: ymin[s])))
    L.append('  Max Y: ' + ' '.join(f'{s}{hm(ymax[s])}' for s in sorted(eu, key=lambda s: ymax[s])))
    odd = [s for s in eu if ymax[s].hour >= 14 or ymax[s].hour < 3]
    if odd: L.append(f'  ⚠ massimo di Y fuori dal mattino: {", ".join(odd)}')
    # raggruppamenti: stazioni con minimo entro 2 minuti
    ts = sorted(ymin.items(), key=lambda kv: kv[1]); groups = []
    for s, t in ts:
        if groups and (t - groups[-1][-1][1]).total_seconds() <= 120: groups[-1].append((s, t))
        else: groups.append([(s, t)])
    big = [g for g in groups if len(g) >= 3]
    L.append('  Gruppi (≥3 stazioni entro 2 min): ' +
             ('; '.join(f'{hm(g[0][1])} ' + ','.join(s for s, _ in g) for g in big) if big else 'nessuno'))
    zmin = {s: D[s].Z.idxmin() for s in eu}
    med = sorted(zmin.values())[len(zmin) // 2]
    L.append('  Min Z: ' + ' '.join(f'{s}{hm(zmin[s])}' for s in sorted(eu, key=lambda s: zmin[s])))
    if 'IZN' in zmin:
        lead = (zmin['IZN'] - med).total_seconds() / 3600
        L.append(f'  IZN minimo Z rispetto alla mediana: {lead:+.1f} h')
        reg['IZN_Zlead_h'] = round(lead, 1)
    reg['Yamp_min'] = round(min(amp.values()), 1); reg['Yamp_max'] = round(max(amp.values()), 1)
    reg['Ymin_groups'] = len(big)
    if big: reg['Ymin_group1'] = hm(big[0][0][1])

def sec_night_evening(D, L, reg, day):
    eu = [s for s in EU if s in D]
    L.append('\n══ 5. NOTTE, OSCILLAZIONI, FASCIA 22 UT ══')
    ns = {s: D[s].between_time('00:00', '03:00').Y.std() for s in eu}
    L.append('  Agitazione notte (std Y 00–03): ' + ' '.join(f'{s}{v:.1f}' for s, v in ns.items()))
    ref = 'NGK' if 'NGK' in D else eu[0]
    y = D[ref].Y.interpolate()
    hp = (y - y.rolling(21, center=True).mean()).rolling(3, center=True).mean()
    pw = hp.pow(2).resample('1h').mean()
    L.append(f'  Potenza oscillazioni {ref} per ora: ' + ' '.join(f'{v:.2f}' for v in pw))
    L.append(f'  Ora più "vibrante": {pw.idxmax().hour:02d} UT ({pw.max():.2f})')
    ext, mean = {}, {}
    for s in eu + [e for e in AURORAL + GLOBAL if e in D]:
        x = D[s].X
        a = x.loc[f'{day} {WIN22[0]}':f'{day} {WIN22[1]}'].mean()
        b = x.loc[f'{day} {WIN22[2]}':f'{day} {WIN22[3]}'] - a
        if b.notna().any():
            ext[s] = b.loc[b.abs().idxmax()]; mean[s] = b.mean()
    em = np.median([ext[s] for s in eu if s in ext])
    L.append(f'  Fascia 22 UT, picco X (criterio fisso): mediana EU {em:+.1f} → segno {"+" if em > 0 else "−"}')
    L.append('    per stazione: ' + ' '.join(f'{s}{v:+.1f}' for s, v in ext.items()))
    split = len({np.sign(ext[s]) for s in eu if s in ext}) > 1
    if split: L.append('    ⚠ stazioni EU con segno discorde (rete divisa)')
    reg['night_std_med'] = round(float(np.median(list(ns.values()))), 2)
    reg['osc_max_hour'] = int(pw.idxmax().hour)
    reg['X22_ext'] = round(float(em), 1)
    reg['X22_split'] = int(split)

def sec_trend(D, L, reg):
    L.append('\n══ 6. MEDIE GIORNALIERE ══')
    L.append('  Z medio: ' + ' '.join(f'{s}{D[s].Z.mean():.1f}' for s in D))
    for s in ('THY', 'NGK', 'IZN', 'LON'):
        if s in D: reg[f'{s}_Zmean'] = round(D[s].Z.mean(), 1)

# ─────────────────────────────────────────────────────────────────────────────
def load_noaa(folder):
    def latest(key):
        recs = []
        for f in sorted(glob.glob(str(folder / f'NOAA_{key}_*.json'))):
            try: recs += json.load(open(f))
            except Exception: pass
        if not recs: return None
        seen, out = set(), []
        for r in recs:                           # unione di più download, senza doppioni
            k = (r.get('time_tag'), r.get('source'))
            if k not in seen: seen.add(k); out.append(r)
        return out
    w, m, e = latest('wind'), latest('mag'), latest('eph')
    if not w:
        return None
    w = pd.DataFrame(w); w['t'] = pd.to_datetime(w.time_tag)
    w = w[w.active == True].set_index('t').sort_index()
    for c in ('proton_speed', 'proton_density'):
        w[c] = pd.to_numeric(w[c], errors='coerce')
    w = w[(w.proton_speed > 200) & (w.proton_speed < 1500) & (w.proton_density > 0.1) & (w.proton_density < 100)]
    w['P'] = 1.6726e-6 * w.proton_density * w.proton_speed ** 2
    mm = None
    if m:
        mm = pd.DataFrame(m); mm['t'] = pd.to_datetime(mm.time_tag)
        mm = mm[mm.active == True].set_index('t').sort_index()
        for c in ('bt', 'bz_gsm'): mm[c] = pd.to_numeric(mm[c], errors='coerce')
    dist = L1_DEFAULT_KM
    if e:
        ee = pd.DataFrame(e)
        xk = [k for k in ee.columns if 'x' in k.lower() and 'gse' in k.lower()]
        if xk:
            v = pd.to_numeric(ee[xk[0]], errors='coerce').dropna()
            if len(v):
                v = float(v.iloc[-1]); dist = v * 6371 if abs(v) < 1000 else abs(v)
    return w, mm, dist

def sec_noaa(folder, L, reg, day):
    L.append('\n══ 7. VENTO SOLARE (NOAA, sorgente attiva, propagato a Terra) ══')
    r = load_noaa(folder)
    if r is None:
        L.append('  file NOAA assenti'); return []
    w, mm, dist = r
    src = w.source.mode().iloc[0] if 'source' in w else '?'
    v = w.proton_speed.resample('1min').mean().interpolate(limit=10)
    lag = (dist / v).rename('lag')
    P = w.P.resample('1min').median().interpolate(limit=5).rolling(5, center=True, min_periods=2).median()
    L.append(f'  sorgente {src}, distanza {dist/1e6:.2f} Mkm, velocità {v.median():.0f} km/s, ritardo tipico {lag.median()/60:.0f} min')
    arrivals = []
    dP = (P.shift(-5) - P.shift(1)) / P.shift(1)
    for a, b, k in runs(dP.abs() > 0.25):
        t = dP.loc[a:b].abs().idxmax(); val = dP[t]
        # persistenza: livello medio 10 min dopo vs 10 min prima
        pre, post = P.loc[t - pd.Timedelta(minutes=12):t - pd.Timedelta(minutes=2)].mean(), P.loc[t + pd.Timedelta(minutes=5):t + pd.Timedelta(minutes=15)].mean()
        if np.isnan(pre) or np.isnan(post) or abs(post / pre - 1) < 0.2:
            continue
        ta = t + pd.Timedelta(seconds=float(lag.get(t, lag.median())))
        arrivals.append((ta, f'P{post/pre-1:+.0%}'))
    if mm is not None:
        bt = mm.bt.resample('1min').mean()
        d3 = bt.diff(3)
        for a, b, k in runs(d3.abs() > 2.5):
            t = d3.loc[a:b].abs().idxmax()
            ta = t + pd.Timedelta(seconds=float(lag.get(t, lag.median())))
            if not any(abs((ta - x).total_seconds()) < 300 for x, _ in arrivals):
                arrivals.append((ta, f'Bt{d3[t]:+.1f}'))
        bz = mm.bz_gsm.resample('1min').mean()
        bzs = bz.rolling(10, min_periods=5).mean()
        turn = bzs.shift(-15) - bzs
        for a_, b_, k_ in runs(turn.abs() > 3.0):
            t = turn.loc[a_:b_].abs().idxmax()
            pre = bzs.loc[t - pd.Timedelta(minutes=10):t]; post = bzs.loc[t:t + pd.Timedelta(minutes=20)]
            if len(pre) and len(post) and ((pre.min() < 0 < post.max()) if turn[t] > 0 else (pre.max() > 0 > post.min())):
                ta = t + pd.Timedelta(seconds=float(lag.get(t, lag.median())))
                if not any(abs((ta - x).total_seconds()) < 600 for x, _ in arrivals):
                    arrivals.append((ta, f'svolta Bz {"→nord" if turn[t] > 0 else "→sud"}'))
        south = [r for r in runs(bz < -5) if r[2] >= 20]
        if south:
            L.append('  Bz < −5 nT per ≥20 min (ora L1): ' + ' '.join(f'{a:%d/%m} {hm(a)}–{hm(b)}' for a, b, _ in south))
    d0 = pd.Timestamp(day); d1 = d0 + pd.Timedelta(days=1)
    arrivals = sorted(a for a in arrivals if d0 <= a[0] < d1)
    L.append(f'  copertura L1 → Terra: {hm(P.index[0]+pd.Timedelta(seconds=float(lag.median())))} del {(P.index[0]+pd.Timedelta(seconds=float(lag.median()))).strftime("%d/%m")} in poi')
    L.append(f'  Discontinuità a L1 (arrivo stimato a Terra): ' +
             (' '.join(f'{hm(t)}[{lab}]' for t, lab in arrivals) if arrivals else 'nessuna netta'))
    L.append('  (pressione su → X su atteso; pressione giù → X giù)')
    reg['L1_disc'] = len(arrivals)
    global NOAA_COVER
    NOAA_COVER = (P.index[0] + pd.Timedelta(seconds=float(lag.median())), P.index[-1] + pd.Timedelta(seconds=float(lag.median())))
    return arrivals

# ─────────────────────────────────────────────────────────────────────────────
def sec_registry(folder, day, reg, L):
    rp = folder.parent / 'SKY_REGISTRO.csv'
    old = pd.read_csv(rp) if rp.exists() else pd.DataFrame()
    if len(old): old = old[old.day != day]
    row = pd.DataFrame([{'day': day, **reg}])
    new = pd.concat([old, row], ignore_index=True).sort_values('day')
    new.to_csv(rp, index=False)
    L.append('\n══ 14. CONFRONTO CON I GIORNI PRECEDENTI (SKY_REGISTRO.csv) ══')
    prev = new[new.day < day].tail(10)
    if prev.empty:
        L.append('  registro vuoto: primo giorno'); return
    if 'X22_ext' in prev:
        seq = ' '.join(('+' if v > 0 else '−') for v in list(prev.X22_ext) + [reg.get('X22_ext', 0)])
        L.append(f'  Segno fascia 22 UT (ultimi giorni → oggi): {seq}')
        alt = all(np.sign(a) != np.sign(b) for a, b in zip(list(prev.X22_ext)[-1:], [reg.get('X22_ext', 0)]))
        L.append(f'  Alternanza rispetto a ieri: {"SÌ" if alt else "NO"}')
    for k in ('THY_Zmean', 'Yamp_max', 'night_std_med', 'sync_1min', 'IZN_Zlead_h'):
        if k in prev and k in reg:
            vals = list(prev[k].dropna())[-5:]
            L.append(f'  {k}: ' + ' → '.join(str(v) for v in vals) + f' → OGGI {reg[k]}')

# ─── ANALISI NON CONVENZIONALI (solo numeri, nessuna etichetta causale) ─────────

def sec_regional(D, L, reg):
    """Scarti dal movimento comune della rete: né locale (1 stazione) né di rete (tutte)."""
    eu = [s for s in EU if s in D]
    L.append('\n══ 9. MOVIMENTI REGIONALI (2–4 stazioni insieme, contro il resto della rete) ══')
    found = 0
    for c in 'XYZ':
        d1 = pd.concat([D[s][c].diff().rename(s) for s in eu], axis=1)
        r = d1.sub(d1.median(axis=1), axis=0)
        thr = r.abs().median() * 1.4826 * 6
        big = r.abs().gt(thr.clip(lower=0.8), axis=1)
        items = []
        for t in big.index[big.sum(axis=1).between(2, 4)]:
            sts = [s for s in eu if big.at[t, s]]
            sg = {np.sign(r.at[t, s]) for s in sts}
            if len(sg) == 1:
                items.append(f'{hm(t)} ' + ','.join(sts) + f'{"+" if sg.pop() > 0 else "−"}')
        found += len(items)
        if items:
            L.append(f'  {c}: ' + ' ; '.join(items[:20]) + (' …' if len(items) > 20 else ''))
    if not found: L.append('  nessuno')
    reg['regional'] = found

def _km(a, b):
    la, lo = COORD[a]; lb, lob = COORD[b]
    x = (lob - lo) * 111.32 * math.cos(math.radians((la + lb) / 2)); y = (lb - la) * 110.57
    return x, y

def sec_direction(D, L):
    """Per ogni evento sincrono: ritardi fra stazioni (correlazione, sotto il minuto) → direzione apparente."""
    eu = [s for s in EU if s in D and s in COORD]
    L.append('\n══ 10. DIREZIONE APPARENTE DEGLI EVENTI SINCRONI (ritardi fra stazioni) ══')
    if not SYNC_LOG: L.append('  nessun evento'); return
    ref = 'THY' if 'THY' in eu else eu[0]
    out = []
    # 1° passaggio: ritardo di ogni stazione su ogni evento → ritardo FISSO (mediana) da sottrarre
    raw = {}
    for e in SYNC_LOG:
        a, b, c = e['t'] - pd.Timedelta(minutes=15), e['end'] + pd.Timedelta(minutes=15), e['comp']
        r0 = D[ref][c].loc[a:b].diff().fillna(0).values
        for s in eu:
            v = D[s][c].loc[a:b].diff().fillna(0).values
            if len(v) != len(r0) or np.std(v) == 0: continue
            cc = np.correlate(v - v.mean(), r0 - r0.mean(), 'full'); k = int(np.argmax(cc)); m0 = len(r0) - 1
            frac = 0
            if 0 < k < len(cc) - 1:
                y0, y1, y2 = cc[k - 1], cc[k], cc[k + 1]; den = y0 - 2 * y1 + y2
                frac = 0.5 * (y0 - y2) / den if den != 0 else 0
            raw.setdefault(s, []).append(k - m0 + frac)
    FIX = {s: float(np.median(v)) for s, v in raw.items() if len(v) >= 3}
    L.append('  ritardo fisso per stazione rispetto a ' + ref + ' (mediana sugli eventi, s): ' +
             ' '.join(f'{s}{v*60:+.0f}' for s, v in FIX.items()))
    L.append('  (un ritardo fisso non è propagazione: orologio o convenzione del minuto — viene sottratto qui sotto)')
    for s, v in FIX.items(): REG_LAGS[s] = round(v * 60, 1)
    for e in SYNC_LOG:
        a, b, c = e['t'] - pd.Timedelta(minutes=15), e['end'] + pd.Timedelta(minutes=15), e['comp']
        r0 = D[ref][c].loc[a:b].diff().fillna(0).values
        xs, ys, lags = [], [], []
        for s in eu:
            v = D[s][c].loc[a:b].diff().fillna(0).values
            if len(v) != len(r0) or np.std(v) == 0: continue
            cc = np.correlate(v - v.mean(), r0 - r0.mean(), 'full'); k = int(np.argmax(cc)); m0 = len(r0) - 1
            if 0 < k < len(cc) - 1:
                y0, y1, y2 = cc[k - 1], cc[k], cc[k + 1]; den = y0 - 2 * y1 + y2
                frac = 0.5 * (y0 - y2) / den if den != 0 else 0
            else: frac = 0
            lag = (k - m0 + frac) - FIX.get(s, 0.0)
            if abs(lag) > 5: continue
            x, y = _km(ref, s); xs.append(x); ys.append(y); lags.append(lag)
        if len(lags) < 5: continue
        A = np.c_[xs, ys, np.ones(len(xs))]
        sol, *_ = np.linalg.lstsq(A, np.array(lags), rcond=None)
        sx, sy = sol[0], sol[1]; smag = math.hypot(sx, sy)
        rng = max(lags) - min(lags)
        if rng < 0.3 or smag == 0:
            out.append(f'{hm(e["t"])}{c}: simultaneo entro ±{rng/2*60:.0f} s (dopo correzione)'); continue
        if False:
            out.append(f'{hm(e["t"])}{c}: simultaneo (±{rng/2*60:.0f} s)')
        else:
            az = (math.degrees(math.atan2(sx, sy)) + 360) % 360   # direzione verso cui si propaga
            v = 1 / smag / 60
            out.append(f'{hm(e["t"])}{c}: verso {az:.0f}°, velocità apparente {v:.0f} km/s (escursione ritardi {rng*60:.0f} s)')
    L += ['  ' + o for o in out] if out else ['  ritardi non determinabili']

def _gmst_deg(t):
    jd = t.to_julian_date()
    return (280.46061837 + 360.98564736629 * (jd - 2451545.0)) % 360

def _alt(ra, dec, lat, lon, t):
    H = math.radians((_gmst_deg(t) + lon - ra) % 360)
    d, f = math.radians(dec), math.radians(lat)
    return math.degrees(math.asin(math.sin(d) * math.sin(f) + math.cos(d) * math.cos(f) * math.cos(H)))

def sec_radiants(folder, D, day, L, reg):
    """Geometria dei radianti degli sciami attivi (dal RIEPILOGO) su SD e sulle stazioni."""
    L.append('\n══ 11. GEOMETRIA DEGLI SCIAMI ATTIVI (altezza del radiante) ══')
    rp = folder / 'RIEPILOGO.txt'
    active = []
    try:
        import sky_catalog as CAT
    except ImportError:
        CAT = None
    sunl = None
    if rp.exists():
        for line in rp.read_text(encoding='utf-8').splitlines():
            if line.startswith('☉ Sole:'):
                try: sunl = float(line.split(':')[1].split('°')[0])
                except Exception: pass
            if line.startswith('☄️') and 'SCIAMI' not in line:
                name = line.replace('☄️', '').strip().split('  ')[0].strip()
                rd = CAT.radiant_radec(name, sunl) if CAT else RADIANTS.get(name)
                if CAT and rd and not name.startswith('Antielio'):
                    for sh in CAT.SHOWERS:
                        if sh['name'] == name and not sh['drift_known']: DRIFT_EST.add(name)
                if rd: RADIANTS[name] = rd; active.append(name); ACTIVE_SH[name] = rd
    if not active: L.append('  nessuno sciame con radiante noto nel RIEPILOGO'); return
    times = pd.date_range(day, periods=24 * 60, freq='1min')
    for name in active:
        ra, dec = RADIANTS[name]
        alt = pd.Series([_alt(ra, dec, *COORD['SD'], t) for t in times], index=times)
        up = alt[alt > 0]
        cul = alt.idxmax()
        rise = [t for t in alt.index[1:] if alt[t] > 0 >= alt[t - pd.Timedelta(minutes=1)]]
        sett = [t for t in alt.index[1:] if alt[t] <= 0 < alt[t - pd.Timedelta(minutes=1)]]
        CULM[name] = cul
        for e in SYNC_LOG:
            off = ((e['t'] - cul).total_seconds() / 60 + 720) % 1440 - 720
            e[f'culm_{name}'] = round(off, 1)
        L.append(f'  {name}{" (deriva stimata)" if name in DRIFT_EST else ""} (AR {ra:.1f}°, Dec {dec:+.1f}°) su SD: sorge {",".join(hm(t) for t in rise) or "—"}, '
                 f'culmina {hm(cul)} ({alt.max():.0f}°), tramonta {",".join(hm(t) for t in sett) or "—"} UT')
        if SYNC_LOG:
            evs = [f'{hm(e["t"])}:{alt[e["t"].floor("min")]:+.0f}°' for e in SYNC_LOG if e['t'].floor('min') in alt.index]
            L.append(f'    altezza del radiante su SD agli eventi sincroni: ' + ' '.join(evs))
            ups = sum(1 for e in SYNC_LOG if alt.get(e['t'].floor('min'), -1) > 0)
            frac_up = len(up) / len(alt)
            L.append(f'    eventi con radiante sopra l\'orizzonte: {ups}/{len(SYNC_LOG)} (atteso per caso: {frac_up*len(SYNC_LOG):.1f})')
            reg[f'rad_up_{name.replace(" ", "")}'] = f'{ups}/{len(SYNC_LOG)}'

def _sun_alt(lat, lon, t):
    jd = t.to_julian_date(); n = jd - 2451545.0
    Ls = (280.460 + 0.9856474 * n) % 360; g = math.radians(357.528 + 0.9856003 * n)
    lam = math.radians(Ls + 1.915 * math.sin(g) + 0.020 * math.sin(2 * g)); e = math.radians(23.44)
    ra = math.degrees(math.atan2(math.cos(e) * math.sin(lam), math.cos(lam))) % 360
    dec = math.degrees(math.asin(math.sin(e) * math.sin(lam)))
    return _alt(ra, dec, lat, lon, t)

def _spearman(a, b):
    a, b = pd.Series(a).rank(), pd.Series(b).rank()
    return float(np.corrcoef(a, b)[0, 1]) if a.std() > 0 and b.std() > 0 else np.nan

def sec_radiant_station(D, day, L, reg):
    """ESPLORATIVO: la stazione che si discosta dal resto della rete è quella con il radiante più alto?
    Confronto con due 'geometrie di controllo': altezza del Sole e longitudine (geografia fissa)."""
    eu = [s for s in EU if s in D and s in COORD]
    L.append('\n══ 15. RADIANTE SOPRA OGNI STAZIONE vs SCARTO DALLA RETE (esplorativo) ══')
    L.append('  (valori esplorativi, non verifiche: correlazione fra stazioni ora per ora; + = scarta di più dove il radiante è più alto)')
    if not ACTIVE_SH or len(eu) < 5:
        L.append('  non calcolabile'); return
    # scarto dalla rete: oscillazioni 2–30 min di ogni stazione meno la mediana della rete, ora per ora
    res = {}
    for s in eu:
        acc = 0
        for c in 'XYZ':
            v = D[s][c].interpolate(limit=5)
            hp = v - v.rolling(31, center=True, min_periods=10).mean()
            acc = acc + hp
        res[s] = acc
    R = pd.DataFrame(res)
    R = R.sub(R.median(axis=1), axis=0).abs().resample('1h').mean()
    hours = R.index
    geos = {'Sole (controllo)': lambda lat, lon, t: _sun_alt(lat, lon, t),
            'Longitudine (controllo)': lambda lat, lon, t: lon,
            'Latitudine (controllo)': lambda lat, lon, t: lat}
    for name, (ra, dec) in ACTIVE_SH.items():
        geos[name] = (lambda ra_, dec_: (lambda lat, lon, t: _alt(ra_, dec_, lat, lon, t)))(ra, dec)
    for gname, f in geos.items():
        cors, up = [], 0
        for h in hours:
            t = h + pd.Timedelta(minutes=30)
            g = [f(*COORD[s], t) for s in eu]
            if '(controllo)' not in gname:
                if max(g) <= 0: continue
                up += 1
            r = [R.at[h, s] for s in eu]
            if np.any(np.isnan(r)): continue
            c = _spearman(g, r)
            if not np.isnan(c): cors.append(c)
        if cors:
            m = np.mean(cors); pos = sum(c > 0.5 for c in cors)
            L.append(f'  {gname}: correlazione media {m:+.2f} su {len(cors)} ore' +
                     (f' (radiante sopra l\'orizzonte in {up} ore)' if up else '') + f'; ore con correlazione > 0.5: {pos}')
            reg[f'geo_{gname.split(" ")[0]}'] = round(float(m), 2)

def sec_orphans(L, reg):
    L.append('\n══ 12. EVENTI "ORFANI" (sincroni, senza discontinuità a L1 entro 15 min) ══')
    orf = [e for e in SYNC_LOG if e['orphan'] == 1]
    unc = [e for e in SYNC_LOG if e['orphan'] == -1]
    if unc: L.append(f'  fuori copertura NOAA (non valutabili): {len(unc)} eventi')
    L.append('  ' + (' '.join(f'{hm(e["t"])}{e["comp"]}({e["n"]}) {e["st"]}{e["val"]:+.1f}' for e in orf) if orf else 'nessuno'))
    L.append('  (se i file NOAA mancano o non coprono l\'orario, tutti risultano orfani: controllare la sezione 7)')
    reg['orphans'] = len(orf)

def sec_recurrence(folder, day, L):
    """Confronta gli eventi di oggi con quelli dei giorni precedenti: stesso orario UT, deriva siderale, deriva lunare."""
    ep = folder.parent / 'SKY_EVENTI.csv'
    old = pd.read_csv(ep, parse_dates=['t']) if ep.exists() else pd.DataFrame(columns=['day', 't', 'comp', 'n', 'st', 'val', 'orphan'])
    old = old[old.day != day]
    today = pd.DataFrame([{'day': day, **{k: v for k, v in e.items() if k != 'end'}} for e in SYNC_LOG])
    pd.concat([old, today], ignore_index=True).to_csv(ep, index=False)
    L.append('\n══ 13. RICORRENZE CON I GIORNI PRECEDENTI (SKY_EVENTI.csv, tolleranza ±3 min) ══')
    if old.empty or today.empty:
        L.append('  storico insufficiente'); return
    d0 = pd.Timestamp(day)
    for label, shift in (('stesso orario UT', 0.0), ('deriva siderale', SIDEREAL_SHIFT_MIN), ('deriva lunare', LUNAR_SHIFT_MIN)):
        hits = []; tested = 0
        for _, e in today.iterrows():
            m_today = e.t.hour * 60 + e.t.minute
            for _, o in old.iterrows():
                dd = (d0 - pd.Timestamp(o.day)).days
                if dd <= 0 or dd > 10: continue
                tested += 1
                m_old = (o.t.hour * 60 + o.t.minute + shift * dd) % 1440
                if abs((m_today - m_old + 720) % 1440 - 720) <= 3 and o.comp == e.comp:
                    hits.append(f'{hm(e.t)}{e.comp}↔{o.day[5:]} {hm(o.t)}')
        exp = tested * 7 / 1440 / 3     # prob. di cadere entro ±3 min, stessa componente (≈1/3)
        L.append(f'  {label}: {len(hits)} coincidenze (attese per caso ≈{exp:.1f})' + (': ' + ' '.join(hits[:12]) if hits else ''))
    # terzo orologio: stessa distanza dalla culminazione del radiante (radiante che si sposta)
    for col in [c for c in today.columns if c.startswith('culm_') and c in old.columns]:
        hits = []; tested = 0
        for _, e in today.iterrows():
            for _, o in old.iterrows():
                dd = (d0 - pd.Timestamp(o.day)).days
                if dd <= 0 or dd > 10 or pd.isna(o[col]) or pd.isna(e[col]): continue
                tested += 1
                if abs(e[col] - o[col]) <= 3 and o.comp == e.comp:
                    hits.append(f'{hm(e.t)}{e.comp}↔{o.day[5:]} {hm(o.t)} ({e[col]:+.0f} min)')
        exp = tested * 7 / 1440 / 3
        L.append(f'  culminazione {col[5:]}: {len(hits)} coincidenze (attese per caso ≈{exp:.1f})' + (': ' + ' '.join(hits[:10]) if hits else ''))

# ─── SINCRONISMI (primi controlli: il cielo) ──────────────────────────────────

SYNC_TOL_MIN = 10       # evento ↔ momento geometrico
REC_TOL_MIN  = 4        # stessa distanza dal momento geometrico in giorni diversi
REC_DAYS     = 30       # quanti giorni indietro cercare le ricorrenze

def _geo_moments(day):
    """sorgere, culminazione, passaggio inferiore, tramonto (UT) visti da SD per Sole, Luna, pianeti, nodo
    e radianti degli sciami attivi (con la loro deriva)"""
    import sky_d, sky_catalog as CAT, datetime as _dt
    d0 = _dt.datetime.fromisoformat(day)
    planets, sun = sky_d.get_planets(d0 + _dt.timedelta(hours=12))
    e = math.radians(23.44)
    def radec(l):
        l = math.radians(l)
        return math.degrees(math.atan2(math.sin(l) * math.cos(e), math.cos(l))) % 360, math.degrees(math.asin(math.sin(l) * math.sin(e)))
    bodies = {n.split(' ', 1)[1]: radec(p['lon']) for n, p in planets.items()}
    yf = d0.year + (d0.timetuple().tm_yday - 1) / 365.25
    for sh in CAT.active_showers(sun, yf):
        rd = CAT.radiant_radec(sh['name'], sun, yf)
        if rd: bodies['☄ ' + sh['name']] = rd
    lat, lon = COORD['SD']
    out = []
    for name, (ra, dec) in bodies.items():
        prev = None
        for m in range(0, 1441):
            t = pd.Timestamp(d0) + pd.Timedelta(minutes=m)
            H = (_gmst_deg(t) + lon - ra) % 360
            a = _alt(ra, dec, lat, lon, t)
            if prev:
                pH, pa = prev
                if pa <= 0 < a: out.append((name, 'sorge', t))
                if pa > 0 >= a: out.append((name, 'tramonta', t))
                if pH > 300 and H < 60: out.append((name, 'culmina', t))
                if pH < 180 <= H: out.append((name, 'passaggio inferiore', t))
            prev = (H, a)
    return out

def sec_sincronismi(folder, day, D, L, reg):
    L.append('\n══ 0. SINCRONISMI — eventi del giorno e momenti del cielo visti da SD (±%d min) ══' % SYNC_TOL_MIN)
    try:
        geo = _geo_moments(day)
    except Exception as ex:
        L.append(f'  calcolo non riuscito: {ex}'); return
    # eventi: sincroni + minimo e massimo di Y (mediana rete EU)
    eu = [s for s in EU if s in D]
    ev = [(e['t'], f"{e['comp']}({e['n']}) {e['st']}{e['val']:+.1f}") for e in SYNC_LOG]
    ymed = pd.concat([D[s].Y - D[s].Y.mean() for s in eu], axis=1).median(axis=1)
    ev += [(ymed.idxmin(), 'minimo Y'), (ymed.idxmax(), 'massimo Y')]
    rows = []
    for t, lab in sorted(ev):
        hits = [(n, k, g) for n, k, g in geo if abs((t - g).total_seconds()) <= SYNC_TOL_MIN * 60]
        txt = '; '.join(f"{n} {k} {g:%H:%M} ({(t - g).total_seconds() / 60:+.0f})" for n, k, g in hits)
        L.append(f'  {t:%H:%M} {lab}: ' + (txt if txt else '—'))
        for n, k, g in hits:
            rows.append({'day': day, 't': t.strftime('%H:%M'), 'evento': lab, 'corpo': n, 'momento': k,
                         'g': g.strftime('%H:%M'), 'scarto': round((t - g).total_seconds() / 60, 1)})
    L.append(f'  (per confronto: in media cadono {len(geo) * 2 * SYNC_TOL_MIN / 1440:.1f} momenti del cielo in una finestra di ±{SYNC_TOL_MIN} min)')
    # archivio e ricorrenze: stesso corpo, stesso momento, stesso scarto (±REC_TOL) in giorni diversi
    sp = folder.parent / 'SKY_SINCRONISMI.csv'
    old = pd.read_csv(sp) if sp.exists() else pd.DataFrame(columns=['day', 't', 'evento', 'corpo', 'momento', 'g', 'scarto'])
    old = old[old.day != day]
    allr = pd.concat([old, pd.DataFrame(rows)], ignore_index=True).sort_values(['day', 't'])
    allr.to_csv(sp, index=False)
    lim = (pd.Timestamp(day) - pd.Timedelta(days=REC_DAYS)).strftime('%Y-%m-%d')
    past = allr[(allr.day < day) & (allr.day >= lim)]
    L.append(f'  RICORRENZE (ultimi {REC_DAYS} giorni: stesso corpo e momento, stesso scarto ±{REC_TOL_MIN} min — l\'orario può slittare con il corpo):')
    found = 0
    for r in rows:
        m = past[(past.corpo == r['corpo']) & (past.momento == r['momento']) & ((past.scarto - r['scarto']).abs() <= REC_TOL_MIN)]
        if len(m):
            found += 1
            serie = ', '.join(f"{a[5:]} {b}" for a, b in zip(m.day, m.t)) + f", OGGI {r['t']}"
            L.append(f"  ★ {r['corpo']} {r['momento']} (scarto {r['scarto']:+.0f} min): {serie}")
    if not found: L.append('  nessuna')
    reg['sincronismi'] = len(rows); reg['ricorrenze_sinc'] = found

# ─── SOLE, NODI, PERCORSI (D2 + D7 nella stessa tabella) ──────────────────────

def _latest(folder, key):
    f = sorted(glob.glob(str(folder / f'SUN_{key}_*')))
    if not f: return None
    txt = open(f[-1], encoding='utf-8', errors='ignore').read()
    if f[-1].endswith('.json'):
        try: return json.loads(txt)
        except Exception: return None
    return txt

def _cls_rank(c):
    try: return {'A': 0, 'B': 1, 'C': 2, 'M': 3, 'X': 4}[c[0].upper()] + float(c[1:] or 0) / 10
    except Exception: return -1

def planets_on_nodes(day, tol=2.0):
    """pianeti (e Luna) entro tol° dal radiante eclittico di uno sciame attivo"""
    try:
        import sky_d, sky_catalog as CAT, datetime as _dt
        dt = _dt.datetime.fromisoformat(day) + _dt.timedelta(hours=12)
        planets, sun = sky_d.get_planets(dt)
        yf = dt.year + (dt.timetuple().tm_yday - 1) / 365.25
        out = []
        for sh in CAT.active_showers(sun, yf):
            for n, p in planets.items():
                if 'Sole' in n: continue
                d = abs((p['lon'] - sh['rad_lon_date'] + 180) % 360 - 180)
                if d <= tol: out.append(f"{n.split(' ', 1)[1]}@{sh['name']} {d:.1f}°")
        return out
    except Exception as e:
        return [f'(calcolo non riuscito: {e})']

def sec_sun(folder, day, L, reg):
    L.append('\n══ 16. SOLE, NODI E PERCORSI (D2 + D7) ══')
    nodes = planets_on_nodes(day)
    L.append('  Pianeti/Luna sul nodo di uno sciame attivo (≤2°): ' + ('; '.join(nodes) if nodes else 'nessuno'))
    fl = []
    for x in (_latest(folder, 'flares') or []):
        c = x.get('max_class') or x.get('current_class') or ''
        t = x.get('max_time') or x.get('begin_time')
        if c and t: fl.append((pd.Timestamp(t).tz_localize(None) if pd.Timestamp(t).tzinfo else pd.Timestamp(t), c))
    for x in (_latest(folder, 'dflr') or []):
        c, t = x.get('classType') or '', x.get('peakTime') or x.get('beginTime')
        if c and t:
            tt = pd.Timestamp(t.replace('Z', '')); 
            if not any(abs((tt - a).total_seconds()) < 900 and b[0] == c[0] for a, b in fl): fl.append((tt, c))
    cmes = []
    for c in (_latest(folder, 'dcme') or []):
        sp = [a.get('speed') for a in c.get('analyses', []) if a.get('speed')]
        if c.get('startTime'): cmes.append((pd.Timestamp(c['startTime'].replace('Z', '')), max(sp) if sp else None))
    hss = [pd.Timestamp(h['eventTime'].replace('Z', '')) for h in (_latest(folder, 'dhss') or []) if h.get('eventTime')]
    d0 = pd.Timestamp(day); d1 = d0 + pd.Timedelta(days=1)
    fday = [(t, c) for t, c in fl if d0 <= t < d1]
    mx = [(t, c) for t, c in fday if c[:1].upper() in 'MX']
    cday = [(t, v) for t, v in cmes if d0 <= t < d1]
    hday = [t for t in hss if d0 <= t < d1]
    newest = max([t for t, _ in fl] + [t for t, _ in cmes] + hss) if (fl or cmes or hss) else None
    if newest is None:
        L.append('  dati solari assenti (file SUN_* non presenti)')
    else:
        if newest < d0 - pd.Timedelta(days=2): L.append(f'  ⚠ dati solari non aggiornati (ultimo evento {newest:%d/%m})')
        L.append(f'  Brillamenti del giorno: {len(fday)} (M/X: {len(mx)})' +
                 (' — ' + ' '.join(f'{c}@{t:%H:%M}' for t, c in sorted(mx)) if mx else ''))
        if fday: L.append(f'  Classe massima: {max((c for _, c in fday), key=_cls_rank)}')
        L.append(f'  CME del giorno: {len(cday)}' + (' — ' + ' '.join(f'{t:%H:%M}' + (f'({v:.0f} km/s)' if v else '') for t, v in cday) if cday else ''))
        L.append(f'  Arrivo di flusso veloce (buco coronale) a Terra: ' + (' '.join(f'{t:%H:%M}' for t in hday) if hday else 'no'))
    idx = _latest(folder, 'indices') or ''
    ssn = f107 = None
    for line in idx.splitlines():
        parts = line.split()
        if len(parts) > 4 and parts[:3] == [f'{d0.year}', f'{d0.month:02d}', f'{d0.day:02d}']:
            try: f107, ssn = float(parts[3]), float(parts[4])
            except Exception: pass
    if f107 is not None: L.append(f'  Flusso radio F10.7: {f107:.0f} — macchie solari: {ssn:.0f}')
    # archivio giornaliero del Sole (per guardare indietro 1–4 giorni)
    sp = folder.parent / 'SKY_SOLE.csv'
    old = pd.read_csv(sp) if sp.exists() else pd.DataFrame(columns=['day'])
    old = old[old.day != day]
    rowd = {'day': day, 'MX': len(mx), 'cme': len(cday), 'hss': len(hday), 'nodi': len(nodes), 'nodi_elenco': ' | '.join(nodes)}
    allrows = pd.concat([old, pd.DataFrame([rowd])], ignore_index=True).sort_values('day')
    allrows.to_csv(sp, index=False)
    # percorsi
    prev = allrows[(allrows.day < day) & (allrows.day >= (d0 - pd.Timedelta(days=4)).strftime('%Y-%m-%d'))]
    chain = prev[((prev.MX > 0) | (prev.cme > 0)) & (prev.nodi > 0)]
    cnt = {'diretto': 0, 'attraverso il Sole': 0, 'standard': 0, 'non valutabile': 0}
    lines = []
    for e in SYNC_LOG:
        if e['orphan'] == -1: k = 'non valutabile'
        elif e['orphan'] == 1: k = 'diretto'
        else: k = 'attraverso il Sole' if len(chain) else 'standard'
        cnt[k] += 1; lines.append(f"{hm(e['t'])}{e['comp']}:{k}")
    L.append('  Percorsi degli eventi sincroni: ' + ', '.join(f'{k} {v}' for k, v in cnt.items()))
    if len(chain): L.append('  (nei 4 giorni prima: eruzioni in giorni con pianeta sul nodo → ' + ', '.join(chain.day.str[5:]) + ')')
    L.append('  "diretto" = senza causa a L1; "attraverso il Sole" = causa a L1 + eruzione 1–4 giorni prima in giorno con nodo; "standard" = causa a L1 senza quel legame')
    reg.update({'sun_MX': len(mx), 'sun_cme': len(cday), 'sun_hss': len(hday), 'sun_f107': f107, 'sun_ssn': ssn,
                'nodi': len(nodes), 'nodi_elenco': ' | '.join(nodes),
                'p_diretto': cnt['diretto'], 'p_sole': cnt['attraverso il Sole'], 'p_standard': cnt['standard']})

# ─────────────────────────────────────────────────────────────────────────────
def main():
    if len(sys.argv) < 2:
        sys.exit('Uso: python3 sky_scan.py <cartella SKY_AAAA-MM-GG>')
    folder = pathlib.Path(sys.argv[1]).expanduser()
    day, D = load_all(folder)
    L = [f'SKY SCAN — {day} — stazioni: {", ".join(D)}',
         'Report automatico: solo numeri e segnalazioni, nessuna attribuzione di causa.']
    reg = {}
    rp = folder / 'RIEPILOGO.txt'
    if rp.exists():
        for line in rp.read_text(encoding='utf-8').splitlines():
            l = line.strip()
            if 'Saturno Δ=' in l and 'NODO' not in l: reg['Saturno_nodo'] = l.split('Δ=')[1].split('°')[0]
            if l.startswith('☽ Luna') and 'illum.' in l: reg['Luna_illum'] = l.split('illum.')[1].split('%')[0].strip()
    sec_integrity(D, L, reg)
    sec_jumps(D, L, reg)
    Lnoaa = []
    arrivals = sec_noaa(folder, Lnoaa, reg, day)
    sec_sync(D, L, reg, arrivals)
    sec_diurnal(D, L, reg)
    sec_night_evening(D, L, reg, day)
    sec_trend(D, L, reg)
    L += Lnoaa
    sec_regional(D, L, reg)
    sec_direction(D, L)
    sec_radiants(folder, D, day, L, reg)
    reg.update({f'lag_{k}': v for k, v in REG_LAGS.items()})
    sec_radiant_station(D, day, L, reg)
    sec_orphans(L, reg)
    sec_recurrence(folder, day, L)
    Ls = []; sec_sincronismi(folder, day, D, Ls, reg); L[2:2] = Ls   # i sincronismi in cima al report
    sec_sun(folder, day, L, reg)
    sec_registry(folder, day, reg, L)
    out = folder / f'SCAN_{day}.txt'
    out.write_text('\n'.join(L) + '\n', encoding='utf-8')
    print('\n'.join(L))
    print(f'\n→ salvato: {out}')

if __name__ == '__main__':
    main()
