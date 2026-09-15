#!/usr/bin/env python3
"""
SKY — Monitor geomagnetico corridoio San Daniele (stateless, GitHub Actions)
Indicatori: Bollinger Bands + Z-score rolling + Rate of Change
Filtro: alert solo se ≥ MIN_STATIONS stazioni simultanee E il check precedente era pulito
Notifiche ntfy SOLO su eventi: ⚡ inizio anomalia, ✅ rientro, 🌅 riassunto 09:00 UT
Nessuno stato su disco: ogni run scarica i dati freschi e ricalcola tutto.
"""

import requests
import numpy as np
import datetime
import os

# ─── CONFIGURAZIONE ───────────────────────────────────────────────────────────

NTFY_TOPIC      = os.environ.get("NTFY_TOPIC")
MORNING_ALREADY = os.environ.get("MORNING_ALREADY", "") == "true"
FORCE_MORNING   = os.environ.get("FORCE_MORNING", "").lower() == "true"

BASE = "https://imag-data.bgs.ac.uk/GIN_V1/hapi"
OBSERVATORIES = {
    'wic': 'Conrad, Austria',
    'lon': 'Lonjsko Polje, Croatia',
    'thy': 'Tihany, Hungary',
}
OBS_ID = {
    'wic': 'wic/best-avail/PT1M/xyzf',
    'lon': 'lon/best-avail/PT1M/xyzf',
    'thy': 'thy/best-avail/PT1M/xyzf',
}

CHECK_STEP_MIN = 15            # = frequenza del cron
LAG_FALLBACKS  = (3, 6, 12, 24)   # lag pubblicazione INTERMAGNET (ore)

MIN_STATIONS  = 2              # con WIC fermo: serve anomalia su LON E THY insieme
BB_WINDOW     = 120
BB_SIGMA      = 2.0
ZSCORE_WINDOW = 360
ZSCORE_THRESH = 2.5
ROC_WINDOW    = 30
ROC_THRESH    = 8.0
COMPONENTS    = ['X', 'Y', 'Z', 'F']

GIORNI_ITA = ['lunedì','martedì','mercoledì','giovedì','venerdì','sabato','domenica']
MESI_ITA   = ['','gennaio','febbraio','marzo','aprile','maggio','giugno',
              'luglio','agosto','settembre','ottobre','novembre','dicembre']

HEADERS = {'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) '
                         'AppleWebKit/537.36 (KHTML, like Gecko) Chrome/126.0 Safari/537.36'}

# ─── NTFY ─────────────────────────────────────────────────────────────────────

def ntfy_send(message, title="SKY", priority="default"):
    if not NTFY_TOPIC:
        print("  [ntfy] topic mancante — messaggio solo nei log:")
        print(message)
        return
    safe_title = title.encode('ascii', 'ignore').decode().strip() or "SKY"
    try:
        r = requests.post(
            f"https://ntfy.sh/{NTFY_TOPIC}",
            data=message.encode('utf-8'),
            headers={'Title': safe_title, 'Priority': priority, 'Tags': 'satellite'},
            timeout=15)
        if r.status_code == 200:
            print("  [ntfy] inviato ✓")
        else:
            print(f"  [ntfy] errore HTTP {r.status_code}: {r.text[:150]}")
    except Exception as e:
        print(f"  [ntfy] errore: {e}")

# ─── UTILITY ──────────────────────────────────────────────────────────────────

def now_utc():
    return datetime.datetime.now(datetime.timezone.utc)

def data_italiana(dt):
    return f"{GIORNI_ITA[dt.weekday()]} {dt.day} {MESI_ITA[dt.month]} {dt.year}"

# ─── FETCH INTERMAGNET ────────────────────────────────────────────────────────

def fetch_range(obs_code, start_dt, end_dt, lag_list=LAG_FALLBACKS):
    """Scarica [start, end] provando a spostare la finestra indietro se il server
    rifiuta ('time outside valid range' = stazione con lag maggiore o ferma)."""
    last_err = ""
    for lag_h in lag_list:
        s = start_dt - datetime.timedelta(hours=lag_h)
        e = end_dt   - datetime.timedelta(hours=lag_h)
        url = (f"{BASE}/data?id={OBS_ID[obs_code]}"
               f"&time.min={s.strftime('%Y-%m-%dT%H:%M:%SZ')}"
               f"&time.max={e.strftime('%Y-%m-%dT%H:%M:%SZ')}&format=json")
        try:
            r = requests.get(url, headers=HEADERS, timeout=30)
        except Exception as ex:
            last_err = str(ex)
            continue
        if r.status_code != 200:
            last_err = f"HTTP {r.status_code} (lag {lag_h}h)"
            if "1405" in r.text:      # fuori range → riprova più indietro
                continue
            continue
        records = r.json().get('data', [])
        if records:
            return parse_records(records)
        last_err = f"0 record (lag {lag_h}h)"
    print(last_err, end=' ')
    return None

def parse_records(records):
    times, X, Y, Z, F = [], [], [], [], []
    FILL = 99999.0
    for rec in records:
        try:
            times.append(datetime.datetime.fromisoformat(rec[0].replace('Z', '+00:00')))
            xyz = rec[1]
            x, y, z, f = float(xyz[0]), float(xyz[1]), float(xyz[2]), float(rec[2])
            X.append(np.nan if abs(x) >= FILL else x)
            Y.append(np.nan if abs(y) >= FILL else y)
            Z.append(np.nan if abs(z) >= FILL else z)
            F.append(np.nan if abs(f) >= FILL else f)
        except Exception:
            continue
    return times, {'X': np.array(X), 'Y': np.array(Y),
                   'Z': np.array(Z), 'F': np.array(F)}

# ─── INDICATORI (logica identica al tuo script originale) ─────────────────────

def bollinger_signal(arr):
    if len(arr) < BB_WINDOW // 2:
        return None
    window = arr[-BB_WINDOW:] if len(arr) >= BB_WINDOW else arr
    mean, std = np.nanmean(window), np.nanstd(window)
    if std < 0.1 or np.isnan(std):
        return None
    current = arr[-1]
    z = (current - mean) / std
    if current > mean + BB_SIGMA * std:   signal = 'BB_SOPRA'
    elif current < mean - BB_SIGMA * std: signal = 'BB_SOTTO'
    else: return None
    return {'type': signal, 'current': round(float(current), 2), 'z': round(float(z), 2)}

def zscore_signal(arr):
    if len(arr) < 30:
        return None
    window = arr[-ZSCORE_WINDOW:] if len(arr) >= ZSCORE_WINDOW else arr
    mean, std = np.nanmean(window), np.nanstd(window)
    if std < 0.1 or np.isnan(std):
        return None
    z = (arr[-1] - mean) / std
    return round(float(z), 2) if abs(z) > ZSCORE_THRESH else None

def roc_signal(arr):
    if len(arr) < ROC_WINDOW:
        return None
    change = arr[-1] - arr[-ROC_WINDOW]
    if np.isnan(change):
        return None
    return round(float(change), 2) if abs(change) > ROC_THRESH else None

def evaluate_station(times, data, cut_dt=None):
    """Indicatori sui dati fino a cut_dt (escluso). cut_dt=None → tutta la serie."""
    alerts = []
    if times is None or len(times) < 30:
        return alerts
    if cut_dt is not None:
        idx = [i for i, t in enumerate(times) if t < cut_dt]
    else:
        idx = list(range(len(times)))
    if len(idx) < 30:
        return alerts
    for comp in COMPONENTS:
        arr = data[comp][idx]
        arr = arr[~np.isnan(arr)]
        if len(arr) < 30:
            continue
        bb = bollinger_signal(arr)
        if bb:
            bb['component'] = comp
            alerts.append(bb)
        zs = zscore_signal(arr)
        if zs is not None:
            alerts.append({'type': 'ZSCORE', 'z': zs,
                           'current': round(float(arr[-1]), 2), 'component': comp})
        roc = roc_signal(arr)
        if roc is not None:
            alerts.append({'type': 'ROC', 'change_nT': roc,
                           'current': round(float(arr[-1]), 2), 'component': comp})
    return alerts

# ─── CHECK ANOMALIA ───────────────────────────────────────────────────────────

def format_alert(stations):
    lines = ["⚡ Anomalia geomagnetica rilevata", ""]
    for obs, alerts in stations.items():
        lines.append(f"{obs}:")
        for a in alerts:
            c = a['component']
            if a['type'].startswith('BB'):
                d = '↑' if 'SOPRA' in a['type'] else '↓'
                lines.append(f"  BB{d} {c} {a['current']:.1f}nT z={a['z']:.1f}σ")
            elif a['type'] == 'ZSCORE':
                lines.append(f"  Z {c} z={a['z']:.1f}σ")
            else:
                s = '+' if a['change_nT'] > 0 else ''
                lines.append(f"  ROC {c} {s}{a['change_nT']:.1f}nT/30'")
    return '\n'.join(lines)

def run_check():
    now   = now_utc()
    start = now - datetime.timedelta(minutes=ZSCORE_WINDOW + 120)

    stations_now, stations_prev, lastz = {}, {}, {}
    for obs in OBSERVATORIES:
        print(f"  {obs.upper()}...", end=' ')
        res = fetch_range(obs, start, now)
        if not res:
            print("→ nessun dato")
            continue
        times, data = res
        lag_min = (now - times[-1]).total_seconds() / 60
        print(f"{len(times)} rec, ultimo {times[-1]:%H:%M} UT (lag {lag_min:.0f}')")
        z = data['Z'][~np.isnan(data['Z'])]
        if len(z):
            lastz[obs.upper()] = float(z[-1])
        a_now  = evaluate_station(times, data)
        a_prev = evaluate_station(times, data,
                                  cut_dt=times[-1] - datetime.timedelta(minutes=CHECK_STEP_MIN))
        if a_now:  stations_now[obs.upper()]  = a_now
        if a_prev: stations_prev[obs.upper()] = a_prev

    alert_now  = len(stations_now)  >= MIN_STATIONS
    alert_prev = len(stations_prev) >= MIN_STATIONS
    print(f"  in anomalia: ora={sorted(stations_now)} | {CHECK_STEP_MIN}' fa={sorted(stations_prev)}")

    if alert_now and not alert_prev:
        ntfy_send(format_alert(stations_now), title="SKY ALERT anomalia in corso", priority="high")
    elif alert_now and alert_prev:
        print("  anomalia in corso (già notificata)")
    elif alert_prev and not alert_now:
        ntfy_send("Anomalia rientrata — situazione tornata normale.", title="SKY rientro", priority="low")
    else:
        print("  OK — " + "  ".join(f"{o} Z={v:.1f}" for o, v in lastz.items()))

# ─── RIASSUNTO MATTUTINO (ricalcolato dai dati grezzi delle 24h) ──────────────

def run_morning():
    now = now_utc()
    day_end   = now.replace(hour=0, minute=0, second=0, microsecond=0)
    day_start = day_end - datetime.timedelta(days=1)
    fetch_start = day_start - datetime.timedelta(minutes=ZSCORE_WINDOW + 60)

    station_data = {}
    zvals = {o: [] for o in OBSERVATORIES}
    for obs in OBSERVATORIES:
        print(f"  {obs.upper()}...", end=' ')
        res = fetch_range(obs, fetch_start, day_end, lag_list=(0,))
        if not res:
            print("→ nessun dato")
            continue
        times, data = res
        print(f"{len(times)} rec")
        station_data[obs] = res
        for tt, v in zip(times, data['Z']):
            if day_start <= tt < day_end and not np.isnan(v):
                zvals[obs].append(float(v))

    # scansione ogni 5' con la STESSA logica del check live
    events = []
    cut = day_start
    while cut < day_end:
        by_station = {}
        for obs, (times, data) in station_data.items():
            a = evaluate_station(times, data, cut_dt=cut)
            if a:
                by_station[obs.upper()] = a
        if len(by_station) >= MIN_STATIONS:
            for obs, alerts in by_station.items():
                for a in alerts:
                    events.append((cut, obs, a))
        cut += datetime.timedelta(minutes=5)

    msg = format_morning(day_start, events, zvals)
    print(msg)
    ntfy_send(msg, title="SKY riassunto mattutino")

def format_morning(day_start, events, zvals):
    lines = [f"🌅 SKY — {data_italiana(day_start)}",
             "Anomalie mezzanotte → mezzanotte UTC", ""]
    if not events:
        lines.append("Nessuna anomalia nelle ultime 24h.")
    else:
        counts = {'BB': 0, 'ZSCORE': 0, 'ROC': 0}
        for _, _, a in events:
            t = a['type']
            counts['BB' if t.startswith('BB') else t] += 1
        lines.append(f"Totale: {len(events)}  BB:{counts['BB']} Z:{counts['ZSCORE']} ROC:{counts['ROC']}")
        lines.append("")
        by_cut = {}
        for cut, obs, a in events:
            by_cut.setdefault(cut, []).append((obs, a))
        shown = 0
        for cut in sorted(by_cut):
            if shown >= 50:
                lines.append(f"... (+{len(events) - shown} eventi)")
                break
            items = by_cut[cut]
            stations = sorted({o for o, _ in items})
            lines.append(f"{cut:%H:%M} UT  [{', '.join(stations)}]")
            for obs, a in items:
                c = a['component']
                if a['type'].startswith('BB'):
                    d = '↑' if 'SOPRA' in a['type'] else '↓'
                    lines.append(f"  {obs} BB{d} {c} {a['current']:.1f}nT z={a['z']:.1f}σ")
                elif a['type'] == 'ZSCORE':
                    lines.append(f"  {obs} Z {c} z={a['z']:.1f}σ")
                else:
                    s = '+' if a['change_nT'] > 0 else ''
                    lines.append(f"  {obs} ROC {c} {s}{a['change_nT']:.1f}nT/30'")
                shown += 1
    lines.append("")
    for obs in OBSERVATORIES:
        v = zvals[obs]
        if v:
            lines.append(f"📡 {obs.upper()} Z: mean={np.mean(v):.2f} "
                         f"min={np.min(v):.2f} max={np.max(v):.2f} nT (n={len(v)})")
    return '\n'.join(lines)

# ─── MAIN ─────────────────────────────────────────────────────────────────────

def main():
    now = now_utc()
    print("=" * 50)
    print(f"  SKY — {now:%Y-%m-%d %H:%M} UT")
    print("=" * 50)

    # riassunto: alle 09:00 UT, oppure forzato; UNA volta al giorno (flag su repo)
    if (FORCE_MORNING or now.hour == 9) and not MORNING_ALREADY:
        print("→ Riassunto mattutino")
        try:
            run_morning()
            with open('.morning_flag', 'w') as f:
                f.write(now.isoformat())
        except Exception as e:
            print(f"  errore riassunto: {e}")

    print("→ Check anomalia")
    try:
        run_check()
    except Exception as e:
        print(f"  errore check: {e}")

if __name__ == "__main__":
    main()
