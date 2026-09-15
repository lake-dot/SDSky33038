#!/usr/bin/env python3
"""
SKY — Monitor geomagnetico corridoio San Daniele (STATELESS, per GitHub Actions)
Indicatori: Bollinger Bands + Z-score rolling + Rate of Change
Filtro: alert solo se ≥2 stazioni simultanee E il check precedente (15' fa) era pulito
Notifiche: ntfy → iPhone
Nessuno stato su disco: ogni run scarica i dati freschi da INTERMAGNET e ricalcola tutto.
"""

import requests
import numpy as np
import datetime
import os

# ─── CONFIGURAZIONE ───────────────────────────────────────────────────────────

NTFY_TOPIC    = os.environ.get("NTFY_TOPIC")
FORCE_MORNING = os.environ.get("FORCE_MORNING", "").lower() == "true"

OBSERVATORIES = {
    'wic': 'Conrad, Austria (220km)',
    'lon': 'Lonjsko Polje, Croatia (290km)',
    'thy': 'Tihany, Hungary (370km)',
}
OBS_ID = {
    'wic': 'wic/best-avail/PT1M/xyzf',
    'lon': 'lon/best-avail/PT1M/xyzf',
    'thy': 'thy/best-avail/PT1M/xyzf',
}

BASE          = "https://imag-data.bgs.ac.uk/GIN_V1/hapi"
STEP_MIN      = 15     # minuti tra un check e il successivo (= cron GitHub)
DATA_LAG_H    = 3      # lag di pubblicazione INTERMAGNET (come nel tuo script)

MIN_STATIONS  = 2
BB_WINDOW     = 120
BB_SIGMA      = 2.0
ZSCORE_WINDOW = 360
ZSCORE_THRESH = 2.5
ROC_WINDOW    = 30
ROC_THRESH    = 8.0

COMPONENTS = ['X', 'Y', 'Z', 'F']

GIORNI_ITA = ['lunedì','martedì','mercoledì','giovedì','venerdì','sabato','domenica']
MESI_ITA   = ['','gennaio','febbraio','marzo','aprile','maggio','giugno',
              'luglio','agosto','settembre','ottobre','novembre','dicembre']

# ─── NTFY ─────────────────────────────────────────────────────────────────────

def ntfy_send(message, title="SKY ⚡", silent=False, tags="zap"):
    if not NTFY_TOPIC:
        print("  [ntfy] ⚠️ topic non configurato — messaggio solo nei log:")
        print(message)
        return
    try:
        requests.post(
            f"https://ntfy.sh/{NTFY_TOPIC}",
            data=message.encode('utf-8'),
            headers={'Title': title,
                     'Priority': 'low' if silent else 'high',
                     'Tags': tags},
            timeout=15)
        print("  [ntfy] inviato ✓")
    except Exception as e:
        print(f"  [ntfy] errore: {e}")

# ─── UTILITY ──────────────────────────────────────────────────────────────────

def now_utc():
    return datetime.datetime.now(datetime.timezone.utc)

def data_italiana(dt):
    return f"{GIORNI_ITA[dt.weekday()]} {dt.day} {MESI_ITA[dt.month]} {dt.year}"

# ─── FETCH INTERMAGNET ────────────────────────────────────────────────────────

def fetch_range(obs_code, start_dt, end_dt):
    """Scarica i dati tra start ed end. Ritorna (times, dict di array) o None."""
    obs_id = OBS_ID.get(obs_code)
    if not obs_id:
        return None
    for id_fmt in [obs_id,
                   obs_id.replace('best-avail', 'adjusted'),
                   obs_id.replace('best-avail', 'reported')]:
        url = (f"{BASE}/data?id={id_fmt}"
               f"&time.min={start_dt.strftime('%Y-%m-%dT%H:%M:%SZ')}"
               f"&time.max={end_dt.strftime('%Y-%m-%dT%H:%M:%SZ')}"
               f"&format=json")
        try:
            r = requests.get(url, timeout=30)
            if r.status_code == 200:
                records = r.json().get('data', [])
                if records:
                    return parse_records(records)
        except Exception:
            continue
    return None

def parse_records(records):
    times, X, Y, Z, F = [], [], [], [], []
    FILL = 99999.0
    for rec in records:
        try:
            ts = rec[0].replace('Z', '+00:00')
            times.append(datetime.datetime.fromisoformat(ts))
            xyz = rec[1]
            x, y, z = float(xyz[0]), float(xyz[1]), float(xyz[2])
            f = float(rec[2])
            X.append(np.nan if abs(x) >= FILL else x)
            Y.append(np.nan if abs(y) >= FILL else y)
            Z.append(np.nan if abs(z) >= FILL else z)
            F.append(np.nan if abs(f) >= FILL else f)
        except Exception:
            continue
    return times, {'X': np.array(X), 'Y': np.array(Y),
                   'Z': np.array(Z), 'F': np.array(F)}

# ─── INDICATORI (identici al tuo script) ──────────────────────────────────────

def bollinger_signal(arr):
    if len(arr) < BB_WINDOW // 2 or np.isnan(arr[-1]):
        return None
    window = arr[-BB_WINDOW:] if len(arr) >= BB_WINDOW else arr
    mean, std = np.nanmean(window), np.nanstd(window)
    if std < 0.1 or np.isnan(std):
        return None
    current = arr[-1]
    z = (current - mean) / std
    if current > mean + BB_SIGMA * std:  signal = 'BB_SOPRA'
    elif current < mean - BB_SIGMA * std: signal = 'BB_SOTTO'
    else: return None
    return {'type': signal, 'current': round(float(current), 2), 'z': round(float(z), 2)}

def zscore_signal(arr):
    if len(arr) < 30 or np.isnan(arr[-1]):
        return None
    window = arr[-ZSCORE_WINDOW:] if len(arr) >= ZSCORE_WINDOW else arr
    mean, std = np.nanmean(window), np.nanstd(window)
    if std < 0.1 or np.isnan(std):
        return None
    z = (arr[-1] - mean) / std
    return round(float(z), 2) if abs(z) > ZSCORE_THRESH else None

def roc_signal(arr):
    if len(arr) < ROC_WINDOW or np.isnan(arr[-1]):
        return None
    change = arr[-1] - arr[-ROC_WINDOW]
    if np.isnan(change):
        return None
    return round(float(change), 2) if abs(change) > ROC_THRESH else None

def evaluate_station(times, data, cut_dt=None):
    """Indicatori sui dati fino a cut_dt (escluso). cut_dt=None → tutti i dati."""
    alerts = []
    if times is None or len(times) < 30:
        return alerts
    if cut_dt is not None:
        keep = [i for i, t in enumerate(times) if t < cut_dt]
        if len(keep) < 30:
            return alerts
        data = {c: data[c][keep] for c in COMPONENTS}
    for comp in COMPONENTS:
        arr = data[comp]
        arr = arr[~np.isnan(arr)]          # come il tuo buffer: solo valori validi
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

# ─── CHECK (run singolo) ──────────────────────────────────────────────────────

def format_alert(alerts_by_station):
    t = (now_utc() - datetime.timedelta(hours=DATA_LAG_H)).strftime('%H:%M UT')
    lines = [f"SKY {t} — {', '.join(alerts_by_station)}"]
    for obs, alerts in alerts_by_station.items():
        for a in alerts:
            c = a['component']
            if a['type'].startswith('BB'):
                d = '↑' if 'SOPRA' in a['type'] else '↓'
                lines.append(f"{obs} {c} {d} z={a['z']:.1f}σ")
            elif a['type'] == 'ZSCORE':
                lines.append(f"{obs} {c} z={a['z']:.1f}σ")
            elif a['type'] == 'ROC':
                sign = '+' if a['change_nT'] > 0 else ''
                lines.append(f"{obs} {c} ROC {sign}{a['change_nT']:.0f}nT/30'")
    return '\n'.join(lines)

def run_check():
    end   = now_utc() - datetime.timedelta(hours=DATA_LAG_H)
    start = end - datetime.timedelta(minutes=ZSCORE_WINDOW + 120)
    cut   = end - datetime.timedelta(minutes=STEP_MIN)

    stations_now, stations_prev, lastz = {}, {}, {}

    for obs in OBSERVATORIES:
        print(f"  {obs.upper()}...", end=' ')
        res = fetch_range(obs, start, end)
        if not res:
            print("ERRORE fetch")
            continue
        times, data = res
        print(f"{len(times)} record")
        z = data['Z'][~np.isnan(data['Z'])]
        if len(z):
            lastz[obs.upper()] = float(z[-1])
        a_now = evaluate_station(times, data)
        a_prev = evaluate_station(times, data, cut_dt=cut)
        if a_now:  stations_now[obs.upper()]  = a_now
        if a_prev: stations_prev[obs.upper()] = a_prev

    alert_now  = len(stations_now)  >= MIN_STATIONS
    alert_prev = len(stations_prev) >= MIN_STATIONS
    print(f"  stazioni in anomalia: ora={len(stations_now)}, {STEP_MIN}' fa={len(stations_prev)}")

    if alert_now and not alert_prev:
        ntfy_send(format_alert(stations_now), title="SKY ⚡ anomalia in corso")
    elif alert_now and alert_prev:
        print("  anomalia in corso (già notificata — nessun invio)")
    elif alert_prev and not alert_now:
        ntfy_send("Anomalia rientrata — situazione tornata normale.",
                  title="SKY ✅ rientro", silent=True, tags="white_check_mark")
    else:
        print("  OK — " + "  ".join(f"{o} Z={v:.1f}" for o, v in lastz.items()))

# ─── RIASSUNTO MATTUTINO (ricalcolato dai dati grezzi, zero memoria) ──────────

def run_morning():
    now = now_utc()
    day_end   = now.replace(hour=0, minute=0, second=0, microsecond=0)
    day_start = day_end - datetime.timedelta(days=1)
    fetch_start = day_start - datetime.timedelta(minutes=ZSCORE_WINDOW + 60)

    print("  scarico 24h + buffer caldo...")
    station_data, zvals = {}, {obs: [] for obs in OBSERVATORIES}
    for obs in OBSERVATORIES:
        res = fetch_range(obs, fetch_start, day_end)
        if res:
            station_data[obs] = res
            for tt, v in zip(res[0], res[1]['Z']):
                if day_start <= tt < day_end and not np.isnan(v):
                    zvals[obs].append(float(v))
        else:
            print(f"  {obs.upper()}: ERRORE fetch")

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
                    a['obs'] = obs
                    a['cut'] = cut.isoformat()
                    events.append(a)
        cut += datetime.timedelta(minutes=5)

    msg = format_morning(day_start, events, zvals)
    print(msg)
    ntfy_send(msg, title="🌅 SKY riassunto mattutino", tags="sunrise")

def format_morning(day_start, events, zvals):
    lines = [f"🌅 SKY — {data_italiana(day_start)}",
             f"Anomalie mezzanotte→mezzanotte UTC", ""]
    if not events:
        lines.append("Nessuna anomalia registrata nelle ultime 24h.")
    else:
        counts = {'BB': 0, 'ZSCORE': 0, 'ROC': 0}
        for a in events:
            t = a['type']
            counts['BB' if t.startswith('BB') else t] += 1
        lines.append(f"Totale: {len(events)}  BB:{counts['BB']} Z:{counts['ZSCORE']} ROC:{counts['ROC']}")
        lines.append("")
        by_cut = {}
        for a in events:
            by_cut.setdefault(a['cut'], {}).setdefault(a['obs'], []).append(a)
        shown = 0
        for cut in sorted(by_cut):
            if shown >= 60:
                lines.append(f"... e altri {len(events) - sum(1 for _ in range(shown))} eventi")
                break
            lines.append(f"{cut[11:16]} UT  [{', '.join(sorted(by_cut[cut]))}]")
            for obs in sorted(by_cut[cut]):
                for a in by_cut[cut][obs]:
                    c = a['component']
                    if a['type'].startswith('BB'):
                        d = '↑' if 'SOPRA' in a['type'] else '↓'
                        lines.append(f"  {obs} BB{d} {c} {a['current']:.1f}nT z={a['z']:.1f}σ")
                    elif a['type'] == 'ZSCORE':
                        lines.append(f"  {obs} Z {c} z={a['z']:.1f}σ")
                    else:
                        sign = '+' if a['change_nT'] > 0 else ''
                        lines.append(f"  {obs} ROC {c} {sign}{a['change_nT']:.1f}nT/30'")
                    shown += 1
    lines.append("")
    for obs in OBSERVATORIES:
        v = zvals[obs]
        if v:
            lines.append(f"📡 {obs.upper()} Z: mean={round(float(np.mean(v)),2)} "
                         f"min={round(float(np.min(v)),2)} max={round(float(np.max(v)),2)} nT (n={len(v)})")
    return '\n'.join(lines)

# ─── MAIN ─────────────────────────────────────────────────────────────────────

def main():
    now = now_utc()
    print("=" * 60)
    print(f"  SKY stateless — {now.strftime('%Y-%m-%d %H:%M UT')}")
    print("=" * 60)
    if not NTFY_TOPIC:
        print("  ⚠️  NTFY_TOPIC non configurato (aggiungi il Secret!)")

    if FORCE_MORNING or (now.hour == 9 and now.minute < STEP_MIN):
        print("→ Riassunto mattutino")
        try:
            run_morning()
        except Exception as e:
            print(f"  errore riassunto: {e}")

    print("→ Check anomalia")
    run_check()

if __name__ == "__main__":
    main()
