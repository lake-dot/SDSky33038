#!/usr/bin/env python3
"""
SKY — Monitor geomagnetico corridoio San Daniele (GitHub Actions)
Stazioni: WIC (AT), LON (HR), THY (HU), CLF (FR)
Modello dati: INTERMAGNET NRT ≈ istantaneo — si chiede SEMPRE fino ad adesso.
Una stazione il cui ultimo dato è più vecchio di FREEZE_AFTER_MIN è in FREEZE:
viene esclusa dal check live e non conta per MIN_STATIONS.
Indicatori: Bollinger Bands + Z-score rolling + Rate of Change
Filtro: alert solo se ≥ MIN_STATIONS stazioni vive simultanee, transizioni gestite
        con stato persistito su repo (.sky_event)
Intensità: 🔴🟠🟡🟢⚪ da z_max equivalente
Effemeridi: pyswisseph (opzionale) nel riassunto mattutino 09:00 UT
Notifiche ntfy SOLO su eventi: ⚡ anomalia, ✅ rientro, 🌅 riassunto, 🔴 errore workflow
"""

import requests
import numpy as np
import datetime
import bisect
import json
import os

try:
    import swisseph as swe
    SWE_OK = True
except ImportError:
    SWE_OK = False
    print("ATTENZIONE: pyswisseph non installato — effemeridi disattivate")

# ─── CONFIGURAZIONE ───────────────────────────────────────────────────────────

NTFY_TOPIC      = os.environ.get("NTFY_TOPIC")
MORNING_ALREADY = os.environ.get("MORNING_ALREADY", "") == "true"
FORCE_MORNING   = os.environ.get("FORCE_MORNING", "").lower() == "true"

BASE = "https://imag-data.bgs.ac.uk/GIN_V1/hapi"
OBSERVATORIES = {
    'wic': 'Conrad, Austria (220km)',
    'lon': 'Lonjsko Polje, Croatia (290km)',
    'thy': 'Tihany, Hungary (370km)',
    'clf': 'Chambon-la-Forêt, France (~950km)',
}
OBS_ID = {
    'wic': 'wic/best-avail/PT1M/xyzf',
    'lon': 'lon/best-avail/PT1M/xyzf',
    'thy': 'thy/best-avail/PT1M/xyzf',
    'clf': 'clf/best-avail/PT1M/xyzf',
}

CHECK_STEP_MIN    = 15     # frequenza del cron
FREEZE_AFTER_MIN  = 30     # ultimo dato più vecchio di questo → stazione in FREEZE
WIDEN_FALLBACKS_H = (3, 6, 12, 24)  # allargamento finestra (solo per stazioni ferme)
EVENT_FILE        = '.sky_event'

MIN_STATIONS  = 2
BB_WINDOW     = 120
BB_SIGMA      = 2.0
ZSCORE_WINDOW = 360
ZSCORE_THRESH = 2.5
ROC_WINDOW    = 30
ROC_THRESH    = 8.0
COMPONENTS    = ['X', 'Y', 'Z', 'F']

SD_LON            = 13.0
ASPECT_TIGHT      = 2.0
SHOWER_PEAK_THRESH   = 20.0
SHOWER_PLANET_THRESH = 8.0

GIORNI_ITA = ['lunedì','martedì','mercoledì','giovedì','venerdì','sabato','domenica']
MESI_ITA   = ['','gennaio','febbraio','marzo','aprile','maggio','giugno',
              'luglio','agosto','settembre','ottobre','novembre','dicembre']

HEADERS = {'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) '
                         'AppleWebKit/537.36 (KHTML, like Gecko) Chrome/126.0 Safari/537.36'}

# ─── ASTRONOMIA ───────────────────────────────────────────────────────────────

PLANET_NAMES = {
    swe.SUN: '☉ Sole', swe.MOON: '☽ Luna', swe.MERCURY: '☿ Mercurio',
    swe.VENUS: '♀ Venere', swe.MARS: '♂ Marte', swe.JUPITER: '♃ Giove',
    swe.SATURN: '♄ Saturno', swe.URANUS: '♅ Urano', swe.NEPTUNE: '♆ Nettuno',
    swe.TRUE_NODE: '☊ Nodo Nord',
} if SWE_OK else {}

SIGN_NAMES = ['Ariete','Toro','Gemelli','Cancro','Leone','Vergine',
              'Bilancia','Scorpione','Sagittario','Capricorno','Acquario','Pesci']

METEOR_SHOWERS = {
    'Perseidi':       {'lambda': 46.0,  'peak_sun': 140.0, 'active': (100, 160), 'peak_date': '12 ago'},
    'Kappa Cignidi':  {'lambda': 321.4, 'peak_sun': 145.0, 'active': (135, 155), 'peak_date': '17 ago'},
    'Tauridi Sud':    {'lambda': 50.0,  'peak_sun': 220.0, 'active': (170, 270), 'peak_date': '5 nov'},
    'Tauridi Nord':   {'lambda': 58.0,  'peak_sun': 230.0, 'active': (185, 285), 'peak_date': '12 nov'},
    'June Tauridi S': {'lambda': 271.0, 'peak_sun':  80.0, 'active': (60,  120), 'peak_date': '27 giu'},
    'June Tauridi N': {'lambda': 279.0, 'peak_sun':  85.0, 'active': (60,  120), 'peak_date': '1 lug'},
    'Delta Aquaridi': {'lambda': 333.0, 'peak_sun': 125.0, 'active': (95,  155), 'peak_date': '30 lug'},
    'Alpha Capric.':  {'lambda': 307.0, 'peak_sun': 127.0, 'active': (95,  165), 'peak_date': '2 ago'},
    'Leonidi inv.':   {'lambda': 152.0, 'peak_sun': 235.0, 'active': (100, 160), 'peak_date': '17 nov'},
    'Encke nodo':     {'lambda': 334.0, 'peak_sun':   0.0, 'active': (0,   360), 'peak_date': '—'},
}

def lon_delta(a, b):
    d = abs(a - b) % 360
    return d if d <= 180 else 360 - d

def get_planets(dt):
    if not SWE_OK:
        return {}, 0.0
    jd = swe.julday(dt.year, dt.month, dt.day,
                    dt.hour + dt.minute/60.0 + dt.second/3600.0)
    planets = {}
    for pid, name in PLANET_NAMES.items():
        lon, retro = None, False
        for flags in (swe.FLG_SWIEPH | swe.FLG_SPEED, swe.FLG_MOSEPH | swe.FLG_SPEED):
            try:
                result, _ = swe.calc_ut(jd, pid, flags)
                lon, retro = result[0], result[3] < 0
                break
            except Exception:
                continue
        if lon is not None:
            planets[name] = {'lon': round(lon, 2),
                             'sign': SIGN_NAMES[int(lon/30) % 12],
                             'deg': round(lon % 30, 2), 'retro': retro}
    sun_lon = planets.get('☉ Sole', {}).get('lon', 0.0)
    return planets, sun_lon

def active_showers_probable(sun_lon, planets):
    result = []
    for sname, shower in METEOR_SHOWERS.items():
        a_min, a_max = shower['active']
        is_active = (a_min <= sun_lon <= a_max) if a_max > a_min else \
                    (sun_lon >= a_min or sun_lon <= a_max)
        if not is_active:
            continue
        dist_peak = lon_delta(sun_lon, shower['peak_sun'])
        planet_hit = any(lon_delta(p['lon'], shower['lambda']) <= SHOWER_PLANET_THRESH
                         for p in planets.values())
        if dist_peak <= SHOWER_PEAK_THRESH or planet_hit:
            result.append((sname, shower, dist_peak))
    return result

def astro_lines(dt):
    if not SWE_OK:
        return ["(effemeridi non disponibili)"]
    planets, sun_lon = get_planets(dt)
    lines = [f"☉ Sole: {sun_lon:.2f}° {SIGN_NAMES[int(sun_lon/30) % 12]}", ""]
    retro = [n for n, p in planets.items()
             if p['retro'] and '☉' not in n and '☽' not in n]
    if retro:
        lines.append("℞ RETROGRADI: " + ", ".join(retro)); lines.append("")
    sd_hits = [f"  {n} Δ={lon_delta(p['lon'], SD_LON):.2f}°"
               for n, p in planets.items()
               if lon_delta(p['lon'], SD_LON) <= ASPECT_TIGHT]
    if sd_hits:
        lines.append(f"📍 NODO SD 13°E (Δ≤{ASPECT_TIGHT}°):")
        lines.extend(sd_hits); lines.append("")
    pl = [(n, p) for n, p in planets.items() if '☽' not in n]
    conj = []
    for i in range(len(pl)):
        for j in range(i + 1, len(pl)):
            d = lon_delta(pl[i][1]['lon'], pl[j][1]['lon'])
            if d <= ASPECT_TIGHT:
                conj.append(f"  {pl[i][0]} ☌ {pl[j][0]}  Δ={d:.2f}°")
            elif abs(d - 180) <= ASPECT_TIGHT:
                conj.append(f"  {pl[i][0]} ☍ {pl[j][0]}  Δ={abs(d-180):.2f}°")
    if conj:
        lines.append(f"⚡ ASPETTI (Δ≤{ASPECT_TIGHT}°):")
        lines.extend(conj); lines.append("")
    showers = active_showers_probable(sun_lon, planets)
    if showers:
        lines.append(f"☄️ SCIAMI PROBABILI (Sole {sun_lon:.1f}°):")
        for sname, shower, dist_peak in showers:
            lines.append(f"  {sname}  λ={shower['lambda']:.0f}°  "
                         f"Δpicco={dist_peak:.1f}°  picco={shower['peak_date']}")
            for pname, p in planets.items():
                d = lon_delta(p['lon'], shower['lambda'])
                if d <= SHOWER_PLANET_THRESH:
                    stars = ' ★★' if d <= 2 else ' ★' if d <= 5 else ''
                    lines.append(f"    → {pname}  Δ={d:.2f}°{stars}")
        lines.append("")
    return lines

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

def intensita_evento(alerts_by_station):
    max_z = 0.0
    for obs, alerts in alerts_by_station.items():
        for a in alerts:
            z = abs(a.get('z', 0.0))
            if a.get('type') == 'ROC':
                z = abs(a.get('change_nT', 0.0)) / 8.0 * 2.0
            max_z = max(max_z, z)
    if max_z >= 7.0:   return '🔴 ECCEZIONALE', max_z
    elif max_z >= 5.0: return '🟠 ESTREMA',     max_z
    elif max_z >= 3.5: return '🟡 FORTE',        max_z
    elif max_z >= 2.5: return '🟢 MEDIA',        max_z
    else:              return '⚪ BASSA',         max_z

# ─── STATO EVENTO (persistito su repo) ────────────────────────────────────────

def load_event_state():
    try:
        with open(EVENT_FILE) as f:
            s = json.load(f)
        if isinstance(s, dict) and isinstance(s.get('active'), bool):
            return s
    except Exception:
        pass
    return {'active': False, 'started': None}

def save_event_state(state):
    try:
        with open(EVENT_FILE, 'w') as f:
            json.dump(state, f)
    except Exception as e:
        print(f"  [stato] errore salvataggio: {e}")

def durata_str(start_iso):
    try:
        t0 = datetime.datetime.fromisoformat(start_iso)
        mins = int((now_utc() - t0).total_seconds() // 60)
        h, m = divmod(mins, 60)
        return f"{h}h {m:02d}m" if h else f"{m}m"
    except Exception:
        return "?"

# ─── FETCH INTERMAGNET ────────────────────────────────────────────────────────

def fetch_range(obs_code, start_dt, end_dt, widen=WIDEN_FALLBACKS_H):
    """Scarica dati tra start ed end. time.max è SEMPRE end_dt (≈ adesso):
    nessun limite preventivo — se la stazione pubblica, arrivano i dati freschi.
    'widen' allarga solo il passato (per stazioni ferme da ore) — mai il futuro."""
    last_err = ""
    for extra_h in (0,) + tuple(widen):
        s = start_dt - datetime.timedelta(hours=extra_h)
        url = (f"{BASE}/data?id={OBS_ID[obs_code]}"
               f"&time.min={s.strftime('%Y-%m-%dT%H:%M:%SZ')}"
               f"&time.max={end_dt.strftime('%Y-%m-%dT%H:%M:%SZ')}&format=json")
        try:
            r = requests.get(url, headers=HEADERS, timeout=30)
        except Exception as ex:
            last_err = str(ex)
            continue
        if r.status_code != 200:
            last_err = f"HTTP {r.status_code} (finestra +{extra_h}h)"
            if "1405" in r.text:      # 'time outside valid range' → allarga al passato
                continue
            continue
        records = r.json().get('data', [])
        if records:
            return parse_records(records)
        last_err = f"0 record (finestra +{extra_h}h)"
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

# ─── INDICATORI (logica identica allo script originale) ───────────────────────

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
    alerts = []
    if times is None or len(times) < 30:
        return alerts
    n = bisect.bisect_left(times, cut_dt) if cut_dt is not None else len(times)
    if n < 30:
        return alerts
    for comp in COMPONENTS:
        arr = data[comp][:n]
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
    t = now_utc().strftime('%H:%M UT')
    intensita, max_z = intensita_evento(stations)
    lines = [f"⚡ SKY {t} — {', '.join(stations)}",
             f"{intensita} (z_max={max_z:.1f}σ)", ""]
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

    stations_alert, lastz, fresh, frozen = {}, {}, [], []
    for obs in OBSERVATORIES:
        print(f"  {obs.upper()}...", end=' ')
        res = fetch_range(obs, start, now)
        if not res:
            print("→ nessun dato")
            continue
        times, data = res
        age_min = (now - times[-1]).total_seconds() / 60
        if age_min > FREEZE_AFTER_MIN:
            frozen.append(obs.upper())
            print(f"{len(times)} rec ma FREEZE — ultimo dato {age_min/60:.1f}h fa → escluso")
            continue
        fresh.append(obs.upper())
        z = data['Z'][~np.isnan(data['Z'])]
        if len(z):
            lastz[obs.upper()] = float(z[-1])
        alerts = evaluate_station(times, data)
        if alerts:
            stations_alert[obs.upper()] = alerts
            print(f"{len(times)} rec (lag {age_min:.0f}') → {len(alerts)} segnali")
        else:
            print(f"{len(times)} rec (lag {age_min:.0f}') → ok")

    alert_now = len(stations_alert) >= MIN_STATIONS
    state = load_event_state()
    print(f"  vive={fresh} frozen={frozen} | anomalia ora={sorted(stations_alert)} | "
          f"stato evento: {'ATTIVO' if state['active'] else 'no'}")

    if alert_now and not state['active']:
        ntfy_send(format_alert(stations_alert),
                  title="SKY ALERT anomalia in corso", priority="high")
        state = {'active': True, 'started': now.isoformat(),
                 'stations': sorted(stations_alert)}
        save_event_state(state)
    elif alert_now and state['active']:
        print("  evento già attivo — nessun invio")
    elif not alert_now and state['active']:
        if not fresh:
            print("  nessuna stazione viva — rientro non dichiarabile, stato invariato")
        else:
            ntfy_send(f"Anomalia rientrata dopo {durata_str(state.get('started'))} "
                      f"— situazione normale.",
                      title="SKY rientro", priority="low")
            state = {'active': False, 'started': None}
            save_event_state(state)
    else:
        print("  OK — " + "  ".join(f"{o} Z={v:.1f}" for o, v in lastz.items()))

# ─── RIASSUNTO MATTUTINO ──────────────────────────────────────────────────────

def run_morning():
    now = now_utc()
    day_end   = now.replace(hour=0, minute=0, second=0, microsecond=0)
    day_start = day_end - datetime.timedelta(days=1)
    fetch_start = day_start - datetime.timedelta(minutes=ZSCORE_WINDOW + 60)

    station_data = {}
    zvals = {o: [] for o in OBSERVATORIES}
    for obs in OBSERVATORIES:
        print(f"  {obs.upper()}...", end=' ')
        # dati STORICI (ieri): finestra esatta, nessun allargamento
        res = fetch_range(obs, fetch_start, day_end, widen=(0,))
        if not res:
            print("→ nessun dato")
            continue
        times, data = res
        print(f"{len(times)} rec")
        station_data[obs] = res
        for tt, v in zip(times, data['Z']):
            if day_start <= tt < day_end and not np.isnan(v):
                zvals[obs].append(float(v))

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

    msg = format_morning(now, day_start, events, zvals)
    print(msg)
    ntfy_send(msg, title="SKY riassunto mattutino")

def format_morning(now, day_start, events, zvals):
    lines = [f"🌅 SKY — {data_italiana(now)}", ""]
    lines.extend(astro_lines(now))

    if events:
        lines.append("─" * 24)
        counts = {'BB': 0, 'ZSCORE': 0, 'ROC': 0}
        for _, _, a in events:
            t = a['type']
            counts['BB' if t.startswith('BB') else t] += 1
        lines.append(f"ANOMALIE ({data_italiana(day_start)}) — "
                     f"Totale: {len(events)}  BB:{counts['BB']} Z:{counts['ZSCORE']} ROC:{counts['ROC']}")
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

    stat_lines = []
    for obs in OBSERVATORIES:
        v = zvals[obs]
        if v:
            stat_lines.append(f"📡 {obs.upper()} Z: mean={np.mean(v):.2f} "
                              f"min={np.min(v):.2f} max={np.max(v):.2f} nT (n={len(v)})")
    if stat_lines:
        lines.append("─" * 24)
        lines.extend(stat_lines)
    return '\n'.join(lines)

# ─── MAIN ─────────────────────────────────────────────────────────────────────

def main():
    now = now_utc()
    print("=" * 50)
    print(f"  SKY — {now:%Y-%m-%d %H:%M} UT")
    print("=" * 50)

    if (FORCE_MORNING or now.hour == 9) and not MORNING_ALREADY:
        print("→ Riassunto mattutino")
        run_morning()
        with open('.morning_flag', 'w') as f:
            f.write(now.isoformat())

    print("→ Check anomalia")
    run_check()

    # raise RuntimeError("test rosso")   # ← decommenta per testare la notifica di errore

if __name__ == "__main__":
    main()
