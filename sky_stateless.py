#!/usr/bin/env python3
"""
SKY — Monitor geomagnetico corridoio San Daniele (GitHub Actions)
Stazioni: WIC (AT), LON (HR), THY (HU), CLF (FR)

Modello dati (empirico): il feed INTERMAGNET/BGS pubblica a LOTTI irregolari,
tipicamente con ~3h di ritardo. Strategia:
  - fetch: si prova SEMPRE prima la finestra fino ad 'adesso'; se il server
    rifiuta ('time outside valid range'), la finestra retrocede a passi.
  - freeze: una stazione è STALE se il suo ultimo dato è >= STALE_GAP_MIN più
    vecchio di quello della stazione più aggiornata della rete (auto-calibrato).
    Le stale sono escluse dal check live.

Indicatori: Bollinger Bands + Z-score rolling + Rate of Change
Filtro: alert solo se >= MIN_STATIONS stazioni VIVE simultanee; transizioni
        gestite con stato persistito su repo (.sky_event)

Notifiche ntfy:
  ⚡ inizio anomalia            → subito, alta priorità
  🔄 evento in corso            → a cambio classe di intensità, o heartbeat ogni UPDATE_EVERY_MIN
                                 (durata, intensità attuale vs precedente, trend, picco)
  ✅ rientro in ciclo normale   → con durata totale e picco dell'evento
  🌅 riassunto mattutino 09:00 UT (con effemeridi se pyswisseph presente)
  🔴 errore workflow            → push (gestito nel workflow, step if: failure())
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

CHECK_STEP_MIN = 15                  # frequenza del cron
STALE_GAP_MIN  = 90                  # stazione STALE se >=90' indietro rispetto
                                     # alla stazione più aggiornata della rete
BACK_STEPS_H   = (0, 3, 6, 12, 24)   # retrocessione finestra fetch (mai il futuro)
EVENT_FILE     = '.sky_event'

UPDATE_EVERY_MIN = 60                # heartbeat durante evento attivo
TREND_DZ         = 0.3               # Δσ per dichiarare l'intensità in aumento/diminuzione

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

def age_str(ts, now):
    mins = max(0, int((now - ts).total_seconds() // 60))
    h, m = divmod(mins, 60)
    return f"{h}h {m:02d}m" if h else f"{m}m"

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
            json.dump(state, f, ensure_ascii=False)
    except Exception as e:
        print(f"  [stato] errore salvataggio: {e}")

def durata_str(start_iso):
    try:
        t0 = datetime.datetime.fromisoformat(start_iso)
        mins = max(0, int((now_utc() - t0).total_seconds() // 60))
        h, m = divmod(mins, 60)
        return f"{h}h {m:02d}m" if h else f"{m}m"
    except Exception:
        return "?"

# ─── FETCH INTERMAGNET ────────────────────────────────────────────────────────

def fetch_range(obs_code, start_dt, end_dt, back_steps=BACK_STEPS_H):
    """Prova PRIMA la finestra ideale fino ad 'adesso'; se il server rifiuta
    ('time outside valid range'), retrocede TUTTA la finestra a passi."""
    last_err = ""
    for off_h in back_steps:
        off = datetime.timedelta(hours=off_h)
        s = start_dt - off
        e = end_dt   - off
        url = (f"{BASE}/data?id={OBS_ID[obs_code]}"
               f"&time.min={s.strftime('%Y-%m-%dT%H:%M:%SZ')}"
               f"&time.max={e.strftime('%Y-%m-%dT%H:%M:%SZ')}&format=json")
        try:
            r = requests.get(url, headers=HEADERS, timeout=30)
        except Exception as ex:
            last_err = str(ex)
            continue
        if r.status_code != 200:
            last_err = f"HTTP {r.status_code} (end −{off_h}h)"
            if "1405" in r.text:
                continue
            continue
        records = r.json().get('data', [])
        if records:
            return parse_records(records)
        last_err = f"0 record (end −{off_h}h)"
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

# ─── FORMATTEZZA ALERT / UPDATE ───────────────────────────────────────────────

def _station_lines(stations):
    lines = []
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
    return lines

def format_alert(stations, epoch_str):
    t = now_utc().strftime('%H:%M UT')
    intensita, max_z = intensita_evento(stations)
    lines = [f"⚡ SKY {t} — {', '.join(stations)}",
             f"{intensita} (z_max={max_z:.1f}σ)",
             f"dati fino alle {epoch_str}", ""]
    lines.extend(_station_lines(stations))
    return '\n'.join(lines)

def format_update(stations, state, epoch_str):
    """Update evento in corso: durata, intensità attuale vs notifica precedente,
    trend, picco finora."""
    classe, max_z = intensita_evento(stations)
    prev_class = state.get('last_class') or classe
    prev_z     = float(state.get('last_zmax', max_z))
    if max_z > prev_z + TREND_DZ:    trend = '↗ in aumento'
    elif max_z < prev_z - TREND_DZ:  trend = '↘ in diminuzione'
    else:                            trend = '→ stabile'
    lines = [f"⚡ SKY — evento in corso da {durata_str(state.get('started'))}",
             f"{', '.join(stations)}",
             f"{classe} (z_max={max_z:.1f}σ) {trend}  [era {prev_class} {prev_z:.1f}σ]",
             f"picco finora: {state.get('peak_class','—')} ({float(state.get('peak_zmax',0)):.1f}σ)",
             f"dati fino alle {epoch_str}", ""]
    lines.extend(_station_lines(stations))
    return '\n'.join(lines)

# ─── CHECK ANOMALIA ───────────────────────────────────────────────────────────

def run_check():
    now   = now_utc()
    start = now - datetime.timedelta(minutes=ZSCORE_WINDOW + 120)

    results, frozen = {}, []
    for obs in OBSERVATORIES:
        print(f"  {obs.upper()}...", end=' ')
        res = fetch_range(obs, start, now)
        if not res or not res[0]:
            print("→ nessun dato disponibile")
            frozen.append(obs.upper())
            continue
        times, data = res
        results[obs] = res
        print(f"{len(times)} rec | ultimo dato {times[-1]:%H:%M} UT "
              f"({age_str(times[-1], now)} fa)")

    # ── il presente della rete = la stazione più aggiornata
    last_ts = {obs: times[-1] for obs, (times, data) in results.items()}
    if not last_ts:
        print("  NESSUNA stazione risponde — rete muta, stato evento invariato")
        return
    net_fresh = max(last_ts.values())
    epoch_str = net_fresh.strftime('%H:%M UT')
    live  = sorted(o for o, t in last_ts.items()
                   if (net_fresh - t).total_seconds() / 60 <= STALE_GAP_MIN)
    stale = sorted(o for o in last_ts if o not in live)
    print(f"  📡 dati disponibili fino alle {epoch_str} | "
          f"vive={sorted(o.upper() for o in live)} "
          f"stale={sorted(o.upper() for o in stale) + frozen}")

    # ── indicatori SOLO sulle stazioni vive
    stations_alert, lastz = {}, {}
    for obs in live:
        times, data = results[obs]
        z = data['Z'][~np.isnan(data['Z'])]
        if len(z):
            lastz[obs.upper()] = (float(z[-1]), times[-1])
        alerts = evaluate_station(times, data)
        if alerts:
            stations_alert[obs.upper()] = alerts

    alert_now = len(stations_alert) >= MIN_STATIONS
    state = load_event_state()
    now_iso = now.isoformat()

    if alert_now:
        classe, max_z = intensita_evento(stations_alert)

    if alert_now and not state['active']:
        # ── INIZIO evento
        ntfy_send(format_alert(stations_alert, epoch_str),
                  title="SKY ALERT anomalia in corso", priority="high")
        state = {'active': True, 'started': now_iso,
                 'last_notify': now_iso,
                 'last_class': classe, 'last_zmax': max_z,
                 'peak_class': classe, 'peak_zmax': max_z,
                 'stations': sorted(stations_alert), 'epoch': epoch_str}
        save_event_state(state)

    elif alert_now and state['active']:
        # ── evento continua: aggiorna il picco, decidi se mandare update
        if max_z > float(state.get('peak_zmax', 0)):
            state['peak_class'], state['peak_zmax'] = classe, max_z
        last_notify = state.get('last_notify')
        try:
            elapsed = (now - datetime.datetime.fromisoformat(last_notify)).total_seconds() / 60
        except Exception:
            elapsed = float('inf')
        class_changed = classe != state.get('last_class')
        if class_changed or elapsed >= UPDATE_EVERY_MIN:
            ntfy_send(format_update(stations_alert, state, epoch_str),
                      title="SKY evento in corso", priority="high")
            state['last_notify'] = now_iso
            state['last_class'], state['last_zmax'] = classe, max_z
            state['stations'] = sorted(stations_alert)
            state['epoch'] = epoch_str
            save_event_state(state)
        else:
            next_in = max(0, UPDATE_EVERY_MIN - elapsed)
            print(f"  evento in corso ({classe} z_max={max_z:.1f}σ) — "
                  f"update tra {next_in:.0f}' o a cambio intensità")

    elif not alert_now and state['active']:
        # ── RIENTRO in ciclo normale
        if not live:
            print("  nessuna stazione viva — rientro non dichiarabile, stato invariato")
        else:
            ntfy_send("✅ Rientro in ciclo normale "
                      f"dopo {durata_str(state.get('started'))}\n"
                      f"picco evento: {state.get('peak_class','—')} "
                      f"(z_max={float(state.get('peak_zmax',0)):.1f}σ)\n"
                      f"dati fino alle {epoch_str}",
                      title="SKY rientro", priority="low")
            state = {'active': False, 'started': None}
            save_event_state(state)

    else:
        print("  OK — " + "  ".join(f"{o} Z={v:.1f} ({t:%H:%M} UT)"
                                    for o, (v, t) in lastz.items()))

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
        res = fetch_range(obs, fetch_start, day_end, back_steps=(0,))
        if not res or not res[0]:
            print("→ nessun dato")
            continue
        times, data = res
        print(f"{len(times)} rec | ultimo {times[-1]:%H:%M} UT")
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
