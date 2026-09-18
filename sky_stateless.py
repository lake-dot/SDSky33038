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

Indicatori: Bollinger Bands + Z-score rolling + Rate of Change
Filtro: alert solo se >= MIN_STATIONS stazioni VIVE simultanee

Notifiche ntfy (macchina a stati su repo, .sky_event):
  ⚡ inizio anomalia                     → subito, alta priorità
  🔄 escalation                          → SOLO se classe intensità SUPERIORE all'ultima
                                          notificata, e al massimo uno ogni UPDATE_MIN_GAP_MIN
  ✅ rientro in ciclo normale            → dopo CLEAN_CHECKS_RIENTRO check puliti consecutivi
                                          (durata evento + picco raggiunto)
  🌅 effemeridi mattutine 09:00 UT       → astro a matematica pura (nessuna dipendenza)
  🔴 errore workflow                     → push (gestito nel workflow, step if: failure())
"""

import math
import requests
import numpy as np
import datetime
import bisect
import json
import os

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

UPDATE_MIN_GAP_MIN    = 180          # max un aggiornamento ogni 3h (e solo escalation)
CLEAN_CHECKS_RIENTRO  = 2            # check puliti consecutivi prima del rientro (2×15')

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

# ─── ASTRONOMIA: MATEMATICA PURA (dal tuo sky.py, nessuna dipendenza) ─────────

DEG = math.pi / 180.0
AU_LIGHT_DAYS = 0.0057755183

UNCERT = {
    'Sole': 0.02, 'Luna': 0.20, 'Mercurio': 0.30,
    'Venere': 0.30, 'Marte': 0.30, 'Giove': 0.30,
    'Saturno': 0.30, 'Urano': 0.30, 'Nettuno': 0.40,
    'Nodo Nord': 0.50,
}

SYMBOLS = {
    'Sole': '☉ Sole', 'Luna': '☽ Luna', 'Mercurio': '☿ Mercurio',
    'Venere': '♀ Venere', 'Marte': '♂ Marte', 'Giove': '♃ Giove',
    'Saturno': '♄ Saturno', 'Urano': '♅ Urano', 'Nettuno': '♆ Nettuno',
    'Nodo Nord': '☊ Nodo Nord',
}

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

def norm360(x):
    return x % 360.0

def wrap180(x):
    return ((x + 180.0) % 360.0) - 180.0

def precession_deg(jd):
    T = (jd - 2451545.0) / 36525.0
    return (5029.0966 * T + 1.11113 * T * T) / 3600.0

def julian_day(dt):
    y, m = dt.year, dt.month
    d = dt.day + (dt.hour + dt.minute / 60.0 + dt.second / 3600.0) / 24.0
    if m <= 2:
        y -= 1
        m += 12
    A = y // 100
    B = 2 - A + A // 4
    return int(365.25 * (y + 4716)) + int(30.6001 * (m + 1)) + d + B - 1524.5

def sun_lon(jd):
    n = jd - 2451545.0
    L = norm360(280.460 + 0.9856474 * n)
    g = (357.528 + 0.9856003 * n) * DEG
    return norm360(L + 1.915 * math.sin(g) + 0.020 * math.sin(2 * g))

def moon_lon(jd):
    T = (jd - 2451545.0) / 36525.0
    Lp = 218.3164477 + 481267.88123421 * T
    D = 297.8501921 + 445267.1114034 * T
    M = 357.5291092 + 35999.0502909 * T
    Mp = 134.9633964 + 477198.8675055 * T
    F = 93.2720950 + 483202.0175233 * T
    lon = (Lp
           + 6.288774 * math.sin(Mp * DEG)
           + 1.274027 * math.sin((2 * D - Mp) * DEG)
           + 0.658314 * math.sin(2 * D * DEG)
           + 0.213618 * math.sin(2 * Mp * DEG)
           - 0.185116 * math.sin(M * DEG)
           - 0.114332 * math.sin(2 * F * DEG)
           + 0.058793 * math.sin((2 * D - 2 * Mp) * DEG)
           + 0.057066 * math.sin((2 * D - M - Mp) * DEG)
           + 0.053322 * math.sin((2 * D + Mp) * DEG)
           + 0.045758 * math.sin((2 * D - M) * DEG)
           - 0.040923 * math.sin((M - Mp) * DEG)
           - 0.034720 * math.sin(D * DEG)
           - 0.030383 * math.sin((M + Mp) * DEG)
           + 0.015327 * math.sin((2 * D - 2 * F) * DEG)
           - 0.012528 * math.sin((Mp - 2 * F) * DEG))
    return norm360(lon)

def moon_node_lon(jd):
    T = (jd - 2451545.0) / 36525.0
    D = 297.8501921 + 445267.1114034 * T
    mean_node = 125.04452 - 1934.136261 * T + 0.0020708 * T * T
    return norm360(mean_node - 1.4979 * math.sin(2 * D * DEG))

PLANET_ELEMS = {
    'Mercurio': dict(a=0.38709927, adot=0.00000037, e=0.20563593, edot=0.00001906,
                     I=7.00497902, Idot=-0.00594749, L=252.25032350, Ldot=149472.67411175,
                     w=77.45779628, wdot=0.16047689, O=48.33076593, Odot=-0.12534081),
    'Venere':   dict(a=0.72333566, adot=0.00000390, e=0.00677672, edot=-0.00004107,
                     I=3.39467605, Idot=-0.00078890, L=181.97909950, Ldot=58517.81538729,
                     w=131.60246718, wdot=0.00268329, O=76.67984255, Odot=-0.27769418),
    'terra':    dict(a=1.00000261, adot=0.00000562, e=0.01671123, edot=-0.00004392,
                     I=-0.00001531, Idot=-0.01294668, L=100.46457166, Ldot=35999.37244981,
                     w=102.93768193, wdot=0.32327364, O=0.0, Odot=0.0),
    'Marte':    dict(a=1.52371034, adot=0.00001847, e=0.09339410, edot=0.00007882,
                     I=1.84969142, Idot=-0.00813131, L=-4.55343205, Ldot=19140.30268499,
                     w=-23.94362959, wdot=0.44441088, O=49.55953891, Odot=-0.29257343),
    'Giove':    dict(a=5.20288700, adot=-0.00011607, e=0.04838624, edot=-0.00013253,
                     I=1.30439695, Idot=-0.00183714, L=34.39644051, Ldot=3034.74612775,
                     w=14.72847983, wdot=0.21252668, O=100.47390909, Odot=0.20469106),
    'Saturno':  dict(a=9.53667594, adot=-0.00125060, e=0.05386179, edot=-0.00050991,
                     I=2.48599187, Idot=0.00193609, L=49.95424423, Ldot=1222.49362201,
                     w=92.59887831, wdot=-0.41897216, O=113.66242448, Odot=-0.28867794),
    'Urano':    dict(a=19.18916464, adot=-0.00196176, e=0.04725744, edot=-0.00004397,
                     I=0.77263783, Idot=-0.00242939, L=313.23810451, Ldot=428.48202785,
                     w=170.95427630, wdot=0.40805281, O=74.01692503, Odot=0.04240589),
    'Nettuno':  dict(a=30.06992276, adot=0.00026291, e=0.00859048, edot=0.00005105,
                     I=1.77004347, Idot=0.00035372, L=-55.12002969, Ldot=218.45945325,
                     w=44.96476227, wdot=-0.32241464, O=131.78422574, Odot=-0.00508664),
}

def kepler_solve(M_rad, e):
    E = M_rad + e * math.sin(M_rad)
    for _ in range(10):
        dE = (E - e * math.sin(E) - M_rad) / (1.0 - e * math.cos(E))
        E -= dE
        if abs(dE) < 1e-9:
            break
    return E

def helio_ecl(el, jd):
    T = (jd - 2451545.0) / 36525.0
    a = el['a'] + el['adot'] * T
    e = el['e'] + el['edot'] * T
    I = (el['I'] + el['Idot'] * T) * DEG
    L = el['L'] + el['Ldot'] * T
    w = el['w'] + el['wdot'] * T
    O = el['O'] + el['Odot'] * T
    M = norm360(L - w) * DEG
    om = (w - O) * DEG
    Om = O * DEG
    E = kepler_solve(M, e)
    xp = a * (math.cos(E) - e)
    yp = a * math.sqrt(1.0 - e * e) * math.sin(E)
    co, so = math.cos(om), math.sin(om)
    ci, si = math.cos(Om), math.sin(Om)
    cI, sI = math.cos(I), math.sin(I)
    x = (co * ci - so * si * cI) * xp + (-so * ci - co * si * cI) * yp
    y = (co * si + so * ci * cI) * xp + (-so * si + co * ci * cI) * yp
    z = (so * sI) * xp + (co * sI) * yp
    return x, y, z

def planet_geo_lon(name, jd, with_precession=True):
    el = PLANET_ELEMS[name]
    tx, ty, tz = helio_ecl(PLANET_ELEMS['terra'], jd)
    px, py, pz = helio_ecl(el, jd)
    gx, gy, gz = px - tx, py - ty, pz - tz
    dist = math.sqrt(gx * gx + gy * gy + gz * gz)
    for _ in range(2):
        jd2 = jd - dist * AU_LIGHT_DAYS
        px, py, pz = helio_ecl(el, jd2)
        gx, gy, gz = px - tx, py - ty, pz - tz
        dist = math.sqrt(gx * gx + gy * gy + gz * gz)
    lon = norm360(math.degrees(math.atan2(gy, gx)))
    if with_precession:
        lon = norm360(lon + precession_deg(jd))
    return lon

def body_lon(name, jd):
    if name == 'Sole':
        return sun_lon(jd)
    if name == 'Luna':
        return moon_lon(jd)
    if name == 'Nodo Nord':
        return moon_node_lon(jd)
    return planet_geo_lon(name, jd)

def get_planets(dt):
    jd = julian_day(dt)
    planets = {}
    bodies = ['Sole', 'Luna', 'Mercurio', 'Venere', 'Marte',
              'Giove', 'Saturno', 'Urano', 'Nettuno', 'Nodo Nord']
    for name in bodies:
        lon = body_lon(name, jd)
        unc = UNCERT.get(name, 0.5)
        speed = None
        if name != 'Nodo Nord':
            l1 = body_lon(name, jd - 1.5)
            l2 = body_lon(name, jd + 1.5)
            speed = wrap180(l2 - l1) / 3.0
        retro = (speed is not None and speed < -0.02)
        stationary = (speed is not None and abs(speed) <= 0.02)
        deg = lon % 30.0
        border = (deg < 2 * unc) or (deg > 30.0 - 2 * unc)
        planets[SYMBOLS[name]] = {'lon': round(lon, 2),
                                  'sign': SIGN_NAMES[int(lon // 30) % 12],
                                  'deg': round(deg, 2),
                                  'retro': retro,
                                  'stationary': stationary,
                                  'speed': speed,
                                  'border': border}
    return planets, planets['☉ Sole']['lon']

def flag_mark(planets, *names):
    for n in names:
        p = planets.get(n)
        if p and (p.get('border') or p.get('stationary')):
            return ' (+-)'
    return ''

def lon_delta(a, b):
    d = abs(a - b) % 360
    return d if d <= 180 else 360 - d

def active_showers_probable(sun_lon_val, planets):
    result = []
    for sname, shower in METEOR_SHOWERS.items():
        a_min, a_max = shower['active']
        is_active = (a_min <= sun_lon_val <= a_max) if a_max > a_min else \
                    (sun_lon_val >= a_min or sun_lon_val <= a_max)
        if not is_active:
            continue
        dist_peak = lon_delta(sun_lon_val, shower['peak_sun'])
        planet_hit = any(lon_delta(pl['lon'], shower['lambda']) <= SHOWER_PLANET_THRESH
                         for pl in planets.values())
        if dist_peak <= SHOWER_PEAK_THRESH or planet_hit:
            result.append((sname, shower, dist_peak))
    return result

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

def durata_str(start_iso):
    try:
        t0 = datetime.datetime.fromisoformat(start_iso)
        mins = max(0, int((now_utc() - t0).total_seconds() // 60))
        h, m = divmod(mins, 60)
        return f"{h}h {m:02d}m" if h else f"{m}m"
    except Exception:
        return "?"

INTENSITY_ORDER = {'⚪ BASSA': 0, '🟢 MEDIA': 1, '🟡 FORTE': 2,
                   '🟠 ESTREMA': 3, '🔴 ECCEZIONALE': 4}

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
    return {'active': False, 'started': None, 'clean_streak': 0}

def save_event_state(state):
    try:
        with open(EVENT_FILE, 'w') as f:
            json.dump(state, f, ensure_ascii=False)
    except Exception as e:
        print(f"  [stato] errore salvataggio: {e}")

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
    """Aggiornamento SOLO in escalation: intensità superiore all'ultima notificata."""
    classe, max_z = intensita_evento(stations)
    lines = [f"⚡ SKY — evento in corso da {durata_str(state.get('started'))}",
             f"{', '.join(stations)} — {classe} (z_max={max_z:.1f}σ) ↗ in aumento",
             f"era: {state.get('last_class')} {float(state.get('last_zmax', 0)):.1f}σ · "
             f"picco finora: {state.get('peak_class')} ({float(state.get('peak_zmax', 0)):.1f}σ)",
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
        state['clean_streak'] = 0

    if alert_now and not state['active']:
        # ── INIZIO evento
        ntfy_send(format_alert(stations_alert, epoch_str),
                  title="SKY ALERT anomalia in corso", priority="high")
        state = {'active': True, 'started': now_iso,
                 'last_notify': now_iso,
                 'last_class': classe, 'last_zmax': max_z,
                 'peak_class': classe, 'peak_zmax': max_z,
                 'clean_streak': 0,
                 'stations': sorted(stations_alert), 'epoch': epoch_str}
        save_event_state(state)

    elif alert_now and state['active']:
        # ── evento continua: aggiorna picco, valuta SOLO escalation
        if max_z > float(state.get('peak_zmax', 0)):
            state['peak_class'], state['peak_zmax'] = classe, max_z
        cur_order  = INTENSITY_ORDER.get(classe, 0)
        prev_order = INTENSITY_ORDER.get(state.get('last_class', '⚪ BASSA'), 0)
        escalated  = cur_order > prev_order
        try:
            gap_min = (now - datetime.datetime.fromisoformat(
                       state.get('last_notify'))).total_seconds() / 60
        except Exception:
            gap_min = float('inf')
        if escalated and gap_min >= UPDATE_MIN_GAP_MIN:
            ntfy_send(format_update(stations_alert, state, epoch_str),
                      title="SKY intensita in aumento", priority="high")
            state['last_notify'] = now_iso
            state['last_class'], state['last_zmax'] = classe, max_z
            state['stations'] = sorted(stations_alert)
            state['epoch'] = epoch_str
            save_event_state(state)
        else:
            if escalated:
                print(f"  escalation {classe} rilevata ma attendsi "
                      f"(min {UPDATE_MIN_GAP_MIN}' tra notifiche)")
            else:
                print(f"  evento in corso ({classe} z_max={max_z:.1f}σ) — "
                      f"nessun invio (solo escalation notificano)")

    elif not alert_now and state['active']:
        # ── RIENTRO: serve quiete CONFERMATA su CHECK puliti consecutivi
        if not live:
            print("  nessuna stazione viva — rientro non dichiarabile, stato invariato")
        else:
            state['clean_streak'] = int(state.get('clean_streak', 0)) + 1
            if state['clean_streak'] >= CLEAN_CHECKS_RIENTRO:
                ntfy_send("✅ Rientro in ciclo normale "
                          f"dopo {durata_str(state.get('started'))}\n"
                          f"picco evento: {state.get('peak_class', '—')} "
                          f"(z_max={float(state.get('peak_zmax', 0)):.1f}σ)\n"
                          f"dati fino alle {epoch_str}",
                          title="SKY rientro", priority="low")
                state = {'active': False, 'started': None, 'clean_streak': 0}
                save_event_state(state)
            else:
                save_event_state(state)
                print(f"  quiete {state['clean_streak']}/{CLEAN_CHECKS_RIENTRO} — "
                      f"attendo conferma prima del rientro")

    else:
        print("  OK — " + "  ".join(f"{o} Z={v:.1f} ({t:%H:%M} UT)"
                                    for o, (v, t) in lastz.items()))

# ─── EFFEMERIDI MATTUTINE (SOLO astro, matematica pura) ───────────────────────

def build_morning(dt):
    date_ita = data_italiana(dt)
    planets, sun_val = get_planets(dt)

    lines = [f"🌅 SKY — {date_ita}",
             f"☉ Sole: {sun_val:.2f}° {SIGN_NAMES[int(sun_val / 30) % 12]}",
             ""]

    retro = [name for name, p in planets.items()
             if p['retro'] and 'Sole' not in name and 'Luna' not in name]
    if retro:
        lines.append("℞ RETROGRADI: " + ", ".join(retro))
        lines.append("")

    sd_hits = []
    for name, p in planets.items():
        d = lon_delta(p['lon'], SD_LON)
        if d <= ASPECT_TIGHT:
            sd_hits.append(f"  {name} Δ={d:.2f}°{flag_mark(planets, name)}")
    if sd_hits:
        lines.append(f"📍 NODO SD 13°E (Δ≤{ASPECT_TIGHT}°):")
        lines.extend(sd_hits)
        lines.append("")

    planet_list = [(name, p) for name, p in planets.items() if 'Luna' not in name]
    conj_lines = []
    for i in range(len(planet_list)):
        for j in range(i + 1, len(planet_list)):
            n1, p1 = planet_list[i]
            n2, p2 = planet_list[j]
            d = lon_delta(p1['lon'], p2['lon'])
            if d <= ASPECT_TIGHT:
                approx = ' (≈)' if abs(d - ASPECT_TIGHT) < 0.5 else ''
                conj_lines.append(f"  {n1} ☌ {n2}  Δ={d:.2f}°{approx}{flag_mark(planets, n1, n2)}")
            elif abs(d - 180) <= ASPECT_TIGHT:
                approx = ' (≈)' if abs(abs(d - 180) - ASPECT_TIGHT) < 0.5 else ''
                conj_lines.append(f"  {n1} ☍ {n2}  Δ={abs(d - 180):.2f}°{approx}{flag_mark(planets, n1, n2)}")
    if conj_lines:
        lines.append(f"⚡ ASPETTI (Δ≤{ASPECT_TIGHT}°):")
        lines.extend(conj_lines)
        lines.append("")

    shower_lines = []
    for sname, shower, dist_peak in active_showers_probable(sun_val, planets):
        shower_lines.append(f"☄️ {sname}  λ={shower['lambda']:.0f}°  "
                            f"Δpicco={dist_peak:.1f}°  picco={shower['peak_date']}")
        for name, p in planets.items():
            d = lon_delta(p['lon'], shower['lambda'])
            if d <= ASPECT_TIGHT:
                shower_lines.append(f"    → {name} Δ={d:.2f}° **{flag_mark(planets, name)}")
    if shower_lines:
        lines.append("☄️ SCIAMI PROBABILI:")
        lines.extend(shower_lines)
        lines.append("")

    if any('(≈)' in l or '(+-)' in l for l in lines):
        lines.append("ℹ️ (≈) aspetto al bordo soglia · (+-) corpo vicino a confine segno o stazionario")
    return '\n'.join(lines)

def run_morning():
    now = now_utc()
    print("→ Effemeridi (matematica pura, nessun fetch)")
    msg = build_morning(now)
    print(msg)
    ntfy_send(msg, title="SKY effemeridi")

# ─── MAIN ─────────────────────────────────────────────────────────────────────

def main():
    now = now_utc()
    print("=" * 50)
    print(f"  SKY — {now:%Y-%m-%d %H:%M} UT")
    print("=" * 50)

    if (FORCE_MORNING or now.hour == 9) and not MORNING_ALREADY:
        print("→ Effemeridi mattutine")
        run_morning()
        with open('.morning_flag', 'w') as f:
            f.write(now.isoformat())

    print("→ Check anomalia")
    run_check()

    # raise RuntimeError("test rosso")   # ← decommenta per testare la notifica di errore

if __name__ == "__main__":
    main()
