#!/usr/bin/env python3
"""
SKY-DUMP — Estrazione dati geomagnetici di una giornata intera + effemeridi
(ossatura di sky.py: stessi OBS_ID, stesso fetch GIN, stesso parsing)
(Astronomia: matematica pura, identica alla versione GitHub Actions — nessuna dipendenza)

Output — cartella ~/Desktop/SKY_YYYY-MM-GG/:
  • un TXT per stazione: risoluzione 1 minuto, X Y Z F — archivio completo
  • RIEPILOGO.TXT: UNICO file di sintesi — effemeridi complete della giornata
    + statistiche stazione per stazione
  • INCROCIO_<step>min.TXT: tutte le stazioni allineate sulla stessa timeline,
    scarti dalla mediana del giorno — il file pensato per essere letto
    INTEGRALMENTE da un'AI in un messaggio

Uso:
    python3 sky_d.py                  → ieri, passo 10 min
    python3 sky_d.py 2026-09-15       → giorno specifico
    python3 sky_d.py 2026-09-15 5     → giorno specifico, passo 5 min

Sola lettura: nessun alert, nessun log, nessuna notifica, non tocca sky.py.
"""

import sys
import os
import time
import math
import requests
import numpy as np
import datetime
import pathlib

# ─── CONFIGURAZIONE ───────────────────────────────────────────────────────────

OBSERVATORIES = {
    # — Corridoio primario —
    'wic': 'Conrad, Austria (220km)',
    'lon': 'Lonjsko Polje, Croatia (290km)',
    'thy': 'Tihany, Hungary (370km)',
    # — 2ª linea: ovest —
    'clf': 'Chambon-la-Forêt, France (~830km)',
    'bfo': 'Black Forest, Germania (~430km)',
    # — 2ª linea: nord (meridiano ~13°E, allineata con SD) —
    'ngk': 'Niemegk, Germania (~660km)',
    # — 2ª linea: est —
    'bel': 'Belsk, Polonia (~850km)',
    'izn': 'Iznik, Turchia (~1500km)',
    # — 2ª linea: sud —
    'dur': 'Duronia, Italia (~520km)',
    # — Aurorali (test sottotempeste) —
    'abk': 'Abisko, Svezia (aurorale ~68°N)',
    'tro': 'Tromsø, Norvegia (aurorale ~70°N)',
    # — Extra-europee (test globale/regionale) —
    'kak': 'Kakioka, Giappone',
    'frd': 'Fredericksburg, USA',
    'ott': 'Ottawa, Canada',
    'brd': 'Brandon, Canada',
}

OBS_ID = {
    'wic': 'wic/best-avail/PT1M/xyzf',
    'lon': 'lon/best-avail/PT1M/xyzf',
    'thy': 'thy/best-avail/PT1M/xyzf',
    'clf': 'clf/best-avail/PT1M/xyzf',
    'bfo': 'bfo/best-avail/PT1M/xyzf',
    'ngk': 'ngk/best-avail/PT1M/xyzf',
    'bel': 'bel/best-avail/PT1M/xyzf',
    'izn': 'izn/best-avail/PT1M/xyzf',
    'dur': 'dur/best-avail/PT1M/xyzf',
    'abk': 'abk/best-avail/PT1M/xyzf',
    'tro': 'tro/best-avail/PT1M/xyzf',
    'kak': 'kak/best-avail/PT1M/xyzf',
    'frd': 'frd/best-avail/PT1M/xyzf',
    'ott': 'ott/best-avail/PT1M/xyzf',
    'brd': 'brd/best-avail/PT1M/xyzf',
}

# DUR pubblica in un sistema ruotato: declinazione IGRF 2026 ≈ 4.36°, dai dati ≈ 0.09°
# → nel file unico si aggiungono anche DURg_X / DURg_Y ruotate nel sistema geografico.
DUR_ROT_DEG = 4.27

# File unico con tutte le stazioni al minuto (CSV, colonne STAZ_X, STAZ_Y, STAZ_Z, STAZ_F)
WRITE_ALL_CSV       = True
WRITE_STATION_FILES = os.environ.get('SKY_STATION_FILES', '1') == '1'   # TXT per stazione (0 = solo file unico)

# Vento solare NOAA (ultime 24h, sorgenti SOLAR1/IMAP/ACE) — salvato solo se si estrae ieri/oggi
NOAA_RTSW = {
    'wind': 'https://services.swpc.noaa.gov/json/rtsw/rtsw_wind_1m.json',
    'mag':  'https://services.swpc.noaa.gov/json/rtsw/rtsw_mag_1m.json',
    'eph':  'https://services.swpc.noaa.gov/json/rtsw/rtsw_ephemerides_1h.json',
}

# Fallback su codice alternativo se quello ufficiale non risponde.
# VUOTO di proposito: un alias silenzioso può sostituire una stazione
# con un'altra a centinaia di km di distanza.
ALIASES = {}

BASE        = "https://imag-data.bgs.ac.uk/GIN_V1/hapi"
CHUNK_HOURS = 6        # giornata scaricata a blocchi di 6h (4 richieste/stazione)
PAUSE_S     = 0.4      # pausa cortese tra le richieste
FILL        = 99999.0  # valore fill INTERMAGNET → campo vuoto nel txt

USE_LOCAL_DAY = False  # False = giornata UT (00:00–23:59 UTC, come i dati grezzi)
LOCAL_UTC_OFFSET = datetime.timedelta(hours=2)   # usata solo se USE_LOCAL_DAY = True

import os
OUT_ROOT = pathlib.Path(os.environ.get('SKY_OUT', pathlib.Path.home() / "Desktop"))   # Cartella output (Desktop del Mac)

CROSS_STEP_DEFAULT = 10          # minuti tra le righe dell'INCROCIO
CROSS_COMPS = ['X', 'Y', 'Z']    # F resta nei file per stazione

# Effemeridi: matematica pura (nessuna dipendenza esterna)
INCLUDE_EPHEMERIDES  = True   # sezione effemeridi nel RIEPILOGO.txt
EPHEMERIDES_REF_HOUR = 12     # ora UT di riferimento del giorno per le posizioni

SD_LON            = 13.0
ASPECT_TIGHT      = 2.0
SHOWER_PEAK_THRESH   = 20.0
SHOWER_PLANET_THRESH = 8.0

# ─── ASTRONOMIA: MATEMATICA PURA (identica al monitor GitHub Actions) ─────────

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

GIORNI_ITA = ['lunedì','martedì','mercoledì','giovedì','venerdì','sabato','domenica']
MESI_ITA   = ['','gennaio','febbraio','marzo','aprile','maggio','giugno',
              'luglio','agosto','settembre','ottobre','novembre','dicembre']

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
    'Alpha Aurigidi': {'lambda': 90.8,  'peak_sun': 158.6, 'active': (148, 168), 'peak_date': '1 set'},
    'Eps Perseidi':   {'lambda': 56.6,  'peak_sun': 166.7, 'active': (160, 176), 'peak_date': '9 set'},
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

def data_italiana(dt):
    return f"{GIORNI_ITA[dt.weekday()]} {dt.day} {MESI_ITA[dt.month]} {dt.year}"

def fase_lunare(d_moon):
    """d_moon = wrap180(Luna − Sole): >0 crescente, <0 calante."""
    e = abs(d_moon)
    crescente = d_moon > 0
    if e <= 22.5:  return 'novilunio'
    if e <= 67.5:  return 'falce crescente' if crescente else 'falce calante'
    if e <= 112.5: return 'primo quarto'    if crescente else 'ultimo quarto'
    if e <= 157.5: return 'gibbosa crescente' if crescente else 'gibbosa calante'
    return 'plenilunio'

# ─── EFFEMERIDI DELLA GIORNATA (sezione del RIEPILOGO) ────────────────────────

def build_ephemeridi(dt):
    """Sezione effemeridi: posizioni, retrogradi, nodo SD, aspetti, sciami."""
    date_ita = data_italiana(dt)
    planets, sun_val = get_planets(dt)
    moon = planets['☽ Luna']
    d_moon = wrap180(moon['lon'] - sun_val)
    elong = abs(d_moon)
    illum = (1.0 - math.cos(elong * DEG)) / 2.0 * 100.0

    lines = ["─" * 70,
             "EFFEMERIDI DELLA GIORNATA",
             f"{date_ita} — riferimento {dt:%Y-%m-%d %H:%M} UT",
             "Posizioni geocentriche in longitudine eclittica — matematica pura (±0.3°)",
             "─" * 70, "",
             f"☉ Sole: {sun_val:.2f}° {SIGN_NAMES[int(sun_val / 30) % 12]}",
             f"☽ Luna: {moon['lon']:.2f}° — {fase_lunare(d_moon)} "
             f"(elong. {elong:.1f}°, illum. {illum:.0f}%)",
             ""]

    lines.append("🪐 POSIZIONI:")
    for name, p in planets.items():
        r = ''
        if p['speed'] is not None:
            if p['retro']:
                r = ' ℞'
            elif p['stationary']:
                r = ' ℞?'
        b = ' !' if p['border'] else ''
        lines.append(f"  {name:14s} {p['lon']:7.2f}°  {p['sign']:12s} "
                     f"{p['deg']:5.2f}°{r}{b}")
    lines.append("")

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
    try:
        import sky_catalog as CAT   # catalogo IMO 2026, λ☉ convertita alla data (vedi nota nel file)
        yfrac = dt.year + (dt.timetuple().tm_yday - 1) / 365.25
        for sh in CAT.active_showers(sun_val, yfrac):
            if sh['dist_peak'] is None:
                shower_lines.append(f"☄️ {sh['name']}  radiante mobile λ={sh['rad_lon_date']:.0f}°  V={sh['v']} km/s  (attivo fino al 20 set)")
            else:
                shower_lines.append(f"☄️ {sh['name']}  λ={sh['rad_lon_date']:.0f}°  AR={sh['ra']:.0f}° Dec={sh['dec']:+.0f}°  "
                                    f"V={sh['v']} km/s  Δpicco={sh['dist_peak']:.1f}°  picco λ☉={sh['peak_date_lon']:.1f}°")
            for name, p in planets.items():
                d = lon_delta(p['lon'], sh['rad_lon_date'])
                if d <= ASPECT_TIGHT:
                    shower_lines.append(f"    → {name} Δ={d:.2f}°{flag_mark(planets, name)}")
        title = "☄️ SCIAMI ATTIVI (catalogo IMO 2026, λ☉ alla data):"
    except ImportError:
        for sname, shower, dist_peak in active_showers_probable(sun_val, planets):
            shower_lines.append(f"☄️ {sname}  λ={shower['lambda']:.0f}°  "
                                f"Δpicco={dist_peak:.1f}°  picco={shower['peak_date']}")
        title = "☄️ SCIAMI PROBABILI (catalogo interno — sky_catalog.py non trovato):"
    if shower_lines:
        lines.append(title)
        lines.extend(shower_lines)
        lines.append("")

    lines.append("ℹ️ ℞ retrogrado · ℞? stazionario · ! bordo segno · "
                 "(≈) al bordo soglia · (+-) incertezza elevata")
    return lines

# ─── UTILITY ──────────────────────────────────────────────────────────────────

def now_utc():
    return datetime.datetime.now(datetime.UTC)

# ─── INTERMAGNET (identico a sky.py, timeout maggiorato: payload più grandi) ──

def fetch_data(obs_code, start_dt, end_dt):
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
                    return records
        except Exception:
            continue
    return None

def parse_records(records):
    times, X, Y, Z, F = [], [], [], [], []
    FILL = 99999.0
    for rec in records:
        try:
            times.append(rec[0])
            xyz = rec[1]
            x, y, z = float(xyz[0]), float(xyz[1]), float(xyz[2])
            f = float(rec[2])
            X.append(np.nan if abs(x) >= FILL else x)
            Y.append(np.nan if abs(y) >= FILL else y)
            Z.append(np.nan if abs(z) >= FILL else z)
            F.append(np.nan if abs(f) >= FILL else f)
        except Exception:
            continue
    return times, np.array(X), np.array(Y), np.array(Z), np.array(F)

def fetch_day(obs_code, start, end):
    """Scarica l'intero intervallo a blocchi di CHUNK_HOURS, deduplicando i timestamp."""
    merged = {}
    cur = start
    while cur < end:
        chunk_end = min(cur + datetime.timedelta(hours=CHUNK_HOURS), end)
        records = fetch_data(obs_code, cur, chunk_end)
        if records:
            for rec in records:
                try:
                    merged[rec[0]] = rec      # dedup per timestamp
                except Exception:
                    pass
        if chunk_end >= end:
            break
        cur = chunk_end
        time.sleep(PAUSE_S)
    return [merged[ts] for ts in sorted(merged.keys())]

# ─── GIORNATA DA ESTRARRE ─────────────────────────────────────────────────────

def day_bounds(date_arg=None):
    """(start, end, label) — da argomento da riga di comando oppure 'ieri'."""
    if date_arg:
        try:
            d = datetime.datetime.strptime(date_arg, '%Y-%m-%d').date()
        except ValueError:
            sys.exit(f"Data non valida: '{date_arg}' — usa AAAA-MM-GG (es. 2026-09-15)")
    else:
        if USE_LOCAL_DAY:
            local = now_utc() + LOCAL_UTC_OFFSET
            d = (local - datetime.timedelta(days=1)).date()
        else:
            d = (now_utc() - datetime.timedelta(days=1)).date()

    if USE_LOCAL_DAY:
        start = datetime.datetime(d.year, d.month, d.day) - LOCAL_UTC_OFFSET
    else:
        start = datetime.datetime(d.year, d.month, d.day)

    end = start + datetime.timedelta(days=1) - datetime.timedelta(seconds=1)
    return start, end, d.isoformat()

# ─── FORMATTAZIONE / STATISTICHE ──────────────────────────────────────────────

def fmt_val(v):
    """nT con 2 decimali; fill/NaN → campo vuoto (import pulito in fogli di calcolo)."""
    try:
        if v is None or np.isnan(v):
            return ''
    except TypeError:
        return ''
    return f"{v:.2f}"

def z_stats(Z):
    Zv = Z[~np.isnan(Z)]
    if len(Zv) == 0:
        return None
    return float(np.mean(Zv)), float(np.min(Zv)), float(np.max(Zv))

def comp_stats(arr):
    v = arr[~np.isnan(arr)]
    if len(v) == 0:
        return None
    return float(np.min(v)), float(np.max(v))

# ─── OUTPUT ───────────────────────────────────────────────────────────────────

def write_station_file(path, used, label, day_label, times, X, Y, Z, F, n, expected):
    """File per singola stazione: 1 minuto, X Y Z F. Prima riga compatibile
    con sky_stats.py (header 'COD — nome')."""
    L = []
    L.append("═" * 70)
    L.append(f"{used.upper()} — {label}")
    L.append(f"Giornata: {day_label} — INTERMAGNET GIN, cadenza 1 minuto, componenti X Y Z F (nT)")
    L.append(f"Record: {n}/{expected}  (blocchi da {CHUNK_HOURS}h, id {used}/best-avail)")
    st = z_stats(Z)
    if st:
        L.append(f"Z (nT):  mean={st[0]:.1f}  min={st[1]:.1f}  max={st[2]:.1f}")
    L.append("Valori mancanti (fill 99999) → campo vuoto")
    L.append("═" * 70)
    L.append(f"{'TIME (UT)':<27s}{'X':>12s}{'Y':>12s}{'Z':>12s}{'F':>12s}")
    for ts, x, y, z, f in zip(times, X, Y, Z, F):
        L.append(f"{ts:<27s}{fmt_val(x):>12s}{fmt_val(y):>12s}"
                 f"{fmt_val(z):>12s}{fmt_val(f):>12s}")
    path.write_text('\n'.join(L) + '\n', encoding='utf-8')

def write_station_nodata(path, code, name, day_label):
    L = ["═" * 70,
         f"{code} — {name}",
         f"Giornata: {day_label}",
         "NESSUN DATO DISPONIBILE per questa giornata",
         "═" * 70]
    path.write_text('\n'.join(L) + '\n', encoding='utf-8')

def build_cross(stations, start, end, step, day_label):
    """Matrice unica: righe = tempo (passo 'step'), colonne = stazione_componente.
    Valori = scarto dalla mediana del giorno (nT). Le stazioni diventano
    direttamente confrontabili (baseline rimosse) e il file resta abbastanza
    piccolo da essere letto intero da un'AI."""
    cols, series, medians = [], {}, {}
    for st in stations:
        if not st.get('used'):
            continue
        code = st['used'].upper()
        for comp in CROSS_COMPS:
            m, vals = {}, []
            for ts, v in zip(st['times'], st[comp]):
                if np.isnan(v):
                    continue
                m[ts[:16]] = v          # chiave 'YYYY-MM-DDTHH:MM' (robusta a Z/secondi)
                vals.append(v)
            if not vals:
                continue
            cols.append((code, comp))
            series[(code, comp)] = m
            medians[(code, comp)] = float(np.median(vals))

    L = [f"SKY INCROCIO — {day_label} — passo {step} min — orario UT",
         "Valori: scarto dalla mediana del giorno, per stazione/componente (nT, 1 decimale).",
         "Cella vuota = dato mancante in quel minuto. F e risoluzione 1 min: nei file stazione.",
         "time," + ",".join(f"{c}_{comp}" for c, comp in cols)]

    t = start
    while t <= end:
        key = t.strftime('%Y-%m-%dT%H:%M')
        row = [t.strftime('%H:%M')]
        for c, comp in cols:
            v = series[(c, comp)].get(key)
            row.append('' if v is None else f"{v - medians[(c, comp)]:+.1f}")
        L.append(','.join(row))
        t += datetime.timedelta(minutes=step)
    return L

def build_riepilogo(stations, day_label, expected, tz_note, eph_lines=None):
    """Documento unico: effemeridi complete + statistiche stazioni."""
    L = ["═" * 70,
         "SKY — RIEPILOGO GIORNATA",
         f"Giornata: {day_label} ({tz_note}) — estrazione {now_utc():%Y-%m-%d %H:%M} UT",
         f"Fonte: INTERMAGNET GIN — cadenza 1 minuto — attesi {expected} record/stazione",
         "═" * 70]

    # Sezione effemeridi (completa)
    if eph_lines:
        L.append("")
        L.extend(eph_lines)
        L.append("")

    # Sezione dati stazioni
    L.append("─" * 70)
    L.append("DATI STAZIONI")
    L.append("─" * 70)
    for st in stations:
        code = (st.get('used') or st['obs']).upper()
        L.append("")
        L.append(f"{code} — {st['label']}")
        if not st.get('used'):
            L.append("  NESSUN DATO DISPONIBILE per questa giornata")
            continue
        L.append(f"  record {st['n']}/{expected} ({st['n']/expected*100:.0f}%)"
                 f" — id {st['used']}/best-avail")
        for comp in ['X', 'Y', 'Z', 'F']:
            cs = comp_stats(st[comp])
            if cs is None:
                L.append(f"  {comp}: non pubblicata")
            else:
                extra = ''
                if comp == 'Z':
                    zs = z_stats(st[comp])
                    extra = f"  mean={zs[0]:.1f}"
                L.append(f"  {comp}: {cs[0]:.1f} … {cs[1]:.1f}  (Δ{cs[1]-cs[0]:.1f}){extra}")
    L.append("")
    L.append("In questa cartella: un TXT per stazione (1 min) + INCROCIO (timeline comune)")
    return L

# ─── FILE UNICO + NOAA ────────────────────────────────────────────────────────

def write_all_csv(path, stations, start, expected):
    import math as _m
    idx = {}
    cols = []
    for st in stations:
        if not st.get('used'):
            continue
        code = st['used'].upper()
        cols.append(code)
        idx[code] = {t[:16]: (st['X'][i], st['Y'][i], st['Z'][i], st['F'][i])
                     for i, t in enumerate(st['times'])}
    has_dur = 'DUR' in idx
    head = ['time_UT'] + [f'{c}_{k}' for c in cols for k in 'XYZF']
    if has_dur:
        head += ['DURg_X', 'DURg_Y']
    th = _m.radians(DUR_ROT_DEG)
    lines = ['# SKY file unico — nT, 1 min, UT — vuoto = dato mancante',
             f'# DURg_X/DURg_Y = DUR ruotata di {DUR_ROT_DEG}° nel sistema geografico',
             ','.join(head)]
    for k in range(expected):
        t = (start + datetime.timedelta(minutes=k)).strftime('%Y-%m-%dT%H:%M')
        row = [t]
        for c in cols:
            v = idx[c].get(t)
            row += [fmt_val(x) for x in v] if v else ['', '', '', '']
        if has_dur:
            v = idx['DUR'].get(t)
            if v and not (np.isnan(v[0]) or np.isnan(v[1])):
                xg = v[0]*_m.cos(th) - v[1]*_m.sin(th)
                yg = v[0]*_m.sin(th) + v[1]*_m.cos(th)
                row += [f'{xg:.2f}', f'{yg:.2f}']
            else:
                row += ['', '']
        lines.append(','.join(row))
    path.write_text('\n'.join(lines) + '\n', encoding='utf-8')

def save_noaa(out_dir):
    stamp = now_utc().strftime('%Y%m%dT%H%MZ')
    for key, url in NOAA_RTSW.items():
        try:
            r = requests.get(url, timeout=30, headers={'User-Agent': 'Mozilla/5.0'})
            if r.status_code == 200:
                (out_dir / f'NOAA_{key}_{stamp}.json').write_bytes(r.content)
                print(f"  NOAA {key}: salvato")
            else:
                print(f"  NOAA {key}: HTTP {r.status_code}")
        except Exception as e:
            print(f"  NOAA {key}: errore ({e})")

# ─── MAIN ─────────────────────────────────────────────────────────────────────

def main():
    args = sys.argv[1:]
    date_arg = args[0] if len(args) >= 1 else None
    try:
        step = int(args[1]) if len(args) >= 2 else CROSS_STEP_DEFAULT
    except ValueError:
        sys.exit(f"Passo non valido: '{args[1]}' — usa i minuti (es. 5 o 10)")
    step = max(1, step)

    start, end, day_label = day_bounds(date_arg)
    expected = int((end - start).total_seconds() // 60) + 1   # ~1440

    out_root = OUT_ROOT if OUT_ROOT.exists() else pathlib.Path.home()
    out_dir = out_root / f"SKY_{day_label}"
    out_dir.mkdir(parents=True, exist_ok=True)

    tz_note = 'ora locale CEST' if USE_LOCAL_DAY else 'UTC'
    print("=" * 60)
    print(f"  SKY-DUMP — Giornata {day_label} ({tz_note})")
    print(f"  Intervallo: {start:%Y-%m-%d %H:%M} → {end:%Y-%m-%d %H:%M} UT")
    print(f"  Stazioni: {', '.join(s.upper() for s in OBSERVATORIES)}")
    print(f"  Attesi {expected} record/minuto per stazione")
    print(f"  INCROCIO: passo {step} min")
    print("=" * 60)

    stations, incomplete = [], []

    for obs, name in OBSERVATORIES.items():
        print(f"  {obs.upper()}...", end=' ', flush=True)

        records = fetch_day(obs, start, end)
        used, label = obs, name
        if not records and obs in ALIASES:
            alt_code, alt_name = ALIASES[obs]
            records = fetch_day(alt_code, start, end)
            if records:
                used, label = alt_code, alt_name

        if not records:
            print("NESSUN DATO")
            if WRITE_STATION_FILES:
                write_station_nodata(out_dir / f"{obs.upper()}.txt", obs.upper(), name, day_label)
            stations.append({'obs': obs, 'used': None, 'label': name, 'n': 0})
            continue

        times, X, Y, Z, F = parse_records(records)
        n = len(times)
        print(f"{n}/{expected} record ({n/expected*100:.0f}%) — servita come {used.upper()}")
        if n < expected:
            incomplete.append(f"{obs.upper()} ({n}/{expected})")

        if WRITE_STATION_FILES:
            write_station_file(out_dir / f"{used.upper()}.txt",
                               used, label, day_label, times, X, Y, Z, F, n, expected)
        stations.append({'obs': obs, 'used': used, 'label': label, 'n': n,
                         'times': times, 'X': X, 'Y': Y, 'Z': Z, 'F': F})

    # ── Effemeridi della giornata (matematica pura, calcolo istantaneo) ──
    eph_lines = None
    if INCLUDE_EPHEMERIDES:
        eph_dt = start + datetime.timedelta(hours=EPHEMERIDES_REF_HOUR)
        eph_lines = build_ephemeridi(eph_dt)

    cross_name = f"INCROCIO_{step}min.txt"
    (out_dir / cross_name).write_text(
        '\n'.join(build_cross(stations, start, end, step, day_label)) + '\n',
        encoding='utf-8')
    (out_dir / "RIEPILOGO.txt").write_text(
        '\n'.join(build_riepilogo(stations, day_label, expected, tz_note, eph_lines)) + '\n',
        encoding='utf-8')

    if WRITE_ALL_CSV:
        write_all_csv(out_dir / f"TUTTE_{day_label}.csv", stations, start, expected)
        print(f"  TUTTE_{day_label}.csv — file unico con tutte le stazioni")
    if (now_utc().date() - datetime.date.fromisoformat(day_label)).days <= 1:
        save_noaa(out_dir)
    else:
        print("  NOAA: saltato (i file coprono solo le ultime 24h)")

    files = ', '.join((s.get('used') or s['obs']).upper() + '.txt' for s in stations)
    print("\n" + "=" * 60)
    print(f"  Cartella: {out_dir}")
    print(f"  Per stazione (1 min): {files}")
    print(f"  {cross_name}  ← da incollare INTERO all'AI per gli incroci")
    print(f"  RIEPILOGO.txt — documento unico: effemeridi complete + dati stazioni")
    if incomplete:
        print(f"  ⚠️ Dati incompleti su: {', '.join(incomplete)}")
        print("    (il GIN pubblica con ~3h di ritardo: rilancia più tardi per completare)")
    print(f"  Apri la cartella:  open \"{out_dir}\"")
    print("=" * 60)

if __name__ == "__main__":
    main()
