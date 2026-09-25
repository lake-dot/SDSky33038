"""
sky_catalog.py — catalogo sciami meteorici condiviso da sky_d.py e sky_scan.py

NOTA DI METODO (da non dimenticare)
───────────────────────────────────
• Orologio comune = LONGITUDINE DEL SOLE (λ☉): la posizione della Terra sull'orbita,
  indipendente dal calendario — come facevano gli antichi.
• Gli script calcolano λ☉ ALLA DATA (equinozio dell'anno in corso, longitudine tropicale).
• IMO pubblica i λ☉ riferiti all'EQUINOZIO 2000 (J2000).
• Qui i valori sono scritti in J2000 (come nella fonte) e convertiti ALLA DATA aggiungendo
  la precessione (~50.29"/anno → ~0.37° nel 2026, ≈ 9 ore di spostamento sui picchi).
  Così RIEPILOGO, SCAN e registro usano un unico sistema.
• Fonte: lista di lavoro IMO 2026 (via elenco Wikipedia "List of meteor showers", 2026).
• Rami di Encke separati (Tauridi Sud / Nord): niente più voce generica "Encke nodo".
• Sorgente dell'Antielio (ANT, include le ex Piscidi del sud, picco storico ~20 set):
  radiante mobile, calcolato come punto antisolare sull'eclittica (approssimazione).
"""
import math, datetime

PRECESSION_ARCSEC_YR = 50.29
OBLIQ_J2000 = 23.4393

def precession_to_date(year_frac):
    return (year_frac - 2000.0) * PRECESSION_ARCSEC_YR / 3600.0

def _sun_lon_of_date(d):
    jd = (d - datetime.datetime(2000, 1, 1, 12)).total_seconds() / 86400 + 2451545.0
    n = jd - 2451545.0
    L = (280.460 + 0.9856474 * n) % 360
    g = math.radians(357.528 + 0.9856003 * n)
    return (L + 1.915 * math.sin(g) + 0.020 * math.sin(2 * g)) % 360

def _date_to_lam_j2000(dm):
    """'17 Jul' (date IMO 2026) → λ☉ J2000"""
    d = datetime.datetime.strptime(dm + ' 2026', '%d %b %Y')
    return (_sun_lon_of_date(d) - precession_to_date(2026.5)) % 360

def _ecl_lon(ra_deg, dec_deg):
    e = math.radians(OBLIQ_J2000); a = math.radians(ra_deg); d = math.radians(dec_deg)
    return math.degrees(math.atan2(math.sin(a) * math.cos(e) + math.tan(d) * math.sin(e), math.cos(a))) % 360

# nome, attivo dal, al, λ picco J2000, AR (ore), Dec, V (km/s), ZHR, corpo progenitore
_RAW = [
    ('Quadrantidi',          '28 Dec', '12 Jan', 283.15, 15.3, +49, 41, '80', '2003 EH1'),
    ('Gamma Ursae Minoridi', '10 Jan', '22 Jan', 298.0,  15.2, +67, 31, '3',  '?'),
    ('Alpha Centauridi',     '31 Jan', '20 Feb', 319.4,  14.1, -58, 58, '6',  '?'),
    ('Liridi',               '14 Apr', '30 Apr', 32.32,  18.1, +34, 49, '18', 'C/1861 G1 Thatcher'),
    ('Pi Puppidi',           '15 Apr', '28 Apr', 33.5,   7.3,  -45, 18, 'var', '26P/Grigg-Skjellerup'),
    ('Eta Aquaridi',         '19 Apr', '28 May', 45.5,   22.5, -1,  66, '50', '1P/Halley'),
    ('Eta Liridi',           '3 May',  '14 May', 50.0,   19.4, +43, 43, '3',  'C/1983 H1'),
    ('Arietidi diurne',      '14 May', '24 Jun', 76.7,   2.9,  +24, 38, '30', '1566 Icarus?'),
    ('Bootidi di giugno',    '22 Jun', '2 Jul',  90.3,   14.7, +48, 18, 'var', '7P/Pons-Winnecke'),
    ('Pegasidi di luglio',   '1 Jul',  '20 Jul', 108.0,  23.1, +11, 63, '3',  'C/1979 Y1'),
    ('Gamma Draconidi lug.', '25 Jul', '31 Jul', 125.13, 18.7, +51, 27, '5',  '?'),
    ('Delta Aquaridi Sud',   '12 Jul', '23 Aug', 128.0,  22.7, -16, 41, '25', 'P/2008 Y12'),
    ('Alpha Capricornidi',   '3 Jul',  '15 Aug', 128.0,  20.5, -10, 23, '5',  '169P/NEAT'),
    ('Eta Eridanidi',        '31 Jul', '19 Aug', 135.0,  2.7,  -11, 64, '3',  'C/1852 K1'),
    ('Perseidi',             '17 Jul', '24 Aug', 140.0,  3.2,  +58, 59, '100', '109P/Swift-Tuttle'),
    ('Kappa Cignidi',        '3 Aug',  '28 Aug', 144.0,  19.1, +59, 23, '3',  '2002 GJ8?'),
    ('Alpha Aurigidi',       '28 Aug', '5 Sep',  158.6,  6.1,  +39, 66, '6',  'C/1911 N1 Kiess'),
    ('Eps Perseidi',         '5 Sep',  '21 Sep', 166.7,  3.2,  +40, 64, '8',  '?'),
    ('Lincidi di settembre', '10 Sep', '8 Oct',  170.0,  7.5,  +56, 60, '3',  '?'),
    ('Sextantidi diurne',    '20 Sep', '6 Oct',  188.0,  10.4, -2,  32, '5',  '2005 UD'),
    ('Camelopardalidi ott.', '5 Oct',  '6 Oct',  192.58, 10.9, +79, 47, '5',  '?'),
    ('Draconidi',            '6 Oct',  '10 Oct', 195.4,  17.5, +54, 20, '5',  '21P/Giacobini-Zinner'),
    ('Epsilon Geminidi',     '14 Oct', '27 Oct', 205.0,  6.8,  +27, 70, '3',  'C/1964 N1'),
    ('Orionidi',             '2 Oct',  '7 Nov',  208.0,  6.3,  +16, 66, '20', '1P/Halley'),
    ('Leonis Minoridi',      '19 Oct', '27 Oct', 211.0,  10.8, +37, 62, '2',  'C/1739 K1'),
    ('Tauridi Sud',          '20 Sep', '20 Nov', 223.0,  3.5,  +15, 27, '7',  '2P/Encke'),
    ('Tauridi Nord',         '20 Oct', '10 Dec', 230.0,  3.9,  +22, 29, '5',  '2004 TG10 (Encke)'),
    ('Leonidi',              '6 Nov',  '30 Nov', 235.27, 10.1, +22, 71, '15', '55P/Tempel-Tuttle'),
    ('Alpha Monocerotidi',   '15 Nov', '25 Nov', 239.32, 7.8,  +1,  65, 'var', '?'),
    ('Orionidi di novembre', '13 Nov', '6 Dec',  246.0,  6.1,  +16, 44, '3',  '?'),
    ('Fenicidi',             '1 Dec',  '5 Dec',  249.5,  0.5,  -27, 15, 'var', '289P/Blanpain'),
    ('Puppidi-Velidi',       '1 Dec',  '15 Dec', 255.0,  8.2,  -45, 44, '10', '?'),
    ('Monocerotidi',         '1 Dec',  '19 Dec', 257.0,  6.7,  +8,  41, '3',  'C/1917 F1'),
    ('Sigma Idridi',         '3 Dec',  '20 Dec', 257.0,  8.3,  +2,  58, '7',  'C/2023 P1?'),
    ('Geminidi',             '4 Dec',  '20 Dec', 262.2,  7.5,  +33, 35, '150', '3200 Phaethon'),
    ('Comae Berenicidi',     '4 Dec',  '30 Jan', 271.0,  10.9, +29, 65, '3',  '?'),
    ('Ursidi',               '17 Dec', '26 Dec', 270.7,  14.5, +76, 33, '10', '8P/Tuttle'),
]

SHOWERS = []
for name, a, b, lp, ra_h, dec, v, zhr, parent in _RAW:
    ra = ra_h * 15.0
    SHOWERS.append({'name': name, 'start': _date_to_lam_j2000(a), 'end': _date_to_lam_j2000(b),
                    'peak': lp, 'ra': ra, 'dec': dec, 'v': v, 'zhr': zhr, 'parent': parent,
                    'rad_lon': _ecl_lon(ra, dec)})
# Antielio: attivo 10 dic → 20 set, radiante mobile (punto antisolare)
ANTIHELION = {'name': 'Antielio (ANT, ex Piscidi)', 'start': _date_to_lam_j2000('10 Dec'),
              'end': _date_to_lam_j2000('20 Sep'), 'v': 30, 'zhr': '4', 'parent': 'vari'}

def _in(lam, s, e):
    return (s <= lam <= e) if e >= s else (lam >= s or lam <= e)

def _d(a, b):
    return abs((a - b + 180) % 360 - 180)

def active_showers(sun_lon_of_date, year_frac):
    """Sciami attivi per λ☉ ALLA DATA. Restituisce valori già convertiti alla data."""
    p = precession_to_date(year_frac)
    lam = (sun_lon_of_date - p) % 360            # λ☉ in J2000 per il confronto col catalogo
    out = []
    for s in SHOWERS:
        if _in(lam, s['start'], s['end']):
            out.append({**s, 'peak_date_lon': (s['peak'] + p) % 360,
                        'rad_lon_date': (s['rad_lon'] + p) % 360,
                        'dist_peak': _d(lam, s['peak'])})
    if _in(lam, ANTIHELION['start'], ANTIHELION['end']):
        ant_lon = (sun_lon_of_date + 180) % 360
        out.append({**ANTIHELION, 'peak': None, 'peak_date_lon': None, 'rad_lon_date': ant_lon,
                    'ra': None, 'dec': None, 'dist_peak': None})
    out.sort(key=lambda x: 999 if x['dist_peak'] is None else x['dist_peak'])
    return out

def radiant_radec(name, sun_lon_of_date=None):
    """AR/Dec (gradi, J2000) del radiante; per l'Antielio: punto antisolare sull'eclittica."""
    for s in SHOWERS:
        if s['name'] == name:
            return s['ra'], s['dec']
    if name.startswith('Antielio') and sun_lon_of_date is not None:
        lon = math.radians((sun_lon_of_date + 180) % 360); e = math.radians(OBLIQ_J2000)
        ra = math.degrees(math.atan2(math.sin(lon) * math.cos(e), math.cos(lon))) % 360
        dec = math.degrees(math.asin(math.sin(lon) * math.sin(e)))
        return ra, dec
    return None
