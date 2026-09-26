"""solar_fetch.py — attività solare giornaliera, compatta (pochi KB).
Lanciato dal workflow sky_noaa.yml insieme a noaa_fetch.py → archivio data/sun/
Fonti: NOAA SWPC (brillamenti, regioni attive, indici) e NASA DONKI (brillamenti, CME, flussi veloci)."""
import json, pathlib, datetime, requests
OUT = pathlib.Path('data/sun'); OUT.mkdir(parents=True, exist_ok=True)
now = datetime.datetime.utcnow(); stamp = now.strftime('%Y%m%dT%H%MZ')
d0 = (now - datetime.timedelta(days=7)).strftime('%Y-%m-%d'); d1 = now.strftime('%Y-%m-%d')
DONKI = 'https://kauai.ccmc.gsfc.nasa.gov/DONKI/WS/get/'
SRC = {
    'flares':  'https://services.swpc.noaa.gov/json/goes/primary/xray-flares-7-day.json',
    'regions': 'https://services.swpc.noaa.gov/json/solar_regions.json',
    'indices': 'https://services.swpc.noaa.gov/text/daily-solar-indices.txt',
    'dflr':    f'{DONKI}FLR?startDate={d0}&endDate={d1}',
    'dcme':    f'{DONKI}CME?startDate={d0}&endDate={d1}',
    'dhss':    f'{DONKI}HSS?startDate={d0}&endDate={d1}',
}
for key, url in SRC.items():
    try:
        r = requests.get(url, timeout=60, headers={'User-Agent': 'Mozilla/5.0'}); r.raise_for_status()
        ext = 'txt' if url.endswith('.txt') else 'json'
        data = r.text
        if ext == 'json':
            js = r.json()
            if key == 'regions':      # tieni solo gli ultimi 3 giorni
                lim = (now - datetime.timedelta(days=3)).strftime('%Y-%m-%d')
                js = [x for x in js if str(x.get('observed_date', '9999')) >= lim]
            if key == 'dcme':         # tieni solo i campi utili
                js = [{'startTime': c.get('startTime'), 'note': (c.get('note') or '')[:200],
                       'analyses': [{'speed': a.get('speed'), 'type': a.get('type'), 'isMostAccurate': a.get('isMostAccurate'),
                                     'lat': a.get('latitude'), 'lon': a.get('longitude'), 'halfAngle': a.get('halfAngle')}
                                    for a in (c.get('cmeAnalyses') or [])]} for c in js]
            data = json.dumps(js, separators=(',', ':'))
        (OUT / f'SUN_{key}_{stamp}.{ext}').write_text(data)
        print(key, 'ok', len(data), 'byte')
    except Exception as e:
        print(key, 'ERRORE', e)
limit = now - datetime.timedelta(days=10)
for f in OUT.glob('SUN_*'):
    try:
        if datetime.datetime.strptime(f.stem.split('_')[-1], '%Y%m%dT%H%MZ') < limit: f.unlink()
    except Exception:
        pass
