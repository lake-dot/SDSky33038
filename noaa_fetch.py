"""noaa_fetch.py — salva il vento solare NOAA (ultime 24h) in forma compatta: solo sorgente attiva, solo i campi utili.
Lanciato 2 volte al giorno dal workflow sky_noaa.yml → archivio in data/noaa/"""
import json, pathlib, datetime, requests
OUT = pathlib.Path('data/noaa'); OUT.mkdir(parents=True, exist_ok=True)
URL = {'wind': ('https://services.swpc.noaa.gov/json/rtsw/rtsw_wind_1m.json', ['proton_speed', 'proton_density']),
       'mag':  ('https://services.swpc.noaa.gov/json/rtsw/rtsw_mag_1m.json',  ['bt', 'bz_gsm', 'by_gsm']),
       'eph':  ('https://services.swpc.noaa.gov/json/rtsw/rtsw_ephemerides_1h.json', None)}
stamp = datetime.datetime.utcnow().strftime('%Y%m%dT%H%MZ')
for key, (url, fields) in URL.items():
    try:
        r = requests.get(url, timeout=60, headers={'User-Agent': 'Mozilla/5.0'}); r.raise_for_status()
        data = r.json()
        if fields:
            data = [{'time_tag': d['time_tag'], 'active': True, 'source': d.get('source'), **{f: d.get(f) for f in fields}}
                    for d in data if d.get('active')]
        (OUT / f'NOAA_{key}_{stamp}.json').write_text(json.dumps(data, separators=(',', ':')))
        print(key, len(data), 'record')
    except Exception as e:
        print(key, 'ERRORE', e)
# pulizia: tieni solo gli ultimi 4 giorni di archivio grezzo (i dati utili finiscono nelle cartelle giornaliere)
limit = datetime.datetime.utcnow() - datetime.timedelta(days=4)
for f in OUT.glob('NOAA_*.json'):
    try:
        t = datetime.datetime.strptime(f.stem.split('_')[-1], '%Y%m%dT%H%MZ')
        if t < limit: f.unlink()
    except Exception:
        pass
