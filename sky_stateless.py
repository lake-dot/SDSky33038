#!/usr/bin/env python3
"""
SKY — TEST v2: notifica di avvio + fetch dati dalle 3 stazioni INTERMAGNET
Fix: emoji nel titolo (ntfy via JSON), User-Agent browser, diagnosi errori HTTP
"""

import requests
import datetime
import os
from urllib.parse import quote

NTFY_TOPIC = os.environ.get("NTFY_TOPIC")

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

DATA_LAG_H = 3   # lag di pubblicazione INTERMAGNET

# ci fingiamo browser: alcuni server rifiutano python-requests con 400/403
HEADERS = {'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) '
                         'AppleWebKit/537.36 (KHTML, like Gecko) Chrome/126.0 Safari/537.36'}


def ntfy_send(message, title="SKY"):
    """Invia notifica via ntfy in modalità JSON: UTF-8 e emoji funzionano ovunque."""
    if not NTFY_TOPIC:
        print("  [ntfy] topic mancante — messaggio solo nei log:")
        print(message)
        return
    try:
        r = requests.post(
            "https://ntfy.sh/",
            json={'topic': NTFY_TOPIC, 'title': title,
                  'message': message, 'priority': 'high', 'tags': ['satellite']},
            timeout=15)
        if r.status_code == 200:
            print("  [ntfy] inviato ✓")
        else:
            print(f"  [ntfy] errore HTTP {r.status_code}: {r.text[:150]}")
    except Exception as e:
        print(f"  [ntfy] errore: {e}")


def fetch_last_z(obs_code, hours_back=2):
    """Scarica le ultime ore di dati e ritorna (n_record, ultimo Z valido, ora).
    Prova più varianti dell'id e mostra la risposta del server in caso di errore."""
    end   = datetime.datetime.now(datetime.timezone.utc) - datetime.timedelta(hours=DATA_LAG_H)
    start = end - datetime.timedelta(hours=hours_back)
    tmin = start.strftime('%Y-%m-%dT%H:%M:%SZ')
    tmax = end.strftime('%Y-%m-%dT%H:%M:%SZ')

    base_id = OBS_ID[obs_code]
    variants = [base_id,
                base_id.replace('best-avail', 'adjusted'),
                base_id.replace('best-avail', 'reported'),
                quote(base_id, safe='')]          # id interamente URL-encoded

    last_err = ""
    for vid in variants:
        url = f"{BASE}/data?id={vid}&time.min={tmin}&time.max={tmax}&format=json"
        try:
            r = requests.get(url, headers=HEADERS, timeout=30)
            if r.status_code == 200:
                records = r.json().get('data', [])
                if not records:
                    last_err = "200 ma 0 record"
                    continue
                for rec in reversed(records):     # ultimo Z valido (no fill 99999)
                    z = float(rec[1][2])
                    if abs(z) < 99999.0:
                        return len(records), z, rec[0][11:16]
                last_err = "nessun Z valido nei record"
                continue
            last_err = f"HTTP {r.status_code}: {r.text[:200]}"
            print(f"\n    [{vid.split('/')[0]}/{vid.split('/')[1]}] {last_err}")
        except Exception as e:
            last_err = str(e)
            print(f"\n    [fetch] {last_err}")
    raise RuntimeError(last_err or "nessuna variante ha funzionato")


def main():
    now = datetime.datetime.now(datetime.timezone.utc)
    print("=" * 50)
    print(f"  SKY test v2 — {now.strftime('%Y-%m-%d %H:%M UT')}")
    print("=" * 50)

    # 1) notifica immediata
    ntfy_send("SKY è partito su GitHub — scarico i dati dalle 3 stazioni...",
              title="SKY 🚀 avvio")

    # 2) fetch dati
    lines = []
    for code in OBSERVATORIES:
        print(f"  {code.upper()}...", end=' ')
        try:
            n, z, t = fetch_last_z(code)
            print(f"OK — {n} record, ultimo Z={z:.1f} nT ({t} UT)")
            lines.append(f"{code.upper()} ({OBSERVATORIES[code]}): Z={z:.1f} nT  [{t} UT, {n} rec]")
        except Exception as e:
            print(f"ERRORE: {e}")
            lines.append(f"{code.upper()}: ERRORE — {str(e)[:100]}")

    # 3) notifica risultati
    ntfy_send("\n".join(lines), title="SKY 📡 letture stazioni")


if __name__ == "__main__":
    main()
