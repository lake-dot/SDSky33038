#!/usr/bin/env python3
"""
SKY — TEST BASE: notifica di avvio + fetch dati dalle 3 stazioni INTERMAGNET
Nessun calcolo, nessun alert: solo verifica che la catena funzioni.
"""

import requests
import datetime
import os

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


def ntfy_send(message, title="SKY"):
    if not NTFY_TOPIC:
        print("  [ntfy] topic mancante — messaggio solo nei log:")
        print(message)
        return
    try:
        r = requests.post(
            f"https://ntfy.sh/{NTFY_TOPIC}",
            data=message.encode('utf-8'),
            headers={'Title': title, 'Priority': 'high', 'Tags': 'satellite'},
            timeout=15)
        if r.status_code == 200:
            print("  [ntfy] inviato ✓")
        else:
            print(f"  [ntfy] errore HTTP {r.status_code}: {r.text[:150]}")
    except Exception as e:
        print(f"  [ntfy] errore: {e}")


def fetch_last_z(obs_code, hours_back=2):
    """Scarica le ultime ore di dati e ritorna (n_record, ultimo Z valido, ora)."""
    end   = datetime.datetime.now(datetime.timezone.utc) - datetime.timedelta(hours=DATA_LAG_H)
    start = end - datetime.timedelta(hours=hours_back)
    url = (f"{BASE}/data?id={OBS_ID[obs_code]}"
           f"&time.min={start.strftime('%Y-%m-%dT%H:%M:%SZ')}"
           f"&time.max={end.strftime('%Y-%m-%dT%H:%M:%SZ')}"
           f"&format=json")
    r = requests.get(url, timeout=30)
    if r.status_code != 200:
        raise RuntimeError(f"HTTP {r.status_code}")
    records = r.json().get('data', [])
    if not records:
        return 0, None, None
    # cerca l'ultimo valore Z valido (esclude i fill value 99999)
    for rec in reversed(records):
        z = float(rec[1][2])
        if abs(z) < 99999.0:
            return len(records), z, rec[0][11:16]
    return len(records), None, None


def main():
    now = datetime.datetime.now(datetime.timezone.utc)
    print("=" * 50)
    print(f"  SKY test — {now.strftime('%Y-%m-%d %H:%M UT')}")
    print("=" * 50)

    # 1) notifica immediata: verifica catena GitHub → ntfy → iPhone
    ntfy_send("SKY è partito su GitHub — scarico i dati dalle 3 stazioni...",
              title="SKY 🚀 avvio")

    # 2) fetch dati
    lines = []
    for code in OBSERVATORIES:
        print(f"  {code.upper()}...", end=' ')
        try:
            n, z, t = fetch_last_z(code)
            if z is not None:
                print(f"OK — {n} record, ultimo Z={z:.1f} nT ({t} UT)")
                lines.append(f"{code.upper()} ({OBSERVATORIES[code]}): Z={z:.1f} nT  [{t} UT, {n} rec]")
            else:
                print(f"{n} record ma nessun Z valido")
                lines.append(f"{code.upper()}: nessun dato valido")
        except Exception as e:
            print(f"ERRORE: {e}")
            lines.append(f"{code.upper()}: ERRORE — {e}")

    # 3) notifica con i risultati
    ntfy_send("\n".join(lines), title="SKY 📡 letture stazioni")


if __name__ == "__main__":
    main()
