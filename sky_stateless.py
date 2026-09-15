#!/usr/bin/env python3
"""
SKY — TEST v3
- ntfy: POST diretto al topic (come il test browser che funzionava), titolo ASCII
- fetch: finestra temporale adattiva (alcune stazioni pubblicano con più lag)
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

HEADERS = {'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) '
                         'AppleWebKit/537.36 (KHTML, like Gecko) Chrome/126.0 Safari/537.36'}


def ntfy_send(message, title="SKY"):
    """Titolo SOLO ASCII (header HTTP); emoji possibili solo nel corpo."""
    if not NTFY_TOPIC:
        print("  [ntfy] topic mancante — messaggio solo nei log:")
        print(message)
        return
    safe_title = title.encode('ascii', 'ignore').decode().strip() or "SKY"
    try:
        r = requests.post(
            f"https://ntfy.sh/{NTFY_TOPIC}",
            data=message.encode('utf-8'),
            headers={'Title': safe_title, 'Priority': 'high', 'Tags': 'satellite'},
            timeout=15)
        if r.status_code == 200:
            print("  [ntfy] inviato ✓")
        else:
            print(f"  [ntfy] errore HTTP {r.status_code}: {r.text[:150]}")
    except Exception as e:
        print(f"  [ntfy] errore: {e}")


def fetch_last_z(obs_code, hours_back=2):
    """Ultimo Z valido. Se il server dice 'time outside valid range',
    riprova con finestre più indietro (stazioni con lag di pubblicazione maggiore)."""
    now = datetime.datetime.now(datetime.timezone.utc)
    last_err = ""

    for lag_h in (3, 6, 12, 24):
        end   = now - datetime.timedelta(hours=lag_h)
        start = end - datetime.timedelta(hours=hours_back)
        url = (f"{BASE}/data?id={OBS_ID[obs_code]}"
               f"&time.min={start.strftime('%Y-%m-%dT%H:%M:%SZ')}"
               f"&time.max={end.strftime('%Y-%m-%dT%H:%M:%SZ')}"
               f"&format=json")
        try:
            r = requests.get(url, headers=HEADERS, timeout=30)
        except Exception as e:
            last_err = str(e)
            continue

        if r.status_code != 200:
            last_err = f"HTTP {r.status_code} (lag {lag_h}h)"
            if "1405" in r.text:      # time outside valid range → prova più indietro
                continue
            raise RuntimeError(f"{last_err}: {r.text[:120]}")

        records = r.json().get('data', [])
        for rec in reversed(records):
            try:
                z = float(rec[1][2])
                if abs(z) < 99999.0:
                    age_h = (now - datetime.datetime.fromisoformat(
                             rec[0].replace('Z', '+00:00'))).total_seconds() / 3600
                    return len(records), z, rec[0][11:16], age_h
            except (IndexError, ValueError, TypeError):
                continue
        last_err = f"200 ma nessun Z valido (lag {lag_h}h)"

    raise RuntimeError(last_err or "nessun dato trovato")


def main():
    now = datetime.datetime.now(datetime.timezone.utc)
    print("=" * 50)
    print(f"  SKY test v3 — {now.strftime('%Y-%m-%d %H:%M UT')}")
    print("=" * 50)

    ntfy_send("SKY e partito su GitHub — scarico i dati dalle 3 stazioni...",
              title="SKY avvio")

    lines = []
    for code in OBSERVATORIES:
        print(f"  {code.upper()}...", end=' ')
        try:
            n, z, t, age = fetch_last_z(code)
            print(f"OK — {n} rec, Z={z:.1f} nT ({t} UT, età {age:.1f}h)")
            lines.append(f"{code.upper()} ({OBSERVATORIES[code]}): Z={z:.1f} nT  "
                         f"[{t} UT, {age:.1f}h fa, {n} rec]")
        except Exception as e:
            print(f"ERRORE: {e}")
            lines.append(f"{code.upper()}: ERRORE — {str(e)[:100]}")

    ntfy_send("\n".join(lines), title="SKY letture stazioni")


if __name__ == "__main__":
    main()
