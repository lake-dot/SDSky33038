"""send_mail.py — invia il report giornaliero con allegati. Secrets: MAIL_USER, MAIL_PASS (password per app), MAIL_TO"""
import os, sys, smtplib, pathlib, mimetypes
from email.message import EmailMessage
if not os.environ.get('MAIL_USER'):
    print('email non configurata: salto'); sys.exit(0)
folder = pathlib.Path(sys.argv[1]); day = folder.name.replace('SKY_', '')
msg = EmailMessage()
msg['Subject'] = f'SKY {day} — dati e SCAN'
msg['From'] = os.environ['MAIL_USER']; msg['To'] = os.environ.get('MAIL_TO', os.environ['MAIL_USER'])
scan = folder / f'SCAN_{day}.txt'
msg.set_content(scan.read_text(encoding='utf-8') if scan.exists() else 'SCAN non disponibile')
for pat in (f'SCAN_{day}.txt', 'RIEPILOGO.txt', f'TUTTE_{day}.csv.gz', 'NOAA_*.json'):
    for f in sorted(folder.glob(pat)):
        ctype = mimetypes.guess_type(f.name)[0] or 'application/octet-stream'
        maint, sub = ctype.split('/', 1)
        msg.add_attachment(f.read_bytes(), maintype=maint, subtype=sub, filename=f.name)
with smtplib.SMTP_SSL(os.environ.get('MAIL_SMTP', 'smtp.gmail.com'), 465) as s:
    s.login(os.environ['MAIL_USER'], os.environ['MAIL_PASS']); s.send_message(msg)
print('email inviata a', msg['To'])
