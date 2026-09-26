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
• Fonte: IAU Meteor Data Center, lista degli sciami STABILITI (113), aggiornamento 21 set 2026
  (file streamestablisheddata2026.txt). Per ogni sciame si usa la soluzione con finestra di
  attività e deriva note e con più meteore (N). Include gli sciami DIURNI (visibili solo ai radar).
  Finestra mancante nella fonte → ±10° attorno al picco (segnalato "finestra stimata").
• Rami di Encke separati (Tauridi Sud / Nord): niente più voce generica "Encke nodo".
• Deriva dei radianti: AR e Dec si spostano ogni giorno (gradi per grado di λ☉ ≈ per giorno).
  Valori IMO/letteratura dove noti; per gli sciami senza dato si usa +1.0°/0.0° (tipico) e il
  report lo segnala con "(deriva stimata)". Con la deriva, il passaggio del radiante può restare
  quasi alla stessa ora solare ogni giorno: per questo l'orologio "radiante" va testato a parte.
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


IAU = [["KSE", "kappa-Serpentids", 6.0, 20.0, 13.2, 230.7, 20.4, 0.66, -0.28, 42.2, 374, ""], ["APS", "Daytime April Piscids", 15.0, 46.0, 29.3, 4.9, 5.5, 1.01, -0.45, 25.7, 3178, ""], ["PPU", "pi-Puppids", 27.0, 39.0, 30.1, 110.1, -46.4, 1.59, -0.43, 15.9, 12, "26P/Grigg-Skjellerup"], ["LYR", "April Lyrids", 26.0, 37.0, 32.3, 272.2, 33.4, 0.71, -0.32, 46.6, 9980, "C/1861 G1 (Thatcher)"], ["HVI", "h-Virginids", 30.0, 47.0, 38.8, 203.0, -10.9, 0.29, -0.25, 18.8, 733, ""], ["ARC", "April rho-Cygnids", 20.0, 47.0, 40.0, 321.6, 47.3, 0.07, 0.02, 40.8, 662, ""], ["BAQ", "beta-Aquariids", 38.0, 62.0, 41.8, 319.6, -1.4, 0.79, 0.22, 68.3, 339, ""], ["SMA", "Southern Daytime May Arietids", 21.0, 65.0, 44.9, 27.1, 6.3, 0.92, 0.34, 25.6, 10276, ""], ["NOC", "Northern Daytime omega-Cetids", 16.0, 58.0, 45.5, 9.0, 17.3, 0.95, 0.36, 36.8, 1256, ""], ["OCE", "Southern Daytime omega-Cetids", 18.0, 62.0, 45.5, 20.5, -6.1, 0.93, 0.44, 36.9, 837, "C/2003 Q1 (SOHO)?"], ["ETA", "eta-Aquariids", 14.0, 100.0, 45.7, 338.8, -0.4, 0.75, 0.36, 66.3, 46032, "1P/Halley"], ["GAQ", "gamma-Aquilids", 42.0, 55.0, 48.7, 304.6, 14.5, 0.91, 0.24, 62.4, 197, "C/1853 G1 (Schweizer)"], ["ELY", "eta-Lyrids", 45.0, 55.0, 50.0, 290.9, 43.7, -0.16, 0.04, 43.9, 990, "C/1983 H1 (IRAS-Araki-Alcock)"], ["JMC", "June mu-Cassiopeiids", 50.0, 85.0, 70.5, 8.3, 52.3, 1.05, 0.57, 43.0, 726, ""], ["ARI", "Daytime Arietids", 64.0, 88.0, 74.5, 41.7, 23.6, 0.6, 0.19, 39.1, 2142, ""], ["ZPE", "Daytime zeta-Perseids", 58.0, 88.0, 74.5, 57.4, 23.4, 1.0, 0.2, 26.4, 843, ""], ["TAH", "tau-Herculids", 72.0, 80.0, 77.9, 231.0, 36.4, 0.57, -0.24, 15.2, 2370, "73P/Schwassmann-Wachmann 3"], ["AVB", "alpha-Virginids", 60.0, 111.0, 83.4, 229.3, 47.1, 0.29, 1.13, 14.2, 416, ""], ["JRC", "June rho-Cygnids", 82.0, 87.0, 84.2, 320.8, 44.4, 1.12, 0.23, 49.9, 229, ""], ["SSG", "Southern mu-Sagittariids", 64.0, 110.0, 85.0, 272.5, -29.6, 0.71, -0.03, 25.9, 2528, ""], ["DLT", "Daytime lambda-Taurids", 70.0, 98.0, 85.5, 56.7, 11.5, 0.82, 0.27, 36.4, 406, ""], ["JBO", "June Bootids", 91.8, 93.2, 92.3, 223.5, 47.2, 0.84, -0.22, 13.8, 55, "7P/Pons-Winnecke"], ["BTA", "Daytime beta-Taurids", 90.0, 100.0, 93.5, 82.0, 20.0, 0.89, 0.04, 27.4, 288, ""], ["JIP", "June iota-Pegasids", 93.0, 96.0, 94.0, 331.8, 29.2, 1.04, 0.31, 58.8, 244, ""], ["COR", "Corvids", 93.0, 100.0, 94.9, 191.6, -19.2, None, None, 14.1, 7470, ""], ["SZC", "Southern June Aquilids", 64.0, 124.0, 101.8, 317.8, -27.7, 0.93, 0.27, 39.4, 3085, ""], ["NZC", "Northern June Aquilids", 73.0, 148.0, 104.1, 312.7, -2.9, 0.81, 0.26, 38.6, 10321, "C/2009 U10 (SOHO)"], ["ALA", "alpha-Lacertids", 100.0, 115.0, 107.0, 343.0, 53.8, 1.13, 0.33, 37.5, 78, ""], ["ZCS", "zeta-Cassiopeiids", 98.0, 120.0, 107.0, 6.1, 50.5, 1.41, 0.29, 57.2, 1740, "109P/Swift-Tuttle"], ["CAN", "c-Andromedids", 90.0, 136.0, 108.2, 29.9, 47.6, 1.21, 0.36, 57.2, 2927, ""], ["FAN", "49-Andromedids", 78.0, 130.0, 109.9, 17.8, 45.0, 0.87, 0.35, 59.8, 2232, "C/2001 W2 (BATTERS)"], ["JXA", "July xi-Arietids", 84.0, 133.0, 111.8, 35.4, 8.6, 0.75, 0.22, 68.9, 1216, "C/1964 N1 (Ikeya)"], ["JPE", "July Pegasids", 98.0, 146.0, 113.1, 351.7, 12.0, 0.85, 0.29, 63.9, 1821, "C?1979 Y1 (Bradfield)"], ["PCA", "psi-Cassiopeiids", 110.0, 118.0, 115.3, 26.0, 70.0, 1.02, 0.37, 42.2, 176, ""], ["PPS", "phi-Piscids", 81.0, 152.0, 117.1, 28.1, 30.5, 0.85, 0.37, 66.5, 10512, "1P/Halley"], ["GDR", "July gamma-Draconids", 117.0, 131.0, 125.1, 279.9, 50.5, 0.01, 0.07, 27.4, 1147, "C/1919 q2 (Metcalf)"], ["CAP", "alpha-Capricornids", 99.0, 146.0, 125.8, 304.8, -9.5, 0.61, 0.24, 22.7, 13532, "169P/NEAT"], ["SDA", "Southern delta-Aquariids", 113.0, 157.0, 127.6, 340.6, -16.2, 0.81, 0.21, 40.6, 74507, "342P/SOHO"], ["PAU", "Piscis Austrinids", 121.0, 151.0, 134.1, 352.5, -21.5, 0.76, 0.28, 44.2, 1178, ""], ["ERI", "eta-Eridanids", 120.0, 149.0, 134.4, 41.4, -13.3, 0.85, 0.28, 64.4, 4272, "C/1852 K1 (Chacornac)"], ["PER", "Perseids", 115.0, 175.0, 140.1, 47.9, 57.9, 1.45, 0.2, 58.8, 141256, "109P/Swift-Tuttle"], ["KCG", "kappa-Cygnids", 108.0, 168.0, 140.8, 285.6, 51.4, 0.6, 0.79, 22.7, 2754, "2021 HK12"], ["BHY", "beta-Hydrusids", 138.0, 144.0, 142.3, 23.2, -80.9, 2.15, 0.84, 21.8, 5, ""], ["NDA", "Northern delta-Aquariids", 122.0, 180.0, 145.4, 351.6, 4.1, 0.79, 0.34, 38.4, 9033, "96P/Machholz"], ["AUD", "August Draconids", 127.0, 168.0, 148.5, 269.8, 59.3, -0.87, 0.25, 21.4, 4502, ""], ["OMG", "omicron-Geminids", 135.0, 166.0, 148.9, 99.4, 43.2, 1.23, -0.06, 56.4, 619, ""], ["AGC", "August gamma-Cepheids", 147.0, 161.0, 154.6, 358.9, 76.7, 0.55, 0.2, 43.7, 744, ""], ["AUR", "Aurigids", 153.0, 166.0, 159.0, 91.8, 39.0, 1.0, -0.01, 66.77, 1700, ""], ["PSO", "pi6-Orionids", 138.0, 183.0, 162.0, 69.3, -1.2, 0.9, 0.16, 66.1, 5421, "1P/Halley"], ["NPI", "Northern delta-Piscids", 153.0, 174.0, 162.3, 1.5, 5.8, 0.79, 0.39, 27.4, 235, ""], ["NUE", "nu-Eridanids", 150.0, 180.0, 166.9, 67.1, 0.7, 0.79, 0.37, 66.1, 5674, "1P/Halley"], ["SPE", "September epsilon-Perseids", 160.8, 199.71, 170.36, 51.7, 39.5, 1.19, 0.05, 64.6, 299, ""], ["CCY", "chi-Cygnids", 120.0, 182.0, 172.2, 300.5, 32.4, -0.44, 1.09, 14.9, 1639, ""], ["NIA", "Northern iota-Aquariids", 126.0, 248.0, 174.8, 7.5, 7.5, 0.74, 0.29, 26.9, 6984, ""], ["KLE", "Daytime kappa-Leonids", 164.0, 200.0, 183.0, 162.3, 14.9, 0.62, -0.3, 43.3, 1366, ""], ["DSX", "Daytime Sextantids", 174.0, 197.0, 186.0, 154.3, -1.0, 0.56, -0.54, 31.3, 1292, "2005 UD"], ["OCT", "October Camelopardalids", 192.0, 192.9, 192.5, 169.1, 78.6, 0.51, -0.48, 45.9, 333, ""], ["DRA", "October Draconids", 194.7, 195.8, 195.4, 262.8, 55.9, 0.27, 0.04, 20.8, 1058, "21P/Giacobini-Zinner"], ["OCU", "October Ursae Majorids", 195.0, 210.0, 202.4, 145.4, 64.1, 1.82, -0.62, 55.5, 935, ""], ["XDR", "xi-Draconids", 181.0, 230.0, 204.7, 161.1, 75.8, 1.48, -0.3, 38.0, 148, ""], ["EGE", "epsilon-Geminids", 188.0, 222.0, 205.2, 101.1, 27.9, 0.85, -0.1, 69.2, 2321, ""], ["TCA", "tau-Cancrids", 171.0, 227.0, 206.3, 136.6, 29.7, 0.9, 0.0, 67.2, 1479, ""], ["STA", "Southern Taurids", 155.0, 270.0, 207.0, 40.1, 10.7, 0.74, 0.16, 27.8, 35213, "2P/Encke"], ["LMI", "Leonis Minorids", 202.0, 224.0, 209.0, 160.1, 36.7, 1.1, -0.31, 61.4, 2054, "C/1739 K1 (Zanotti)"], ["ORI", "Orionids", 198.0, 236.0, 209.5, 96.3, 15.7, 0.78, 0.03, 66.1, 72660, "1P/Halley"], ["LUM", "lambda-Ursae Majorids", 211.0, 218.0, 214.4, 157.9, 49.3, 1.0, -0.47, 60.8, 283, "C/1975 T2 (Suzuki-Saigusa-Mori)"], ["SLD", "Southern lambda-Draconids", 218.0, 224.0, 221.4, 161.3, 68.2, 0.52, -0.5, 49.3, 140, ""], ["CTA", "chi-Taurids", 186.0, 248.0, 222.0, 63.7, 25.8, 0.95, 0.1, 42.6, 2173, ""], ["KUM", "kappa-Ursae Majorids", 220.0, 237.0, 223.7, 145.5, 45.2, 1.18, -0.36, 65.2, 735, ""], ["RPU", "rho-Puppids", 210.0, 237.0, 225.3, 124.8, -25.4, 1.17, -0.07, 57.4, 517, "C/1879 M1 (Swift)"], ["NTA", "Northern Taurids", 150.0, 270.0, 227.3, 56.6, 22.3, 0.83, 0.2, 28.2, 29599, "2P/Encke"], ["OER", "omicron-Eridanids", 186.0, 311.0, 228.1, 56.1, -2.2, 0.6, -0.11, 27.8, 5946, "2015 KK"], ["AND", "Andromedids", 210.0, 252.0, 230.3, 22.2, 31.7, 0.18, 0.76, 17.4, 2479, "3D/Biela"], ["LEO", "Leonids", 213.0, 274.0, 236.3, 154.4, 21.4, 0.63, -0.38, 70.3, 13815, "55P/Tempel-Tuttle"], ["THA", "November theta-Aurigids", 233.0, 239.0, 237.0, 89.0, 34.7, 1.49, 0.14, 33.8, 1180, ""], ["AMO", "alpha-Monocerotids", 230.0, 243.0, 239.6, 117.3, 0.8, 0.85, -0.1, 62.1, 370, ""], ["ACA", "alpha-Canis Majorids", 221.0, 260.0, 242.7, 98.0, -18.9, 0.75, 0.17, 43.7, 1644, ""], ["NOO", "November Orionids", 218.0, 270.0, 247.0, 90.7, 15.4, 0.74, -0.03, 42.9, 8909, ""], ["ORS", "Southern chi-Orionids", 224.0, 270.0, 249.3, 78.4, 17.9, 0.85, 0.04, 27.2, 4395, "2010 LU108"], ["PHO", "Phoenicids", 249.0, 250.0, 249.6, 7.4, -27.8, 0.08, -1.07, 8.9, 5, ""], ["DKD", "December kappa-Draconids", 244.0, 261.0, 251.2, 187.1, 69.9, 1.05, -0.5, 43.5, 2387, ""], ["DPC", "December phi-Cassiopeiids", 250.0, 254.0, 252.0, 19.7, 57.1, -0.6, 1.3, 16.1, 261, "3D/Biela"], ["PSU", "psi-Ursae Majorids", 247.0, 263.0, 252.7, 170.0, 43.2, 1.02, -0.61, 61.7, 457, ""], ["DAD", "December alpha-Draconids", 231.0, 270.0, 256.0, 210.5, 58.3, 0.53, -0.23, 41.5, 2871, "2009 WN25"], ["HYD", "sigma-Hydrids", 240.0, 298.0, 256.6, 125.6, 2.4, 0.85, -0.29, 58.9, 15659, "C/2023 P1 (Nishimura)"], ["DRV", "December rho-Virginids", 243.0, 275.0, 257.1, 188.6, 12.7, 0.89, -0.11, 68.4, 541, "C/1961 T1 (Seki)"], ["MON", "December Monocerotids", 242.0, 275.0, 258.6, 100.8, 8.1, 0.66, -0.13, 41.4, 5710, "C/1917 F1 (Mellish)"], ["EHY", "eta-Hydrids", 240.0, 282.0, 260.6, 135.9, 1.5, 0.84, -0.22, 62.4, 2196, ""], ["GEM", "Geminids", 225.0, 273.0, 261.6, 113.0, 32.4, 1.04, -0.16, 33.8, 104490, "3200 Phaethon"], ["XVI", "December chi-Virginids", 245.0, 287.0, 262.8, 191.4, -10.5, 0.67, -0.35, 68.5, 1481, ""], ["URS", "Ursids", 255.0, 286.0, 270.5, 218.5, 75.9, 2.43, -0.4, 32.7, 5606, "8P/Tuttle"], ["COM", "Comae Berenicids", 246.0, 335.0, 271.3, 164.3, 29.2, 0.86, -0.42, 63.2, 13811, ""], ["OSE", "omega-Serpentids", 256.0, 291.0, 273.9, 241.2, -1.0, 0.74, -0.04, 38.7, 5350, ""], ["SSE", "sigma-Serpentids", 255.0, 291.0, 275.0, 242.4, -0.1, 0.64, 0.03, 42.3, 1075, ""], ["DSV", "December sigma-Virginids", 243.0, 306.0, 275.1, 211.2, 3.6, 0.83, -0.18, 66.4, 1687, "C/1846 J1 (Brorsen)"], ["AHY", "alpha-Hydrids", 251.0, 305.0, 277.1, 122.4, -9.2, 0.67, 0.06, 43.8, 3231, ""], ["SCC", "Southern delta-Cancrids", 270.0, 300.0, 280.9, 112.2, 17.1, 0.85, -0.14, 27.2, 581, "2P/Encke"], ["JLE", "January Leonids", 279.0, 287.0, 282.0, 148.2, 23.7, 0.7, -0.13, 52.3, 1160, ""], ["QUA", "Quadrantids", 271.0, 294.0, 283.1, 229.6, 49.8, 0.86, -0.38, 40.4, 13594, "2003 EH1"], ["XCB", "xi-Coronae Borealids", 287.0, 304.0, 295.0, 247.0, 30.3, 0.39, 0.11, 44.8, 2621, ""], ["TCB", "theta-Coronae Borealids", 287.0, 304.0, 296.0, 233.6, 34.4, 0.3, 0.16, 37.7, 3560, ""], ["XUM", "January xi-Ursae Majorids", 292.0, 305.0, 298.5, 169.6, 32.8, 0.5, -0.42, 40.8, 492, ""], ["LBO", "lambda-Bootids", 279.0, 330.0, 298.8, 232.2, 39.9, 0.85, -0.27, 47.2, 2932, ""], ["NCC", "Northern delta-Cancrids", 290.0, 308.0, 299.0, 131.4, 17.6, 0.4, -0.2, 27.73, 900, ""], ["GUM", "gamma-Ursae Minorids", 294.0, 304.0, 299.0, 231.8, 66.8, 0.7, -0.57, 31.8, 694, ""], ["ECV", "eta-Corvids", 267.0, 359.0, 300.5, 190.3, -18.0, 0.83, -0.16, 67.8, 1881, ""], ["ACB", "alpha-Coronae Borealids", 307.0, 316.0, 308.0, 231.5, 26.0, 1.8, -1.0, 56.71, 500, ""], ["OHY", "omicron-Hydrids", 291.0, 325.0, 309.6, 177.5, -34.2, 0.95, -0.35, 59.2, 1195, ""], ["FEV", "February epsilon-Virginids", 300.0, 335.0, 310.7, 197.8, 12.5, 0.82, -0.38, 62.7, 953, ""], ["AAN", "alpha-Antliids", 295.0, 332.0, 312.0, 160.7, -12.3, 0.745, -0.36, 43.2, 1228, ""], ["FED", "February eta-Draconids", 313.0, 317.0, 315.2, 239.6, 62.0, 0.69, 0.33, 35.2, 139, ""], ["XHE", "x-Herculids", 343.0, 353.0, 351.2, 253.7, 49.2, 0.91, -0.28, 34.7, 157, ""], ["EVI", "eta-Virginids", 347.0, 4.0, 357.0, 185.1, 3.7, 0.63, -0.2, 26.9, 955, "D/1766 G1 (Helfnzrieder)"]]   # code, nome, λ inizio, λ fine, λ picco (J2000), AR, Dec, dAR, dDec, Vg, N, progenitore

SHOWERS = []
for code, name, a, b, lp, ra, dec, dra, ddec, vg, n, parent in IAU:
    est = a is None or b is None
    SHOWERS.append({'code': code, 'name': name,
                    'start': (lp - 10) % 360 if est else a, 'end': (lp + 10) % 360 if est else b,
                    'window_est': est, 'peak': lp, 'ra': ra, 'dec': dec, 'v': vg, 'zhr': '?',
                    'n': n, 'parent': parent, 'rad_lon': _ecl_lon(ra, dec),
                    'drift': (dra if dra is not None else 1.0, ddec if ddec is not None else 0.0),
                    'drift_known': dra is not None,
                    'daytime': name.lower().startswith('daytime')})

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

def radiant_radec(name, sun_lon_of_date=None, year_frac=2026.7):
    """AR/Dec (gradi) del radiante ALLA DATA, con la deriva giornaliera; Antielio: punto antisolare."""
    for s in SHOWERS:
        if s['name'] == name:
            if sun_lon_of_date is None:
                return s['ra'], s['dec']
            lam = (sun_lon_of_date - precession_to_date(year_frac)) % 360
            dl = (lam - s['peak'] + 180) % 360 - 180
            dra, ddec = s['drift']
            return (s['ra'] + dra * dl) % 360, s['dec'] + ddec * dl
    if name.startswith('Antielio') and sun_lon_of_date is not None:
        lon = math.radians((sun_lon_of_date + 180) % 360); e = math.radians(OBLIQ_J2000)
        ra = math.degrees(math.atan2(math.sin(lon) * math.cos(e), math.cos(lon))) % 360
        dec = math.degrees(math.asin(math.sin(lon) * math.sin(e)))
        return ra, dec
    return None
