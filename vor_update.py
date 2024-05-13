from EQdetect.utils.update_db import upd_db
from EQdetect.utils.vorstat import voronoi_sta
from datetime import datetime
import time

vorcel = voronoi_sta()
upd_db()
vorcel.update_trig()
t0 = datetime.utcnow()
while True:
    diff = datetime.utcnow()-t0
    if diff.total_seconds() > 5*60:
        a=time.time()
        upd_db()
        try:
            vorcel.update_trig()
            print(t0,"Update Database",time.time()-a)
        except Exception as e:
            print(e)
            print(a, "Gagal update voronoi cell")
        t0 = datetime.utcnow()