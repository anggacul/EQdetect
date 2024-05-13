import websocket, os, time
from datetime import datetime
import threading
from collections import deque
from EQdetect.utils.vorstat import voronoi_sta
from EQdetect.utils.update_db import upd_db
from EQdetect.utils.config import Config
from EQdetect.utils.report import report
from EQdetect.core.sourcecal import Phase, EQsrc
import numpy as np


def EwMsgEval(d):
    tl = bytearray(d)[0:1]
    t = ''.join([chr(i) for i in bytearray(d)[1:tl[0]+1]])
    data = 0
    if t == 'TYPE_EEW':
        data = ''.join([chr(i) for i in bytearray(d)[1+tl[0]:]])
    return data


def on_message(ws, message):
    global received_data
    pick = EwMsgEval(message)
    if pick != 0:
        received_data.append(pick.split())


def on_error(ws, error):
    print("Error:", error)


def on_close(ws, close_status_code, close_msg):
    print("Closed:", close_status_code, close_msg)
    ws.close()


def on_open(ws):
    print("WebSocket opened")

def update_database_status(vorcel):
    upd_db()
    vorcel.update_trig()

def process_stat(eqp_1, vorcel):
    # global eqp
    eqp_1.eq_stat(vorcel)
    if eqp_1.status == 1:
        if eqp_1.update <= 50:
            eqp_1.cal_hypo()
            eqp_1.endtime = datetime.utcnow()
            eqp_1.update = eqp_1.update + 1
    if eqp_1.update == 0:
        eqp_1.cal_hypo()
    # for i in range(len(eqp)):
    #     if eqp[i].noeq == eqp_1.noeq:
    #         eqp[i] = eqp_1
    #         break
    #out_queue.put(eqp_1)
    return eqp_1

def pending_eq(eqdata):
    eqdata.cal_hypo()
    if eqdata.optsrc[3] > 10.0 or eqdata.optsrc[3] <= 0.0:
        return False, eqdata
    eqdata.endtime = datetime.utcnow()
    eqdata.update = eqdata.update + 1
    print(
        f"Pending Earthquake No {eqdata.noeq} with paramater {eqdata.starttime}\n")
    file = open("output.txt", "a")
    file.write(
        f"Pending Earthquake No {eqdata.noeq} with paramater {eqdata.starttime}\n")
    file.write(
        f"Longitude : {round(eqdata.optsrc[0],3)}+-{round(eqdata.optsrc[5],3)}\nLatitude : {round(eqdata.optsrc[1],3)}+-{round(eqdata.optsrc[6],3)}\nDepth : {round(eqdata.optsrc[2],5)}\nOT : {round(eqdata.optsrc[4],2)}\nMag : {round(eqdata.optsrc[3],2)}\n")
    file.write(
        f"First Trigger : {eqdata.first.sta} {eqdata.triggroup}\n")
    file.write(
        f"First Trigger prior : {eqdata.first.weight}\n\n")
    file.close()
    return True, eqdata


def main():
    print(os.getcwd())
    # global variable of receiving data from webscoket
    global received_data
    received_data = deque(maxlen=100)
    # WebSocket URL
    ws_url = "ws://localhost:8080/"

    # Create a WebSocket instance
    ws = websocket.WebSocketApp(ws_url,
                                on_message=on_message,
                                on_error=on_error,
                                on_close=on_close)

    def run_websocket():
        ws.run_forever()

    websocket_thread = threading.Thread(target=run_websocket)
    websocket_thread.start()

    last_pick = []
    cfg = Config()
    vorcel = voronoi_sta()
    upd_db()
    vorcel.update_trig()
    t0 = datetime.utcnow()
    picking = []
    noeq = 0
    eqp=[]
    t11 = time.time()
    while True:
        if time.time()-t11>=1:
            t11 = time.time()
            print(time.time())
        diff = datetime.utcnow()-t0
        if diff.total_seconds() > 15*60:
            upd_db()
            vorcel.update_trig()
            t0 = datetime.utcnow()

        # if len(eqp) != 0:
        #     threads = [(eqp_l, vorcel) for eqp_l in eqp]
        #     with Pool(processes=5) as pool:
        #         results = pool.map(process_stat, threads)
        #     eqp = [result for result in results]
        #     eqp.sort(key=lambda x: x.noeq)
    
        neqp = len(eqp)
        if len(eqp) != 0:
            k = 0
            while k < len(eqp):
                eqp[k].eq_stat(vorcel)
                if len(eqp[k].rmphase) != 0:
                    for rmpick in eqp[k].rmphase:
                        picking.append(rmpick)
                    eqp[k].rmphase = []
                if eqp[k].status == -1:
                    print("Remove of pending eartquake ", eqp[k].starttime)
                    eqp.pop(k)
                    if len(eqp) == 0:
                        k = neqp + 1
                    else:
                        k = k
                elif eqp[k].status == 1:
                    if eqp[k].update<=50:
                        eqp[k].cal_hypo()
                        eqp[k].endtime = datetime.utcnow()
                        eqp[k].update = eqp[k].update + 1
                        if eqp[k].update > 1:
                            file = open("output1.txt", "a")
                            file.write(
                                f"EQ detected No {eqp[k].noeq} Update {eqp[k].update} with paramater {eqp[k].starttime} {datetime.utcnow()} \n")
                            file.write(
                                f"Longitude : {round(eqp[k].optsrc[0],3)}+-{round(eqp[k].optsrc[5]*111.1,3)}\nLatitude : {round(eqp[k].optsrc[1],3)}+-{round(eqp[k].optsrc[6]*111.1,3)}\nDepth : {round(eqp[k].optsrc[2],5)}+-{round(eqp[k].optsrc[7],3)}\nOT : {round(eqp[k].optsrc[4],2)}\nMag : {round(eqp[k].optsrc[3],2)}\n")
                            file.write(f"First Trigger : {eqp[k].first.sta}\n")
                            file.write(f"Count Phase : {len(eqp[k].phase)}\n")
                            file.write(
                                f"Resampling : {eqp[k].need_resample}\n\n")
                            if eqp[k].update > 6:
                                diffs = abs(np.diff(np.array(eqp[k].var), axis=0))
                                converged = np.sqrt(
                                    diffs[-1][0]**2+diffs[-1][1]**2)*111.1
                                if converged < 1.0:
                                    if eqp[k].update >= 50:
                                        eqp[k].status = 2
                                    file.write(
                                        "The values in the array converge.\n\n")
                                else:
                                    file.write(
                                        "The values in the array do not converge.\n\n")
                            if len(eqp[k].phase) >= 3:
                                if eqp[k].update >= 4 and eqp[k].report is False:
                                    eqp[k].report = True
                                if eqp[k].report is True:
                                    eqp[k].report = 1
                                    report(eqp[k], './eew_report', False)
                                    
                                else:
                                    report(eqp[k], './eew_report', True)
                            file.close()
                    k = k + 1
                else:
                    k = k+1

        if len(list(received_data)) == 0:
            # time.sleep(0.1)
            continue
        data_pick = list(received_data.copy())

        for pick in data_pick:
            if pick in last_pick or float(pick[13]) > 4:
                continue
            print(pick)
            picking.append(Phase(pick))
        last_pick = data_pick

        npick = len(picking)
        if npick != 0:
            j = 0
            while j < len(picking):
                if len(eqp) == 0:
                    if picking[j].set_first(vorcel) is False:
                        picking.pop(j)
                    else:
                        noeq += 1
                        allow_eq, eqpl = pending_eq(EQsrc(picking[j], cfg, noeq))
                        if allow_eq: 
                            eqp.append(eqpl)
                        picking.pop(j)
                    if len(picking) == 0:
                        j = npick+1
                    else:
                        j = j
                else:
                    for k in range(len(eqp)):
                        statpick = 0
                        if eqp[k].status == 0:
                            if eqp[k].check_trgg(picking[j]) is True:
                                picking.pop(j)
                                statpick = 1
                                break
                            else:
                                print(picking[j].sta, "Not triggroup with",
                                      eqp[k].first.sta, picking[j].reserr, eqp[k].endtime)
                        elif eqp[k].status >= 1:
                            if eqp[k].check_estg(picking[j]) is True:
                                if eqp[k].status == 2:
                                    print(picking[j].sta, "Late Pick with",
                                          eqp[k].first.sta, picking[j].reserr)
                                picking.pop(j)
                                statpick = 1
                                break
                            else:
                                print(picking[j].sta, "Not estgroup with",
                                      eqp[k].first.sta, picking[j].reserr)
                    if statpick == 0:
                        if picking[j].set_first(vorcel) is False:
                            picking.pop(j)
                        else:
                            noeq += 1
                            allow_eq, eqpl = pending_eq(EQsrc(picking[j], cfg, noeq))
                            if allow_eq: 
                                eqp.append(eqpl)
                            picking.pop(j)
                    if len(picking) == 0:
                        j = npick+1
                    else:
                        j = j

if __name__ == "__main__":
    main()
