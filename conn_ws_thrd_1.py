import websocket, os, time
from datetime import datetime
import threading
from collections import deque
from EQdetect.utils.vorstat import voronoi_sta
from EQdetect.utils.update_db import upd_db
from EQdetect.utils.config import Config
from EQdetect.utils.report import report
from EQdetect.core.sourcecal import Phase, EQsrc, pending_eq
import numpy as np
from multiprocessing import Pool, Process
from itertools import repeat


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
        pick = pick.split()
        if float(pick[13]) < 4:
            received_data.append(pick)

def on_error(ws, error):
    print("Error:", error)
    ws.close()


def on_close(ws, close_status_code, close_msg):
    print("Closed:", close_status_code, close_msg)
    ws.close()

def on_open(ws):
    print("WebSocket opened")

def process_stat(eqp_1, vorcel):
    if eqp_1.status == 1:
        eqp_1.cal_hypo(vorcel)
    return eqp_1

def check_eq(eq1, eq2):
    src1 = eq1.optsrc
    src2 = eq2.optsrc
    if ((abs(src1[0]-src2[0])<=0.5 and abs(src1[1]-src2[1])<=0.5 and (abs(src1[4]-src2[4])<=10))) or (abs(src1[0]-src2[0])<=0.4 and abs(src1[1]-src2[1])<=0.4 and (abs(src1[4]-src2[4])<=60)):
        return True
    elif eq1.status == 0 and eq1.check_trgg(eq2.first):
        return True
    elif eq1.status >= 1 and eq1.check_estg(eq2.first):
        return True
    else:
        return False

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
    vorcel.update_trig()
    picking = []
    noeq = 0
    eqp = []
    eqonl = []
    last_report_time = time.time()
    while True:
        #Pending Earthquake
        tmp_eqp = []
        for eqp_l in eqp:
            eqp_l.eq_stat(vorcel)
            if len(eqp_l.rmphase) != 0:
                for rmpick in eqp_l.rmphase:
                    picking.append(rmpick)
                eqp_l.rmphase = []
            if eqp_l.status == -1:
                print("Remove of Pending eartquake ", eqp_l.starttime)
                continue
            elif eqp_l.status >= 1:
                eqonl.append(eqp_l)
                continue
            if eqp_l.status==0:                             
                tmp_eqp.append(eqp_l)
        eqp = tmp_eqp
        
        #Merge Event
        for i in range(len(eqonl)):
            if eqonl[i].status == -1:
                continue
            for j in range(len(eqonl)):
                if j>i:
                    if eqonl[j].status == -1:
                        continue
                    if eqonl[i].noeq == eqonl[j].noeq:
                        continue
                    status = check_eq(eqonl[i], eqonl[j])
                    if status:
                        print(eqonl[i].noeq, eqonl[j].noeq, i, j)
                        eqonl[j].status = -1
                        eqonl[j].rmphase = []
                        for k in range(len(eqonl[j].phase)):
                            eqonl[i].phase.append(eqonl[j].phase[k])                
            
        #On Going Earthquake
        # if len(eqonl) > 1:
        #     with Pool(processes=5) as p:
        #         result = p.starmap(process_stat, zip(eqonl, repeat(vorcel)), chunksize=2)
        #     eqonl = [eqon for eqon in result]
        #     eqonl.sort(key=lambda x: x.noeq)
        # else:
        #     for i in range(len(eqonl)):
        #         if eqonl[i].status == 1:
        #             eqonl[i].cal_hypo(vorcel)

        tmp_eqon = []            
        for eqon in eqonl:
            if eqon.status == -1:
                print("Remove of On Going eartquake ", eqon.starttime)
                continue
            eqon.eq_stat(vorcel)
            if len(eqon.rmphase) != 0:
                # for rmpick in eqon.rmphase:
                #     picking.append(rmpick)
                eqon.rmphase = []
            if eqon.status == -1:
                print("Remove of On Going eartquake ", eqon.starttime)
                continue
            elif eqon.status == 1:
                if eqon.update<=50:
                    eqon.cal_hypo(vorcel)
                    eqon.endtime = datetime.utcnow()
                    eqon.update = eqon.update + 1
                    if eqon.update > 1:
                        try:
                            file = open("output1.txt", "a")
                            file.write(
                                f"EQ detected No {eqon.noeq} Update {eqon.update} with paramater {eqon.starttime} {datetime.utcnow()} \n")
                            file.write(
                                f"Longitude : {round(eqon.optsrc[0],3)}+-{round(eqon.optsrc[5]*111.1,3)}\nLatitude : {round(eqon.optsrc[1],3)}+-{round(eqon.optsrc[6]*111.1,3)}\nDepth : {round(eqon.optsrc[2],5)}+-{round(eqon.optsrc[7],3)}\nOT : {round(eqon.optsrc[4],2)}\nMag : {round(eqon.optsrc[3],2)}\n")
                            file.write(f"First Trigger : {eqon.first.sta}\n")
                            file.write(f"Count Phase : {len(eqon.phase)}\n")
                            file.write(
                                f"Resampling : {eqon.need_resample}\n\n")
                        except Exception as e:
                            print(eqon.optsrc, eqon.noeq)
                            raise ValueError(e)
                        
                        if eqon.update > 6:
                            diffs = abs(np.diff(np.array(eqon.var), axis=0))
                            converged = np.sqrt(diffs[-1][0]**2+diffs[-1][1]**2)*111.1
                            if converged < 1.0:
                                if eqon.update >= 50:
                                    eqon.status = 2
                                    file.write(
                                        "The values in the array converge.\n\n")
                                else:
                                    file.write(
                                        "The values in the array do not converge.\n\n")
                        if len(eqon.phase) >= 3:
                            if eqon.update >= 4 and eqon.report is False:
                                eqon.report = True
                                
                            if eqon.report is True:
                                eqon.report = 1
                                if time.time()-last_report_time>60:
                                    report(eqon, './eew_report', False)
                                    last_report_time = time.time()  
                                else:
                                    tele_report = Process(target=report, args=(eqon, './eew_report', True,))
                                    tele_report.daemon = True
                                    tele_report.start()
                                    # report(eqon, './eew_report', True) 
                            else:
                                report(eqon, './eew_report', True)
                        file.close()   
            if eqon.status >= 1:                                 
                tmp_eqon.append(eqon)
        eqonl = tmp_eqon
        
        if len(list(received_data)) == 0:
            # time.sleep(0.1)
            continue
        data_pick = list(received_data.copy())

        #Picking Checking
        for pick in data_pick:
            if pick in last_pick or float(pick[13]) > 4:
                continue
            # print(pick)
            picking.append(Phase(pick))
        last_pick = data_pick

        #Picking Processing
        for pick in picking:
            if not eqp and not eqonl:
                if pick.set_first(vorcel):
                    noeq += 1
                    allow_eq, eqpl = pending_eq(EQsrc(pick, cfg, noeq),vorcel)
                    if allow_eq:
                        eqp.append(eqpl)
                continue
            
            statpick = False
            for eqon in eqonl:
                # pick.calerror(eqon.optsrc)
                if eqon.status >= 1:
                    if eqon.check_estg(pick) is True:
                        statpick = True
                        if eqon.status == 2:
                            print(pick.sta, "Late Pick with", eqon.first.sta)
                        break
                    elif pick.sta in eqp_l.estgroup: 
                        print(pick.sta, "estgroup with", eqon.first.sta, "But Big error")
                    else:
                        print(pick.sta, "Not estgroup with", eqon.first.sta)
                            
            if statpick:
                continue
             
            for eqp_l in eqp:
                if eqp_l.status == 0 and not statpick:
                    # pick.calerror(eqp_l.optsrc)
                    if eqp_l.check_trgg(pick):
                        statpick = True
                        break
                    elif pick.sta in eqp_l.triggroup: 
                        print(pick.sta, "Not triggroup with", eqp_l.first.sta, pick.reserr, eqp_l.endtime, eqp_l.starttime)

            if statpick:
                continue
            
            if pick.set_first(vorcel):
                noeq += 1
                allow_eq, eqpl = pending_eq(EQsrc(pick, cfg, noeq),vorcel)
                if allow_eq:
                    if len(eqp) >= 30:
                        eqp.pop(0)  
                    eqp.append(eqpl)
        picking = []
if __name__ == "__main__":
    main()