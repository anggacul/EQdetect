from datetime import datetime
import telegram, json
import asyncio, requests
asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())

async def send_msg(pesan):
    bot = telegram.Bot(token='')
    await bot.send_message(chat_id=-1001775156336, text=pesan)

async def pushsimap(data):
    payload = {
        "originTime": data[0],
        "sent": data[1],
        "epicenterLon": data[2],
        "epicenterLat": data[3],
        "depth": data[4],
        "magnitude": data[5],
        "identifier": data[6],
        "senderName": "AculCorp",
        "references": "BMKG-EEW"
    }
    json_payload = json.dumps(payload)
    headers = {
        "Content-Type": "application/json"
    }
    try:
        response = requests.post("https://event-eews.seismon.my.id/postevent", data=json_payload, headers=headers, timeout=5.0)
        try:
            res = response.json()
        except:
            try:
                res = response.content.decode()
            except:
                res = "Kesalahan Decode"
        if response.status_code == 200:
            print(response.status_code, res['err'])
        else:
            print(response.status_code, res)                
    except requests.exceptions.RequestException as e:
        print("Error", e)
def report(eqsrc, path, status):
    trp = datetime.utcnow()
    trp = trp.strftime("%Y%m%d%H%M%S")
    outf = f"{path}/{eqsrc.noeq}_{trp}_{eqsrc.update}.dat"
    otp = datetime.utcfromtimestamp(round(eqsrc.optsrc[4],2))
    file = open(outf,"w")
    file.write(
        f"Earthquake No {eqsrc.noeq} report on {eqsrc.endtime}\nOptimal Hypocenter\n")
    file.write("OriginTime    Long    Lat    Depth    Mag(Mpd)    Varlong    Varlat    Vardepth RMSE\n")
    file.write(
        f"{otp.strftime('%Y-%m-%d %H:%M:%S')} {round(eqsrc.optsrc[0],3)} {round(eqsrc.optsrc[1],3)} {round(eqsrc.optsrc[2],5)} {round(eqsrc.optsrc[3],2)} {round(eqsrc.optsrc[5],3)} {round(eqsrc.optsrc[6],3)} {round(eqsrc.optsrc[7],3)} {round(eqsrc.rmse,3)}\n"
    )
    file.write(
        f"First Trigger : {eqsrc.first.sta} {eqsrc.triggroup}\n")
    file.write(
        f"Trigger Group : {eqsrc.triggroup}\n"
    )
    file.write(
        f"Estimat Group : {eqsrc.estgroup}\n"
    )
    file.write("Sta    Lon   Lat   Parr   Perr   pa    pv    pd   wei\n")
    for pick in eqsrc.phase:
        file.write(
            f"{pick.sta} {pick.longitude} {pick.latitude} {pick.picktime} {pick.reserr} {pick.pa} {pick.pv} {pick.pd} {pick.weight}\n")
    file.write("OriginTime    Long    Lat    Depth    Mag   Weight\n")
    for pf in eqsrc.particle:
        file.write(
            f"{pf[4]} {round(pf[0],5)} {round(pf[1],5)} {round(pf[2],5)} {pf[3]} {pf[5]}\n"
        )
    file.close()
    pesan = f"Peringatan Dini Gempa Bumi (MulEW Method)\nOrigin Time : {otp.strftime('%Y-%m-%d %H:%M:%S')} UTC\nLongitude : {round(eqsrc.optsrc[0],3)} +- {round(eqsrc.optsrc[5],3)}\nLatitude : {round(eqsrc.optsrc[1],3)} +- {round(eqsrc.optsrc[6],3)}\nDepth : {round(eqsrc.optsrc[2],5)} km\nMag(Mpd) : {round(eqsrc.optsrc[3],2)}\nHasil merupakan update ke {eqsrc.update}"
    if status is False:
        print("Send TELE")
        asyncio.run(send_msg(pesan))
        #data = [otp.strftime('%Y-%m-%d %H:%M:%S'), trp, round(eqsrc.optsrc[0],3), round(eqsrc.optsrc[1],3), round(eqsrc.optsrc[2],5), round(eqsrc.optsrc[3],2), eqsrc.noeq]
        #asyncio.run(pushsimap(data))
