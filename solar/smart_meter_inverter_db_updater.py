__author__ = 'mjohnpayne'


import glob
from time import sleep as sl
import datetime as dt

import datetime

##make time dict with days and times within them day[11/12/15] = time[], day[xxx][11-30] = [solar gen,export value, import value]

in_db = "/Users/mjohnpayne/Documents/start_meter_data/2017_update/database/pv_database.csv"

update_pv = open("/Users/mjohnpayne/Documents/start_meter_data/2017_update/latest_files/oxley_solar_export_full.csv",'r').readlines()

update_meter = open("/Users/mjohnpayne/Documents/start_meter_data/2017_update/latest_files/Victorian Energy Compare Data.csv",'r').readlines()



#outdb structure date   timewindow-30mins   pv_generation   power_import    power_export    consumption(import + pvgen - export)   pv_consumption(pvgen - export)  (all in kwh)


# print today.strftime('We are the %d, %b %Y')
# 'We are the 22, Nov 2008'
# All the letter after a "%" represent a format for something :
#
# %d is the day number
# %m is the month number
# %b is the month abbreviation
# %y is the year last two digits
# %Y is the all year

def read_db(dbpath):
    db = {}
    existdb = open(dbpath,"r").readlines()
    for i in existdb[1:]:
        print i
        col = i.strip().split(',')
        datet = dt.datetime.strptime(col[0]+col[1], '%d/%m/%Y%H:%M')
        db[datet] = col
    return db

def pv_update(pvfile,db):
    power = 0
    for i in pvfile[1:]:
        col = i.strip().split(',')
        datet = dt.datetime.strptime(col[0]+col[1], '%Y%m%d%H:%M')
        #print datet
        if datet.minute == 30 or datet.minute == 00:
            power += float(col[2])
            if datet not in db:
                db[datet] = [datet.strftime('%d/%m/%Y'),datet.strftime('%H:%M'),str(power/12000),"0",'0','0','0','0','0','0','0','0','0','0',datet.strftime('%Y-%m-%d %H:%M:%S')] #12 5min intervals per hour and convert to kWh
            #print "total: " + str(power/12) + "Wh"
            power = 0
        else:
            power += int(col[2])
    return db


def meter_update(meterfile,db,db_old):
    colnames = meterfile[0].strip().split(',')
    for i in meterfile[1:]:
        col = i.strip().split(',')
        date = col[3]
        for j in range(5,len(col)):
            time = colnames[j]
            time = time.split(" - ")[1][:-1]
            datet = dt.datetime.strptime(date+time, '%d/%m/%Y%H:%M')
            if datet in db:
                #print datet
                if col[2] == "Consumption":
                    if col[j] == "":
                        del db[datet]
                    else:
                        db[datet][3] = col[j]
                    #print col[j],"Cons"
                elif col[2] == "Generation":
                    db[datet][4] = col[j]
                    #print db[datet]
                    consumption = str(float(db[datet][3]) + float(db[datet][2]) - float(db[datet][4]))
                    db[datet][5] = consumption
                    pv_consumption = str(float(db[datet][2]) - float(db[datet][4]))
                    db[datet][6] = pv_consumption
                    cost_frm_import = float(db[datet][3])*0.283
                    db[datet][7] = str(cost_frm_import)
                    saving_frm_export = float(db[datet][4])*0.06
                    db[datet][8] = str(saving_frm_export)
                    saving_frm_pv_usage = float(db[datet][6])*0.2838
                    db[datet][9] = str(saving_frm_pv_usage)
                    total_savings = saving_frm_export + saving_frm_pv_usage
                    db[datet][10] = str(total_savings)
                    daily_charge = 0.021770833333333
                    db[datet][11] = str(daily_charge)
                    total_charges = cost_frm_import+daily_charge-saving_frm_export
                    db[datet][12] = str(total_charges)
                    inout = float(db[datet][4])-float(db[datet][3])
                    db[datet][13] = str(inout)
    return db

def write_db(db,outf):
    out = open(outf,"w")
    out.write("Date,Time,pv_gen,import,export,consumption,pv consumption,import cost,saving from export,saving from pv consumption,total savings,daily_charge,total charge,in-out\n")
    for i in sorted(db.keys()):
        print i,db[i]
        out.write(",".join(db[i])+'\n')

                    #print col[j],"Gen"
                #sl(0.2)


###outdb structure date   timewindow-30mins   pv_generation   power_import    power_export    consumption    pv consumed (all in kwh)

###consumption = import + pvgen - export

#
datab = read_db(in_db)


datac = pv_update(update_pv,datab)

datad = meter_update(update_meter,datac,datab)

write_db(datad,"/Users/mjohnpayne/Documents/start_meter_data/2017_update/database/pv_database.csv")
#

#### Daily totals db

#daily totals db structure:     structure date   pv_generation   power_import    power_export    consumption(import + pvgen - export)   pv_consumption(pvgen - export)  cost from import    savings from export    savings from PV usage    total savings(all in kwh)

update_pv_days = open("/Users/mjohnpayne/Documents/start_meter_data/2017_update/latest_files/oxley_solar_export_days.csv",'r').readlines()
update_meter2 = open("/Users/mjohnpayne/Documents/start_meter_data/2017_update/latest_files/Victorian Energy Compare Data.csv",'r').readlines()
#update_meter2 = open("/Users/mjohnpayne/Documents/start_meter_data/electricity-data.csv",'r').read().split('\r')

dailydb = "/Users/mjohnpayne/Documents/start_meter_data/2017_update/database/pv_daily_database.csv"

def read_daily_db(dbpath):
    db = {}
    existdb = open(dbpath,"r").readlines()
    for i in existdb[1:]:
        col = i.strip().split(',')
        datet = dt.datetime.strptime(col[0], '%d/%m/%Y')
        db[datet] = col
    return db

def meter_days_update(meterfile,db):
    colnames = meterfile[0].strip().split(',')
    for i in meterfile[1:]:
        col = i.strip('\r\n').split(',')
        date = col[3]
        datet = dt.datetime.strptime(date, '%d/%m/%Y')
        tot = 0.0
        for j in range(5,len(col)):
            if col[j] == '':
                continue
            else:
                tot+= float(col[j])
        if datet not in db:
            db[datet] = [datet.strftime('%d/%m/%Y'),'0','0','0','0','0','0','0','0','0','0','0','0']
        if col[2] == "Generation":
            db[datet][3] = str(tot)

        elif col[2] == "Consumption":
            db[datet][2] = str(tot)
    return db

def pv_days_update(pvdaysupdate,db):
    colnames = pvdaysupdate[0].strip().split(',')
    for i in pvdaysupdate[1:]:
        col = i.strip().split(',')
        date = col[0]
        datet = dt.datetime.strptime(date, '%Y%m%d')
        if datet not in db:
            db[datet] = [datet.strftime('%d/%m/%Y'),'0','0','0','0','0','0','0','0','0','0','0','0']
        db[datet][1] = str(float(col[1])/1000)
    return db

def add_extras(db):
    for datet in db:
        consumption = str(float(db[datet][2]) + float(db[datet][1]) - float(db[datet][3]))
        db[datet][4] = consumption
        pv_consumption = str(float(db[datet][1]) - float(db[datet][3]))
        db[datet][5] = pv_consumption
        cost_frm_import = float(db[datet][2])*0.283
        db[datet][6] = str(cost_frm_import)
        saving_frm_export = float(db[datet][3])*0.06
        db[datet][7] = str(saving_frm_export)
        saving_frm_pv_usage = float(db[datet][5])*0.2838
        db[datet][8] = str(saving_frm_pv_usage)
        total_savings = saving_frm_export + saving_frm_pv_usage
        db[datet][9] = str(total_savings)
        db[datet][10] = str(1.045)
        db[datet][11] = str(1.045 + cost_frm_import - saving_frm_export)
    return db

def write_daily_db(db,outf):
    out = open(outf,"w")
    out.write("Date,pv_gen,import,export,consumption,pv consumption,import cost,saving from export,saving from pv consumption,total savings,daily charge,total charges\n")
    for i in sorted(db.keys()):
        out.write(",".join(db[i])+'\n')
    out.close()

# dailydb1 = read_daily_db(dailydb)
#
#
# dailydb2 = meter_days_update(update_meter2,dailydb1)
#
# # for i in sorted(dailydb2.keys()):
# #     print i,dailydb2[i]
# #     sl(0.1)
#
# dailydb3 = pv_days_update(update_pv_days,dailydb2)
#
# dailydb4 = add_extras(dailydb3)
#
# write_daily_db(dailydb4,"/Users/mjohnpayne/Documents/start_meter_data/2017_update/database/pv_daily_database.csv")