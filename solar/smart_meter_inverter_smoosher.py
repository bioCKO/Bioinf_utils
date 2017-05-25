__author__ = 'mjohnpayne'


import glob
from time import sleep as sl

invertlis = glob.glob("/Users/mjohnpayne/Documents/start_meter_data/inverter_data/*.csv")
meter_data = open("/Users/mjohnpayne/Documents/start_meter_data/electricity-data.csv","r").read().split('\r')
mushed_solar_data = open("/Users/mjohnpayne/Documents/start_meter_data/mushed_electricity_data.txt","w")


data = {}

for i in invertlis:
    if len(i.split('/')[-1]) < 25:
        mnth = open(i,"r").read().replace('\x00','').split('\r')[9:]
        for j in mnth:
            col = j.split(';')
            date = str(col[0].strip('\n'))
            if len(col[0]) > 1 and col[2] == "---":
                data[date] = [0.0,0.0,0.0]
            elif len(col[0]) > 1:
                data[date] = [float(col[2]),0.0,0.0]

for i in meter_data[1:]:
    i = i.replace(',,','')
    col = i.strip('\n').split(',')
    date = str(col[3])
    if len(date) < 8:
        date = "0" + date
    date = '/'.join(date.split('/')[:2]) + "/20" + date.split('/')[2]
    if col[2] == "Generation":
        gen = sum(map(float,[x for x in col[5:] if len(x) > 0]))
        # print gen
        data[date][1] = gen
    elif col[2] == "Consumption":
        con = sum(map(float,[x for x in col[5:] if len(x) > 0]))
        # print con
        data[date][2] = con

mushed_solar_data.write("Date\tSolar_generation\tSolar_export\tGrid_import\tSolar_used\n")

for i in sorted(data.keys()):
    mushed_solar_data.write('\t'.join(map(str,[i,data[i][0],data[i][1],data[i][2],data[i][0]-data[i][1]]))+'\n')
