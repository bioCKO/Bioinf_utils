import time

t0= time.time()

#print time.time()

def time_out():
    t= time.time()
    diff = t - t0
    timer = float(diff)
    timer = timer/60
    print '%1.2f'%timer + ' minutes'

time_out()

time.sleep(6)

time_out()

#print time.time()
