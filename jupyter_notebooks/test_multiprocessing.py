
import multiprocessing
import time
import random
from multiprocessing import Lock
import pandas as pd
import numpy as np

start = time.perf_counter()


pid = np.zeros(10)

val = np.zeros(10)

results = pd.DataFrame(data=None)

results['pid'] = pid
results['val'] = val

lock = Lock()


def do_something(i,lock):
    # print(f'Sleeping {seconds} second(s)...')
    pid = str(i+1)
    r = str(random.random())
    
    lock.acquire()
    
    results.loc[i,'pid'] = pid
    results.loc[i,'val'] = r
    
    lock.release()
    
    print(f'random number from process {pid} = {r}\n')
    


    
# do_something()
# do_something()


# p1 = multiprocessing.Process(target=do_something)
# p2 = multiprocessing.Process(target=do_something)

# p1.start()
# p2.start()

# p1.join()
# p2.join()

processes = []

for i in range(10):
    
    p = multiprocessing.Process(target=do_something,args = (i,lock))
    p.start()
    
    processes.append(p)
    
for process in processes:
    process.join()
    





finish = time.perf_counter()


print(f'Finished in {round(finish-start, 2)} seconds(s)')