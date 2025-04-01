import configparser
import os
import sys
import time
import multiprocessing as mp
from multiprocessing import Manager
from multiprocessing import pool



tstart = time.time()


from skymodel.skymodel_continuum import runSkyModel
from skymodel.skymodel_continuum import runCoadd

config = configparser.ConfigParser()
config.read('inis/SDC3a/SDC3_continuum_v4_1.ini')
n_cores = int(config.getfloat("pipeline", "n_cores"))

mp.get_context("fork")
with mp.Manager() as manager:
    pool = mp.Pool(n_cores)
    for i in range(n_cores):
        pool.apply_async(
            runSkyModel(config,i+1,n_cores))
    pool.close()
    pool.join()

if (n_cores >1):
    runCoadd(config,i+1,n_cores)

    
tend = time.time()
print("pipeline finished in {0} seconds.".format(tend - tstart))


