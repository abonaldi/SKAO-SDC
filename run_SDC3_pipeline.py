import configparser
import os
import sys
import time
import multiprocessing
from multiprocessing import Manager
from multiprocessing import pool


doHI = False
docontinuum = True
doobserve = False
tstart = time.time()


if doHI:
    config = configparser.ConfigParser()
    config.read('inis/skymodel/HI/SDC2_HI_'+cver+'.ini')
    from skymodel.skymodel_HI import runSkyModel
    runSkyModel(config)
if docontinuum:
    from skymodel.skymodel_continuum import runSkyModel
    from skymodel.skymodel_continuum import runCoadd

    config = configparser.ConfigParser()
    config.read('inis/SDC3a/SDC3_continuum_v4_1.ini')
    n_cores = int(config.getfloat("pipeline", "n_cores"))

    multiprocessing.get_context("fork")
    with Manager() as manager:
        pool = multiprocessing.Pool(n_cores)
        for i in range(n_cores):
            pool.apply_async(
                runSkyModel(config,1,n_cores))
            pool.close()
            pool.join()

    if (n_cores >1):
        runCoadd(config,i+1,n_cores)

    
if doobserve:
    config = configparser.ConfigParser()
    config.read('inis/observe/SDC2_observe_'+cver+'.ini')
    from observe.observe import run_observe
    run_observe(config)

tend = time.time()
print("pipeline finished in {0} seconds.".format(tend - tstart))


