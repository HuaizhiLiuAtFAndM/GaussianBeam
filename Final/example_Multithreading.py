import BeamProject as bp
import LaguerreGaussianBeam as LGB
import queue
import multiprocessing
import time
import numpy as np
import concurrent.futures

#WaveLength, ApertureRadius, RefractionIndex, PoAngle, p, l, zNod = 0, RayleighRange = 0,BeamWaist = 0
LGB1 = LGB.GaussianBeam(700e-9,0.001,1,0,0,BeamWaist=0.015)
Screen = bp.Screen(3e-2,3e-2,z=0,resolution=256)
ncores = multiprocessing.cpu_count()  - 2
NumRun = 11
Zval = np.linspace(300e-9,10*300e-9,NumRun)
screens = queue.Queue()
screenNames = queue.Queue()

def CheckSave():
    i = 0
    while (i < NumRun):
        print('checking ' + str(i))
        if screens.empty():
            time.sleep(2)
        else:
            bp.Export(screens.get(), VF = 'ellipse', save = True, show = False, name = screenNames.get())
            i = i+1

def simulate(x):
    if x == -1:
        CheckSave()
    else:
        print("start" + str(x))
        Screen = bp.Screen(3e-2,3e-2,z=Zval[x])
        bp.Render(LGB1,Screen,0,0)
        Screen.z = 0
        screens.put(Screen)
        screenNames.put(str(Zval[x]))
        print("finish" + str(x))
    

if __name__ == "__main__":
    with concurrent.futures.ThreadPoolExecutor(max_workers=ncores) as executor:
        executor.map(simulate, range(-1,NumRun))