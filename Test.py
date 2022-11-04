import GaussianBeamProject as GBP
import numpy as np
import concurrent.futures
import multiprocessing
import time
import threading


ncore = multiprocessing.cpu_count()  - 1

NumRun = 31

GB1 = GBP.GaussianBeam(700e-9,0.001,1,0,0,BeamWaist=0.005)
GB2 = GBP.GaussianBeam(700e-9,0.001,1,0,1,BeamWaist=0.005)
Zval = np.linspace(300e-9,10*300e-9,NumRun)
Screens = np.empty(NumRun,dtype=GBP.Screen)

def CheckSave():
    x = 0
    while (x < NumRun):
        print('checking' + str(x))
        if Screens[x] is None:
            time.sleep(2)
        else:
            GBP.Save(Screens[x],str(Zval[x]))
            x = x+1

def simulate(x):
    if x == -1:
        CheckSave()
    else:
        print("start" + str(x))
        Screen = GBP.Screen(3e-2,3e-2,z=Zval[x])
        GBP.Render(GB1,Screen,0,0)
        Screen.z = 0
        GBP.Render(GB2,Screen,0,0)
        #GBP.SaveArray(GBP.GetAngle(Screen),str(Zval[x]))
        Screens[x] = Screen
        print("finish" + str(x))
    



if __name__ == "__main__":
    with concurrent.futures.ThreadPoolExecutor(max_workers=ncore) as executor:
        executor.map(simulate, range(-1,NumRun))



