import LGBeamProject as GBP
import numpy as np
import concurrent.futures
import multiprocessing
import time
import threading
import queue

GB1 = GBP.LGaussianBeam(700e-9,0.001,1,np.pi/2,0,0,BeamWaist=0.005)
Screen = GBP.Screen(3e-2,3e-2,z=0,resolution=100)


GBP.Render(GB1,Screen,0,0)
GBP.Draw(Screen)