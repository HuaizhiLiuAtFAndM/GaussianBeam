import numpy as np
import polarity as polar
import math
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from scipy.special import genlaguerre
from typing import Callable
#plt.use('SVG')

class Screen:
    def __init__(self, x_range, y_range, resolution = 256, z = 0):
        self.x_range = x_range
        self.y_range = y_range
        self.resolution = resolution
        self.screen = [[polar.vector(0,0)]*self.resolution]*self.resolution

        self.x = np.linspace(-x_range,x_range,self.resolution)
        self.y = np.linspace(-y_range,y_range,self.resolution)

        self.z = z
        self.Beams = np.array(Beam)

class Beam:
    def __init__(self, parameter, expression:callable):
        self.expression = expression
        self.BasicParameter = parameter


def Render(TheBeam: Beam, Screen: Screen, x, y):
    z = Screen.z
    Screen.Beams = np.append(Screen.Beams, TheBeam)
    Theta = np.zeros((Screen.resolution, Screen.resolution))
    Rho = np.zeros((Screen.resolution, Screen.resolution))
    for i in range (Screen.resolution):
        #print("Preparing: " + str(i/Screen.resolution*100) +"%")
        for j in range (Screen.resolution):
            Theta[i,j] = np.arctan2(Screen.x[Screen.resolution-i-1]-x,Screen.y[j]-y)
            Rho[Screen.resolution-i-1,j] = Dist(Screen.x[i], Screen.y[j],x,y)
    X = np.linspace(-Screen.x_range,Screen.x_range,Screen.resolution)
    Y = np.linspace(-Screen.y_range,Screen.y_range,Screen.resolution)
    for i in range(Screen.resolution):
        #print("Rendering: " + str(i/Screen.resolution*100) +"%")
        b = [polar.vector(0,0)]*Screen.resolution
        for j in range(Screen.resolution):
            #print(str(i)+str(j))
            a = Screen.screen[i][j] + TheBeam.expression(TheBeam.BasicParameter, X[j], Y[Screen.resolution-i-1], Screen.z, Rho[i,j], Theta[i,j])
            b[j] = a
        Screen.screen[i] = b
    
def Export(Screen: Screen, VF = "Linear" , save = False , show = True, name: str = "Beam"):
    fig = plt.subplot()
    ABS = [[0]*Screen.resolution]*Screen.resolution
    for i in range(Screen.resolution):
        b = [0]*Screen.resolution
        for j in range (Screen.resolution):
            a = Screen.screen[i][j].Mag()
            b[j] = a
        ABS[i] = b
    extent = [-Screen.x_range,Screen.x_range,-Screen.y_range,Screen.y_range]
    fig.imshow(ABS, cmap='gray', extent= extent)

    if VF == "linear":
        if Screen.resolution < 32:
            VFRes = Screen.resolution
            CoordPick = np.linspace(0,Screen.resolution-1,Screen.resolution,dtype=int)
        else:
            VFRes = 20
            CoordPick = np.linspace(0,Screen.resolution-1,20,dtype=int)

        x = np.linspace(-Screen.x_range,Screen.x_range,VFRes, dtype = float)
        y = np.linspace(-Screen.y_range,Screen.y_range,VFRes, dtype = float)
        #print(x)
        #print(y)
        #print(np.float64(pic[0][0].x))
        u = np.zeros((x.size,y.size))
        v = np.zeros((x.size,y.size))
        for i in range(x.size):
           for j in range(y.size):
               #print(str(i) + str(j))
               u[i,j] = np.float64(np.abs(Screen.screen[CoordPick[i]][CoordPick[j]].x))
               v[i,j] = np.float64(np.abs(Screen.screen[CoordPick[i]][CoordPick[j]].y))
               mag = np.sqrt(u[i,j]**2 + v[i,j]**2)
               u[i,j] = u[i,j]/mag
               v[i,j] = v[i,j]/mag
        fig.quiver(x, y, u, v, color='g')
    elif VF == "ellipse":
        if Screen.resolution < 32:
            VFRes = Screen.resolution
            CoordPick = np.linspace(0,Screen.resolution-1,Screen.resolution,dtype=int)
        else:
            VFRes = 21
            CoordPick = np.linspace(0,Screen.resolution-1,21,dtype=int)
        x = np.linspace(-Screen.x_range,Screen.x_range,VFRes, dtype = float)
        y = np.linspace(-Screen.y_range,Screen.y_range,VFRes, dtype = float)
        I = np.zeros((x.size,y.size))
        Q = np.zeros((x.size,y.size))
        U = np.zeros((x.size,y.size))
        V = np.zeros((x.size,y.size))
        Psi = np.zeros((x.size,y.size))## rotation angle
        Chi = np.zeros((x.size,y.size))
        Hor_Rad = np.zeros((x.size,y.size))
        Ver_Rad = np.zeros((x.size,y.size))
        mag = Screen.x_range/VFRes
        for i in range (x.size):
            for j in range (y.size):
                I[i,j] = np.abs(Screen.screen[CoordPick[i]][CoordPick[j]].x)**2 + np.abs(Screen.screen[CoordPick[i]][CoordPick[j]].y)**2
                Q[i,j] = np.abs(Screen.screen[CoordPick[i]][CoordPick[j]].x)**2 - np.abs(Screen.screen[CoordPick[i]][CoordPick[j]].y)**2
                U[i,j] = 2*np.real(Screen.screen[CoordPick[i]][CoordPick[j]].x*np.conj(Screen.screen[CoordPick[i]][CoordPick[j]].y))
                V[i,j] = -2*np.imag(Screen.screen[CoordPick[i]][CoordPick[j]].x*np.conj(Screen.screen[CoordPick[i]][CoordPick[j]].y))
                Psi[i,j] = np.arctan2(U[i,j],Q[i,j])/2
                Chi[i,j] = np.arctan2(V[i,j],np.sqrt(Q[i,j]**2 + U[i,j]**2))/2
                Hor_Rad[i,j] = np.cos(Chi[i,j])*mag
                Ver_Rad[i,j] = np.sin(Chi[i,j])*mag
        #print(Q)
        #print(U)
        print(Psi)
        for j in range (y.size):
            for i in range (x.size):
                if Chi[j,i] > 0:
                    color = 'r'
                elif Chi[j,i] < 0:
                    color = 'b'
                else:
                    color = 'g'
                ellipse = Ellipse(xy=(x[i], y[y.size-j-1]), width=Hor_Rad[j,i], height=Ver_Rad[j,i], angle=Psi[j,i]*180/np.pi, edgecolor=color, fc='None', lw=0.5)
                fig.add_patch(ellipse)
    
    plt.title('Number of Beams : ' + str(Screen.Beams.size) +"\nRed = Left, Blue = Right, Green = Linear")
    
    if save:
        plt.savefig("Output/" + name + ".png")
    if show:
        plt.show()



def GetAngle(Screen:Screen):
    angle = np.angle(Screen.screen)
    return angle

def Dist(x,y,xb,yb):
    return np.sqrt((x-xb)**2 + (y-yb)**2)

def Angle(x,y,xb,yb): 
    if (x-xb)>=0:
        return np.arctan((y-yb)/(x-xb))
    else:
        return np.arctan((y-yb)/(x-xb)) + np.pi