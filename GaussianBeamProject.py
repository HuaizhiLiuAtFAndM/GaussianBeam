import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.special import genlaguerre
#plt.use('SVG')



class Screen:
    def __init__(self, x_range, y_range, resolution = 256, z = 0):
        self.x_range = x_range
        self.y_range = y_range
        self.resolution = resolution
        self.screen = np.zeros((self.resolution, self.resolution), dtype='complex')
        self.x = np.linspace(-x_range,x_range,resolution)
        self.y = np.linspace(-y_range,y_range,resolution)
        self.z = z
        self.Beams = np.array(GaussianBeam)


class GaussianBeam:
    def __init__(self, WaveLength, ApertureRadius, ReflectionIndex, p, l, zNod = 0, RayleighRange = 0,BeamWaist = 0):
        
        self.p = p
        self.l = l
        
        self.wNod = BeamWaist
        self.Zr = RayleighRange
        self.Lambda = WaveLength
        self.Ar = ApertureRadius
        self.n = ReflectionIndex
        #should be refractive
        self.k = 2*np.pi/self.Lambda
        self.z0 = zNod
        
        if self.wNod == 0 & self.Zr == 0:
            raise Exception("Please define RayleighRange or Beam Waist")
        elif self.wNod != 0 & self.Zr == 0:
            self.Zr = np.pi*self.wNod**2*self.n/self.Lambda
        elif self.wNod == 0 & self.Zr != 0:
            self.wNod = np.sqrt(self.Zr * self.Lambda / np.pi)

        self.theta = self.Lambda/(np.pi * self.wNod)


def RenderAmplitude(GaussianBeam: GaussianBeam, Screen: Screen, z):
    GausBeam_A = np.zeros((Screen.resolution, Screen.resolution))
    for i in range(Screen.resolution):
        print("Rendering Amplitude: " + str(i/Screen.resolution*100) +"%")
        for j in range(Screen.resolution):
            GausBeam_A[i,j] = Amp(GaussianBeam.p,GaussianBeam.l,Screen.Rho[i,j],z,GaussianBeam.wNod,GaussianBeam.Zr)
    return GausBeam_A

def Render(GaussianBeam: GaussianBeam, Screen: Screen, x, y):
    z = Screen.z
    #Screen.Beams = np.append(Screen.Beams, GaussianBeam)
    Theta = np.zeros((Screen.resolution, Screen.resolution))
    Rho = np.zeros((Screen.resolution, Screen.resolution))
    for i in range (Screen.resolution):
        #print("Preparing: " + str(i/Screen.resolution*100) +"%")
        for j in range (Screen.resolution):
            Theta[Screen.resolution-i-1,j] = Angle(Screen.x[i],Screen.y[j],x,y)
            Rho[Screen.resolution-i-1,j] = Dist(Screen.x[i], Screen.y[j],x,y)  
    for i in range(Screen.resolution):
        #print("Rendering: " + str(i/Screen.resolution*100) +"%")
        for j in range(Screen.resolution):
            Screen.screen[i,j] = Screen.screen[i,j] + Beam(GaussianBeam.p,GaussianBeam.l,Rho[i,j],Theta[i,j],z,GaussianBeam.k,GaussianBeam.wNod,GaussianBeam.Zr,GaussianBeam.Lambda)

def Draw(Screen: Screen):
    ABS = np.abs(Screen.screen)
    extent = [-Screen.x_range,Screen.x_range,-Screen.y_range,Screen.y_range]
    plt.imshow(ABS, cmap='gray', extent= extent)
    plt.title('Gaussian Beam with Cylindrical Coordinates\nNumber of Beams : ' + str(Screen.Beams.size))
    plt.show()

def Save(Screen: Screen, name):
    ABS = np.abs(Screen.screen)
    extent = [-Screen.x_range,Screen.x_range,-Screen.y_range,Screen.y_range]
    plt.imshow(ABS, cmap='gray', extent= extent)
    plt.title('Gaussian Beam with Cylindrical Coordinates\nNumber of Beams : ' + str(Screen.Beams.size))
    plt.savefig("Output/" + name + ".png")


#def RenderAngle(GaussianBeam: GaussianBeam, Screen: Screen, x, y, z):
#    if (np.any(Screen.screen) != 0):
#        raise Exception("Screen may not be rendered with angle of multiple beams")
#    Theta = np.zeros((Screen.resolution, Screen.resolution))
#    Rho = np.zeros((Screen.resolution, Screen.resolution))
#    for i in range (Screen.resolution):
#        print("Preparing: " + str(i/Screen.resolution*100) +"%")
#        for j in range (Screen.resolution):
#            Theta[i,j] = Angle(Screen.x[i],Screen.y[j],x,y)
#            Rho[i,j] = Dist(Screen.x[i], Screen.y[j],x,y)  
#    for i in range(Screen.resolution):
#        print("Rendering: " + str(i/Screen.resolution*100) +"%")
#        for j in range(Screen.resolution):
#            Screen.screen[i,j] = np.angle(U(GaussianBeam.p,GaussianBeam.l,Rho[i,j],Theta[i,j],z,GaussianBeam.k,GaussianBeam.wNod,GaussianBeam.Zr))

def GetAngle(Screen:Screen):
    angle = np.angle(Screen.screen)
    return angle

def DrawArray(Array):
    plt.imshow(Array, cmap = 'gray')
    plt.show()


def SaveArray(Array, name):
    plt.imshow(Array, cmap = 'gray')
    plt.savefig("Output/" + name + ".png")


def w(wNod,z,Zr):
    return wNod*np.sqrt(1+z**2/Zr**2)
def A(p,l):
    return math.factorial(p)*(2/(np.pi*math.factorial(p)*math.factorial((np.abs(l)+p))))**0.5
def R(z,Zr):
    return z*(1+Zr**2/z**2)
def Phi(z,Zr): 
    return np.arctan(z/Zr)
def Amp(p,l,r,z,wNod,Zr):
    # Amplitude of an LG mode
    return (((np.sqrt(2)*r)/(w(wNod,z,Zr)))**(np.abs(l))*genlaguerre(p,np.abs(l))((2*r**2)/w(wNod,z,Zr)**2)*np.exp((-r**2)/(w(wNod,z,Zr)**2)))

def U(p,l,r,theta,z,k,wNod,Zr):
    #High Order Gaussian Beams in Cylindrical Coodinates
    if z == 0:
        return A(p,l)/w(wNod,z,Zr)*Amp(p,l,r,z,wNod,Zr)*np.exp(1j*l*theta)
    else:
        return (A(p,l)/w(wNod,z,Zr)*Amp(p,l,r,z,wNod,Zr)
        *np.exp((1j*k*r**2)/(2*R(z,Zr)))
        *np.exp(1j*l*theta)
        *np.exp(-1j*Phi(z,Zr)))
def Beam(p,l,r,theta,z,k,wNod,Zr,Lambda):
    return U(p,l,r,theta,z,k,wNod,Zr)*np.exp(1j*(2*np.pi/Lambda)*z)
def Dist(x,y,xb,yb):
    return np.sqrt((x-xb)**2 + (y-yb)**2)
def Angle(x,y,xb,yb): 
    if (x-xb)>=0:
        return np.arctan((y-yb)/(x-xb))
    else:
        return np.arctan((y-yb)/(x-xb)) + np.pi
    
