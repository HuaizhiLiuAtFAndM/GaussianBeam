import BeamProject as bp
import numpy as np
import math
import polarity as polar
from scipy.special import genlaguerre


#Following Code are implication of Lagere Gaussian Beam
def LGw(wNod,z,Zr):
    return wNod*np.sqrt(1+z**2/Zr**2)
def LGA(p,l):
    return math.factorial(p)*(2/(np.pi*math.factorial(p)*math.factorial((np.abs(l)+p))))**0.5
def LGR(z,Zr):
    return z*(1+Zr**2/z**2)
def LGPhi(z,Zr): 
    return np.arctan(z/Zr)
def LGAmp(p,l,r,z,wNod,Zr):
    # Amplitude of an LG mode
    return (((np.sqrt(2)*r)/(LGw(wNod,z,Zr)))**(np.abs(l))*genlaguerre(p,np.abs(l))((2*r**2)/LGw(wNod,z,Zr)**2)*np.exp((-r**2)/(LGw(wNod,z,Zr)**2)))
def LG_U(p,l,r,theta,z,k,wNod,Zr):
    #High Order Gaussian Beams in Cylindrical Coodinates
    if z == 0:
        return LGA(p,l)/LGw(wNod,z,Zr)*LGAmp(p,l,r,z,wNod,Zr)*np.exp(1j*l*theta)
    else:
        return (LGA(p,l)/LGw(wNod,z,Zr)*LGAmp(p,l,r,z,wNod,Zr)
        *np.exp((1j*k*r**2)/(2*LGR(z,Zr)))
        *np.exp(1j*l*theta)
        *np.exp(-1j*LGPhi(z,Zr)))

def LGBeamPh(p,l,r,theta,z,k,wNod,Zr,Lambda):
    return LG_U(p,l,r,theta,z,k,wNod,Zr)*np.exp(1j*(2*np.pi/Lambda)*z)

def LGBeamPhPo(Parameter, x, y, z, Rho, Theta):
#self.expression = [self.p,self.l,self.wNod, self.Zr, self.Lambda, self.Ar, self.n, self.k, self.z0, self.theta]
    one = polar.e_R*1/np.sqrt(2)*LGBeamPh(0,1,Rho,Theta,z,Parameter[7],Parameter[2],Parameter[3],Parameter[4])
    zero = polar.e_L*1/np.sqrt(2)*LGBeamPh(0,0,Rho,Theta,z,Parameter[7],Parameter[2],Parameter[3],Parameter[4])

    return one + zero

class GaussianBeam(bp.Beam):
    def __init__(self, WaveLength, ApertureRadius, RefractionIndex, p, l, zNod = 0, RayleighRange = 0,BeamWaist = 0):
        self.p = p
        self.l = l
        self.wNod = BeamWaist
        self.Zr = RayleighRange
        self.Lambda = WaveLength
        self.Ar = ApertureRadius
        self.n = RefractionIndex
        self.k = 2*np.pi/self.Lambda
        self.z0 = zNod
        if self.wNod == 0 & self.Zr == 0:
            raise Exception("Please define RayleighRange or Beam Waist")
        elif self.wNod != 0 & self.Zr == 0:
            self.Zr = np.pi*self.wNod**2*self.n/self.Lambda
        elif self.wNod == 0 & self.Zr != 0:
            self.wNod = np.sqrt(self.Zr * self.Lambda / np.pi)
        self.theta = self.Lambda/(np.pi * self.wNod)


        self.BasicParameter = [self.p,self.l,self.wNod, self.Zr, self.Lambda, self.Ar, self.n, self.k, self.z0, self.theta]
        self.expression = LGBeamPhPo

LGB1 = GaussianBeam(700e-9,0.001,1,0,0,BeamWaist=0.015)
#p and l are hard coded, so the constructring will have no impact
Screen = bp.Screen(3e-2,3e-2,z=0,resolution=256)
bp.Render(LGB1,Screen,0,0)
bp.Export(Screen, VF="ellipse")
