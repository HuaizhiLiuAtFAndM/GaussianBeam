import matplotlib.pyplot as plt
import numpy as np
class vector:
    def __init__(self,x=0,y=0):
        self.x = x
        self.y = y

    def __add__(self, v):
        x = self.x + v.x
        y = self.y + v.y
        return vector(x,y)

    def __sub__(self, v):
        x = self.x - v.x
        y = self.y - v.y
        return vector(x,y)

    def __mul__(self, c):
        x = self.x * c
        y = self.y * c
        return vector(x,y)


    def __truediv__(self, c):
        x = self.x / c
        y = self.y / c
        return vector(x,y)

    def Angle(self):
        if self.x >= 0:
            return np.arctan(self.y/self.x)
        else:
            return np.arctan(self.y/self.x) + np.pi

    def Mag(self):
        abs_x = np.abs(self.x)
        abs_y = np.abs(self.y)
        return np.sqrt(abs_x**2 + abs_y**2)

    def __repr__(self):
        return "{" + str(self.x) +", " + str(self.y) + "}"

e_x = vector(1,0)
e_y = vector(0,1)

e_D = (e_x + e_y)/np.sqrt(2)
e_A = (e_x * (-1) + e_y)/np.sqrt(2)
e_R = (e_x - e_y * 1j)/np.sqrt(2)
e_L = (e_x + e_y * 1j)/np.sqrt(2)

def e_phi(Phi):
    return e_x * np.cos(Phi) + e_y * np.sin(Phi)

def e_theta_chi(Theta, Chi):
    return e_R * np.exp(1j * Theta) * np.cos(Chi) + e_L * np.exp(-1j * Theta) * np.sin(Chi)