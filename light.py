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
        return np.sqrt(self.x**2 + self.y**2)

    def __repr__(self):
        return "[" + str(self.x) +", " + str(self.y) + "]"


class Screen:
    def __init__(self, x_range, y_range, resolution = 256, z = 0):
        self.x_range = x_range
        self.y_range = y_range
        self.resolution = resolution
        self.screen = [[vector(0,0)]*x_range]*y_range
        self.x = np.linspace(-x_range,x_range,resolution)
        self.y = np.linspace(-y_range,y_range,resolution)
        self.z = z


a = Screen(2,2)
print (a.screen)

e_x = vector(1,0)
e_y = vector(0,1)

def e_phi(Phi):
    return e_x * np.cos(Phi) + e_y * np.sin(Phi)

e_D = (e_x + e_y)/np.sqrt(2)
e_A = (e_x * (-1) + e_y)/np.sqrt(2)
e_R = (e_x - e_y * 1j)/np.sqrt(2)

print (e_R)