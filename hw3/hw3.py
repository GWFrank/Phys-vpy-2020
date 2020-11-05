from vpython import *

G_CONSTANT = 6.673E-11
mass = {'earth': 5.97E24, 'moon': 7.36E22, 'sun':1.99E30}
radius = {'earth': 6.371E6*10, 'moon': 1.317E6*10, 'sun':6.95E8*10}
earth_orbit = {'r': 1.495E11, 'v': 2.9783E4}
moon_orbit = {'r': 3.84E8, 'v': 1.022E3}
theta = 5.145*pi/180.0


class Solar_Sys_Obj(sphere):
    def __init__(self, mass, velocity):
        self.m = mass
        self.v = velocity
        self.grav_field = G_CONSTANT*self.m
    
    def calc_pos(self, other_objs):
        for i in other_objs:
            if i is not self:
                self.v += (i.grav_field*self.m) / (mag(self.pos-i.calc_pos))**2

