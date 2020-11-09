from vpython import *

t = 0
dt = 10**8

G_CONSTANT = 6.673E-11
mass = {'earth': 5.97E24, 'moon': 7.36E22, 'sun':1.99E30}
radius = {'earth': 6.371E6*10, 'moon': 1.317E6*10, 'sun':6.95E8*10}
earth_orbit = {'r': 1.495E11, 'v': 2.9783E4}
moon_orbit = {'r': 3.84E8, 'v': 1.022E3}
theta = 5.145*pi/180.0

my_canvas = canvas(width=700, height=400)
my_canvas.forward = vec(0, -1, 0)

class Solar_Sys_Obj(sphere):
    def __init__(self, name, mass, rad, velocity=vec(0, 0, 0)):
        self.m = mass
        self.v = velocity
        self.radius = rad
        self.grav_field = G_CONSTANT*self.m
        self.name = name
    
    def calc_pos(self, other_objs):
        for i in other_objs:
            if i is not self:
                self.v += dt * (i.grav_field) / (mag(self.pos-i.calc_pos))**2

objs = ["earth", "moon"]
solar_sys = dict.fromkeys(objs)
for i in solar_sys.keys():
    solar_sys[i] = Solar_Sys_Obj(i, mass[i], radius[i])
solar_sys["earth"].pos = vec(0, 0, 0)
solar_sys["moon"].pos = vec(moon_orbit["r"], 0, 0)
solar_sys["moon"].pos = vec(0, 0, moon_orbit["v"])


while True:
    t += dt
    
