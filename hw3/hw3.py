from vpython import *

t = 0
dt = 6*10**2

G_CONSTANT = 6.673E-11
mass = {'earth': 5.97E24, 'moon': 7.36E22, 'sun':1.99E30}
radius = {'earth': 6.371E6*10, 'moon': 1.317E6*10, 'sun':6.95E8*10}
earth_orbit = {'r': 1.495E11, 'v': 2.9783E4}
moon_orbit = {'r': 3.84E8, 'v': 1.022E3}
theta = 5.145*pi/180.0

my_canvas = canvas(width=500, height=500, align="left")

class Solar_Sys_Obj(sphere):
    def __init__(self, name, mass, rad, velocity=vec(0, 0, 0)):
        super().__init__(radius=rad)
        self.m = mass
        self.v = velocity
        # self.radius = rad
        self.grav_field = G_CONSTANT*self.m
        self.name = name
    
    def calc_pos(self, other_objs):
        for i in other_objs:
            if i is not self:
                self.v += (i.grav_field) / (mag(self.pos-i.pos))**2\
                           *norm(i.pos-self.pos)*dt
        self.pos += self.v*dt


objs = ["earth", "moon", "sun"]
# objs = ["earth", "moon"]
solar_sys = dict.fromkeys(objs)
for i in solar_sys.keys():
    solar_sys[i] = Solar_Sys_Obj(i, mass[i], radius[i])

EM_CONST = mass["moon"] / (mass["moon"]+mass["earth"])
EM_POS =  moon_orbit["r"]*EM_CONST
EM_V = moon_orbit["v"]*EM_CONST

solar_sys["earth"].pos = vec((-EM_POS)*cos(-theta) + earth_orbit["r"],
                             (-EM_POS)*sin(-theta),
                             0)
solar_sys["moon"].pos = vec((moon_orbit["r"]-EM_POS)*cos(-theta) + earth_orbit["r"],
                            (moon_orbit["r"]-EM_POS)*sin(-theta),
                            0)
solar_sys["sun"].pos = vec(0, 0, 0)

solar_sys["earth"].v = vec(0, 0, -earth_orbit["v"]+EM_V)
solar_sys["moon"].v = vec(0, 0, -moon_orbit["v"]-earth_orbit["v"]+EM_V)

solar_sys["earth"].texture = textures.earth
solar_sys["moon"].texture = textures.rough
solar_sys["sun"].color = color.yellow

my_canvas.center = solar_sys["earth"].pos
my_canvas.range = 2.2*earth_orbit["r"]

moon_orbit_vec = arrow(pos = solar_sys["moon"].pos,
                       axis = (10**8)*norm(cross(solar_sys["moon"].pos-solar_sys["earth"].pos, solar_sys["moon"].v-solar_sys["earth"].v)),
                       color = color.gray(0.3)
                       )

orbit_vec_canvas = canvas(width=500, height=500, align="right")
moon_orbit_vec_iso = arrow(axis = norm(moon_orbit_vec.axis),
                           color = color.gray(0.3))
precession_trail = sphere(pos = moon_orbit_vec_iso.pos+moon_orbit_vec_iso.axis,
                          make_trail = True, trail_radius = 0.002,
                          trail_type = "curve", interval = 100, opacity = 0)
orbit_vec_canvas.forward = vec(0, -1, 0)
orbit_vec_canvas.range = 1


checkpoint1 = False
checkpoint2 = False
while True:
    t += dt
    for i in solar_sys.keys():
        solar_sys[i].calc_pos(solar_sys.values())
    
    my_canvas.center = solar_sys["earth"].pos
    
    moon_orbit_vec.pos = solar_sys["moon"].pos
    moon_orbit_vec.axis = (10**8)*norm(cross((solar_sys["moon"].pos)-solar_sys["earth"].pos, (solar_sys["moon"].v - solar_sys["earth"].v)))
    moon_orbit_vec_iso.axis = norm(moon_orbit_vec.axis)
    
    precession_trail.pos = moon_orbit_vec_iso.pos+moon_orbit_vec_iso.axis
    
    cot = abs(precession_trail.pos.x/precession_trail.pos.z)

    if cot < 0.001 and precession_trail.pos.z < 0:
        checkpoint1 = True
    if cot < 0.001 and precession_trail.pos.z > 0:
        checkpoint2 = True
    if checkpoint1 and checkpoint2 and cot >= 10**5:
        period_year = t / (60*60*24*365.24)
        print(f"precession period : {period_year:.2f} years")
        break