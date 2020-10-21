from vpython import *

t = 0
dt = 0.0001
k = 150000
g = vector(0, -9.8, 0)
N = 2
mass = 1
radius = 0.2
rope_len = 2

zero_grav_pot_h = -0.5*rope_len
my_grey = vec(0.6, 0.6, 0.6)

ke_ug_graph = graph(width=400, height=400, align="right",
                    title="<i>U<sub>g</sub></i> & <i>K.E.</i> - <i>t</i>", xtitle="<i>t</i>",
                    foreground=color.black, background=vec(0.05, 0.05, 0.05),
                    ymin=0, ymax=1.2)
ke_curve = gcurve(graph=ke_ug_graph, color=color.blue,
                  width=4, label="<i>K.E.</i>")
ug_curve = gcurve(graph=ke_ug_graph, color=color.green,
                  width=4, label="<i>U</i><sub>g</sub>")

avg_ke_ug_graph = graph(width=400, height=400, align="right",
                        title="average <i>U</i><sub>g</sub></i> & <i>K.E.</i> - <i>t</i>", xtitle="<i>t</i>",
                        foreground=color.black, background=vec(0.05, 0.05, 0.05),
                        ymin=0, ymax=1.2)
avg_ke_curve = gcurve(graph=avg_ke_ug_graph, color=color.blue,
                      width=4, label="average <i>K.E.</i>")
avg_ug_curve = gcurve(graph=avg_ke_ug_graph, color=color.green,
                      width=4, label="average <i>U</i><sub>g</sub>")


my_canvas = canvas(width=500, height=400, align="left")


class Pendulum():
    def __init__(self, x_cord, ball_rad, mass, rope_len):
        self.pivot = sphere(pos=vec(x_cord, 0.5*rope_len, 0), radius=0.2*ball_rad,
                            color=my_grey, texture=textures.metal)
        self.ball = sphere(pos=vec(x_cord, -0.5*rope_len, 0), radius=ball_rad,
                           color=my_grey, texture=textures.metal)
        self.rope = cylinder(pos=self.pivot.pos, axis=self.ball.pos - self.pivot.pos, radius=0.05*ball_rad,
                             color=my_grey, texture=textures.metal)
        self.mass = mass
        self.rope_len = rope_len
        self.velocity = vector(0, 0, 0)

    def cal_pos(self, dt):
        self.velocity += vec(g + (-k)*(mag(self.rope.axis) -
                                       rope_len)*norm(self.rope.axis))*dt
        self.ball.pos += self.velocity*dt
        self.rope.axis = self.ball.pos - self.pivot.pos

    def pull_pend(self, height):
        self.ball.pos += vec(-sqrt(self.rope_len**2 -
                                   (self.rope_len-height)**2), height, 0)
        self.rope.axis = self.ball.pos - self.pivot.pos
        # print(mag(self.rope.axis))


def collide(pend1, pend2):
    if mag(pend1.ball.pos - pend2.ball.pos) < (pend1.ball.radius + pend2.ball.radius):
        collide_perp = norm(pend1.ball.pos - pend2.ball.pos)
        v1_perp = proj(pend1.velocity, collide_perp)
        v2_perp = proj(pend2.velocity, collide_perp)
        pend1.velocity += -v1_perp + v2_perp
        pend2.velocity += -v2_perp + v1_perp


cradle = [Pendulum(2*radius*(i-2), radius, mass, rope_len) for i in range(5)]

for i in range(N):
    cradle[i].pull_pend(0.05)

avg_ke = 0
avg_ug = 0
for i in cradle:
    avg_ke += 0.5*i.mass*(mag(i.velocity)**2)
    avg_ug += i.mass*mag(g)*(i.ball.pos.y - zero_grav_pot_h)

while True:
    rate(5000)
    
    t += dt
    for i in cradle:
        i.cal_pos(dt)
    
    for i in range(5-1):
        collide(cradle[i], cradle[i+1])

    cradle_ke = 0
    cradle_ug = 0
    for i in cradle:
        cradle_ke += 0.5*i.mass*(mag(i.velocity)**2)
        cradle_ug += i.mass*mag(g)*(i.ball.pos.y - zero_grav_pot_h)
    avg_ke = (t*avg_ke + dt*cradle_ke)/(t+dt)
    avg_ug = (t*avg_ug + dt*cradle_ug)/(t+dt)

    ke_curve.plot(t, cradle_ke)
    ug_curve.plot(t, cradle_ug)
    avg_ke_curve.plot(t, avg_ke)
    avg_ug_curve.plot(t, avg_ug)
