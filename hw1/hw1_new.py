# Please note that this code is using f-string, which is a feature only available in
# Python 3.6 or later versions
# For the best user experience, please make sure your monitor is at least 1300 px wide 
# and 400 px high
from vpython import *

playspeed = 1
fresh_rate = 2*10**3
t = 0
dt = 1/fresh_rate

arrow_tweak = 0.4
bounce_count = 0
bounce_limit = 3
travel_dist = 0.0
displacement = 0.0
max_h = 0.0

g = 9.8
ball_size = 0.25
C_drag = 0.9
theta = pi/4
grav = vec(0, -g, 0)

v_graph = graph(width = 300, height = 400, align = "right", title = "<i>s</i><sub>ball</sub> - <i>t</i>"
                , xtitle = "<i>t</i>"
                , foreground = color.black, background = vec(0.05, 0.05, 0.05)
                #, xmax = 5.0
                , ymin = 0.0, ymax = 30.0)
v_over_t = gcurve(graph = v_graph, color=color.cyan, width=4, label = "speed of the ball")

x_graph = graph(width = 300, height = 400, align = "right", title = "distance & displacement - <i>t</i>"
                , xtitle = "<i>t</i>"
                , foreground = color.black, background = vec(0.05, 0.05, 0.05)
                # , xmax = 5.0
                , ymin = 0.0, ymax = 30.0)
x_over_t = gcurve(graph = x_graph, color=color.blue, width=4, label = "displacement")
d_over_t = gcurve(graph = x_graph, color=color.green, width=4, label = "distance traveled")


my_canvas = canvas(width = 700, height = 400, align = "left", background = vec(0.05, 0.05, 0.05))
my_canvas.camera.pos = vector(0,1,1)


ground = box(length = 50, width = 7.5, height = 0.2, pos = vec(0,-0.1,0), texture = textures.wood_old)

ball = sphere(pos = vec(-15, ball_size, 0), radius = ball_size, color = color.cyan, make_trail = True)
ball_v = vec(20*cos(theta), 20*sin(theta), 0)

v_arrow = arrow(pos = ball.pos, axis = ball_v*arrow_tweak, shaftwidth = 0.1)


while True:
    rate(playspeed*fresh_rate)
    
    ball.pos += ball_v*dt
    ball_v += (grav + -C_drag*ball_v)*dt
    travel_dist += mag(ball_v)*dt
    displacement = mag(ball.pos - vec(-15, ball_size, 0))
    t += dt

    if ball.pos.y > max_h:
        max_h = ball.pos.y

    v_arrow.pos = ball.pos
    v_arrow.axis = ball_v*arrow_tweak

    v_over_t.plot(pos = (t, mag(ball_v)))
    x_over_t.plot(pos = (t, displacement))
    d_over_t.plot(pos = (t, travel_dist))

    if ball.pos.y <= ball_size:
        ball_v.y *= -1
        bounce_count += 1
        if bounce_count >= bounce_limit:
            displacement = round(displacement, 2)
            travel_dist = round(travel_dist, 2)
            max_h = round(max_h, 2)
            break

msg = text(text = f"displacement : {displacement:.2f}\ndistance traveled : {travel_dist:.2f}\nlargest height : {max_h:.2f}"
        , depth = 0.1, height = 1.4, pos = vec(0, -4, 0), align = "center")
