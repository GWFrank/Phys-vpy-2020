from vpython import *
import numpy as np

from histogram import *


N = 200
m, size = 4E-3/6E23, 31E-12*10
L = ((24.4E-3/(6E23))*N)**(1/3.0)/2 + size
k, T = 1.38E-23, 298.0

t, dt = 0, 3E-13

vrms = (2*k*1.5*T/m)**0.5
stage = 0
atoms = []


deltav = 50.
vdist = graph(x=800, y=0, ymax=N*deltav/1000., width=500,
              height=300, xtitle='v', ytitle='dN', align='left')

theory_low_T = gcurve(color=color.cyan)
dv = 10.

for v in arange(0., 4201.+dv, dv):
    theory_low_T.plot(pos=(v, (deltav/dv)*N*4.*pi *
                           ((m/(2.*pi*k*T))**1.5)*exp((-0.5*m*v**2)/(k*T))*(v**2)*dv))
observation = ghistogram(graph=vdist, bins=arange(
    0., 4200., deltav), color=color.red)


scene = canvas(width=500, height=500,
               background=vector(0.2, 0.2, 0), align='left')
container = box(length=2*L, height=2*L, width=2 *
                L, opacity=0.2, color=color.yellow)
p_a, v_a = np.zeros((N, 3)), np.zeros((N, 3))
for i in range(N):
    p_a[i] = [2 * L*random() - L, 2 * L*random() - L, 2 * L*random() - L]
    if i == N-1:
        atom = sphere(pos=vector(p_a[i, 0], p_a[i, 1], p_a[i, 2]),
                      radius=size, color=color.yellow, make_trail=True, retain=50)
    else:
        atom = sphere(pos=vector(p_a[i, 0], p_a[i, 1], p_a[i, 2]), radius=size, color=vector(
            random(), random(), random()))
    ra = pi*random()
    rb = 2*pi*random()
    v_a[i] = [vrms*sin(ra)*cos(rb), vrms*sin(ra)*sin(rb), vrms*cos(ra)]
    atoms.append(atom)


def vcollision(a1p, a2p, a1v, a2v):
    v1prime = a1v - (a1p - a2p) * sum((a1v-a2v)*(a1p-a2p)) / sum((a1p-a2p)**2)
    v2prime = a2v - (a2p - a1p) * sum((a2v-a1v)*(a2p-a1p)) / sum((a2p-a1p)**2)
    return v1prime, v2prime

def chgmod(evt):
    # print("press")
    global stage, scene, observation, theory_low_T
    k = evt.key
    if k == "n" and stage <= 2:
        stage += 1
    

scene.bind('keydown', chgmod)

momentum = 0
gamma = 5/3
drawn = False
v_W = 0
while True:
    t += dt
    rate(10000)

    p_a += v_a*dt
    for i in range(N):
        atoms[i].pos = vector(p_a[i, 0], p_a[i, 1], p_a[i, 2])
    if stage != 1:
        observation.plot(data=np.sqrt(np.sum(np.square(v_a), -1)))
    
    # collision between atoms
    r_array = p_a-p_a[:, np.newaxis]
    rmag = np.sqrt(np.sum(np.square(r_array), -1))
    hit = np.less_equal(rmag, 2*size)-np.identity(N)
    hitlist = np.sort(np.nonzero(hit.flat)[0]).tolist()
    for ij in hitlist:
        i, j = divmod(ij, N)
        hitlist.remove(j*N+i)
        if sum((p_a[i]-p_a[j])*(v_a[i]-v_a[j])) < 0:
            v_a[i], v_a[j] = vcollision(p_a[i], p_a[j], v_a[i], v_a[j])
    
    # collision between atoms and walls
    for i in range(N):
        if abs(p_a[i][0]) >= container.length/2 - size and p_a[i][0]*v_a[i][0] > 0:
            delt_v = abs(v_a[i][0]) + (v_W)
            v_a[i][0] = -np.sign(p_a[i][0])*(v_W + delt_v)
            momentum += 2*abs(delt_v*m)
        if abs(p_a[i][1]) >= container.height/2 - size and p_a[i][1]*v_a[i][1] > 0:
            v_a[i][1] = - v_a[i][1]
            momentum += 2*abs(v_a[i][1]*m)
        if abs(p_a[i][2]) >= container.width/2 - size and p_a[i][2]*v_a[i][2] > 0:
            v_a[i][2] = - v_a[i][2]
            momentum += 2*abs(v_a[i][2]*m)


    if stage == 1:
        v_W = L/(20000.0*dt)
        container.length -= v_W*2*dt
        if container.length <= L:
            stage += 1

    if stage == 2 and not drawn:
        drawn = True
        v_W = 0
        observation = ghistogram(graph=vdist, bins=arange(
                                 0., 4200., deltav), color=color.blue)
        theory_low_T = gcurve(color=color.cyan)
        dv = 10.
        for v in arange(0., 4201.+dv, dv):
            theory_low_T.plot(pos=(v, (deltav/dv)*N*4.*pi *
                                ((m/(2.*pi*k*T))**1.5)*exp((-0.5*m*v**2)/(k*T))*(v**2)*dv))

    if stage == 3:
        container.length = 2*L

    if (t//dt)%1000 == 0:
        KE = 0.5*m*np.sum(np.square(v_a))
        T = KE/((3/2)*N*k)
        V = container.height*container.width*container.length
        surface_area = 2*(container.height*container.width
                          + container.height*container.length
                          + container.length*container.width)
        p = momentum/(1000*dt*surface_area)
        print("============= Gas Data Report =============")
        print(f"Stage       : {stage}")
        print(f"Temperature : {T}")
        print(f"Pressure    : {p}")
        print(f"Volume      : {V}")
        print(f"P*V         : {p*V}")
        print(f"N*k*T       : {N*k*T}")
        print(f"P*V^gamma   : {p*(V**gamma)}")
        print("===========================================")
        momentum = 0

