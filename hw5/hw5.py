from vpython import *
from diatomic import *


# 20 molecules
N = 20
# 2L is the length of the cubic container box, the number is made up
L = ((24.4E-3/(6E23))*N)**(1/3.0)/50
# average mass of O and C
m = 14E-3/6E23
# some constants to set up the initial speed
k, T = 1.38E-23, 298.0
# some constant
initial_v = (3*k*T/m)**0.5


scene = canvas(width=400, height=400, align='left',
               background=vec(0.1, 0.1, 0.1))
container = box(length=2*L, height=2*L, width=2*L
                , opacity=0.4, color=color.yellow)

energies = graph(width=600, align='left', ymin=0)
c_avg_com_K = gcurve(color=color.green)
c_avg_v_P = gcurve(color=color.red)
c_avg_v_K = gcurve(color=color.purple)
c_avg_r_K = gcurve(color=color.blue)


COs = []

for i in range(N):
    # initialize the 20 CO molecules

    # random() yields a random number between 0 and 1
    O_pos = vec(random()-0.5, random()-0.5, random()-0.5)*L
    # generate one CO molecule
    CO = CO_molecule(pos=O_pos, axis=vector(1.0*d, 0, 0))
    # set up the initial velocity of C randomly
    CO.C.v = vector(initial_v*random(), initial_v*random(), initial_v*random())
    # set up the initial velocity of O randomly
    CO.O.v = vector(initial_v*random(), initial_v*random(), initial_v*random())

    COs.append(CO)


# number of loops that has been run
times = 0

dt = 5E-16
t = 0

avg_com_K, avg_v_K, avg_v_P, avg_r_K = 0, 0, 0, 0

while True:
    rate(3000)

    for CO in COs:
        CO.time_lapse(dt)

    for i in range(N-1):
        # the first N-1 molecules
        for j in range(i+1, N):
            # from i+1 to the last molecules, to avoid double checking
            ## change this to check and handle the collisions between
            ## the atoms of different molecules
            atoms_i = [COs[i].C, COs[i].O]
            atoms_j = [COs[j].C, COs[j].O]
            for a_i in atoms_i:
                for a_j in atoms_j:
                    if (mag(a_i.pos-a_j.pos) <= 2*size
                            and dot(a_i.v-a_j.v, a_i.pos-a_j.pos) < 0):
                        a_i.v, a_j.v = collision(a_i, a_j)
            
    for CO in COs:
        ## change this to check and handle the
        ## collision of the atoms of all molecules on all 6 walls
        atoms = [CO.C, CO.O]
        for a in atoms:
            if ((a.pos.x+size >= L and a.v.x > 0) 
                    or (a.pos.x-size <= -L and a.v.x < 0)):
                a.v.x = -a.v.x
            
            if ((a.pos.y+size >= L and a.v.y > 0) 
                    or (a.pos.y-size <= -L and a.v.y < 0)):
                a.v.y = -a.v.y
            
            if ((a.pos.z+size >= L and a.v.z > 0) 
                    or (a.pos.z-size <= -L and a.v.z < 0)):
                a.v.z = -a.v.z
    
    ## sum com_K, v_K, v_P, and r_Kforall molecules, respectively, to get
    ## total_com_K, total_v_K, total_v_P, total_r_K at the current moment
    total_com_K, total_v_K, total_v_P, total_r_K = 0, 0, 0, 0
    for CO in COs:
        total_com_K += CO.com_K()
        total_v_K += CO.v_K()
        total_v_P += CO.v_P()
        total_r_K += CO.r_K()

    ## calculate avg_com_K to be the time average of total_com_K
    ## since the beginning of the simulation,
    ## and do the same for others.
    avg_com_K = (avg_com_K*t + total_com_K*dt) / (t + dt)
    avg_v_K = (avg_v_K*t + total_v_K*dt) / (t + dt)
    avg_v_P = (avg_v_P*t + total_v_P*dt) / (t + dt)
    avg_r_K = (avg_r_K*t + total_r_K*dt) / (t + dt)

    t += dt
    times += 1

    ## plot avg_com_K, avg_v_K, avg_v_P, and avg_r_K
    # c_avg_com_K.plot(t, avg_com_K)
    # c_avg_v_K.plot(t, avg_v_K)
    # c_avg_v_P.plot(t, avg_v_P)
    # c_avg_r_K.plot(t, avg_r_K)
    
    c_avg_com_K.plot(times, avg_com_K)
    c_avg_v_K.plot(times, avg_v_K)
    c_avg_v_P.plot(times, avg_v_P)
    c_avg_r_K.plot(times, avg_r_K)
