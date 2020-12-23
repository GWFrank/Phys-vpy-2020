import numpy as np
from vpython import *

A, N, omega = 0.10, 50, 2*pi/1.0
size, m, k, d = 0.06, 0.1, 10.0, 0.4

def getAngFreq(n, periods):
    #==============================
    Unit_K = 2 * pi/(N*d)
    Wavevector = n * Unit_K
    phase = Wavevector * np.arange(N) * d

    ball_pos = np.arange(N)*d + A*np.sin(phase)
    ball_orig = np.arange(N)*d
    ball_v = np.zeros(N)
    spring_len = np.ones(N)*d
    #==============================

    pass_orig = False
    pass_times = -1

    t, dt = 0, 0.0003
    while pass_times <= 2*periods:
        t += dt
        
        spring_len[:-1] = ball_pos[1:] - ball_pos[:-1]
        spring_len[-1] = (ball_pos[0]+N*d) - ball_pos[-1]
        
        ball_v[1:] +=  ((k/m)*(spring_len[1:]-d) - (k/m)*(spring_len[:-1]-d))*dt 
        ball_v[0] += ((k/m)*(spring_len[0]-d) - (k/m)*(spring_len[-1]-d))*dt

        ball_pos +=  ball_v*dt
        ball_disp = ball_pos - ball_orig
        
        highest = max(ball_disp)
        if not pass_orig and highest <= 10**-2:
            pass_orig = True
        
        if pass_orig and highest >= 9*10**-2:
            pass_orig = False
            pass_times += 1
            if pass_times == 0:
                t = 0
    
    T = t/periods
    ang_freq = 2*pi/T
    return ang_freq

ang_freq_graph = graph(width = 800, height = 400, align = "left"
                       , title = "Angular Frequency - wavevector"
                       , xtitle = "wavevector", ytitle = "angular frequency"
                       , foreground = color.black, background = vec(0.05, 0.05, 0.05)
                       , xmin = 0, xmax = 25
                       )

ang_freq_curve = gcurve(graph = ang_freq_graph, color=color.cyan, width=4)

for i in range(1, N//2):
    ang_freq_curve.plot(i, getAngFreq(i, 5))
# ang_freqs = [getAngFreq(i, 1) for i in range(1, N//2)]



