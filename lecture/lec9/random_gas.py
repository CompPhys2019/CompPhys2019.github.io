#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#####################################################################
### This very simple script aims in emulating the ###################
###  expansion of a gas in a random process. I of the ###############
#####################################################################

################################################################################
##### Drawing setup ############################################################
################################################################################
fig, axs = plt.subplots()
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

axs.set_aspect('equal')
plt.axis('off')

### parameters ######################################################
L = 400  ## Box linear size
N = 300 ## Number of particles
max_time = N * 2

### This routine draws the box ############################
def draw_initial():
    
    linestyle = '-'
    linewidth = 2.0
    line_color = 'dodgerblue'

    axs.plot([0.0, 2*L], [0.0, 0.0], linestyle = linestyle, \
    linewidth = linewidth, color = line_color)
    axs.plot([0.0, 2*L], [L, L], linestyle = linestyle, \
    linewidth = linewidth, color = line_color)

    axs.plot([0.0, 0.0], [0, L], linestyle = linestyle, \
    linewidth = linewidth, color = line_color)
    axs.plot([2*L, 2*L], [0, L], linestyle = linestyle, \
    linewidth = linewidth, color = line_color)

    axs.plot([L, L], [0, L/2. - 0.05*L], linestyle = linestyle, \
    linewidth = linewidth, color = line_color)

    axs.plot([L, L], [L/2. + 0.05*L, L], linestyle = linestyle, \
    linewidth = linewidth, color = line_color)
    
# === This routine initialize a random configuration in the left side ==== #
def initialize(N, L):
    
    x = L*np.random.rand(N)    ### N random numbers in [0,L)
    y = L*np.random.rand(N)
    
    nleft = N
    
    return x, y, nleft

# === This routine changes particle positions in the box ================= #
def step(N, L):
    
    global nleft
    
    rnd_particle = np.random.randint(N)
    
    if (x[rnd_particle] < L):
        nleft-=1   ### move to right
        x[rnd_particle] = L*np.random.rand() + L   ### random number in [L,2L)
        y[rnd_particle] = L*np.random.rand()
    else:
        nleft+=1   ### move to right
        x[rnd_particle] = L*np.random.rand()   ### random number in [0,L)
        y[rnd_particle] = L*np.random.rand()
        
    return x, y, nleft

### Here starts the fun ############
draw_initial()
time_range = range(0, max_time)

x, y, nleft = initialize(N, L)

particles, = axs.plot(x, y, linestyle = 'none', marker = 'o', color = 'crimson', markersize = 4)

text_template = r'$n_{\rm left}/N = %.2f$'
nleft_text = axs.text(0.50,1.1, '', \
    transform=axs.transAxes, fontsize=16, fontweight='bold', color = 'black',\
        horizontalalignment='center',)

def init():
    """initialize animation"""
    
    global nleft

    nleft_text.set_text(r'$n_{\rm left}/N = %.2f$' %(float(nleft)/N))
    
    return particles, nleft_text
    
def animate(time_range):
    
    #print x, y
    
    x, y, nleft = step(N, L)
    
    particles.set_data(x, y)
    
    nleft_text.set_text(r'$n_{\rm left}/N = %.2f$' %(float(nleft)/N))
    
    return particles, nleft_text

#### Updating animation
anim = animation.FuncAnimation(fig, animate, time_range,\
    interval=20, repeat=False, blit=False, init_func=init)

## Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=18000)

anim.save('N%dL%d.mp4' %(N, L), fps=15, dpi=200)

#plt.show()
