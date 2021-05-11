import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from scipy.stats import norm
import matplotlib.mlab as mlab

kB=1.380649e-23

class particle:
  def __init__(self, position, velocity, S_ID):
    self.position = np.array(position)
    self.velocity = np.array(velocity)
    self.S_ID = S_ID


def readparticles(filename):
    with open(filename) as file:
        lines = [line.rstrip('\n') for line in file]

    parts = []
    
    for index, line in enumerate(lines):
        cols = line.split()
        if cols[0] != '%':
            parts.append(particle([np.float(cols[1]), np.float(cols[2]), np.float(cols[3])], [np.float(cols[4]), np.float(cols[5]), np.float(cols[6])], int(cols[9])))
    print('Read', len(parts), 'particles from file.')
    return parts


def plot_timestep(parts, vectors = False, vscale = 1e-4):
    fig = plt.figure()
    s1 = 0
    s2 = 0
    for part in parts:
        if part.S_ID == 1:
            s1 = s1 + 1
        else:
            s2 = s2 + 1
    colors = ['r' if part.S_ID == 2 else 'b' for part in parts]        
    plt.scatter([part.position[0] for part in parts], [part.position[1] for part in parts], color = colors, s=2.0)
    if vectors:
        for part in parts:
            plt.arrow(part.position[0], part.position[1], part.velocity[0]*vscale, part.velocity[1]*vscale)#, head_width=0.05, head_length=0.1)
    #plt.xlim([-0.12, 0.12])
    #plt.ylim([-0.12, 0.12])
    #plt.plot([-0.5, -0.5, 0.5, 0.5, -0.5], [-0.5, 0.5, 0.5, -0.5, -0.5], '--')
    plt.axis('equal')
    plt.show(block=False)
    

def plot_movie(timesteps, start, finish):
    fig = plt.figure()
    
    
    for parts in timesteps[start:finish]:
        plt.cla()
        colors = []
        s1 = 0
        s2 = 0
        colors = ['r' if part.S_ID == 2 else 'b' for part in parts]        
        plt.scatter([part.position[0] for part in parts], [part.position[1] for part in parts], color = colors, s=2.0)
        #plt.xlim([-0.12, 0.12])
        #plt.ylim([-0.12, 0.12])
        #plt.plot([-0.5, -0.5, 0.5, 0.5, -0.5], [-0.5, 0.5, 0.5, -0.5, -0.5], '--')
        plt.show(block=False)
        plt.axis('equal')
        plt.pause(0.01)

def plot_vdf(timestep):
    plt.figure()
    plt.subplot(211)
    datos = [parts.velocity[0] for parts in timestep]
    n, bins, patches = plt.hist(datos, bins=100, density=1)
    (mu, sigma) = norm.fit(datos)
    y = norm.pdf( bins, mu, sigma)
    l = plt.plot(bins, y, 'r--', linewidth=2)
    plt.xlabel('$v_x$')
    
    plt.subplot(212)
    datos = [parts.velocity[1] for parts in timestep]
    n, bins, patches = plt.hist(datos, bins=100, density=1)
    (mu, sigma) = norm.fit(datos)
    y = norm.pdf( bins, mu, sigma)
    l = plt.plot(bins, y, 'r--', linewidth=2)
    plt.xlabel('$v_y$')

    plt.tight_layout()
    plt.show(block=False)


def plot_vdf_2d(timestep):

    x = [parts.velocity[0] for parts in timestep]
    y = [parts.velocity[1] for parts in timestep]
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.005


    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, 0.2]
    rect_histy = [left + width + spacing, bottom, 0.2, height]

    # start with a rectangular Figure
    plt.figure(figsize=(8, 8))

    ax_scatter = plt.axes(rect_scatter)
    ax_scatter.tick_params(direction='in', top=True, right=True)
    ax_histx = plt.axes(rect_histx)
    ax_histx.tick_params(direction='in', labelbottom=False)
    ax_histy = plt.axes(rect_histy)
    ax_histy.tick_params(direction='in', labelleft=False)

    # the scatter plot:
    ax_scatter.plot(x, y, ls='', marker=',', color = 'k')
    ax_scatter.axis('equal')

    # now determine nice limits by hand:
    binwidth = 500
    lim = np.ceil(np.abs([x, y]).max() / binwidth) * binwidth
    ax_scatter.set_xlim((-lim, lim))
    ax_scatter.set_ylim((-lim, lim))

    bins = np.arange(-lim, lim + binwidth, binwidth)
    ax_histx.hist(x, bins=bins, histtype='step', color='k', density=1)
    ax_histy.hist(y, bins=bins, orientation='horizontal', histtype='step', color='k', density=1)

    ax_histx.set_xlim(ax_scatter.get_xlim())
    ax_histy.set_ylim(ax_scatter.get_ylim())

    plt.show(block=False)

def readfiles(firststep, laststep, every=1, nproc=1):
    timesteps = []
    for time in range(firststep, laststep, every):
        timestep=[]
        for proc in range(nproc):
            name = 'proc_{:05d}_time_{:08d}'.format(proc, time)
            fn = name
            timestep.extend(readparticles(fn))
        timesteps.append(timestep)
    return timesteps

def compute_Tx(timestep, m):
    datos = [parts.velocity[0] for parts in timestep]
    sigma = np.std(datos)
    return m/kB*sigma**2

def compute_Ty(timestep, m):
    datos = [parts.velocity[1] for parts in timestep]
    sigma = np.std(datos)
    return m/kB*sigma**2
    
def compute_Tz(timestep, m):
    datos = [parts.velocity[2] for parts in timestep]
    sigma = np.std(datos)
    return m/kB*sigma**2

def compute_ux(timestep):
    datos = [parts.velocity[0] for parts in timestep]
    vmean = np.mean(datos)
    return vmean

def plot_energy(timesteps, start, finish):
    totens = []
    for parts in timesteps[start:finish]:
        toten = 0
        for part in parts:
            toten += np.linalg.norm(part.velocity)**2
        totens.append(toten)
    plt.figure()
    plt.plot(totens)
    plt.show(block=False)

def subset_species(timesteps, sid):
    return [[particle for particle in timestep if (particle.S_ID == sid)] for timestep in timesteps]

def position_within(particle, bounds):
    return (particle.position[0] >= bounds[0] and particle.position[0] <= bounds[1] and
            particle.position[1] >= bounds[2] and particle.position[1] <= bounds[3] and
            particle.position[2] >= bounds[4] and particle.position[2] <= bounds[5])

def subset_configuration(timesteps, bounds):
    return [[particle for particle in timestep if position_within(particle, bounds)] for timestep in timesteps]

        
##########################################################################
##########################################################################
##########################################################################
##########################################################################
# End of the function definitions.
# Write here what you want the script to do.

# Example 1: load timesteps 0...99 run on only 1 processor, and plot the first time step:
timesteps = readfiles(0, 100, every=1, nproc = 1)
plot_timestep(timesteps[0], vectors = False)



### Example 2: compute x-velocity and temperatures only for particles
### of species 3, that are within the specified rectangle:
##timesteps = readfiles(200000, 200001, every=1, nproc = 16)
##onlye = subset_species(timesteps, 3)
##onlyewithin = subset_configuration(onlye, [0.89, 0.9, 0, 0.2, -3.14, 3.14])
##plot_vdf_2d(onlyewithin[0])
##print(compute_ux(onlyewithin[0]))
##print(compute_Tx(onlyewithin[0], me))
##print(compute_Ty(onlyewithin[0], me))
##print(compute_Tz(onlyewithin[0], me))

