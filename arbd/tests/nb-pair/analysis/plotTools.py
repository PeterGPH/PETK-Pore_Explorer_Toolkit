import numpy as np
import matplotlib as mpl
import matplotlib.ticker as plticker
import sys

def hsv2rgb(h,s,v):
    # color=hsv_to_rgb( np.reshape( np.array([1./3, 0.5, 1]), [1,1,3] ) )[0,0,:],
    hsv = np.reshape( [h,s,v], [1,1,3] )
    return mpl.colors.hsv_to_rgb( hsv )[0,0,:]

# c0  = hsv2rgb(1/3, 0.5, 0.9)
# c1a = hsv2rgb(1/2, 0.25, 0.9)
# c1b = hsv2rgb(1/2, 0.5, 0.7)
# c1c = hsv2rgb(1/2, 0.75, 0.5)
# c2  = hsv2rgb(2/3, 0.5, 0.9)
# c3  = hsv2rgb(1/6, 0.5, 0.9)
# c3  = hsv2rgb(1/2, 0.5, 0.9)


def hsv2rgb256(h,s,v):
    return hsv2rgb(h/360,s/255,v/255)

# http://mathematica.stackexchange.com/questions/20851/mathematica-color-schemes-for-the-colorblind
# c0  = np.array([51,   34, 136])/255
# c1a = np.array([136, 204, 238])/255
# c1b = np.array([68,  170, 153])/255
# c1c = np.array([17,  119,  51])/255
# c2  = np.array([153, 153,  51])/255
# # c3  = np.array([136,  34,  85])/255
# c3  = np.array([204, 121, 167])/255

# c0  = np.array([153,  34, 136])/255
# c1a = np.array([51,  102, 170])/255
# c1b = np.array([17,  170, 153])/255
# c1c = np.array([102, 170,  85])/255
# c2  = np.array([238,  51,  51])/255
# c3  = np.array([238, 119,  34])/255

# c0  = np.array([136,  46, 114])/255
# c1a = np.array([ 25, 101, 176])/255
# c1b = np.array([ 82, 137, 199])/255
# c1c = np.array([123, 175, 222])/255
# # c2  = np.array([144, 201, 135])/255
# # c2  = np.array([246, 193,  65])/255
# c2  = np.array([241, 147,  45])/255
# c3  = np.array([ 78, 178, 101])/255

# c0  = hsv2rgb256(314, 169, 136)
c0  = hsv2rgb256(280, 226, 94)
# c1a = hsv2rgb256(209, 219, 176)
# c1b = hsv2rgb256(211, 150, 199)
# c1c = hsv2rgb256(208, 114, 222)
c1a = hsv2rgb256(209, 219, 120)
c1b = hsv2rgb256(211, 150, 180)
c1c = hsv2rgb256(208, 114, 230)
c2  = hsv2rgb256( 31, 208, 241)
c3  = hsv2rgb256(133, 143, 178)
c1  = c1b

m0 = 's'
m1 = '<'
m2 = 'o'
m3 = '>'

# ----------------------------------------------------------------- #
# Make colormap based on Paul Tol's best visibility gradients. See  #
# <http://www.sron.nl/~pault/> for more info on these colors. Also  #
# see <http://matplotlib.sourceforge.net/api/colors_api.html> and   #
# <http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps> on some #
# matplotlib examples                                               #
# ----------------------------------------------------------------- #
# Deviation around zero colormap (blue--red)
cols = []
for x in np.linspace(0,1, 256):
    rcol = 0.237 - 2.13*x + 26.92*x**2 - 65.5*x**3 + 63.5*x**4 - 22.36*x**5
    gcol = ((0.572 + 1.524*x - 1.811*x**2)/(1 - 0.291*x + 0.1574*x**2))**2
    bcol = 1/(1.579 - 4.03*x + 12.92*x**2 - 31.4*x**3 + 48.6*x**4 - 23.36*x**5)
    cols.append((rcol, gcol, bcol))
cm_plusmin = mpl.colors.LinearSegmentedColormap.from_list("PaulT_plusmin", cols)

# Linear colormap (white--red)
from scipy.special import erf
cols = []
for x in np.linspace(0,1, 256):
    rcol = (1 - 0.392*(1 + erf((x - 0.869)/ 0.255)))
    gcol = (1.021 - 0.456*(1 + erf((x - 0.527)/ 0.376)))
    bcol = (1 - 0.493*(1 + erf((x - 0.272)/ 0.309)))
    cols.append((rcol, gcol, bcol))
cm_linear = mpl.colors.LinearSegmentedColormap.from_list("PaulT_linear", cols)

# Linear colormap (rainbow)
cols = [(0,0,0)]
for x in np.linspace(0,1, 254):
    rcol = (0.472-0.567*x+4.05*x**2)/(1.+8.72*x-19.17*x**2+14.1*x**3)
    gcol = 0.108932-1.22635*x+27.284*x**2-98.577*x**3+163.3*x**4-131.395*x**5+40.634*x**6
    bcol = 1./(1.97+3.54*x-68.5*x**2+243*x**3-297*x**4+125*x**5)
    cols.append((rcol, gcol, bcol))
cols.append((1,1,1))
cm_rainbow = mpl.colors.LinearSegmentedColormap.from_list("PaulT_rainbow", cols)
cm_default = cm_plusmin


def info(*obj):
    print('INFO: ',*obj , file=sys.stderr)

def blkavg(a,blklength):
    """ 
    a is (n*N+extra)x2 array
    blklength (units of x-axis)
    returns Nx2
    """
    t = a[:,0]
    dt = t[1]-t[0]
    y = a[:,1]

    # n is number of elements in a block, N is the number of blocks
    n = int(np.floor(blklength/dt))
    N = int(np.floor(len(t)/n))
    tb = np.reshape(t[:n*N],(N,n)).mean(axis=1)
    yb = np.reshape(y[:n*N],(N,n)).mean(axis=1)
    r = np.vstack((tb,yb)).transpose()
    # v = np.gradient(yb,dt*n)
    # s = v.std()/np.sqrt(N)
    # return r,v,s
    return r

# def blockAvg(data, pointsPerBlock=0, nblocks=100): 
#     ## http://stackoverflow.com/questions/14229029/block-mean-of-numpy-2d-array
#     if pointsPerBlock > 0:
#         nblocks = len(data)/pointsPerBlock
    
#     l = data.shape[0]
#     pointsPerBlock = int(l/nblocks)
    
#     data = data[
#     data.reshape( (l*nblocks)
    
#     info(nblocks)
#     info(data.shape)

#     blocks = data.reshape((data.shape[0], -1, nblocks))
#     info(blocks)
#     blockaverage = np.mean(blocks, axis=1) 
#     # blockstd = np.std(blocks]) 
#     # return blockaverage, blockstd
#     return blockaverage

def getTickDivisions(xmin,xmax):
    dx = np.abs(xmax-xmin)

    minticksize = dx / 2
    magnitude = np.log(minticksize)/np.log(10)
    tmp = magnitude
    magnitude = 10**np.floor(magnitude)
    residual = minticksize / magnitude 

    if np.abs(np.round(residual*2) - residual*2) < 1e-6 and magnitude > 0.01:
        tick = residual * magnitude
    elif residual > 5:
        tick = 10*magnitude
    elif residual > 2:
        tick =  5*magnitude
    elif residual > 1:
        tick =  2*magnitude
    else:
        tick =  magnitude
    
    # info('range:',dx)
    # info('minticksize:',minticksize)
    # info('tmp:',tmp)
    # info('magnitude:',magnitude)
    # info('residual:',residual)
    # info('tick:',tick)

    # info(dx,magnitude,residual,tick)
    return tick
# for f in np.linspace(0,1,11)[1:]:
#     getTickDivisions(0,f)

def setRange(axes, minmax):
    xmin,xmax,ymin,ymax = minmax
    axes.axis(minmax)

    dx = getTickDivisions(xmin,xmax)
    dy = getTickDivisions(ymin,ymax)

    setTicks(axes,dx,dy)


# def setTicks(axes, dx, dy, numXminor=1, numYminor=1, xoffset=0, yoffset=0):
def setTicks(axes, dx, dy, numXminor=1, numYminor=1):
    axes.xaxis.set_major_locator( plticker.MultipleLocator( base=dx ))
    axes.xaxis.set_minor_locator( plticker.MultipleLocator( base=dx/(numXminor+1) ))
    axes.yaxis.set_major_locator( plticker.MultipleLocator( base=dy ))
    axes.yaxis.set_minor_locator( plticker.MultipleLocator( base=dy/(numXminor+1) ))

