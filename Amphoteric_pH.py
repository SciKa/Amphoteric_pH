import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sympy as syp
from mayavi import mlab
c = 0.05

x = syp.symbols("x")


Ka1_xspace = np.linspace(0,10**-1,10000)
Ka2_yspace = np.linspace(0,10**-1,10000)
Kw_zspace = np.linspace(0,10**-1,10000)

#Ka1_xspace = np.linspace(0,10**-5,10000)
#Ka2_yspace = np.linspace(0,10**-5,10000)
#Kw_zspace = np.linspace(0,10**-5,10000)

#choose one up there.

size = 10000
Ka1_xspace_sap = np.random.choice(Ka1_xspace,size)
Ka2_yspace_sap = np.random.choice(Ka2_yspace,size)
Kw_zspace_sap = np.random.choice(Kw_zspace,size)
delta_dspace = np.array([])

for i in range(size):
    Ka1 = Ka1_xspace_sap[i]
    Ka2 = Ka2_yspace_sap[i]
    Kw = Kw_zspace_sap[i]
    precise_eq = syp.Eq(0, c * x ** 3 / (x + Ka2) + Ka1 * x ** 2 - Ka1 * Ka2 * c * x / (x + Ka2) - Kw * Ka1)
    p_result = syp.solve(precise_eq,x)
    p_result = list(map(lambda x : float(x.coeff('I',0)),p_result))
    tempres = None
    for res in p_result:
        if res > 0 and tempres == None:
            tempres = res
        elif res > 0 and tempres != None:
            raise(Exception,"p_result double positive error.")

    p_result = tempres

    try:
        a_result = np.sqrt(Ka1 * Ka2)
        delta = np.abs(a_result - p_result)
    except TypeError:
        a_result = np.NaN
        delta = np.NaN

    delta_dspace = np.append(delta_dspace,[delta])

    print(i)


#Kw_zspace_sap = Kw_zspace_sap * 10**1
#Ka2_yspace_sap = Ka2_yspace_sap * 10**-2
#Ka1_xspace_sap = Ka1_xspace_sap * 10**-2
#use for wrong scaling correction

figure = mlab.figure('DensityPlot')
pts = mlab.points3d(Ka1_xspace_sap, Ka2_yspace_sap, Kw_zspace_sap,delta_dspace , scale_mode='none', scale_factor=0.0005)
#scale_factor = 0.0000001 or 0.001
mlab.xlabel('Ka1,unit : 10^2')
mlab.ylabel('Ka2,unit : 10^2')
mlab.zlabel('Kw,unit : 10^-1')

mlab.axes(None, extent=[0, 0.1, 0, 0.1, 0, 0.1])
#use for wrong scaling correction too
#mlab.axes()

mlab.show()


