# -*- coding: utf-8 -*-
'''
accpy.simulate.const
Felix Kramer (felix.kramer@physik.hu-berlin.de)
'''

# ratio of circumference to diameter of circle
pi = 3.141592653589793115997963468544185161590576171875
# speed of light----------------------------------------/ (m/s)
cl = 299792458.
# elementary charge-------------------------------------/ (e)=(As)=(C)
qe = 1.602176565e-19
# electron mass-----------------------------------------/ (kg)
me = 9.10938291e-31
# proton mass-------------------------------------------/ (kg)
mp = 1.672621777e-27
# muon mass---------------------------------------------/ (kg)
mu = 1.883531475e-28
# electron restenergy-----------------------------------/ (J)=(Nm)=(Ws)
Ee = me*cl**2
# proton restenergy-------------------------------------/ (J)=(Nm)=(Ws)
Ep = mp*cl**2
# muon restenergy---------------------------------------/ (J)=(Nm)=(Ws)
Eu = mu*cl**2
# classical radius of electron--------------------------/ (m)
re = qe**2/(me*1e7)
# classical radius of proton----------------------------/ (m)
rp = qe**2/(mp*1e7)
# classical radius of muon------------------------------/ (m)
ru = qe**2/(mu*1e7)
# vacuum permeability / magnetic field contant----------/ (N/A^2)=(Henry/m)
u0 = 4*pi*1E-7
# vacuum permittivity / elektrical field const----------/ (As/Vm)
e0 = 1/(u0*cl**2)
# Planck constant---------------------------------------/ (Js)
hp = 6.62606957e-34
# reduced Planck constant-------------------------------/ (Js)
hb = hp/2/pi
# Boltzmann constant------------------------------------/ (J/K)
kb = 1.3806488e-23
# Avogadro number---------------------------------------/ (1/mole)
NA = 6.02214129e23
# gas constant------------------------------------------/ (J/Kmole)
RG = kb*NA
# gravitational constant--------------------------------/ (Nm^2/kg^2)
Gr = 6.67384e-11
# gravitational acceleration Berlin---------------------/ (m/s^2)
ge = 9.812753
