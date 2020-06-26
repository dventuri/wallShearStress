"""
From Prof. H. Kudela:
http://fluid.itcmp.pwr.wroc.pl/~znmp/dydaktyka/fundam_FM/Lecture13.pdf, and
http://www.itcmp.pwr.wroc.pl/~znmp/dydaktyka/fundam_FM/Lecture_no3_Turbulent_flow_Modelling.pdf.


In a fully-developed flow, there is a balance between the pressure gradient
and viscous forces:
P1*pi*r^2 - P2*pi*r^2 - 2*pi*r*L*tau = 0,
or
-Delta{P}/L = (2*tau)/r.

In order for 'tau' to be linearly dependent on the radius and to remove 'r' dependency
on the equation, let
tau = C*r,
since at r=0 -> tau=0.

As in r = D/2 -> tau = tau_w:
C = 2*tau_w/D,
and
-Delta{P}/L = dp/dx = (4*tau_w)/D

From the boundary layer theory, using a turbulent flow solution:
tau_w = (rho*f*Uavg^2)/8,
where the friction factor may be calculated from Prandtl's or Blasius'
expressions.
"""
# TODO: #1 Needs double-check
from scipy.optimize import newton
import numpy as np

U = 7.89
D = 30.5/1000
rho = 1.2
mu = 1.84469e-5

Re = rho*U*D/mu


def fBlasius(Re):
    f = (100*Re)**(-0.25)
    return f


def fPrandtl(Re,f):
    func = 2*np.log10(Re*(f**0.5))-0.8+1/(f**0.5)
    return func


root1 = fBlasius(Re)
print(root1)
root2 = newton(fPrandtl,fBlasius(Re),args=(Re,))
print(root2)

tauW1 = (rho*root1*U**2)/8
print(tauW1)
tauW2 = (rho*root2*U**2)/8
print(tauW2)

dpdx1 = 4*tauW1/D
print(dpdx1)
dpdx2 = 4*tauW2/D
print(dpdx2)
