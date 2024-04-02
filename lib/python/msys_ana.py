########################################################################
# MechSys - Open Library for Mechanical Systems                        #
# Copyright (C) 2009 Dorival M. Pedroso                                #
#                                                                      #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# any later version.                                                   #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with this program. If not, see <http://www.gnu.org/licenses/>  #
########################################################################

from msys_invs import *
from msys_fig import *

# CamClay: closed-form solution: ev, ed
# =====================================
# Proportional loading: k = dq/dp (path)
# p and q: Octahedral invariants
def CamClay_ev_ed_oct (prms, inis, p0, q0, p, q):

    # parameters and initial values
    lam = prms['lam']
    kap = prms['kap']
    nu  = prms['nu']
    v0  = inis['v0']
    M   = phi_calc_M (prms['phi'], 'oct')
    chi = (kap-lam)/v0

    # auxiliary variables
    alp    = (3.0*(1.0-2.0*nu))/(2.0*(nu+1.0))
    dp, dq = p-p0, q-q0
    k      = dq/dp
    r      = q/(M*p)
    r0     = q0/(M*p0)

    # check failure condition
    if r>1.0: raise Exception('[1;31mCamClay_ev_ed_oct: Invalid stress state (failure condition). r=q/(M*p)=%g must be smaller than or equal to 1[0m'%r)

    # elastic strains
    eve = -(kap*log(p/p0))/v0
    ede = sqrt(3.0)*k*kap*log(p/p0)/(2.0*alp*v0)
    if dp<0 or dq<0: return eve, ede

    # elastoplastic strains
    evp = chi*(log((r**2.0+1.0)/(r0**2.0+1.0))+log(p/p0))
    edp = 2.0*chi*( k*M*log(p/p0)/(k**2.0-M**2.0) - k*log((r+1.0)/(r0+1.0))/(2.0*k+2.0*M) + k*log((r-1.0)/(r0-1.0))/(2.0*k-2.0*M) + arctan(r) - arctan(r0) )/(sqrt(3.0)*M)
    ev  = eve + evp
    ed  = ede + edp
    return ev, ed

def CamClay_ev_ed (prms, inis, kcam, pcam0, qcam0, dpcam, npts):
    dp     = dpcam*sqrt(3.0)
    k      = kcam*sqrt(2.0)/3.0
    p0, q0 = pcam0*sqrt(3.0), qcam0*sqrt(2.0/3.0)
    M      = phi_calc_M (prms['phi'], 'oct')
    if k>M:
        dpmax = (M*p0-q0)/(k-M)
        if dp>dpmax: dp = 0.999*dpmax
    pf, qf = p0+dp, q0+k*dp
    p      = linspace (p0, pf, npts)
    q      = q0 + k*(p-p0)
    ev, ed = zeros(npts), zeros(npts)
    for i in range(npts-1): ev[1+i], ed[1+i] = CamClay_ev_ed_oct (prms, inis, p0, q0, p[1+i], q[1+i])
    pcam = p/sqrt(3.0)
    qcam = q*sqrt(3.0/2.0)
    return pcam, qcam, ev, ed
