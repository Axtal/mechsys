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

from numpy import sqrt, matrix, arcsin, sin, pi, array, diag
import numpy.linalg as npyla

# Calculate ev and ed
# ===================
# eps is a column matrix in Mandel's basis
def eps_calc_ev_ed(eps,Type='cam'):
    if Type=='oct':
        ev = (eps[0,0]+eps[1,0]+eps[2,0])/sqrt(3.0)
        ed = sqrt((eps[0,0]-eps[1,0])**2.0 + (eps[1,0]-eps[2,0])**2.0 + (eps[2,0]-eps[0,0])**2.0 + 3.0*(eps[3,0]**2.0))/sqrt(3.0)
    elif Type=='cam':
        ev = eps[0,0]+eps[1,0]+eps[2,0]
        ed = sqrt((eps[0,0]-eps[1,0])**2.0 + (eps[1,0]-eps[2,0])**2.0 + (eps[2,0]-eps[0,0])**2.0 + 3.0*(eps[3,0]**2.0))*(sqrt(2.0)/3.0)
    else: raise Exception('eps_calc_ev_ed: Method not available for invariant Type==%s'%Type)
    return ev, ed


# Calculate p
# ===========
# sig is a column matrix in Mandel's basis
def sig_calc_p(sig,Type='oct'):
    if Type=='oct':
        p = -(sig[0,0]+sig[1,0]+sig[2,0])/sqrt(3.0)
    elif Type=='cam':
        p = -(sig[0,0]+sig[1,0]+sig[2,0])/3.0
    return p


# Calculate p and q
# =================
# sig is a column matrix in Mandel's basis
def sig_calc_p_q(sig,Type='oct',with_derivs=False,qtol=1.0e-8):
    if Type=='oct':
        p = -(sig[0,0]+sig[1,0]+sig[2,0])/sqrt(3.0)
        q = sqrt((sig[0,0]-sig[1,0])**2.0 + (sig[1,0]-sig[2,0])**2.0 + (sig[2,0]-sig[0,0])**2.0 + 3.0*(sig[3,0]**2.0))/sqrt(3.0)
    elif Type=='cam':
        p = -(sig[0,0]+sig[1,0]+sig[2,0])/3.0
        q = sqrt((sig[0,0]-sig[1,0])**2.0 + (sig[1,0]-sig[2,0])**2.0 + (sig[2,0]-sig[0,0])**2.0 + 3.0*(sig[3,0]**2.0))/sqrt(2.0)
    else: raise Exception('sig_calc_p_q: Method not available for invariant Type==%s'%Type)
    if with_derivs:
        if Type=='oct':
            I      = matrix([[1.0],[1.0],[1.0],[0.0]])
            dpdsig = (-1.0/sqrt(3.0))*I
            if q>qtol:
                s      = sig - ((sig[0,0]+sig[1,0]+sig[2,0])/3.0)*I
                dqdsig = (1.0/q)*s
            else: dqdsig = matrix([[0.0],[0.0],[0.0],[0.0]])
        else: raise Exception('sig_calc_p_q: Derivatives not available for invariant Type==%s'%Type)
        return p, q, dpdsig, dqdsig
    else: return p, q


# Calculate p and q
# =================
def s123_calc_p_q(s123,Type='oct',with_derivs=False,qtol=1.0e-8):
    if Type=='oct':
        p = -(s123[0]+s123[1]+s123[2])/sqrt(3.0)
        q = sqrt((s123[0]-s123[1])**2.0 + (s123[1]-s123[2])**2.0 + (s123[2]-s123[0])**2.0)/sqrt(3.0)
    elif Type=='cam':
        p = -(s123[0]+s123[1]+s123[2])/3.0
        q = sqrt((s123[0]-s123[1])**2.0 + (s123[1]-s123[2])**2.0 + (s123[2]-s123[0])**2.0)/sqrt(2.0)
    else: raise Exception('s123_calc_p_q: Method not available for invariant Type==%s'%Type)
    if with_derivs:
        dev_s123 = array([2.0*s123[0]-s123[1]-s123[2], 2.0*s123[1]-s123[2]-s123[0], 2.0*s123[2]-s123[0]-s123[1]])/3.0
        dpds123  = zeros(3)
        dqds123  = zeros(3)
        if Type=='oct':
            for k in range(3):
                dpds123[k] = -sqrt(3.0)/3.0
                if q>qtol: dqds123[k] = dev_s123[k]/q
        elif Type=='cam':
            for k in range(3):
                dpds123[k] = -1.0/3.0
                if q>qtol: dqds123[k] = 1.5*dev_s123[k]/q
        return p, q, dpds123, dqds123
    else: return p, q


# Invariants: p, q, t
# ===================
def sig_calc_pqt(sig,Type='oct',with_derivs=False,qtol=1.0e-8):
    if Type=='oct':
        p = -(sig[0,0]+sig[1,0]+sig[2,0])/sqrt(3.0)
        q = sqrt((sig[0,0]-sig[1,0])**2.0 + (sig[1,0]-sig[2,0])**2.0 + (sig[2,0]-sig[0,0])**2.0 + 3.0*(sig[3,0]**2.0))/sqrt(3.0)
        r = q
    elif Type=='cam':
        p = -(sig[0,0]+sig[1,0]+sig[2,0])/3.0
        q = sqrt((sig[0,0]-sig[1,0])**2.0 + (sig[1,0]-sig[2,0])**2.0 + (sig[2,0]-sig[0,0])**2.0 + 3.0*(sig[3,0]**2.0))/sqrt(2.0)
        r = sqrt(2.0/3.0)*q
    I = matrix([[1.0],[1.0],[1.0],[0.0]])
    s = sig - ((sig[0,0]+sig[1,0]+sig[2,0])/3.0)*I
    t = -3.0*sqrt(6.0)*calc_determinant(s)/(r**3.0) if q>qtol else 0.0
    if with_derivs:
        if Type=='oct':
            dpdsig = (-1.0/sqrt(3.0))*I
            dqdsig = (1.0/q)*s if q>qtol else matrix([[0.0],[0.0],[0.0],[0.0]])
        elif Type=='cam':
            dpdsig = (-1.0/3.0)*I
            dqdsig = (1.5/q)*s if q>qtol else matrix([[0.0],[0.0],[0.0],[0.0]])
        if q>qtol:
            ss     = calc_squared(s)
            devss  = ss - ((ss[0,0]+ss[1,0]+ss[2,0])/3.0)*I
            dtdsig = (-3.0*t/(q**2.0))*s - (3.0*sqrt(6.0)/(q**3.0))*devss
        else: dtdsig = matrix([[0.0],[0.0],[0.0],[0.0]])
        return p, q, t, dpdsig, dqdsig, dtdsig
    else: return p, q, t


# Invariants: p, q, t
# ===================
def s123_calc_pqt(s123,Type='oct',with_derivs=False,qtol=1.0e-8):
    s1, s2, s3 = s123[0], s123[1], s123[2]
    if Type=='oct':
        p = -(s1+s2+s3)/sqrt(3.0)
        q = sqrt((s1-s2)**2.0 + (s2-s3)**2.0 + (s3-s1)**2.0)/sqrt(3.0)
        r = q
    elif Type=='cam':
        p = -(s1+s2+s3)/3.0
        q = sqrt((s1-s2)**2.0 + (s2-s3)**2.0 + (s3-s1)**2.0)/sqrt(2.0)
        r = sqrt(2.0/3.0)*q
    I = array([1.0,1.0,1.0])
    s = s123 - ((s1+s2+s3)/3.0)*I
    t = -3.0*sqrt(6.0)*s[0]*s[1]*s[2]/(r**3.0) if q>qtol else 0.0
    if t<=-1.0: t = -1.0
    if t>= 1.0: t =  1.0
    if with_derivs:
        if Type=='oct':
            dpds123 = (-1.0/sqrt(3.0))*I
            dqds123 = (1.0/q)*s if q>qtol else zeros(3)
        elif Type=='cam':
            dpds123 = (-1.0/3.0)*I
            dqds123 = (1.5/q)*s if q>qtol else zeros(3)
        if q>qtol:
            l       = (s1-s2)*(s2-s3)*(s3-s1)
            B       = array([s3-s2, s1-s3, s2-s1])
            dtds123 = (-sqrt(6.0)*l/(q**5.0))*B
        else: dtds123 = array([0.0,0.0,0.0])
        return p, q, t, dpds123, dqds123, dtds123
    else: return p, q, t


# Euclidian norm of a tensor
# ==========================
# A is a column matrix in Mandel's basis
def calc_norm(A):
    return sqrt((A.T*A)[0,0])


# Determinant of a tensor
# =======================
# A is a column matrix in Mandel's basis
def calc_determinant(A):
    return A[0,0]*A[1,0]*A[2,0] - A[2,0]*A[3,0]*A[3,0]/2.0


# Squared of a tensor
# ===================
# A is a column matrix in Mandel's basis
def calc_squared(A):
    return matrix([[A[0,0]*A[0,0]+A[3,0]*A[3,0]/2.0], [A[1,0]*A[1,0]+A[3,0]*A[3,0]/2.0], [A[2,0]*A[2,0]], [A[0,0]*A[3,0]+A[1,0]*A[3,0]]])


# Characteristic invariants
# =========================
# A is a column matrix in Mandel's basis
def char_invs(A, with_derivs=False):
    I1 = A[0,0] + A[1,0] + A[2,0]
    I2 = A[0,0]*A[1,0] + A[1,0]*A[2,0] + A[2,0]*A[0,0] - A[3,0]*A[3,0]/2.0
    I3 = A[0,0]*A[1,0]*A[2,0] - A[2,0]*A[3,0]*A[3,0]/2.0
    if with_derivs: # calculate derivatives
        I     = matrix([[1.0],[1.0],[1.0],[0.0]])
        dI1dA = I.copy()
        dI2dA = I1*I - A
        dI3dA = calc_squared(A) - I1*A + I2*I
        return I1, I2, I3, dI1dA, dI2dA, dI3dA
    else: return I1, I2, I3


# Calculate t
# ===========
# t = sin(3*theta)
def sig_calc_t(sig,with_derivs=False,qtol=1.0e-8):
    I = matrix([[1.0],[1.0],[1.0],[0.0]])
    s = sig - ((sig[0,0]+sig[1,0]+sig[2,0])/3.0)*I
    q = calc_norm(s)
    if q>qtol:  t = -3.0*sqrt(6.0)*calc_determinant(s)/(q**3.0)
    else:       t = -1.0
    if t<=-1.0: t = -1.0
    if t>= 1.0: t =  1.0
    if with_derivs:
        if q>qtol:
            ss     = calc_squared(s)
            devss  = ss - ((ss[0,0]+ss[1,0]+ss[2,0])/3.0)*I
            dtdsig = (-3.0*t/(q**2.0))*s - (3.0*sqrt(6.0)/(q**3.0))*devss
        else: dtdsig = matrix([[0.0],[0.0],[0.0],[0.0]])
        return t, dtdsig
    return t


# Calculate s123
# ==============
def pqt_calc_s123(p,q,t,Type='oct'):
    if Type=='cam':
        p = p*sqrt(3.0)
        q = q*sqrt(2.0/3.0)
    elif Type=='oct': pass # p=poct, q=qoct
    else: raise Exception('pqt_calc_s123: Method not available for invariant Type==%s'%Type)
    #if t<=-1.0: t = -1.0
    #if t>= 1.0: t =  1.0
    th  = arcsin(t)/3.0
    s1  = -p/sqrt(3.0) + 2.0*q*sin(th-2.0*pi/3.0)/sqrt(6.0)
    s2  = -p/sqrt(3.0) + 2.0*q*sin(th)           /sqrt(6.0)
    s3  = -p/sqrt(3.0) + 2.0*q*sin(th+2.0*pi/3.0)/sqrt(6.0)
    return array([s1, s2, s3])


# Convert p and q
# ===============
def convert_p_q(p,q,t,TypeOri,TypeDest):
    s123 = pqt_calc_s123(p,q,t,TypeOri)
    return s123_calc_p_q(s123,TypeDest)


# Calculate principal stresses
# ============================
# sig is a column matrix in Mandel's basis
def sig_calc_s123(sig,with_projs=False,do_sort=False):
    sq2  = sqrt(2.0)
    smat = matrix([[ sig[0,0]     , sig[3,0]/sq2 ,     0.0  ],
                   [ sig[3,0]/sq2 , sig[1,0]     ,     0.0  ],
                   [     0.0      ,     0.0      , sig[2,0] ]])
    if with_projs: # with eigenprojectors
        if do_sort: raise Exception('sig_calc_s123: do_sort cannot be used when with_projs==True')
        s123, evecs = npyla.eig(smat)
        projs       = [] # eigenprojectors
        for v in evecs:
            Pmat = matrix(v.reshape(3,1)) * matrix(v)
            projs.append(matrix([[Pmat[0,0]],[Pmat[1,1]],[Pmat[2,2]],[Pmat[0,1]*sq2]]))
        return s123, projs
    else:
        s123 = npyla.eigvalsh(smat)
        if do_sort: s123.sort()#reverse=True)
        return s123

# Calculate principal stresses and rotation matrix
# ================================================
# sig is a column matrix (4x1) in Mandel's basis
def sig_calc_rot(sig):
    sq2  = sqrt(2.0)
    smat = matrix([[ sig[0,0]     , sig[3,0]/sq2 ,     0.0  ],
                   [ sig[3,0]/sq2 , sig[1,0]     ,     0.0  ],
                   [     0.0      ,     0.0      , sig[2,0] ]])
    return npyla.eig(smat)

# Calculate principal strains
# ===========================
# eps is a column matrix in Mandel's basis
def eps_calc_e123(eps, do_sort=False):
    sq2  = sqrt(2.0)
    emat = matrix([[ eps[0,0]     , eps[3,0]/sq2 ,     0.0  ],
                   [ eps[3,0]/sq2 , eps[1,0]     ,     0.0  ],
                   [     0.0      ,     0.0      , eps[2,0] ]])
    e123 = npyla.eigvalsh(emat)
    if do_sort: e123.sort()
    return e123


# Octahedral coordinates: sa, sb, sc
# ==================================
def s123_calc_oct(s123):
    sa = sqrt(2.0)*(s123[1]-s123[2])/2.0
    sb = (s123[2]+s123[1]-2.0*s123[0])/sqrt(6.0)
    sc = -(s123[0]+s123[1]+s123[2])/sqrt(3.0)
    return sa, sb, sc


# Octahedral coordinates: sa, sb, sc
# ==================================
def sxyz_calc_oct(sxyz):
    sa = (sxyz[1]-sxyz[0])/sqrt(2.0)
    sb = (sxyz[1]+sxyz[0]-2.0*sxyz[2])/sqrt(6.0)
    sc = -(sxyz[0]+sxyz[1]+sxyz[2])/sqrt(3.0)
    return sa, sb, sc


# Principal stresses: s1, s2, s3
# ==============================
def oct_calc_s123(sa,sb,sc):
    s1 =               - 2.0*sb/sqrt(6.0) - sc/sqrt(3.0)
    s2 =  sa/sqrt(2.0) +     sb/sqrt(6.0) - sc/sqrt(3.0)
    s3 = -sa/sqrt(2.0) +     sb/sqrt(6.0) - sc/sqrt(3.0)
    return array([s1,s2,s3])


# Principal stresses: sx, sy, sz
# ==============================
def oct_calc_sxyz(sa,sb,sc):
    sz =               - 2.0*sb/sqrt(6.0) - sc/sqrt(3.0)
    sy =  sa/sqrt(2.0) +     sb/sqrt(6.0) - sc/sqrt(3.0)
    sx = -sa/sqrt(2.0) +     sb/sqrt(6.0) - sc/sqrt(3.0)
    return array([sx,sy,sz])


# Calculate M
# ===========
# phi: friction angle at compression (degrees)
# M:   max q/p at compression
def phi_calc_M(phi,Type='oct'):
    sphi = sin(phi*pi/180.0)
    if   Type=='oct': M = 2.0*sqrt(2.0)*sphi/(3.0-sphi)
    elif Type=='cam': M = 6.0*sphi/(3.0-sphi)
    elif Type=='smp':
            sq2 = sqrt(2.0)
            sq3 = sqrt(3.0)
            sq6 = sqrt(6.0)
            eta = 2.0*sq2*sphi/(3.0-sphi)
            c   = sqrt((2.0+sq2*eta-2.0*eta*eta)/(3.0*sq3*(sq2*eta+2.0)))
            a   = sqrt((2.0*eta+sq2)/sq6)
            b   = sqrt((sq2-eta)/sq6)
            M   = sqrt((eta*eta+1.0)/(c*c*(a+2.0*b)**2.0)-1.0)
    else: raise Exception('phi_calc_M: Method not available for invariant Type==%s'%Type)
    return M


# Calculate phi
# =============
# M:   max q_cam/p_cam at compression
# phi: friction angle at compression (degrees)
def M_calc_phi(M,Type='oct',Mode='comp'):
    if Mode=='comp':
        if   Type=='oct': sphi = 3.0*M/(M+2.0*sqrt(2.0))
        elif Type=='cam': sphi = 3.0*M/(M+6.0)
        else: raise Exception('M_calc_phi: Method not available for invariant Type==%s'%Type)
    elif Mode=='ext':
        if   Type=='oct': sphi = 3.0*M/(-M+2.0*sqrt(2.0))
        elif Type=='cam': sphi = 3.0*M/(-M+6.0)
        else: raise Exception('M_calc_phi: Method not available for invariant Type==%s'%Type)
    else: raise Exception('M_calc_phi: Mode %s not available, use either Mode=comp or Mode=ext')
    return arcsin(sphi)*180.0/pi


# Calculate cu
# ============
# Calculate undrained cohesion (cu) for given q_failure
def qf_calc_cu(qf, qType="oct", psa=False):
    coef = sqrt(3.0) if psa else 2.0
    if   qType=="oct": return qf/(sqrt(2.0/3.0)*coef)
    elif qType=="cam": return qf/coef
    else: raise Exception("qf_calc_cu: Method is not available for invariant qType==%s"%qType)


# Calculate pf and qf for a proportional load path
# ================================================
# p and q are 'cam' invariants
# k = dq/dp
# M = qf/pf
def prop_calc_pqf(p0,q0,M,k):
    pf = (k*p0-q0)/(k-M)
    qf = M*pf
    return pf, qf


# Output data
# ===========
# sig is a column matrix in Mandel's basis
# eps is a column matrix in Mandel's basis
def output_data(sig,eps,Type='oct'):
    p,  q      = sig_calc_p_q   (sig,Type)
    ev, ed     = eps_calc_ev_ed (eps)
    dat        = {}
    dat['sx']  = sig[0,0]
    dat['sy']  = sig[1,0]
    dat['sz']  = sig[2,0]
    dat['sxy'] = sig[3,0]/sqrt(2.0)
    dat['p']   = p
    dat['q']   = q
    dat['ex']  = eps[0,0]
    dat['ey']  = eps[1,0]
    dat['ez']  = eps[2,0]
    dat['exy'] = eps[3,0]/sqrt(2.0)
    dat['ev']  = ev
    dat['ed']  = ed
    return dat


# Check invariants
# ================
if __name__=='__main__':
    #sig         = matrix([[0.375], [-0.25], [0.375], [0.0*sqrt(2.0)]])
    #sig         = matrix([[-90.0], [-120.0], [-60.0], [30.0*sqrt(2.0)]])
    #sig         = matrix([[-90.0], [-90.0], [-90.0], [0.0*sqrt(2.0)]])
    sig         = matrix([[-90.0], [-20.0], [-20.0], [0.0*sqrt(2.0)]])
    I           = matrix([[1.0],[1.0],[1.0],[0.0]])
    dev_sig     = sig - ((sig[0,0]+sig[1,0]+sig[2,0])/3.0)*I
    s123        = sig_calc_s123 (sig)
    s1,s2,s3    = s123[0], s123[1], s123[2]
    p, q        = sig_calc_p_q  (sig)
    t           = sig_calc_t    (sig)
    p1,q1       = s123_calc_p_q (s123)
    I123        = array([1.0,1.0,1.0])
    dev_s123    = s123 - ((s1+s2+s3)/3.0)*I123
    q2          = sqrt(dot(dev_s123,dev_s123))
    q3          = calc_norm(dev_sig)
    t1          = -3.0*sqrt(6.0)*dev_s123[0]*dev_s123[1]*dev_s123[2]/(q2**3.0) if q2>0.0 else -1.0
    I1,I2,I3, dI1,dI2,dI3 = char_invs(sig, with_derivs=True)
    I1b         = s1+s2+s3
    I2b         = s1*s2+s2*s3+s3*s1
    I3b         = s1*s2*s3
    s1b,s2b,s3b = pqt_calc_s123 (p,q,t)
    sa, sb, sc  = s123_calc_oct (s123)
    s123c       = oct_calc_s123 (sa,sb,sc)
    s1c,s2c,s3c = s123c[0], s123c[1], s123c[2]
    print 'sig                = [%g, %g, %g, %g]  '%(sig[0,0],sig[1,0],sig[2,0],sig[3,0])
    print 'dev_sig            = [%g, %g, %g, %g]  '%(dev_sig[0,0],dev_sig[1,0],dev_sig[2,0],dev_sig[3,0])
    print 'theta              = %g                '%(60.0*arcsin(t)/pi)
    print 'p,  q,  t          = %g, %g, %g        '%(p,  q,  t)
    print 'p1, q1, t1, q2, q3 = %g, %g, %g, %g, %g'%(p1, q1, t1,   q2, q3)
    print 'sa,  sb,  sc       = %g, %g, %g        '%(sa,  sb,  sc)
    print 's1,  s2,  s3       = %g, %g, %g        '%(s1,  s2,  s3)
    print 's1b, s2b, s3b      = %g, %g, %g        '%(s1b, s2b, s3b)
    print 's1c, s2c, s3c      = %g, %g, %g        '%(s1c, s2c, s3c)
    print 'I1 , I2,  I3       = %g, %g, %g        '%(I1 , I2 , I3)
    print 'I1b, I2b, I3b      = %g, %g, %g        '%(I1b, I2b, I3b)
    print 'dI1dsig            = [%g, %g, %g, %g]  '%(dI1[0],dI1[1],dI1[2],dI1[3])
    print 'dI2dsig            = [%g, %g, %g, %g]  '%(dI2[0],dI2[1],dI2[2],dI2[3])
    print 'dI3dsig            = [%g, %g, %g, %g]  '%(dI3[0],dI3[1],dI3[2],dI3[3])
    print 'ERROR(p,q,t)       = %g, %g, %g, %g, %g'%(p-p1, q-q1, t-t1, q-q2, q-q3)
    print 'ERROR(sk-skb)      = %g, %g, %g        '%(s1-s1b, s2-s2b, s3-s3b)
    print 'ERROR(sk-skc)      = %g, %g, %g        '%(s1-s1c, s2-s2c, s3-s3c)
    print 'ERROR(Ik)          = %g, %g, %g        '%(I1-I1b, I2-I2b, I3-I3b)
