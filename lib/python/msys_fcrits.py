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

from mechsys   import *
from msys_invs import *
from msys_fig  import *
from numpy     import ogrid

class FCrits:
    # Constructor
    # ===========
    def __init__(self):
        # data
        self.c       = 0.0                   # cohesion
        self.phi     = 30.0                  # friction angle for FC
        self.kY      = None                  # von Mises coefficient
        self.b       = None                  # anisotropic fc: b coefficient
        self.alp     = None                  # anisotropic fc: alpha coefficient
        self.a       = None                  # anisotropic fc: bedding planes normal
        self.sstar   = None                  # anisotropic fc: sig to pass fc through
        self.R       = None                  # anisotropic fc: R coefficient
        self.obliq1  = False                 # anisotropic fc: oblique projection type 1
        self.obliq2  = False                 # anisotropic fc: oblique projection type 2
        self.Imat    = matrix(diag(ones(3))) # identity matrix 3x3
        self.t       = 1.0                   # t=sin(3th) for p-q plane
        self.pqty    = 'cam'                 # pq type for p-q plane
        self.r       = None                  # radius in Pi plane
        self.sc      = 1.0                   # mean pressure => distance of cross-section to the origin of Haigh-Westergaard space
        self.samin   = None                  # min sa in Pi plane
        self.samax   = None                  # max sa in Pi plane
        self.sbmin   = None                  # min sb in Pi plane
        self.sbmax   = None                  # max sb in Pi plane
        self.scmin   = None                  # min sc in Pi plane (3D only)
        self.scmax   = None                  # max sc in Pi plane (3D only)
        self.sxyz    = True                  # consider sx(right), sy(left), sz(up) instead of s1(up), s2(right), s3(left)
        self.rst_txt = []                    # rosette text
        self.leg_plt = []                    # legend plots (items appended in plot)
        self.leg_txt = []                    # legend texts (items appended in plot)
        self.smpinvs = SMPInvs()             # SMP invariants
        self.pmin    = None                  # min p for pq plane
        self.pmax    = None                  # max p for pq plane
        self.qmin    = None                  # min q for pq plane
        self.qmax    = None                  # max q for pq plane


    # Intersection and friction angle
    # ===============================
    def inter (self, sig_star=2, do_plot=True):
        sphi = sin(self.phi*pi/180.0)
        A    = (1.0+sphi)/(1.0-sphi)
        if   sig_star==0: l = array([-A , -1., -1.])
        elif sig_star==1: l = array([-1., -A , -1.])
        elif sig_star==2: l = array([-1., -1., -A ])
        dev_sig_tr = array([2.*l[0]-l[1]-l[2],
                            2.*l[1]-l[2]-l[0],
                            2.*l[2]-l[0]-l[1]])/3.0
        sig   = array([-1.,-1.,-1.])/sqrt(3.)
        a     = 0.5
        sig_a = sig + a*dev_sig_tr
        if do_plot:
            sa, sb, sc = sxyz_calc_oct (sig_a)
            print sa, sb, sc
            plot ([sa], [sb], 'go')


    # Set sc
    # ======
    def set_sc (self, pcam=1.0/sqrt(3.0)):
        self.sc = pcam*sqrt(3.0)
        print 'sc    = ', self.sc


    # Set constants
    # =============
    def set_ctes (self, phi_deg=30.0, c=None, sY=None, b=None, pTol=0.1, pRef=100.0, psa=False):
        self.c     = 0.0 if c==None else c
        self.phi   = phi_deg
        self.sphi  = sin(self.phi*pi/180.0)
        self.tgphi = tan(self.phi*pi/180.0)
        self.kDP   = 2.0*sqrt(2.0)*self.sphi/(3.0-self.sphi)
        self.kMN   = 9.0 + 8.0*self.tgphi**2.0
        self.kLD   = ((3.0-self.sphi)**3.0)/((1.0+self.sphi)*((1.0-self.sphi)**2))
        self.cbar = sqrt(3.0)*self.c/self.tgphi
        if sY==None:
            if c==None: self.kVM = self.sc*self.kDP
            else:       self.kVM = sqrt(2.0)*self.c if psa else 2.0*sqrt(2.0/3.0)*self.c
        else: self.kVM = sqrt(2.0/3.0)*sY
        if not b==None: self.smpinvs.b = b
        A         = (1.0+self.sphi)/(1.0-self.sphi)
        sig       = [-A, -1., -1., 0., 0., 0.]
        sp, sq    = self.smpinvs.Calc (sig)
        self.kGE  = sq/sp
        self.pTol = pTol
        self.pRef = pRef
        print 'c     = ', self.c
        print 'phi   = ', self.phi
        print 'kDP   = ', self.kDP
        print 'kMN   = ', self.kMN
        print 'kLD   = ', self.kLD
        print 'kVM   = ', self.kVM
        print 'kGE   = ', self.kGE
        print 'pTol  = ', self.pTol
        print 'pRef  = ', self.pRef


    # Generic FC constants
    # ====================
    def calc_arc (self, M):
        den = 1.0 + M*M
        pc  = self.pTol/(1.0-M/sqrt(den))
        pm  = pc/den
        r2  =((pc*M)**2.0)/den
        return pc, pm, r2


    # Set constants of anisotropic criterion
    # ======================================
    def aniso (self, phi_deg=30.0, b=None, alpha=0.1, a=[0.,0.,1.], sig_star=2):

        # set constants
        a = matrix(a).T
        self.phi   = phi_deg
        self.b     = b
        self.alp   = alpha
        self.a     = a / norm(a)
        self.sstar = sig_star
        self.set_ctes (phi_deg=phi_deg)
        print 'phi   = ', self.phi
        print 'b     = ', self.b
        print 'alp   = ', self.alp
        print 'a     = ', self.a[0,0], self.a[1,0], self.a[2,0]
        print 'sstar = ', self.sstar

        # assemble l, vector with principal values
        sphi = sin(self.phi*pi/180.0)
        A    = (1.0+sphi)/(1.0-sphi)
        if   sig_star==0: l = array([-A , -1., -1.])
        elif sig_star==1: l = array([-1., -A , -1.])
        elif sig_star==2: l = array([-1., -1., -A ])
        else: raise Exception('set_anisocrit: sig_star must be 0, 1, or 2. %d is invalid'%sig_star)

        # invariants
        nvec, namp, sig, tau = self.get_nvec_namp_sig_tau (l)
        self.R               = tau/sig


    # Anisotropic criterion invariants
    # ================================
    def get_nvec_namp_sig_tau(self, l, Q=None):

        # calculate nvec
        if self.b==None:
            I2   = l[0]*l[1] + l[1]*l[2] + l[2]*l[0]
            I3   = l[0]*l[1]*l[2]
            nvec = -matrix(sqrt(I3/(I2*l))).T
        else:
            nnew = -matrix([[abs(l[0])**(-self.b)], [abs(l[1])**(-self.b)], [abs(l[2])**(-self.b)]])
            nvec = nnew / norm(nnew)

        # calculate namp
        if Q==None: namp = nvec + self.alp *       self.a
        else:       namp = nvec + self.alp * Q.T * self.a
        namp = namp / norm(namp)

        # calculate sig and tau
        if   self.obliq1: Pamp = (nvec * namp.T) / (nvec.T * namp)[0]
        elif self.obliq2: Pamp = (namp * nvec.T) / (namp.T * nvec)[0]
        else:             Pamp =  namp * namp.T
        Qamp = self.Imat - Pamp
        tamp = diag(l) * namp
        pamp = Pamp*tamp
        qamp = Qamp*tamp
        sp   = norm(pamp)
        sq   = norm(qamp)

        # return invariants
        return nvec, namp, sp, sq


    # Failure criteria (or yield surface)
    # ===================================
    def func (self, sig, typ):

        # yield surface
        ysurf = False
        if typ[0]=='y':
            ysurf = True
            typ   = typ[1:]
            pc    = self.sc/2.0

        # invariants
        p, q = sig_calc_p_q (sig)

        # von Mises
        if typ=='VM':
            if ysurf: raise Exception('func: ysurf is not available with VM')
            else: f = q - self.kVM

        # Drucker/Prager
        elif typ=='DP':
            if ysurf: f = ((p-pc)/pc)**2.0 + (q/(self.kDP*pc))**2.0 - 1.0
            else:     f = q - (p + self.cbar)*self.kDP

        # Mohr/Coulomb
        elif typ=='MC':
            t  = sig_calc_t (sig)
            th = arcsin(t)/3.0
            g  = sqrt(2.0)*self.sphi/(sqrt(3.0)*cos(th)-self.sphi*sin(th))
            if ysurf: raise Exception('func: ysurf is not available with MC')
            else: f = q - (p + self.cbar)*g

        # rounded Mohr/Coulomb
        elif typ=='MCr':
            t  = sig_calc_t (sig)
            th = arcsin(t)/3.0
            g  = sqrt(2.0)*self.sphi/(sqrt(3.0)*cos(th)-self.sphi*sin(th))
            if ysurf: raise Exception('func: ysurf is not available with MC')
            pc, pm, r2 = self.calc_arc (g)
            if (p<pm): f = q*q/(self.pRef**2.0) +((p-pc)**2.0)/(self.pRef**2.0) - r2/(self.pRef**2.0)
            else:      f = q/self.pRef - g*p/self.pRef

        # Nakai/Matsuoka
        elif typ=='MN':
            l = sig_calc_s123(sig)
            if l[0]>0.0 or l[1]>0.0 or l[2]>0.0: return -1.0e+8
            I1,I2,I3 = char_invs(sig)
            if ysurf: raise Exception('func: ysurf is not available with MN')
            else: f = I1*I2 - self.kMN*I3

        # Lade/Duncan
        elif typ=='LD':
            sig0 = self.cbar*matrix([[1.0],[1.0],[1.0],[0.0]])
            sig_ = sig+sig0
            l    = sig_calc_s123(sig_)
            if l[0]>0.0 or l[1]>0.0 or l[2]>0.0: return -1.0e+8
            I1,I2,I3 = char_invs(sig_)
            if ysurf: raise Exception('func: ysurf is not available with LD')
            else: f = I1**3.0 - self.kLD*I3

        # Argyris/Sheng
        elif typ=='AS':
            p, q, t = sig_calc_pqt (sig, 'cam')
            Mcs     = 6.0*self.sphi/(3.0-self.sphi)
            om      = ((3.0-self.sphi)/(3.0+self.sphi))**4.0
            M       = Mcs*(2.0*om/(1.0+om-(1.0-om)*t))**0.25;
            if ysurf: raise Exception('func: ysurf is not available with AS')
            else: f = q/p - M

        # Anisotropic
        elif typ=='AMP' or typ=='AMPb' or typ=='AMPba':
            l, Q = sig_calc_rot (sig)
            if l[0]>0.0 or l[1]>0.0 or l[2]>0.0: return 1.0e+8
            nvec, namp, sp, sq = self.get_nvec_namp_sig_tau (l, Q)
            if ysurf:
                z0 = 2.0*pc/sqrt(3.0)
                f  = log(sp/z0) + (1.0/self.bet)*(sq/(self.R*sp))**self.bet
                #f  = ((sp-pc)/pc)**2.0 + (sq/(self.R*pc))**2.0 - 1.0
            else: f = sq - self.R * sp

        # General
        elif typ=='GE':
            sp, sq     = self.smpinvs.Calc ([sig[0,0],sig[1,0],sig[2,0],sig[3,0]])
            pc, pm, r2 = self.calc_arc (self.kGE)
            if ysurf: raise Exception('func: ysurf is not available with GE')
            if (sp<pm): f = sq*sq/(self.pRef**2.0) +((sp-pc)**2.0)/(self.pRef**2.0) - r2/(self.pRef**2.0)
            else:       f = sq/self.pRef - self.kGE*sp/self.pRef

        # error
        else: raise Exception('func: typ==%s is invalid' % typ)

        # return function value
        return f


    # Failure criteria names
    # ======================
    def names (self, typ):
        if   typ=='VM':    return 'von Mises'
        elif typ=='DP':    return 'Drucker/Prager'
        elif typ=='MC':    return 'Mohr/Coulomb'
        elif typ=='MCr':   return 'Mohr/Coulomb(rounded)'
        elif typ=='MN':    return 'Matsuoka/Nakai'
        elif typ=='MNnl':  return 'Matsuoka/Nakai(non-linear)'
        elif typ=='LD':    return 'Lade/Duncan'
        elif typ=='AS':    return 'Argyris/Sheng'
        elif typ=='AMP':   return 'Anisotropic'
        elif typ=='AMPb':  return r'$b=%g$'%self.b
        elif typ=='AMPba': return r'$b=%g,\,\alpha=%g$'%(self.b,self.alp)
        elif typ=='GE':    return 'Generic(b=%g)'%self.smpinvs.b
        else: raise Exception('failure_crit_names: typ==%s is invalid' % typ)


    # Rosette
    # =======
    # th   : show theta angles
    # ref  : reference lines
    # pos  : positive values
    # full : 360 degrees
    # fsz  : font size
    def rst (self, th=False, ref=False, full=True, pos=False, fsz=10):

        # radius
        r = self.sc*phi_calc_M(self.phi,'oct') if self.r==None else self.r

        # constants
        cr = 1.0
        cf = 0.2 if not full else cr
        l1 = (             0.0  , cr*r            ) # line: 1 end points
        l2 = (-cf*r*cos(pi/6.0) ,-cf*r*sin(pi/6.0)) # line: 2 end points
        l3 = ( cf*r*cos(pi/6.0) ,-cf*r*sin(pi/6.0)) # line: 3 end points
        l4 = (-cr*r*cos(pi/6.0) , cr*r*sin(pi/6.0)) # line: 4 = neg 1 end points
        lo = (-cr*r*cos(pi/3.0) , cr*r*sin(pi/3.0)) # line: origin of cylindrical system

        # main lines
        plot ([0.0,l1[0]],[0.0,l1[1]],'k-', color='grey', zorder=0)
        plot ([0.0,l2[0]],[0.0,l2[1]],'k-', color='grey', zorder=0)
        plot ([0.0,l3[0]],[0.0,l3[1]],'k-', color='grey', zorder=0)

        # reference
        plot ([0.0, l4[0]],[0.0, l4[1]],'--', color='grey', zorder=-1)
        if full:
            plot ([0.0,-l4[0]],[0.0, l4[1]],'--', color='grey', zorder=-1)
            plot ([0.0,   0.0],[0.0,-l1[1]],'--', color='grey', zorder=-1)
        if ref:
            plot ([0.0, lo[0]],[0.0, lo[1]],'--', color='grey', zorder=-1)
            if full:
                plot ([0.0,-lo[0]],[0.0, lo[1]],'--', color='grey', zorder=-1)
                plot ([-cr*r,cr*r],[0.0,0.0],   '--', color='grey', zorder=-1)

        # text
        if self.sxyz: k1,k2,k3 = 'z','y','x'
        else:         k1,k2,k3 = '1','3','2'
        if pos:
            if th: t1 = text(l1[0],l1[1],r'$\sigma_%s,\theta=+30^\circ$'%k1, ha='center', fontsize=fsz)
            else:  t1 = text(l1[0],l1[1],r'$\sigma_%s$'%k1,                  ha='center', fontsize=fsz)
            t2 = text(l2[0],l2[1],r'$\sigma_%s$'%k2,  ha='right',  fontsize=fsz)
            t3 = text(l3[0],l3[1],r'$\sigma_%s$'%k3,  ha='left',   fontsize=fsz)
        else:
            if th: t1 = text(l1[0],l1[1],r'$-\sigma_%s,\theta=+30^\circ$'%k1, ha='center', fontsize=fsz)
            else:  t1 = text(l1[0],l1[1],r'$-\sigma_%s$'%k1,                  ha='center', fontsize=fsz)
            t2 = text(l2[0],l2[1],r'$-\sigma_%s$'%k2,  ha='right',  fontsize=fsz)
            t3 = text(l3[0],l3[1],r'$-\sigma_%s$'%k3,  ha='left',   fontsize=fsz)
        self.rst_txt.append (t1)
        self.rst_txt.append (t2)
        self.rst_txt.append (t3)
        if th:
            t4 = text(lo[0],lo[1],r'$\theta=0^\circ$',   ha='center', fontsize=fsz)
            t5 = text(l4[0],l4[1],r'$\theta=-30^\circ$', ha='center', fontsize=fsz)
            self.rst_txt.append (t4)
            self.rst_txt.append (t5)


    # Get boundaries of active axes
    # ==============================
    def get_bounds (self, sf=0.05):
            xmin, xmax = gca().get_xbound()
            ymin, ymax = gca().get_ybound()
            Dx  , Dy   = xmax-xmin, ymax-ymin
            xmin, xmax = xmin-sf*Dx, xmax+sf*Dx
            ymin, ymax = ymin-sf*Dy, ymax+sf*Dy
            return xmin, xmax, ymin, ymax


    # Plot failure criteria
    # =====================
    def plot (self, typ='MN', clr='red', lst='-', lwd=1, label=None, np=40, leg_set=True, show_phi=False, fsz=10,
              plane='oct', pct_scale=0.1):

        # data
        r = self.sc*phi_calc_M(self.phi,'oct') if self.r==None else self.r
        x, y, f = zeros((np,np)), zeros((np,np)), zeros((np,np))

        # data for octahedral plane
        if plane=='oct':
            samin = -(1.0+pct_scale)*r if self.samin==None else self.samin
            samax =  (1.0+pct_scale)*r if self.samax==None else self.samax
            sbmin = -(1.0+pct_scale)*r if self.sbmin==None else self.sbmin
            sbmax =  (1.0+pct_scale)*r if self.sbmax==None else self.sbmax
            #samax += 0.3*(samax-samin)
            dsa = (samax-samin)/np
            dsb = (sbmax-sbmin)/np
            for i in range(np):
                for j in range(np):
                    x[i,j] = samin + i*dsa
                    y[i,j] = sbmin + j*dsb
                    if self.sxyz: s = oct_calc_sxyz (x[i,j], y[i,j], self.sc)
                    else:         s = oct_calc_s123 (x[i,j], y[i,j], self.sc)
                    sig    = matrix([[s[0]],[s[1]],[s[2]],[0.0]])
                    f[i,j] = self.func (sig, typ)
            if show_phi:
                fmt = '%g'
                s   = '$\phi_{comp}=' + fmt + '^\circ$'
                self.rst_txt.append (text(0.0, 3.0*sbmin, s%self.phi, fontsize=fsz, ha='center', va='top'))
                #self.rst_txt.append (text(0.5, 0.0, s%self.phi, fontsize=fsz, ha='center', va='top', transform=gca().transAxes))

        # data for p-q plane
        elif plane=='pq':
            pmin, pmax, qmin, qmax = self.get_bounds ()
            if not self.pmin==None: pmin = self.pmin
            if not self.pmax==None: pmax = self.pmax
            if not self.qmin==None: qmin = self.qmin
            if not self.qmax==None: qmax = self.qmax
            dp = (pmax-pmin)/np
            dq = (qmax-qmin)/np
            for i in range(np):
                for j in range(np):
                    x[i,j] = pmin + i*dp
                    y[i,j] = qmin + j*dq
                    s123   = pqt_calc_s123 (x[i,j], y[i,j], self.t, self.pqty)
                    sig    = matrix([[s123[0]],[s123[1]],[s123[2]],[0.0]])
                    f[i,j] = self.func (sig, typ)

        # data for s3-s1, s3-s2 plane
        elif plane=='dev':
            cte = -sqrt(3.0)*self.sc
            xmin, xmax, ymin, ymax = self.get_bounds (0.1)
            dx = (xmax-xmin)/np
            dy = (ymax-ymin)/np
            for i in range(np):
                for j in range(np):
                    x[i,j] = xmin + i*dx
                    y[i,j] = ymin + j*dy
                    s1     = (y[i,j] - 2.0*x[i,j] + cte)/3.0
                    s2     = (x[i,j] - 2.0*y[i,j] + cte)/3.0
                    s3     = (x[i,j] +     y[i,j] + cte)/3.0
                    sig    = matrix([[s1], [s2], [s3], [0.0]])
                    f[i,j] = self.func (sig, typ)

        # data for sxyz plane
        elif plane=='xy':
            xmin, xmax, ymin, ymax = self.get_bounds ()
            dx  = (xmax-xmin)/np
            dy  = (ymax-ymin)/np
            sq2 = sqrt(2.0)
            plot ([xmin,xmax], [xmin/sq2,xmax/sq2], 'k-') # hydrostatic line
            for i in range(np):
                for j in range(np):
                    x[i,j] = xmin + i*dx
                    y[i,j] = ymin + j*dy
                    sig    = matrix([[-x[i,j]/sq2], [-x[i,j]/sq2], [-y[i,j]], [0.0]])
                    f[i,j] = self.func (sig, typ)

        # wrong plane
        else: raise Exception('FCrits::plot: plane=%s is invalid'%plane)

        # contour
        contour (x,y,f, [0.0], colors=clr, linestyles=lst, linewidths=lwd)
        #gca().set_xticks([])
        #gca().set_yticks([])
        #gca().set_frame_on (False)
        axis ('equal')

        # legend
        if leg_set:
            self.leg_plt.append (plot ([0],[0], color=clr, linestyle=lst, linewidth=lwd))
            self.leg_txt.append (self.names (typ))


    # Legend
    # ======
    def leg (self, fsz=8, ncol=2):
        #return legend (self.leg_plt, self.leg_txt, bbox_to_anchor=(0,0,1,-1), loc=3, ncol=ncol, mode='expand',
                   #borderaxespad=0., handlelength=3, prop={'size':fsz})
        l = legend (self.leg_plt, self.leg_txt, bbox_to_anchor=(0.,1.02,1.,.102), loc=3, ncol=ncol, mode='expand',
                    borderaxespad=0., handlelength=3, prop={'size':fsz})
        if len(self.rst_txt)>0: return [l.legendPatch] + self.rst_txt
        else:                   return [l.legendPatch]


    # Plot 3D failure criteria
    # ========================
    def plot3d (self, typs=['MC','MN'], clr='red', np=40):

        # radius
        r = 1.*phi_calc_M(self.phi,'oct') if self.r==None else self.r

        # contour
        F = []
        for typ in typs: F.append (zeros((np,np,np)))
        sa    = zeros ((np,np,np))
        sb    = zeros ((np,np,np))
        sc    = zeros ((np,np,np))
        samin = -1.1*r       if self.samin==None else self.samin
        samax =  1.1*r       if self.samax==None else self.samax
        sbmin = -1.1*r       if self.sbmin==None else self.sbmin
        sbmax =  1.1*r       if self.sbmax==None else self.sbmax
        scmin =  0.0         if self.scmin==None else self.scmin
        scmax =  1.1*self.sc if self.scmax==None else self.scmax
        dsa   = (samax-samin)/np
        dsb   = (sbmax-sbmin)/np
        dsc   = (scmax-scmin)/np
        for i in range(np):
            for j in range(np):
                for k in range(np):
                    sa[i,j,k] = samin + i*dsa
                    sb[i,j,k] = sbmin + j*dsb
                    sc[i,j,k] = scmin + k*dsc
                    if self.sxyz: s = oct_calc_sxyz (sa[i,j,k], sb[i,j,k], sc[i,j,k])
                    else:         s = oct_calc_s123 (sa[i,j,k], sb[i,j,k], sc[i,j,k])
                    sig      = matrix([[s[0]],[s[1]],[s[2]],[0.0]])
                    for m, typ in enumerate(typs): F[m][i,j,k] = self.func (sig, typ)

        from enthought.mayavi.mlab import contour3d
        from enthought.mayavi.mlab import show as mlab_show

        for m, typ in enumerate(typs): contour3d (sa, sb, sc, F[m], contours=[0.0])
        mlab_show ()
