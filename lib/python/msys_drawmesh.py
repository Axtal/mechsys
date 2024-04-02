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

from numpy import array, sqrt, linspace, pi, cos, sin, arctan2
from pylab import figure, text, show, axis, gca, gcf
from pylab import matplotlib as MPL
from msys_fig import GetLightClr, GetClr

class DrawMesh:

    # Constructor
    # ============
    #
    # V: vertices
    # C: cells
    #
    # V=[[ 0, -1,   0.0,  0.0],
    #    [ 1, -2,   1.0,  0.0],
    #    [ 2, -4,   0.0,  1.0],
    #    [ 3, -3,   1.0,  1.0]]
    #
    # C=[[ 0, -1, [0,1,3,2], {0:-10,1:-20,2:-30,3:-40}, {}]]
    #
    # pct: percentage of drawing limits to use for icons
    def __init__(self, V,C, Pins={}, pct=0.001, fsz1=8, fsz2=6,
                 yidpct=0.002, icf=0.04, acf=0.06, lcf=0.05, ndl=5, ndf=5, ndc=5,
                 rainbow=True):
        # mesh
        self.V      = V
        self.C      = C
        self.Pins   = Pins
        self.ndim   = len(V[0])-2
        if self.ndim>2: raise Exception("Drawing: NDim=%d is invalid"%self.ndim)

        # constants
        self.pct  = pct
        self.fsz1 = fsz1
        self.fsz2 = fsz2
        self.lw   = 0.5

        # use different colors according to cell tags?
        self.rainbow = rainbow

        # constants for arrows
        self.icf = icf # icons coefficient
        self.acf = acf # arrow coefficient 
        self.lcf = lcf # load coefficient
        self.ndl = ndl # num divisions for load icon (arrows)
        self.ndf = ndf # num divisions for flux icon (arrows)
        self.ndc = ndc # num divisions for convection icon (little S)

        # colors
        self.pink    = (250/255.0,204/255.0,228/255.0)
        self.lblue   = (217/255.0,228/255.0,255/255.0)
        self.lgreen  = (100/255.0,241/255.0,193/255.0)
        self.dblue   = ( 45/255.0,  0/255.0,160/255.0)
        self.orange  = (241/255.0,125/255.0,  0/255.0)
        self.lyellow = (234/255.0,228/255.0,179/255.0)
        #self.purple  = '#9b8de3'
        self.purple  = '#c5a9f3'
        self.dred    = '#b30000'
        self.mblue   = '#84bcdc'

        # assign colors
        self.celledgeclr = self.dblue
        self.lineedgeclr = self.dred

        # drawing limits (bounding box)
        allx      = [v[2] for v in self.V]
        if self.ndim==1: ally = [0.0, 1.0]
        else:            ally = [v[3] for v in self.V]
        self.lims = array([min(allx),max(allx),min(ally),max(ally)])
        self.diag = sqrt((self.lims[1]-self.lims[0])**2.0+(self.lims[3]-self.lims[2])**2.0)

        # noise to move tags and ids
        self.yidnoise = yidpct * self.diag

        # lengths
        self.il = self.icf * self.diag # icon's length
        self.al = self.acf * self.diag # arrows's length

        # matplotlib's structures
        self.PH = MPL.path.Path
        self.PP = MPL.patches
        self.PC = MPL.patches.PathPatch

    # Draw mesh
    #==========
    def draw(self, with_ids=True, with_tags=True,
             edgescells={}, cellscells=[], vertscells=[], vertsverts=[],
             p_vertscells=[], p_verts=[], only_lin_cells=False,
             with_grid=True, rotateIds=False, jointsR=None, lineedgeLws={}):
        # get figure
        fig = gcf()
        ax  = fig.add_subplot(111)
        if with_grid:
            ax.grid(color='gray')
            ax.set_axisbelow(True)

        # draw points at bounding box
        dlim  = sqrt((self.lims[1]-self.lims[0])**2.0+(self.lims[3]-self.lims[2])**2.0)
        limsx = self.lims + array([-self.pct*dlim,+self.pct*dlim,-self.pct*dlim,+self.pct*dlim]) # extended limits
        ax.plot(limsx[:2],limsx[2:],'o',marker='None')

        # draw solid cells
        if not only_lin_cells:
            for c in self.C:
                dat  = []
                con  = c[2] # connectivity
                nnod = len(con)
                if nnod>2:
                    x0 = self.V[con[0]][2]
                    y0 = self.V[con[0]][3]
                    dat.append((self.PH.MOVETO, (x0,y0)))
                    # edges for 3,4
                    if nnod<=4:
                        for j in range(1,nnod):
                            xj = self.V[con[j]][2]
                            yj = self.V[con[j]][3]
                            dat.append((self.PH.LINETO, (xj,yj)))
                    # edges for 6,8
                    if nnod==6:
                        dat.append((self.PH.LINETO, (self.V[con[3]][2], self.V[con[3]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[1]][2], self.V[con[1]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[4]][2], self.V[con[4]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[2]][2], self.V[con[2]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[5]][2], self.V[con[5]][3])))
                    if nnod==8 or nnod==9:
                        dat.append((self.PH.LINETO, (self.V[con[4]][2], self.V[con[4]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[1]][2], self.V[con[1]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[5]][2], self.V[con[5]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[2]][2], self.V[con[2]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[6]][2], self.V[con[6]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[3]][2], self.V[con[3]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[7]][2], self.V[con[7]][3])))
                    if nnod==15:
                        dat.append((self.PH.LINETO, (self.V[con[ 6]][2], self.V[con[ 6]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[ 3]][2], self.V[con[ 3]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[ 7]][2], self.V[con[ 7]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[ 1]][2], self.V[con[ 1]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[ 8]][2], self.V[con[ 8]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[ 4]][2], self.V[con[ 4]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[ 9]][2], self.V[con[ 9]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[ 2]][2], self.V[con[ 2]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[10]][2], self.V[con[10]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[ 5]][2], self.V[con[ 5]][3])))
                        dat.append((self.PH.LINETO, (self.V[con[11]][2], self.V[con[11]][3])))
                    dat.append((self.PH.CLOSEPOLY, (0,0)))
                if len(dat)>0:
                    cmd,vert = zip(*dat)
                    ph0 = self.PH (vert, cmd)
                    if self.rainbow:
                        clr = GetLightClr(abs(c[1]))
                        pc0 = self.PC (ph0, facecolor=clr, edgecolor='black', linewidth=self.lw, clip_on=False)
                    else:
                        pc0 = self.PC (ph0, facecolor=self.lblue, edgecolor=self.celledgeclr, linewidth=self.lw, clip_on=False)
                    ax.add_patch(pc0)

        # draw linear cells
        if self.ndim==1:
            for c in self.C:
                con = c[2] # connectivity
                if len(con)==2:
                    x0 = self.V[con[0]][2]
                    y0 = 0.0
                    x1 = self.V[con[1]][2]
                    y1 = 0.0
                    XY = array([[x0,y0],[x1,y1]])
                    if self.rainbow:
                        clr = GetClr(abs(c[1]))
                        ax.add_patch (MPL.patches.Polygon(XY, closed=False, edgecolor=clr, lw=2))
                    else:
                        ax.add_patch (MPL.patches.Polygon(XY, closed=False, edgecolor=self.celledgeclr, lw=2))
        else:
            for c in self.C:
                con = c[2] # connectivity
                if len(con)==2:
                    lw = 3
                    if c[1] in lineedgeLws: lw = lineedgeLws[c[1]]
                    x0 = self.V[con[0]][2]
                    y0 = self.V[con[0]][3]
                    x1 = self.V[con[1]][2]
                    y1 = self.V[con[1]][3]
                    XY = array([[x0,y0],[x1,y1]])
                    ax.add_patch (MPL.patches.Polygon(XY, closed=False, edgecolor=self.lineedgeclr, lw=lw))

        # text
        if with_ids or with_tags:
            for c in self.C:
                # centre
                if only_lin_cells and len(c[2])>2: continue
                cf = 0.5 if only_lin_cells else 0.3
                xc, yc, alp = self.get_cell_centre(c, cf)
                if not rotateIds: alp = 0.0
                if with_ids:
                    if with_tags: txt = '%d (%d)' % (c[0], c[1])
                    else:         txt = '%d'      %  c[0]
                else:             txt = '%d'      %        c[1]
                if self.rainbow:
                    ax.text(xc,yc, txt, rotation=alp, va='center', ha='center', backgroundcolor='white', fontsize=self.fsz2)
                else:
                    bclr = GetLightClr(abs(c[1]+1))
                    ax.text(xc,yc, txt, rotation=alp, va='center', ha='center', backgroundcolor=bclr, fontsize=self.fsz2)

        # edge tags
        if with_tags and self.ndim>1:
            for c in self.C:
                if len(c)>3: self.edge_tags(ax, c, only_tag=True)

        # draw nodes
        if with_ids:
            for v in self.V:
                s = '%d' % v[0]
                yval = v[3]+self.yidnoise if self.ndim>1 else self.yidnoise
                text(v[2], yval, s, va='bottom', ha='right', color='black', backgroundcolor=self.lyellow, fontsize=self.fsz1)

        # draw nodes
        if with_tags:
            for v in self.V:
                tag = v[1]
                if tag<0 and with_tags:
                    xval = v[2]+self.yidnoise
                    yval = v[3]-self.yidnoise if self.ndim>1 else -self.yidnoise
                    text(xval, yval, '%d'%tag, va='top', ha='left', color='black', backgroundcolor=self.orange, fontsize=self.fsz2)

        # vertscells
        for iv, cells in enumerate(vertscells):
            s = ''
            for k, ic in enumerate(cells):
                s += '%d' % ic
                if k < len(cells)-1: s += ' '
            xval = self.V[iv][2]-2.*self.yidnoise
            yval = self.V[iv][3]-self.yidnoise if self.ndim>1 else -self.yidnoise
            text(xval, yval, s, va='top', ha='right', color='black', backgroundcolor='white', fontsize=self.fsz2)

        # vertsverts
        for iv, verts in enumerate(vertsverts):
            s = ''
            for k, otherv in enumerate(verts):
                s += '%d' % otherv
                if k < len(verts)-1: s += ' '
            xval = self.V[iv][2]+3.*self.yidnoise
            yval = self.V[iv][3]+2.*self.yidnoise if self.ndim>1 else +2.*self.yidnoise
            text(xval, yval, s, va='bottom', ha='left', color='black', backgroundcolor='yellow', fontsize=self.fsz2)

        # cellscells
        for ic, cells in enumerate(cellscells):
            s = ''
            for k, otherc in enumerate(cells):
                s += '%d' % otherc
                if k < len(cells)-1: s += ' '
            xc, yc, _ = self.get_cell_centre(self.C[ic], 0.)
            xc -= 2.*self.yidnoise
            yc -= 2.*self.yidnoise
            text(xc, yc, s, va='top', ha='right', color='black', backgroundcolor='white', fontsize=self.fsz2)

        # edgescells
        for ed, cells in edgescells.items():
            xc = (self.V[ed[0]][2]+self.V[ed[1]][2])/2.
            yc = (self.V[ed[0]][3]+self.V[ed[1]][3])/2. if self.ndim>1 else -3.*self.yidnoise
            s  = ''
            for k, ic in enumerate(cells):
                s += '%d' % ic
                if k < len(cells)-1: s += ' '
            text(xc, yc, s, va='center', ha='center', color='black', backgroundcolor='white', fontsize=self.fsz2)

        # patch information
        for iv, cells in enumerate(p_vertscells):              # for each 'driver' vertex and cells around it
            if len(cells) == 0: continue                       # skip boundary vertex (without cells)
            s = ''
            for k, ic in enumerate(cells):
                s += '%d' % ic
                if k < len(cells)-1:
                    s += ' '
                    if k%3 == 2: s += '\n'
            s += '\n---\n'
            for k, iva in enumerate(p_verts[iv]):
                s += '%d' % iva
                if k < len(p_verts[iv])-1:
                    s += ' '
                    if k%3 == 2: s += '\n'
            xval = self.V[iv][2]+2.*self.yidnoise
            yval = self.V[iv][3]-self.yidnoise if self.ndim>1 else -self.yidnoise
            text(xval, yval, s, va='bottom', ha='left', color='black', backgroundcolor='lightgray', fontsize=self.fsz2)

        # joints
        if jointsR!=None:
            for v in self.V:
                gca().add_patch(self.PP.Circle((v[2],v[3]), jointsR, clip_on=False, edgecolor='black', facecolor='white', lw=1, zorder=10))

        # pins
        for key, val in self.Pins.items():
            # ids
            x = self.V[key][2]
            y = self.V[key][3]
            s = '(%d,' % key
            for i in range(len(val)):
                if i==len(val)-1: s += '%d)' % val[i]
                else:             s += '%d,' % val[i]
            text(x, y, s, va='bottom', ha='left', backgroundcolor=self.lyellow, fontsize=self.fsz1)
            # tag
            tag = self.V[key][1]
            if tag<0 and with_tags:
                text(x, y, '%d'%tag, va='top', color='black', backgroundcolor='none', fontsize=self.fsz2)

    # Draw edge tags
    # ==============
    def edge_tags(self, ax, c, ebcs={}, only_tag=False):
        con  = c[2]     # connectivity
        nnod = len(con) # number of nodes per cell
        for side, tag in c[3].items():
            if tag<0: # has tag
                if nnod==3:
                    if   side==0: na, nb = con[0], con[1]
                    elif side==1: na, nb = con[1], con[2]
                    elif side==2: na, nb = con[2], con[0]
                elif nnod==6:
                    if   side==0: na, nb,  nc, nd = con[0], con[3],  con[3], con[1]
                    elif side==1: na, nb,  nc, nd = con[1], con[4],  con[4], con[2]
                    elif side==2: na, nb,  nc, nd = con[2], con[5],  con[5], con[0]
                elif nnod==4:
                    if   side==0: na, nb = con[0], con[1]
                    elif side==1: na, nb = con[1], con[2]
                    elif side==2: na, nb = con[2], con[3]
                    elif side==3: na, nb = con[3], con[0]
                elif nnod==8 or nnod==9:
                    if   side==0: na, nb,  nc, nd = con[0], con[4],  con[4], con[1]
                    elif side==1: na, nb,  nc, nd = con[1], con[5],  con[5], con[2]
                    elif side==2: na, nb,  nc, nd = con[2], con[6],  con[6], con[3]
                    elif side==3: na, nb,  nc, nd = con[3], con[7],  con[7], con[0]
                elif nnod==15:
                    if   side==0: na,nb, nc,nd, ne,nf, ng,nh = con[0],con[ 6], con[ 6],con[3], con[3],con[ 7], con[ 7],con[1]
                    elif side==1: na,nb, nc,nd, ne,nf, ng,nh = con[1],con[ 8], con[ 8],con[4], con[4],con[ 9], con[ 9],con[2]
                    elif side==2: na,nb, nc,nd, ne,nf, ng,nh = con[2],con[10], con[10],con[5], con[5],con[11], con[11],con[0]
                xa = self.V[na][2]
                ya = self.V[na][3]
                xb = self.V[nb][2]
                yb = self.V[nb][3]
                xc = (xa+xb)/2.0
                yc = (ya+yb)/2.0
                if only_tag:
                    ax.text(xc,yc, '(%d)'%tag, ha='center', va='center', fontsize=self.fsz2, backgroundcolor=self.pink)
                    if nnod==6 or nnod==8 or nnod==9 or nnod==15:
                        xa = self.V[nc][2]
                        ya = self.V[nc][3]
                        xb = self.V[nd][2]
                        yb = self.V[nd][3]
                        xc = (xa+xb)/2.0
                        yc = (ya+yb)/2.0
                        ax.text(xc,yc, '(%d)'%tag, ha='center', va='center', fontsize=self.fsz2, backgroundcolor=self.pink)
                    if nnod==15:
                        xa = self.V[ne][2]
                        ya = self.V[ne][3]
                        xb = self.V[nf][2]
                        yb = self.V[nf][3]
                        xc = (xa+xb)/2.0
                        yc = (ya+yb)/2.0
                        ax.text(xc,yc, '(%d)'%tag, ha='center', va='center', fontsize=self.fsz2, backgroundcolor=self.pink)
                        xa = self.V[ng][2]
                        ya = self.V[ng][3]
                        xb = self.V[nh][2]
                        yb = self.V[nh][3]
                        xc = (xa+xb)/2.0
                        yc = (ya+yb)/2.0
                        ax.text(xc,yc, '(%d)'%tag, ha='center', va='center', fontsize=self.fsz2, backgroundcolor=self.pink)
                else:
                    pa, pb = array((xa,ya)), array((xb,yb))
                    for k, v in ebcs[tag].items():
                        self.draw_bc_over_line(ax, pa, pb, k, v)

    # Get cell centre
    # ===============
    def get_cell_centre(self, c, cf=0.5):
        con  = c[2] # connectivity
        nnod = len(con)
        if nnod==2:
            if self.ndim==1:
                x0 = self.V[con[0]][2]
                y0 = 0.0
                x1 = self.V[con[1]][2]
                y1 = 0.0
                xc = x0 + 0.5*(x1-x0)
                yc = y0 + 0.5*(y1-y0)
            else:
                x0 = self.V[con[0]][2]
                y0 = self.V[con[0]][3]
                x1 = self.V[con[1]][2]
                y1 = self.V[con[1]][3]
                xc = x0 + cf*(x1-x0)
                yc = y0 + cf*(y1-y0)
        else:
            xc = self.V[con[0]][2]
            yc = self.V[con[0]][3]
            for j in range(1,nnod):
                xc += self.V[con[j]][2]
                yc += self.V[con[j]][3]
            xc = xc/nnod
            yc = yc/nnod
        # angle between first and second nodes
        alp = 0.0
        if self.ndim > 1:
            x0, y0 = self.V[con[0]][2], self.V[con[0]][3]
            x1, y1 = self.V[con[1]][2], self.V[con[1]][3]
            dx, dy = x1-x0, y1-y0
            alp = arctan2(dy, dx) * 180./pi
            #if   x1 > x0: alp =  arccos((x1-x0)/l01) * 180./pi
            #elif x0 > x1: alp = -arccos((x0-x1)/l01) * 180./pi
        return xc, yc, alp

    # Show figure
    # ===========
    def show(self):
        axis('equal')
        show()

    # Draw node boundary conditions
    #==============================
    def node_bcs(self, vb, rotate=False, usetip=False, zorder=None, skiptags=[], tag2normal=None):
        ax = gca()
        for v in self.V:
            tag = v[1]
            if tag>=0: continue
            if tag in skiptags: continue
            B = vb[tag]
            # bcs
            T  = False
            ux = False
            uy = False
            wz = False
            fx = False
            fy = False
            sfx, sfy = 0.0, 0.0
            for key in B:
                if key=='T':  T  = True
                if key=='ux': ux = True
                if key=='uy': uy = True
                if key=='wz': wz = True
                if key=='fx': fx, sfx = True, self.sgn(B[key])
                if key=='fy': fy, sfy = True, self.sgn(B[key])
                if key=='T_func':
                    text(v[2],v[3],'f',fontsize=18,color=self.orange)
            # draw icon
            if T:
                n = array([0.0,1.0])
                ax.add_patch(self.PP.CirclePolygon(array([v[2],v[3]]),self.il/2.0,facecolor='none',edgecolor=self.orange,linewidth=2,zorder=zorder))
            if ux or uy:
                # direction
                dpx = abs(v[2]-self.lims[0])
                dpy = abs(v[3]-self.lims[2])
                dqx = abs(v[2]-self.lims[1])
                dqy = abs(v[3]-self.lims[3])
                fix = False
                if ux and not uy:                      # only ux
                    if dpx<0.01: n = array([ 1.0,0.0]) # left
                    else:        n = array([-1.0,0.0]) # ?
                elif uy and not ux:                    # only uy
                    if dpy<0.01: n = array([0.0, 1.0]) # bottom
                    else:        n = array([0.0,-1.0]) # ?
                else:                                  # both ux and uy
                    n = array([ 0.0, 1.0]) # ?
                    if rotate:
                        if   dpx<0.01 and dpy<0.01: n = array([ 0.0, 1.0]) # left-bottom
                        elif dpx<0.01 and dqy<0.01: n = array([ 0.0,-1.0]) # left-top
                        elif dqx<0.01 and dpy<0.01: n = array([ 0.0, 1.0]) # right-bottom
                        elif dqx<0.01 and dqy<0.01: n = array([ 0.0,-1.0]) # right-top
                        elif dpx<0.01:              n = array([ 1.0, 0.0]) # left
                        elif dqx<0.01:              n = array([-1.0, 0.0]) # right
                        elif dpy<0.01:              n = array([ 0.0, 1.0]) # bottom
                        elif dqy<0.01:              n = array([ 0.0,-1.0]) # top
                    fix = True
                if tag2normal!=None: n = tag2normal[tag]
                if wz: self.fixed_rotation (ax,array([v[2],v[3]]),self.il,n)
                else:  self.little_triangle(ax,array([v[2],v[3]]),self.il,n,fix)
            if fx or fy:
                dX = sfx * self.al
                dY = sfy * self.al
                if fx:
                    dx, dy = dX, 0.0
                    ax.add_patch(self.PP.Arrow(v[2],v[3],dx,dy,facecolor='black',edgecolor='none',width=0.5*dx,zorder=zorder))
                if fy:
                    dx, dy, yc = 0.0, dY, v[3]
                    if usetip: yc += abs(dy)
                    ax.add_patch(self.PP.Arrow(v[2],yc,dx,dy,facecolor='black',edgecolor='none',width=0.5*dy,zorder=zorder))

    # Draw edge boundary conditions
    #==============================
    def edge_bcs(self, eb, zorder=None):
        for c in self.C:
            self.edge_tags(gca(), c, eb, only_tag=False)

    # Draw parameters
    # ===============
    def params(self, params, fs=10, rotateIds=True, mn={}, sf=None, zorder=None):
        # mn : normal multiplier to displace text: map ids of cells to mn
        # find qnmax
        qnmax = 0.
        for _, P in params.items():
            for key, vals in P.items():
                if key=='qnqt':
                    qnmax = max(qnmax, abs(vals[0]), abs(vals[1]))
        # calc scale factor
        if sf==None:
            if qnmax > 0.0:
                Dx = self.lims[1]-self.lims[0]
                Dy = self.lims[3]-self.lims[2]
                dd = max([Dx,Dy])
                sf = 0.1 * dd / qnmax
            else: sf = 1.0
        # draw
        ax = gca()
        for c in self.C:
            l    = ''
            P    = params[c[1]]
            keys = P.keys()
            i0, i1 = c[2][0], c[2][1]
            pa, pb = array((self.V[i0][2],self.V[i0][3])), array((self.V[i1][2],self.V[i1][3]))
            for i, key in enumerate(keys):
                if key=='qnqt':
                    self.draw_bc_over_line(gca(), pa,pb,key,P[key],len(c[2])==2,sf,zorder)
                else:
                    #if P[key]>1.0e6: l += r'$%s=%g$' % (key,float(P[key]))
                    #else:            l += r'$%s=%s$' % (key,P[key])
                    l += r'$%s=%g$' % (key,P[key])
                    if i!=len(keys)-1: l+='\n'
            xc, yc, alp = self.get_cell_centre(c)
            if not rotateIds: alp = 0.0
            if c[0] in mn:
                dp  = pb-pa
                dL  = sqrt(dp[0]**2.0+dp[1]**2.0) # edge length
                n   = array([-dp[1],dp[0]])/dL    # unit normal (of beam)
                xc += mn[c[0]] * n[0]
                yc += mn[c[0]] * n[1]
            text(xc, yc, l, rotation=alp, fontsize=fs, ha='center',va='center')#, backgroundcolor='white')

    # Draw bry cond over line
    #========================
    def draw_bc_over_line(self, ax, pa,pb,key,val,beam=False,sf=1,zorder=None):
        dp  = pb-pa
        pm  = (pa+pb)/2.0
        dL  = sqrt(dp[0]**2.0+dp[1]**2.0) # edge length
        n   = array([dp[1],-dp[0]])/dL    # unit normal
        t   = array([ -n[1], n[0]])       # unit tangent
        ddl = 1.0/self.ndl
        ndl = int(dL/ddl)
        if ndl==0: ndl=2
        if beam: n, t = -n, -t
        if key=='qnqt':
            p = pa.copy()
            if beam: qnl, qnr, qt = val[0], val[1], val[2]
            else:    qnl, qnr, qt = val[0], val[0], val[1]
            if abs(qnl) > 0.0 or abs(qnr) > 0.0: #'qn'
                for i in range(ndl+1):
                    xc, yc = p[0], p[1]
                    qn = qnl + sqrt((p[0]-pa[0])**2.+(p[1]-pa[1])**2.) * (qnr-qnl) / dL
                    dx = sf * qn * n[0]
                    dy = sf * qn * n[1]
                    if self.sgn(qnl) < 0 or self.sgn(qnr) < 0:
                        ll = sqrt(dx**2.+dy**2.)
                        xc += n[0]*ll
                        yc += n[1]*ll
                    if zorder==None:
                        ax.add_patch(self.PP.Arrow(xc,yc,dx,dy,width=0.4*self.al,lw=1,facecolor=self.mblue,edgecolor='None'))
                    else:
                        ax.add_patch(self.PP.Arrow(xc,yc,dx,dy,width=0.4*self.al,lw=1,facecolor=self.mblue,edgecolor='None',zorder=zorder))
                    p += dp/ndl
            if abs(qt) > 0.0: #'qt'
                dx = self.lcf*self.sgn(qt)*self.il*t[0]
                dy = self.lcf*self.sgn(qt)*self.il*t[1]
        if key=='conv':
            p  = pa.copy()
            for i in range(self.ndc+1):
                self.little_S(ax,p,self.il,n)
                p += dp/self.ndc
        if key=='flux':
            if fabs(v)<1.0e-7: # insulated
                p0 = pa+(0.3*self.il)*n
                p1 = pb+(0.3*self.il)*n
                xy = array([[p0[0],p0[1]],[p1[0],p1[1]]])
                ax.add_patch(self.PP.Polygon(xy,edgecolor=self.orange,linewidth=4))
                p0 += (0.3*self.il)*n
                p1 += (0.3*self.il)*n
                xy  = array([[p0[0],p0[1]],[p1[0],p1[1]]])
                ax.add_patch(self.PP.Polygon(xy,edgecolor=self.orange,linewidth=4))
            else:
                p = pa.copy()
                for i in range(self.ndf+1):
                    if self.sgn(v)>0 and i==0: p += self.lcf*self.il*n
                    dx = -self.lcf*self.sgn(v)*self.il*n[0]
                    dy = -self.lcf*self.sgn(v)*self.il*n[1]
                    ax.add_patch(self.PP.Arrow(p[0],p[1],dx,dy,facecolor=(241/255.0,125/255.0,0/255.0),edgecolor='none',width=0.7*self.l))
                    p += dp/self.ndf

    # Draw little S
    # =============
    def little_S(self, ax, p0,l,n):
        hh  = l*sqrt(3.0)/4.0     # half-height
        t   = array([-n[1],n[0]]) # unit tangent
        pm  = p0 + hh*n
        p1  = pm + (l/2.0)*t
        p2  = pm - (l/2.0)*t
        p3  = p0 + 2.0*hh*n
        dat = [(self.PH.MOVETO,(p0[0],p0[1])), (self.PH.CURVE4,(p1[0],p1[1])), (self.PH.CURVE4,(p2[0],p2[1])), (self.PH.CURVE4,(p3[0],p3[1]))]
        c,v = zip(*dat)
        ax.add_patch(self.PC(self.PH(v,c), facecolor='none', edgecolor=(241/255.0,125/255.0,0/255.0), linewidth=4))

    # Draw fixed rotation
    # ===================
    def fixed_rotation(self, ax, p0,l,n):
        h   = 0.0                 # height
        t   = array([-n[1],n[0]]) # unit tangent
        pm  = p0 - h*n
        p1  = pm + (l/2.0)*t
        p2  = pm - (l/2.0)*t
        #dat = [(self.PH.MOVETO,(p0[0],p0[1]))]#, (self.PH.LINETO,(p1[0],p1[1])), (self.PH.LINETO,(p2[0],p2[1])), (self.PH.CLOSEPOLY,(0,0))]
        dat = [(self.PH.MOVETO,(p1[0],p1[1])), (self.PH.LINETO,(p2[0],p2[1]))]
        c,v = zip(*dat)
        ax.add_patch(self.PC(self.PH(v,c), edgecolor='black', linewidth=2))
        dat = []
        # ticks
        ld  = l/3.0 # length of details
        nd = 7
        dx = l/nd
        p  = p1
        q  = p - ld*n
        dat.append((self.PH.MOVETO,(p[0],p[1])))
        dat.append((self.PH.LINETO,(q[0],q[1])))
        for i in range(nd):
            p -= dx*t
            q  = p - (l/3.0)*n
            dat.append((self.PH.MOVETO,(p[0],p[1])))
            dat.append((self.PH.LINETO,(q[0],q[1])))
        # create patch
        c,v = zip(*dat)
        ax.add_patch(self.PC(self.PH(v,c), edgecolor='black', linewidth=1))

    # Draw little triangle
    # ====================
    # p0: tip position   (array)
    # l:  side length    (array)
    # n:  unit direction (array)
    def little_triangle(self, ax, p0,l,n,fixed=False):
        h   = l*sqrt(3.0)/2.0     # height
        t   = array([-n[1],n[0]]) # unit tangent
        pm  = p0 - h*n
        p1  = pm + (l/2.0)*t
        p2  = pm - (l/2.0)*t
        dat = [(self.PH.MOVETO,(p0[0],p0[1])), (self.PH.LINETO,(p1[0],p1[1])), (self.PH.LINETO,(p2[0],p2[1])), (self.PH.CLOSEPOLY,(0,0))]
        ld  = l/3.0 # length of details
        if fixed: # fixed support
            nd = 7
            dx = l/nd
            p  = p1
            q  = p - ld*n
            dat.append((self.PH.MOVETO,(p[0],p[1])))
            dat.append((self.PH.LINETO,(q[0],q[1])))
            for i in range(nd):
                p -= dx*t
                q  = p - (l/3.0)*n
                dat.append((self.PH.MOVETO,(p[0],p[1])))
                dat.append((self.PH.LINETO,(q[0],q[1])))
        else:
            pc = pm-(ld/2.0)*n
            p3 = p1 - ld*n
            p4 = p2 - ld*n
            xy = array([[p3[0],p3[1]],[p4[0],p4[1]]])
            ax.add_patch(self.PP.CirclePolygon(pc,ld/2.0,facecolor='none'))
            ax.add_patch(self.PP.Polygon(xy,linewidth=2))
        # create patch
        c,v = zip(*dat)
        ax.add_patch(self.PC(self.PH(v,c), facecolor=(250/255.0,204/255.0,95/255.0), edgecolor='black', linewidth=1))

    # Sign funcion
    # ============
    def sgn(self, val):
        if val<0: return -1
        else:     return  1
