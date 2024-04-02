from msys_plt import *
#import mechsys 

voronoi = 0 # 0 => spheres

if voronoi: # voronoi
    qcamf      = 9.35
    pcamf      = 5.0
    p          = Plotter()
    p.show_k   = False
    p.fc_p     = pcamf*sqrt(3.0)
    p.fc_phi   = ms.M2Phi (qcamf/pcamf, "cam")
    p.fc_cu    = ms.qf2cu (qcamf,       "cam")
    p.fc_c     = 0.0
    p.fc_np    = 40
    p.log_p    = False
    p.div_by_p = True
    p.fc_ty    = ['MN', 'LD', 'MC', 'VM']
    p.pcte     = True
    p.set_eps  = True
    p.fc_np    = 100
    #p.maxed    = 9.0
    #p.maxev    = 4.0
    #p.maxidx   = 10
    p.justone = 3
    p.proport = 1.0 if p.justone==5 else 0.75
    #p.plot ("vorodense/ttt_c_walls.res",   label=r'$e = 0.4$',  marker='s',    markevery=15, clr='blue',   txtmax=True,  txtlst=True,  draw_fl=True,  draw_ros=True)
    if   p.justone==0: legend (loc='lower right', ncol=2)
    elif p.justone==1: legend (loc='lower right', ncol=2)
    elif p.justone==2: legend (loc='lower right', ncol=1)
    elif p.justone==3: legend (loc='upper left',  ncol=1)
    elif p.justone==4: legend (loc='upper left',  ncol=1)
    elif p.justone==5: p.fc_leg()
    elif p.justone==6: legend (loc='lower right', ncol=2)
    savefig("voro_%d.eps"%p.justone)
    print 'friction angle phi = ', p.fc_phi
    p.show()

else:
    qcamf      = 8.31
    pcamf      = 5.0
    p          = Plotter()
    p.show_k   = False
    p.fc_p     = pcamf*sqrt(3.0)
    p.fc_phi   = ms.M2Phi (qcamf/pcamf, "cam")
    p.fc_cu    = ms.qf2cu (qcamf,       "cam")
    p.fc_c     = 0.0
    p.fc_np    = 40
    p.log_p    = False
    p.div_by_p = True
    p.fc_ty    = ['MN', 'LD', 'MC', 'VM']
    p.pcte     = True
    p.justone  = -1
    p.fc_np    = 20
    p.plot ("ttt_c_walls.res",   clr='blue',   txtmax=True,  txtlst=True,  draw_fl=True,  draw_ros=True)
    #p.plot ("th_30/th_30_c_walls.res",   clr='blue',   txtmax=True,  txtlst=True,  draw_fl=True,  draw_ros=True)
    #p.plot ("th_15/th_15_c_walls.res",   clr='orange', txtmax=False, txtlst=False, draw_fl=False, draw_ros=False)
    #p.plot ("th_00/th_00_c_walls.res",   clr='black',  txtmax=False, txtlst=False, draw_fl=False, draw_ros=False)
    #p.plot ("th_-15/th_-15_c_walls.res", clr='red',    txtmax=False, txtlst=False, draw_fl=False, draw_ros=False)
    #p.plot ("th_-30/th_-30_c_walls.res", clr='cyan',   txtmax=False, txtlst=False, draw_fl=False, draw_ros=False)
    p.show ()
