// Moving cube in x-axis.


// MechSys
#include <mechsys/lbmmpm/Domain.h>
#include <math.h>
#include <iostream>
#include <fstream>

struct UserData
{
    double                         * Vel;
    double                           rho;
    double                          rhof;
    double                            Tf;
    double                            Ly;
    double                            Lz;
    double                            nu;
    double                            CD;
    double                            Re;
    double                            Fb;
    Vec3_t                    fluidforce;
    //Array<MPM::Particle *>    EndBeamPar;
    //Array<MPM::Particle *>      ForcePar;
    Array<double >                    x0;
    Array<double >                    x2; 
    std::ofstream                oss_ss1;

};

void Setup (LBMMPM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    ///*
    /*#pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    for (size_t ix=0; ix<dom.LBMDOM.Ndim(0); ++ix)
    for (size_t iy=0; iy<dom.LBMDOM.Ndim(1); ++iy)
    for (size_t iz=0; iz<dom.LBMDOM.Ndim(2); ++iz)
    {
        dom.LBMDOM.BForce[0][ix][iy][iz] = dom.LBMDOM.Rho[0][ix][iy][iz]*dat.fluidforce;
    }*/   
    //#pragma omp parallel for schedule(static) num_threads(dom.Nproc)   
    //for (size_t ip=0; ip < dom.MPMDOM.Particles.Size(); ip++)
    //{ 
        //dom.MPMDOM.Particles[ip]->vf = Vec3_t(0.01,0.0,0.0); 
        //dom.MPMDOM.Particles[ip]->vf = 0.0001;
    //} 

    //std::cout << dom.MPMDOM.Particles[0]->x << std::endl;

}

void Report (LBMMPM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    //std::cout << dom.MPMDOM.Particles[3]->h(0) << " " << dom.MPMDOM.Particles[3]-> h(1) << " " << dom.MPMDOM.Particles[3]->h(2) << std::endl;
    //std::cout << dom.MPMDOM.Particles[0]->h(0) << " " << dom.MPMDOM.Particles[0]-> h(1) << " " << dom.MPMDOM.Particles[0]->h(2) << std::endl;
    if (dom.idx_out==0)
    {
        String fs1;
        fs1.Printf("tlbmmpm02.res");
        dat.oss_ss1.open(fs1.CStr());
        //dat.oss_ss1 << Util::_10_6 << "Time" << Util::_8s << "Re" << Util::_8s << "CD" << Util::_8s << "Dx/b" << Util::_8s << "Dz/b" << "\n";
        dat.oss_ss1 << Util::_10_6 << "Time" << Util::_8s << "Re" << Util::_8s << "CD" << "\n";
    }
    else 
    {
        //double Ly    = 10.0;
        double x     = 0.0;
        double z     = 0.0;
        double xf    = 0.0; // Dx/b
        double zf    = 0.0; // Dz/b
        double vel   = 0.0;
        double force = 0.0;
        int    i     = 0;
        int    j     = 0;

        // calculate Re 
        //for (size_t ix=0; ix<dom.LBMDOM.Ndim(0); ++ix)
        //for (size_t iy=0; iy<dom.LBMDOM.Ndim(1); ++iy)
        //for (size_t iz=0; iz<dom.LBMDOM.Ndim(2); ++iz)
        //{
            //vel += dom.LBMDOM.Vel[0][ix][iy][iz](0);
            //i = i+1;
        //}
        //vel /= i;
        vel = dom.MPMDOM.Particles[0]->v(0);
        std::cout << "vel1 is" << " " << vel << std::endl;
        double Re   = 1.2407*vel*dat.Ly/dat.nu;
        //double Re   = dat.rhof*vel*dat.Ly/dat.nu;
        std::cout << "Re is" << " " << Re << std::endl;

        // calculate CD
        for (size_t ip=0; ip < dom.MPMDOM.Particles.Size(); ip++)
        {
            force = force + dom.MPMDOM.Particles[ip]->h(0);
            //std::cout << "force is" << " " << force << std::endl;

        }
        
        double CD = force/(0.5*dat.rhof*vel*vel*dat.Ly*dat.Lz*1.2090);
        //std::cout << "vel is" << " " << vel << std::endl;
        
        
        std::cout << "CD is" << " " << CD << std::endl;

        dat.oss_ss1 << Util::_10_6 << dom.Time << Util::_8s << Re << Util::_8s << CD << "\n";
    }
}

int main(int argc, char **argv) try
{
    //Number of cores
    /*if (argc<2) throw new Fatal("This program must be called with one argument: the name of the data input file without the '.inp' suffix.\nExample:\t %s filekey\n",argv[0]);
    size_t Nproc = 1; 
    if (argc>=3) Nproc=atoi(argv[2]);
    String filekey  (argv[1]);
    String filename (filekey+".inp");
    if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());
    ifstream infile(filename.CStr());
    double bforce;
    {
        infile >> bforce;           infile.ignore(200,'\n');
    }*/

    size_t Nproc = 0.75*omp_get_max_threads();
    if (argc>1) Nproc = atoi(argv[1]);
    
    double bforce = 0.0;
    //Properties of LBM
    size_t nx = 241;
    size_t ny = 61;
    size_t nz = 61;
    //size_t nx = 420; //number of cells in x-axis for benchmark
    //size_t ny = 320;
    //size_t nz = 340;
    //double dx = 0.5; //benchmark value
    double dx = 0.1;
    double rhof = 1000.0; // density of fluid

    //Properties of MPM
    //size_t ndiv = 15; //number of divisions per x lenght (meshes number in x-axis)
    size_t ndiv = 1; //number of divisions per x lenght
    double K    = 10.0e4; //Bulk modulus
    double Nu   = 0.3; //Poisson ratio
    double  E   = (3.0*(1.0-2.0*Nu))*K; //Young modulus
    double  G   = E/(2.0*(1.0+Nu)); //Shear modulus
    double rho  = 694.95; //density of solid
    //double Lx   = 30.0;   
    //double Ly   = 6.0;
    //double Lz   = 6.0;
    double Lx   = 14.5*dx; //length of the beam in x   
    double Ly   = 14.5*dx;
    double Lz   = 14.5*dx;
    double Cs   = sqrt(E/rho); //Speed of sound
    double h    = Ly/ndiv; 
    double h2   = Lz/ndiv; //length per mesh
    double dt   = 0.2;
    //double Dx   = 2.0*Lx/ndiv;
    double Dx   = 2.0*Lz/ndiv;
    double Bc   = Dx;
    double nu   = dx*dx/(dt*2400.0); //viscosity
    
    

    LBMMPM::Domain dom(D3Q15,nu,iVec3_t(nx,ny,nz),dx,dt);
    UserData dat;
    dom.UserData = &dat;

    dat.Ly = Ly;
    dat.Lz = Lz;
    dat.nu = nu;
    dat.rhof = rhof;

    double vel = 0.001;
    

    //Add beam to the middle of LBM cells
    //dom.MPMDOM.AddRectangularBeamMesh(-1, Vec3_t(0.5*nx*dx-0.5*Lx,0.5*ny*dx-0.5*Ly,0.0), Vec3_t(0.5*nx*dx+0.5*Lx,0.5*ny*dx+0.5*Ly,Lz), rho, ndiv);
    //dom.MPMDOM.AddRectangularBeamMesh(-1, Vec3_t(0.5*nx*dx-0.5*Lx,0.5*ny*dx-0.5*Ly,0.5*nz*dx-0.5*Lz), Vec3_t(0.5*nx*dx+0.5*Lx,0.5*ny*dx+0.5*Ly,0.5*nz*dx+0.5*Lz), rho, ndiv);
    dom.MPMDOM.AddRectangularBeamMesh(-1, Vec3_t(0.3*nx*dx-0.5*Lx,0.5*ny*dx-0.5*Ly,0.5*nz*dx-0.5*Lz), Vec3_t(0.3*nx*dx+0.5*Lx,0.5*ny*dx+0.5*Ly,0.5*nz*dx+0.5*Lz), rho, ndiv);
    //dom.MPMDOM.Particles.Push(new MPM::Particle(-1,Vec3_t(0.5*dx*nx,0.5*dx*ny,0.5*dx*nz), OrthoSys::O, 3000.0*pow(4.0*dx,3.0), pow(4.0*dx,3.0)));
    //dom.MPMDOM.Particles[0]->FixVeloc(0.1*dx/dt,0.0,0.0);
    //dom.MPMDOM.ResizeDomain(Vec3_t(0.5*nx*dx-5.0*Lx,0.5*ny*dx-5.0*Ly,-2.0*Lz), Vec3_t(0.5*nx*dx+5.0*Lx,0.5*ny*dx+5.0*Ly,5.0*Lz),Dx);
    //dom.MPMDOM.ResizeDomain(Vec3_t(0.5*nx*dx-5.0*Lx,0.5*ny*dx-5.0*Ly,0.5*nz*dx-5.0*Lz), Vec3_t(0.5*nx*dx+5.0*Lx,0.5*ny*dx+5.0*Ly,0.5*nz*dx+5.0*Lz),Dx);
    //dom.MPMDOM.ResizeDomain(Vec3_t(0.3*nx*dx-8.0*Lx,0.5*ny*dx-8.0*Ly,0.5*nz*dx-8.0*Lz), Vec3_t(0.3*nx*dx+8.0*Lx,0.5*ny*dx+8.0*Ly,0.5*nz*dx+8.0*Lz),Dx);
    dom.MPMDOM.ResizeDomain(Vec3_t(0.0,0.0,0.0), Vec3_t(nx*dx,ny*dx,nz*dx),Dx);
    
    //Setting properties for the material points
    //double v0 = 1.0;
    for (size_t ip=0; ip < dom.MPMDOM.Particles.Size(); ip++)
    {
        dom.MPMDOM.Particles[ip]->K = K;
        dom.MPMDOM.Particles[ip]->G = G;
        dom.MPMDOM.Particles[ip]->Tag = -2;
        dom.MPMDOM.Particles[ip]->FixVeloc(vel,0.0,0.0);
    }
    dom.MPMDOM.Gn = 0.0e0;
    //force /= dom.MPMDOM.Particles.Size();
    //std::cout << "force is" << force << std::endl;

    //Setting properties for nodes
    for (size_t in=0; in < dom.MPMDOM.Nnodes; in++)
    {
        Vec3_t xn;
        dom.MPMDOM.NodePosition(in,xn);
        dom.MPMDOM.Nodes[in].FixVeloc(vel,0.0,0.0);
    }

    //Setting intial conditions of fluid
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        //Vec3_t v(0.0,0.0,0.0);
        Vec3_t v(0.0,0.0,0.0);
        iVec3_t idx(ix,iy,iz);
        dom.LBMDOM.Initialize(0,idx,1000.0/*rho*/,v);
        if (ix==0||ix==nx-1||iy==0||iy==ny-1||iz==0||iz==nz-1)
        {
            dom.LBMDOM.IsSolid[0][ix][iy][iz] = true;
        }
        /*if (iz==0||iz==nz-1)
        {
            dom.LBMDOM.IsSolid[0][ix][iy][iz] = true;
        }*/
    }  

    dat.Vel = new double[ny];
    //double umax    = 1.0e-3*dx/dt;
    dat.rho = 1000.0;
    
    double Tf = 1.0e5*dt; // Final Time
    //double Tf = dt; // Final Time
    dat.Tf = Tf;
    //dat.fluidforce = Vec3_t(6.0e-2*dx/dt/Tf,0.0,0.0);
    dat.fluidforce = Vec3_t(bforce*dx/dt/Tf,0.0,0.0);
    dat.Fb = 0.01; 
    
    //dom.MPMDOM.BoundaryMesh();
    //dom.Reset();
    //dom.ImprintLattice();
    dom.Solve(Tf,Tf/100,Setup,Report,"tlbmmpm03",true,Nproc);
}
MECHSYS_CATCH
