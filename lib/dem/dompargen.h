#ifndef MECHSYS_DEM_DOMAIN_PARGEN_H
#define MECHSYS_DEM_DOMAIN_PARGEN_H

// Particle generation

inline void Domain::GenSpheres (int Tag, double L, size_t N, double rho,char const * Type, size_t Randomseed, double fraction, double RminFraction)
{
    // find radius from the edge's length
    Util::Stopwatch stopwatch;
    printf("\n%s--- Generating packing of spheres -----------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    srand(Randomseed);
    double R = L/(2.0*N);
    if (strcmp(Type,"Normal")==0)
    {
        for (size_t n=0; n<N*N*N; ++n) 
        {
            Vec3_t pos(-L/2.0+R, -L/2.0+R, -L/2.0+R);
            size_t i = (n%N);
            size_t j = (n/N)%N;
            size_t k = (n/(N*N));
            pos += Vec3_t(2.0*i*R, 2.0*j*R, 2.0*k*R);
            if (rand()<fraction*RAND_MAX) AddSphere (Tag,pos,R*RminFraction/(1.0-(1.0*rand())/RAND_MAX*(1.0-RminFraction)),rho);
        }
    }
    else if (strcmp(Type,"HCP")==0)
    {
        size_t nx = N;
        size_t ny = int(L/(sqrt(3.0)*R));
        size_t nz = int(L/(sqrt(8.0/3.0)*R));
        for (size_t k = 0; k < nz; k++)
        {
            for (size_t j = 0; j < ny; j++)
            {
                Vec3_t X;
                if (k%2==0) X = Vec3_t(-2*R-L/2.0,R-L/2.0,2*R-L/2.0+k*sqrt(8.0/3.0)*R);
                else X = Vec3_t(-R-L/2.0,R+sqrt(1.0/3.0)*R-L/2.0,2*R-L/2.0+k*sqrt(8.0/3.0)*R);
                if (j%2==0) X += Vec3_t(R,j*sqrt(3.0)*R,0.0);
                else X += Vec3_t(0.0,j*sqrt(3.0)*R,0.0);
                for (size_t i = 0; i < nx; i++)
                {
                    X += Vec3_t(2*R,0.0,0.0);
                    if (rand()<fraction*RAND_MAX) AddSphere(Tag,X,R*RminFraction/(1.0-(1.0*rand())/RAND_MAX*(1.0-RminFraction)),rho);
                }
            }
        }
    }
    else throw new Fatal ("Right now there are only two possible packings available the Normal and the HCP, packing %s is not implemented yet",Type);
    printf("%s  Num of particles   = %zd%s\n",TERM_CLR2,Particles.Size(),TERM_RST);
}

inline void Domain::GenSpheresBox (int Tag, Vec3_t const & X0, Vec3_t const & X1, double R, double rho, char const * Type, size_t Randomseed, double fraction, double RminFraction)
{
    // find radius from the edge's length
    Util::Stopwatch stopwatch;
    printf("\n%s--- Generating packing of spheres -----------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    srand(Randomseed);
    if (strcmp(Type,"Normal")==0)
    {
        size_t nx = 0.5*(X1(0)-X0(0))/R;
        size_t ny = 0.5*(X1(1)-X0(1))/R;
        size_t nz = 0.5*(X1(2)-X0(2))/R;
        for (size_t i = 0; i < nx; i++)
        for (size_t j = 0; j < ny; j++)
        for (size_t k = 0; k < nz; k++)
        {
            //Vec3_t pos(-(X1(0)-X0(0))/2.0+R, -(X1(1)-X0(1))/2.0+R, -(X1(2)-X0(2))/2.0+R);
            Vec3_t pos(X0(0)+R,X0(1)+R,X0(2)+R);
            pos += Vec3_t(2.0*i*R, 2.0*j*R, 2.0*k*R);
            if (rand()<fraction*RAND_MAX) AddSphere (Tag,pos,R*RminFraction/(1.0-(1.0*rand())/RAND_MAX*(1.0-RminFraction)),rho);
        }
    }

    else if (strcmp(Type,"HCP")==0)
    {
        size_t nx = 0.5*(X1(0)-X0(0))/R-1;
        size_t ny = int((X1(1)-X0(1))/(sqrt(3.0)*R));
        size_t nz = int((X1(2)-X0(2))/(sqrt(8.0/3.0)*R));
        for (size_t k = 0; k < nz; k++)
        {
            for (size_t j = 0; j < ny; j++)
            {
                Vec3_t X;
                if (k%2==0) X = Vec3_t(-R,R,2*R+k*sqrt(8.0/3.0)*R) + X0;
                else X = Vec3_t(0.0,R+sqrt(1.0/3.0)*R,2*R+k*sqrt(8.0/3.0)*R) + X0;
                if (j%2==0) X += Vec3_t(R,j*sqrt(3.0)*R,0.0);
                else X += Vec3_t(0.0,j*sqrt(3.0)*R,0.0);
                for (size_t i = 0; i < nx; i++)
                {
                    X += Vec3_t(2*R,0.0,0.0);
                    //std::cout << X << X0 << X1 << std::endl;
                    if ((X(0)<X0(0))||(X(0)>X1(0))) continue;
                    if ((X(1)<X0(1))||(X(1)>X1(1))) continue;
                    if ((X(2)<X0(2))||(X(2)>X1(2))) continue;
                    //if (rand()<fraction*RAND_MAX) AddSphere (Tag,X,R*RminFraction+(1.0*rand())/RAND_MAX*(R-R*RminFraction),rho);
                    if (rand()<fraction*RAND_MAX) AddSphere (Tag,X,R*RminFraction/(1.0-(1.0*rand())/RAND_MAX*(1.0-RminFraction)),rho);
                }
            }
        }
    }
    
    else throw new Fatal ("Right now there are only two possible packings available the Normal and the HCP, packing %s is not implemented yet",Type);
    printf("%s  Num of particles   = %zd%s\n",TERM_CLR2,Particles.Size(),TERM_RST);
}

inline void Domain::GenRice (int Tag, double L, size_t N, double R, double rho, size_t Randomseed, double fraction)
{
    Util::Stopwatch stopwatch;
    printf("\n%s--- Generating packing of 'rices' -----------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    srand(Randomseed);
    double dL = L/N;
    for (size_t n=0; n<N*N*N; ++n) 
    {
        Vec3_t pos(-L/2.0+dL, -L/2.0+dL, -L/2.0+dL);
        size_t i = (n%N);
        size_t j = (n/N)%N;
        size_t k = (n/(N*N));
        pos += Vec3_t(2.0*i*dL, 2.0*j*dL, 2.0*k*dL);
        if (rand()<fraction*RAND_MAX) AddRice (Tag, pos, R, dL-2*R, rho);
    }
    printf("%s  Num of particles   = %zd%s\n",TERM_CLR2,Particles.Size(),TERM_RST);
}

inline void Domain::GenBox (int InitialTag, double Lx, double Ly, double Lz, double R, double Cf, bool Cohesion)
{
    /*                         +----------------+
     *                       ,'|              ,'|
     *                     ,'  |  ___       ,'  |
     *     z             ,'    |,'4,'  [1],'    |
     *     |           ,'      |~~~     ,'      |
     *    ,+--y      +'===============+'  ,'|   |
     *  x'           |   ,'|   |      |   |2|   |
     *               |   |3|   |      |   |,'   |
     *               |   |,'   +- - - | +- - - -+
     *               |       ,'       |       ,'
     *               |     ,' [0]  ___|     ,'
     *               |   ,'      ,'5,'|   ,'
     *               | ,'        ~~~  | ,'
     *               +----------------+'
     */

    
    // add faces of box
    Vec3_t axis0(OrthoSys::e0); // rotation of face
    Vec3_t axis1(OrthoSys::e1); // rotation of face
    size_t IIndex = Particles.Size();  // First face index
    AddPlane (InitialTag,   Vec3_t(Lx/2.0,0.0,0.0),  R, Cf*Lz, Cf*Ly, 1.0, M_PI/2.0, &axis1);
    Particles[Particles.Size()-1]->Initialize(Particles.Size()-1);
    AddPlane (InitialTag-1, Vec3_t(-Lx/2.0,0.0,0.0), R, Cf*Lz, Cf*Ly, 1.0, 3.0*M_PI/2.0, &axis1);
    Particles[Particles.Size()-1]->Initialize(Particles.Size()-1);
    AddPlane (InitialTag-2, Vec3_t(0.0,Ly/2.0,0.0),  R, Cf*Lx, Cf*Lz, 1.0, 3.0*M_PI/2.0, &axis0);
    Particles[Particles.Size()-1]->Initialize(Particles.Size()-1);
    AddPlane (InitialTag-3, Vec3_t(0.0,-Ly/2.0,0.0), R, Cf*Lx, Cf*Lz, 1.0, M_PI/2.0, &axis0);
    Particles[Particles.Size()-1]->Initialize(Particles.Size()-1);
    AddPlane (InitialTag-4, Vec3_t(0.0,0.0,Lz/2.0),  R, Cf*Lx, Cf*Ly, 1.0);
    Particles[Particles.Size()-1]->Initialize(Particles.Size()-1);
    AddPlane (InitialTag-5, Vec3_t(0.0,0.0,-Lz/2.0), R, Cf*Lx, Cf*Ly, 1.0, M_PI, &axis0);
    Particles[Particles.Size()-1]->Initialize(Particles.Size()-1);

    // define some tolerance for comparissions
    if (Cohesion)
    {
        double tol1 = 1.0e-8;
        double tol2 = 1.0e-3;
        for (size_t i=0;i<IIndex;i++)
        {
            Particle * P1 = Particles[i];
            for (size_t j=IIndex;j<Particles.Size();j++)
            {
                Particle * P2 = Particles[j];
                for (size_t k=0;k<P1->Faces.Size();k++)
                {
                    Face * F1 = P1->Faces[k];
                    Vec3_t n1,c1;
                    F1->Normal  (n1);
                    F1->Centroid(c1);
                    Face * F2 = P2->Faces[0];
                    Vec3_t n2,c2;
                    F2->Normal  (n2);
                    F2->Centroid(c2);
                    Vec3_t n = 0.5*(n1-n2);
                    n/=norm(n);
                    if ((fabs(dot(n1,n2)+1.0)<tol1)
                       &&(fabs(dot(c2-c1,n)-2*R)<tol2))
                    {
                        BInteractons.Push(new BInteracton(P1,P2,k,1));
                        break;
                    }
                }
            }        
        }
    }
}

inline void Domain::GenBox (int InitialTag, Vec3_t & Xmin, Vec3_t & Xmax, double R, double Cf, bool Cohesion)
{
    double Lx = Xmax(0) - Xmin(0);
    double Ly = Xmax(1) - Xmin(1);
    double Lz = Xmax(2) - Xmin(2);
    GenBox(InitialTag,Lx,Ly,Lz,R,Cf,Cohesion);
    Vec3_t trans = 0.5*(Xmax+Xmin);
    Particles[Particles.Size()-1]->Translate(trans);
    Particles[Particles.Size()-2]->Translate(trans);
    Particles[Particles.Size()-3]->Translate(trans);
    Particles[Particles.Size()-4]->Translate(trans);
    Particles[Particles.Size()-5]->Translate(trans);
    Particles[Particles.Size()-6]->Translate(trans);
}

inline void Domain::GenOpenBox (int InitialTag, double Lx, double Ly, double Lz, double R, double Cf)
{
    /*                         +----------------+
     *                       ,'|              ,'|
     *                     ,'  |  ___       ,'  |
     *     z             ,'    |,'N,'  [1],'    |
     *     |           ,'      |~~~     ,'      |
     *    ,+--y      +'===============+'  ,'|   |
     *  x'           |   ,'|   |      |   |2|   |
     *               |   |3|   |      |   |,'   |
     *               |   |,'   +- - - | +- - - -+
     *               |       ,'       |       ,'
     *               |     ,' [0]  ___|     ,'
     *               |   ,'      ,'4,'|   ,'
     *               | ,'        ~~~  | ,'
     *               +----------------+'
     */


    // Creates an open box without the top lid, acts as a container
    
    // add faces of box
    Vec3_t axis0(OrthoSys::e0); // rotation of face
    Vec3_t axis1(OrthoSys::e1); // rotation of face
    AddPlane (InitialTag,   Vec3_t(Lx/2.0,0.0,0.0),  R, Cf*Lz, Cf*Ly, 1.0, M_PI/2.0, &axis1);
    Particles[Particles.Size()-1]->Initialize(Particles.Size()-1);
    AddPlane (InitialTag-1, Vec3_t(-Lx/2.0,0.0,0.0), R, Cf*Lz, Cf*Ly, 1.0, 3.0*M_PI/2.0, &axis1);
    Particles[Particles.Size()-1]->Initialize(Particles.Size()-1);
    AddPlane (InitialTag-2, Vec3_t(0.0,Ly/2.0,0.0),  R, Cf*Lx, Cf*Lz, 1.0, 3.0*M_PI/2.0, &axis0);
    Particles[Particles.Size()-1]->Initialize(Particles.Size()-1);
    AddPlane (InitialTag-3, Vec3_t(0.0,-Ly/2.0,0.0), R, Cf*Lx, Cf*Lz, 1.0, M_PI/2.0, &axis0);
    Particles[Particles.Size()-1]->Initialize(Particles.Size()-1);
    AddPlane (InitialTag-4, Vec3_t(0.0,0.0,-Lz/2.0), R, Cf*Lx, Cf*Ly, 1.0, M_PI, &axis0);
    Particles[Particles.Size()-1]->Initialize(Particles.Size()-1);
}

inline void Domain::GenBoundingBox (int InitialTag, double R, double Cf,bool Cohesion)
{
    Center();
    Vec3_t minX,maxX;
    BoundingBox(minX,maxX);
    GenBox(InitialTag, maxX(0)-minX(0)+2*R, maxX(1)-minX(1)+2*R, maxX(2)-minX(2)+2*R, R, Cf,Cohesion);
}

inline void Domain::GenBoundingPlane (int InitialTag, int GroupTag, double R, double Cf,bool Cohesion)
{
    Vec3_t minX,maxX;
    BoundingBoxTag(minX,maxX,GroupTag);
    Vec3_t axis0(OrthoSys::e0); // rotation of face
    Vec3_t axis1(OrthoSys::e1); // rotation of face
    size_t IIndex = Particles.Size();  // First face index
    double Lx = maxX(0)-minX(0)+2*R;
    double Ly = maxX(1)-minX(1)+2*R;
    double Lz = maxX(2)-minX(2)+2*R;
    AddPlane (InitialTag  , Vec3_t(0.5*(minX(0)+maxX(0)),maxX(1)+R,0.5*(minX(2)+maxX(2))),  R, Cf*Lx, Cf*Lz, 1.0, 3.0*M_PI/2.0, &axis0);
    Particles[Particles.Size()-1]->Initialize(Particles.Size()-1);
    AddPlane (InitialTag-1, Vec3_t(0.5*(minX(0)+maxX(0)),minX(1)-R,0.5*(minX(2)+maxX(2))), R, Cf*Lx, Cf*Lz, 1.0, M_PI/2.0, &axis0);
    Particles[Particles.Size()-1]->Initialize(Particles.Size()-1);

    // define some tolerance for comparissions
    if (Cohesion)
    {
        double tol1 = 1.0e-8;
        double tol2 = 1.0e-3;
        for (size_t i=0;i<IIndex;i++)
        {
            Particle * P1 = Particles[i];
            for (size_t j=IIndex;j<Particles.Size();j++)
            {
                Particle * P2 = Particles[j];
                for (size_t k=0;k<P1->Faces.Size();k++)
                {
                    Face * F1 = P1->Faces[k];
                    Vec3_t n1,c1;
                    F1->Normal  (n1);
                    F1->Centroid(c1);
                    Face * F2 = P2->Faces[0];
                    Vec3_t n2,c2;
                    F2->Normal  (n2);
                    F2->Centroid(c2);
                    Vec3_t n = 0.5*(n1-n2);
                    n/=norm(n);
                    if ((fabs(dot(n1,n2)+1.0)<tol1)
                       &&(fabs(dot(c2-c1,n)-2*R)<tol2)
                       &&(fabs(Distance(c1,*F2)-2.0*R)<tol2))
                    {
                        BInteractons.Push(new BInteracton(P1,P2,k,1));
                        break;
                    }
                }
            }        
        }
    }
}

inline void Domain::GenFromMesh (Mesh::Generic & M, double R, double rho, bool Cohesion, bool MC, double thickness)
{
    // info
    Util::Stopwatch stopwatch;
    printf("\n%s--- Generating particles from mesh ----------------------------------------------%s\n",TERM_CLR1,TERM_RST);

    size_t IIndex = Particles.Size();

    for (size_t i=0; i<M.Cells.Size(); ++i)
    {

        Array<Vec3_t> V;             // Array of vertices
        Array<Array <int> > E;       // Array of edges
        Array<Array <int> > F;       // array of faces
        if (M.NDim==3)
        {
            if (thickness > 0.0) throw new Fatal("Domain::GenFromMesh: Thickness should not be used in a 3D mesh");
            Array<Mesh::Vertex*> const & verts = M.Cells[i]->V;
            size_t nverts = verts.Size();

            // verts
            V.Resize(nverts);
            for (size_t j=0; j<nverts; ++j)
            {
                V[j] = verts[j]->C;
            }

            // edges
            size_t nedges = Mesh::NVertsToNEdges3D[nverts];
            E.Resize(nedges);
            for (size_t j=0; j<nedges; ++j)
            {
                E[j].Push (Mesh::NVertsToEdge3D[nverts][j][0]);
                E[j].Push (Mesh::NVertsToEdge3D[nverts][j][1]);
            }

            size_t nfaces = Mesh::NVertsToNFaces3D[nverts];
            size_t nvperf = Mesh::NVertsToNVertsPerFace3D[nverts];
            F.Resize(nfaces);
            for (size_t j=0; j<nfaces; ++j)
            {
                for (size_t k=0; k<nvperf; ++k)
                {
                    // TODO: check if face is planar or not
                    F[j].Push(Mesh::NVertsToFace3D[nverts][j][k]);
                }
            }
        }
        else if (M.NDim==2)
        {
            if (thickness <= 0.0) throw new Fatal("Domain::GenFromMesh: Thickness should be positive in a 2D mesh");
            Array<Mesh::Vertex*> const & verts = M.Cells[i]->V;
            size_t nverts = verts.Size();
            V.Resize(2*nverts);
            for (size_t j=0; j<nverts; ++j)
            {
                V[j] = verts[j]->C;
                V[j+nverts] = verts[j]->C + Vec3_t(0.0,0.0,thickness);
            }
            size_t nedges = 3*nverts;
            E.Resize(nedges);
            for (size_t j=0; j<nverts; ++j)
            {
                E[j].Push (Mesh::NVertsToEdge2D[nverts][j][0]);
                E[j].Push (Mesh::NVertsToEdge2D[nverts][j][1]);
                E[j+nverts].Push (Mesh::NVertsToEdge2D[nverts][j][0]+nverts);
                E[j+nverts].Push (Mesh::NVertsToEdge2D[nverts][j][1]+nverts);
                E[j+2*nverts].Push(j);
                E[j+2*nverts].Push(j+nverts);
            }
            size_t nfaces = nverts+2;
            F.Resize(nfaces);
            for (size_t j=0; j<nverts; ++j)
            {
                F[j].Push (Mesh::NVertsToEdge2D[nverts][j][0]);
                F[j].Push (Mesh::NVertsToEdge2D[nverts][j][1]);
                F[j].Push (Mesh::NVertsToEdge2D[nverts][j][1]+nverts);
                F[j].Push (Mesh::NVertsToEdge2D[nverts][j][0]+nverts);
                F[nverts].Push(nverts-1-j);
                F[nverts+1].Push(j+nverts);
            }
        }

        double vol; // volume of the polyhedron
        Vec3_t CM;  // Center of mass of the polyhedron
        Mat3_t It;  // Inertia tensor of the polyhedron
        PolyhedraMP(V,F,vol,CM,It);
        Erosion(V,E,F,R);

        // add particle
        Particles.Push (new Particle(M.Cells[i]->Tag, V,E,F,OrthoSys::O,OrthoSys::O,R,rho));
        Particles[Particles.Size()-1]->Eroded = true;
        Particles[Particles.Size()-1]->Index = Particles.Size()-1;
        if (!MC)
        {
            Particles[Particles.Size()-1]->x       = CM;
            Particles[Particles.Size()-1]->Props.V = vol;
            Particles[Particles.Size()-1]->Props.m = vol*rho;
            Vec3_t I;
            Quaternion_t Q;
            Vec3_t xp,yp,zp;
            Eig(It,I,xp,yp,zp);
            CheckDestroGiro(xp,yp,zp);
            I *= rho;
            Q(0) = 0.5*sqrt(1+xp(0)+yp(1)+zp(2));
            Q(1) = (yp(2)-zp(1))/(4*Q(0));
            Q(2) = (zp(0)-xp(2))/(4*Q(0));
            Q(3) = (xp(1)-yp(0))/(4*Q(0));
            Q = Q/norm(Q);
            Particles[Particles.Size()-1]->I     = I;
            Particles[Particles.Size()-1]->Q     = Q;
            double Dmax = Distance(CM,V[0])+R;
            for (size_t i=1; i<V.Size(); ++i)
            {
                if (Distance(CM,V[i])+R > Dmax) Dmax = Distance(CM,V[i])+R;
            }
            Particles[Particles.Size()-1]->Ekin = 0.0;
            Particles[Particles.Size()-1]->Erot = 0.0;
            Particles[Particles.Size()-1]->Dmax  = Dmax;
            Particles[Particles.Size()-1]->PropsReady = true;
        }
    }

    Array<Array <int> > Neigh(Particles.Size()-IIndex);
    Array<Array <int> > FNeigh(Particles.Size()-IIndex);
    if(Cohesion)
    {
        M.FindNeigh();
        //std::cout << M;
        for (size_t i=0; i<M.Cells.Size(); ++i) 
        {
            for (Mesh::Neighs_t::const_iterator p=M.Cells[i]->Neighs.begin(); p!=M.Cells[i]->Neighs.end(); ++p)
            {
                Neigh[i].Push(p->second.second->ID);
                FNeigh[i].Push(p->second.first);
            }           
        }
        for (size_t i=0; i<Neigh.Size(); ++i)
        {
            for (size_t j=0; j<Neigh[i].Size(); ++j)
            {
                size_t index = Neigh[Neigh[i][j]].Find(i);
                if ((size_t)Neigh[i][j]>i) BInteractons.Push(new BInteracton(Particles[i+IIndex],Particles[Neigh[i][j]+IIndex],FNeigh[i][j],FNeigh[Neigh[i][j]][index]));
            }
        }
        
    }

    // info
    printf("%s  Num of particles   = %zd%s\n",TERM_CLR2,Particles.Size(),TERM_RST);
}

inline void Domain::AddVoroPack (int Tag, double R, double Lx, double Ly, double Lz, size_t nx, size_t ny, size_t nz, double rho
                                 , bool Cohesion, bVec3_t Periodic,size_t Randomseed, double fraction, Vec3_t qin)
{
    // info
    Util::Stopwatch stopwatch;
    printf("\n%s--- Adding Voronoi particles packing --------------------------------------------%s\n",TERM_CLR1,TERM_RST);

    srand(Randomseed);
    const double x_min=-(nx/2.0), x_max=nx/2.0;
    const double y_min=-(ny/2.0), y_max=ny/2.0;
    const double z_min=-(nz/2.0), z_max=nz/2.0;
    voro::container con(x_min,x_max,y_min,y_max,z_min,z_max,nx,ny,nz, Periodic(0),Periodic(1),Periodic(2),8);
    int n = 0;
    for (size_t i=0; i<nx; i++)
    {
        for (size_t j=0; j<ny; j++)
        {
            for (size_t k=0; k<nz; k++)
            {
                double x = x_min+(i+0.5*qin(0)+(1-qin(0))*double(rand())/RAND_MAX)*(x_max-x_min)/nx;
                double y = y_min+(j+0.5*qin(1)+(1-qin(1))*double(rand())/RAND_MAX)*(y_max-y_min)/ny;
                double z = z_min+(k+0.5*qin(2)+(1-qin(2))*double(rand())/RAND_MAX)*(z_max-z_min)/nz;
                con.put (n,x,y,z);
                n++;
            }
        }
    }

    Array<Array <size_t> > ListBpairs(n);
    size_t IIndex = Particles.Size();
    voro::c_loop_all vl(con);
    voro::voronoicell c;
    if(vl.start()) do if(con.compute_cell(c,vl)) 
    {
        {
            if (rand()<fraction*RAND_MAX)
            {
                AddVoroCell(Tag,c,R,rho,true,Vec3_t(Lx/nx,Ly/ny,Lz/nz));
                double *pp = con.p[vl.ijk]+3*vl.q;
                Vec3_t trans(Lx*pp[0]/nx,Ly*pp[1]/ny,Lz*pp[2]/nz);
                Particle * P = Particles[Particles.Size()-1];
                P->Translate(trans);
            }
        }
    } while(vl.inc());


    // info
    printf("%s  Num of particles   = %zd%s\n",TERM_CLR2,Particles.Size(),TERM_RST);


    if (Cohesion)
    {
        //if (fraction<1.0) throw new Fatal("Domain::AddVoroPack: With the Cohesion all particles should be considered, plese change the fraction to 1.0");

        // define some tolerance for comparissions
        double tol1 = 1.0e-8;
        double tol2 = 1.0e-3;
        for (size_t i=IIndex;i<Particles.Size()-1;i++)
        {
            Particle * P1 = Particles[i];
            for (size_t j=i+1;j<Particles.Size();j++)
            {
                Particle * P2 = Particles[j];
                if (Distance(P1->x,P2->x)<P1->Dmax+P2->Dmax)
                {
                    for (size_t k=0;k<P1->Faces.Size();k++)
                    {
                        Face * F1 = P1->Faces[k];
                        Vec3_t n1,c1;
                        F1->Normal  (n1);
                        F1->Centroid(c1);
                        bool found = false;
                        for (size_t l=0;l<P2->Faces.Size();l++)
                        {
                            Face * F2 = P2->Faces[l];
                            Vec3_t n2,c2;
                            F2->Normal  (n2);
                            F2->Centroid(c2);
                            Vec3_t n = 0.5*(n1-n2);
                            n/=norm(n);
                            if ((fabs(dot(n1,n2)+1.0)<tol1)
                               &&(fabs(Distance(c1,*F2)-2*R)<tol2)
                               &&(fabs(Distance(c2,*F1)-2*R)<tol2))
                            {
                                BInteractons.Push(new BInteracton(P1,P2,k,l));
                                found = true;
                                break;
                            }
                        }
                        if (found) break;
                    }
                }
            }
        }
    }
}

inline void Domain::AddVoroPack (int Tag, double R, double Lx, double Ly, double Lz, size_t nx, size_t ny, size_t nz, double rho
                                 , bool Cohesion, bool Periodic,size_t Randomseed, double fraction, Vec3_t qin)
{
    AddVoroPack(Tag,R,Lx,Ly,Lz,nx,ny,nz,rho,Cohesion,bVec3_t(Periodic,Periodic,Periodic),Randomseed,fraction,qin);
}
// Sihgle particle addition

inline void Domain::AddSphere (int Tag,Vec3_t const & X, double R, double rho)
{
    // vertices
    Array<Vec3_t> V(1);
    V[0] = X;

    // edges
    Array<Array <int> > E(0); // no edges

    // faces
    Array<Array <int> > F(0); // no faces

    // add particle
    Particles.Push (new Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));

    Quaternion_t q;
    q(0) = 1.0;
    q(1) = 0.0;
    q(2) = 0.0;
    q(3) = 0.0;
    q = q/norm(q);

    Particles[Particles.Size()-1]->Q          = q;
    Particles[Particles.Size()-1]->Props.V    = (4.0/3.0)*M_PI*R*R*R;
    Particles[Particles.Size()-1]->Props.m    = rho*(4.0/3.0)*M_PI*R*R*R;
    Particles[Particles.Size()-1]->I          = (2.0/5.0)*Particles[Particles.Size()-1]->Props.m*R*R;
    Particles[Particles.Size()-1]->x          = X;
    Particles[Particles.Size()-1]->Ekin       = 0.0;
    Particles[Particles.Size()-1]->Erot       = 0.0;
    Particles[Particles.Size()-1]->Dmax       = R;
    Particles[Particles.Size()-1]->PropsReady = true;
    Particles[Particles.Size()-1]->Index      = Particles.Size()-1;

}

inline void Domain::AddCube (int Tag, Vec3_t const & X, double R, double L, double rho, double Angle, Vec3_t * Axis)
{
    // vertices
    Array<Vec3_t> V(8);
    double l = L/2.0;
    V[0] = -l, -l, -l;
    V[1] =  l, -l, -l;
    V[2] =  l,  l, -l;
    V[3] = -l,  l, -l;
    V[4] = -l, -l,  l;
    V[5] =  l, -l,  l;
    V[6] =  l,  l,  l;
    V[7] = -l,  l,  l;

    // edges
    Array<Array <int> > E(12);
    for (size_t i=0; i<12; ++i) E[i].Resize(2);
    E[ 0] = 0, 1;
    E[ 1] = 1, 2;
    E[ 2] = 2, 3;
    E[ 3] = 3, 0;
    E[ 4] = 4, 5;
    E[ 5] = 5, 6;
    E[ 6] = 6, 7;
    E[ 7] = 7, 4;
    E[ 8] = 0, 4;
    E[ 9] = 1, 5;
    E[10] = 2, 6;
    E[11] = 3, 7;

    // faces
    Array<Array <int> > F(6);
    for (size_t i=0; i<6; i++) F[i].Resize(4);
    F[0] = 4, 7, 3, 0;
    F[1] = 1, 2, 6, 5;
    F[2] = 0, 1, 5, 4;
    F[3] = 2, 3, 7, 6;
    F[4] = 0, 3, 2, 1;
    F[5] = 4, 5, 6, 7;

    // calculate the rotation
    bool ThereisanAxis = true;
    if (Axis==NULL)
    {
        Angle   = (1.0*rand())/RAND_MAX*2*M_PI;
        Axis = new Vec3_t((1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX);
        ThereisanAxis = false;
    }
    Quaternion_t q;
    NormalizeRotation (Angle,(*Axis),q);
    for (size_t i=0; i<V.Size(); i++)
    {
        Vec3_t t;
        Rotation (V[i],q,t);
        V[i] = t+X;
    }

    // add particle
    Particles.Push (new Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));

    // clean up
    if (!ThereisanAxis) delete Axis;
    q(0) = 1.0;
    q(1) = 0.0;
    q(2) = 0.0;
    q(3) = 0.0;
    q = q/norm(q);

    Particles[Particles.Size()-1]->Q          = q;
    Particles[Particles.Size()-1]->Props.V    = L*L*L;
    Particles[Particles.Size()-1]->Props.m    = rho*L*L*L;
    Particles[Particles.Size()-1]->I          = L*L, L*L, L*L;
    Particles[Particles.Size()-1]->I         *= Particles[Particles.Size()-1]->Props.m/6.0;
    Particles[Particles.Size()-1]->x          = X;
    Particles[Particles.Size()-1]->Ekin       = 0.0;
    Particles[Particles.Size()-1]->Erot       = 0.0;
    Particles[Particles.Size()-1]->Dmax       = sqrt(3.0*L*L/4.0)+R;
    Particles[Particles.Size()-1]->PropsReady = true;
    Particles[Particles.Size()-1]->Index      = Particles.Size()-1;
    

}

inline void Domain::AddRecBox (int Tag, Vec3_t const & X, Vec3_t const & L, double R, double rho, double Angle, Vec3_t * Axis)
{
    // vertices
    Array<Vec3_t> V(8);
    double lx = L(0)/2.0;
    double ly = L(1)/2.0;
    double lz = L(2)/2.0;
    V[0] = -lx, -ly, -lz;
    V[1] =  lx, -ly, -lz;
    V[2] =  lx,  ly, -lz;
    V[3] = -lx,  ly, -lz;
    V[4] = -lx, -ly,  lz;
    V[5] =  lx, -ly,  lz;
    V[6] =  lx,  ly,  lz;
    V[7] = -lx,  ly,  lz;

    // edges
    Array<Array <int> > E(12);
    for (size_t i=0; i<12; ++i) E[i].Resize(2);
    E[ 0] = 0, 1;
    E[ 1] = 1, 2;
    E[ 2] = 2, 3;
    E[ 3] = 3, 0;
    E[ 4] = 4, 5;
    E[ 5] = 5, 6;
    E[ 6] = 6, 7;
    E[ 7] = 7, 4;
    E[ 8] = 0, 4;
    E[ 9] = 1, 5;
    E[10] = 2, 6;
    E[11] = 3, 7;

    // faces
    Array<Array <int> > F(6);
    for (size_t i=0; i<6; i++) F[i].Resize(4);
    F[0] = 4, 7, 3, 0;
    F[1] = 1, 2, 6, 5;
    F[2] = 0, 1, 5, 4;
    F[3] = 2, 3, 7, 6;
    F[4] = 0, 3, 2, 1;
    F[5] = 4, 5, 6, 7;

    // calculate the rotation
    bool ThereisanAxis = true;
    if (Axis==NULL)
    {
        Angle   = (1.0*rand())/RAND_MAX*2*M_PI;
        Axis = new Vec3_t((1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX);
        ThereisanAxis = false;
    }
    Quaternion_t q;
    NormalizeRotation (Angle,(*Axis),q);
    for (size_t i=0; i<V.Size(); i++)
    {
        Vec3_t t;
        Rotation (V[i],q,t);
        V[i] = t+X;
    }

    // add particle
    Particles.Push (new Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));

    // clean up
    if (!ThereisanAxis) delete Axis;
    //q(0) = 1.0;
    //q(1) = 0.0;
    //q(2) = 0.0;
    //q(3) = 0.0;
    q = q/norm(q);

    Particles[Particles.Size()-1]->Q          = q;
    Particles[Particles.Size()-1]->Props.V    = L(0)*L(1)*L(2);
    Particles[Particles.Size()-1]->Props.m    = rho*L(0)*L(1)*L(2);
    Particles[Particles.Size()-1]->I          = L(1)*L(1) + L(2)*L(2), L(0)*L(0) + L(2)*L(2), L(0)*L(0) + L(1)*L(1);
    Particles[Particles.Size()-1]->I         *= Particles[Particles.Size()-1]->Props.m/12.0;
    Particles[Particles.Size()-1]->x          = X;
    Particles[Particles.Size()-1]->Ekin       = 0.0;
    Particles[Particles.Size()-1]->Erot       = 0.0;
    Particles[Particles.Size()-1]->Dmax       = 0.5*norm(L)+R;
    Particles[Particles.Size()-1]->PropsReady = true;
    Particles[Particles.Size()-1]->Index      = Particles.Size()-1;
    

}

inline void Domain::AddTetra (int Tag, Vec3_t const & X, double R, double L, double rho, double Angle, Vec3_t * Axis)
{
    // vertices
    double sq8 = sqrt(8.0);
    Array<Vec3_t> V(4);
    V[0] =  L/sq8,  L/sq8, L/sq8;
    V[1] = -L/sq8, -L/sq8, L/sq8;
    V[2] = -L/sq8,  L/sq8,-L/sq8;
    V[3] =  L/sq8, -L/sq8,-L/sq8;

    // edges
    Array<Array <int> > E(6);
    for (size_t i=0; i<6; ++i) E[i].Resize(2);
    E[0] = 0, 1;
    E[1] = 1, 2;
    E[2] = 2, 0;
    E[3] = 0, 3;
    E[4] = 1, 3;
    E[5] = 2, 3;

    // face
    Array<Array <int> > F;
    F.Resize(4);
    for (size_t i=0; i<4; ++i) F[i].Resize(3);
    F[0] = 0, 3, 2;
    F[1] = 0, 1, 3;
    F[2] = 0, 2, 1;
    F[3] = 1, 2, 3;

    // calculate the rotation
    bool ThereisanAxis = true;
    if (Axis==NULL)
    {
        Angle   = (1.0*rand())/RAND_MAX*2*M_PI;
        Axis = new Vec3_t((1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX);
        ThereisanAxis = false;
    }
    Quaternion_t q;
    NormalizeRotation (Angle,(*Axis),q);
    for (size_t i=0; i<V.Size(); i++)
    {
        Vec3_t t;
        Rotation (V[i],q,t);
        V[i] = t+X;
    }

    // add particle
    Particles.Push (new Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));

    // clean up
    if (!ThereisanAxis) delete Axis;
    q(0) = 1.0;
    q(1) = 0.0;
    q(2) = 0.0;
    q(3) = 0.0;
    q = q/norm(q);

    Particles[Particles.Size()-1]->Q          = q;
    Particles[Particles.Size()-1]->Props.V    = sqrt(2.0)*L*L*L/12.0;
    Particles[Particles.Size()-1]->Props.m    = rho*sqrt(2.0)*L*L*L/12.0;
    Particles[Particles.Size()-1]->I          = L*L, L*L, L*L;
    Particles[Particles.Size()-1]->I         *= Particles[Particles.Size()-1]->Props.m/20.0;
    Particles[Particles.Size()-1]->x          = X;
    Particles[Particles.Size()-1]->Ekin       = 0.0;
    Particles[Particles.Size()-1]->Erot       = 0.0;
    Particles[Particles.Size()-1]->Dmax       = sqrt(3.0*L*L/8.0)+R;
    Particles[Particles.Size()-1]->PropsReady = true;
    Particles[Particles.Size()-1]->Index      = Particles.Size()-1;
}

inline void Domain::AddDrill (int Tag, const Vec3_t & X, double R, double Lt, double Ll, double rho)
{
    Array<Vec3_t> V(9);
    V[0] =  Lt/2.0,  Lt/2.0, Ll/2.0;
    V[1] = -Lt/2.0,  Lt/2.0, Ll/2.0;
    V[2] = -Lt/2.0, -Lt/2.0, Ll/2.0;
    V[3] =  Lt/2.0, -Lt/2.0, Ll/2.0;
    V[4] =  Lt/2.0,  Lt/2.0,    0.0;
    V[5] = -Lt/2.0,  Lt/2.0,    0.0;
    V[6] = -Lt/2.0, -Lt/2.0,    0.0;
    V[7] =  Lt/2.0, -Lt/2.0,    0.0;
    V[8] =     0.0,     0.0,-Ll/2.0;
    
    Array<Array <int> > E(12);
    for (size_t i=0; i<12; ++i) E[i].Resize(2);
    for (size_t i=0;i<4;i++)
    {
        E[i]   = i  , (i+1)%4    ;
        E[i+4] = i+4, (i+5)%4 + 4;
        E[i+8] = i+4, 8;
    }

    Array<Array <int> > F;
    F.Resize(9);
    F[0].Resize(4);
    F[0] = 0, 1, 2, 3;
    F[1].Resize(4);
    F[1] = 0, 4, 5, 1;
    F[2].Resize(4);
    F[2] = 1, 5, 6, 2;
    F[3].Resize(4);
    F[3] = 2, 6, 7, 3;
    F[4].Resize(4);
    F[4] = 3, 7, 4, 0;
    F[5].Resize(3);
    F[5] = 4, 8, 5;
    F[6].Resize(3);
    F[6] = 5, 8, 6;
    F[7].Resize(3);
    F[7] = 6, 8, 7;
    F[8].Resize(3);
    F[8] = 7, 8, 4;


    double vol; // volume of the polyhedron
    Vec3_t CM;  // Center of mass of the polyhedron
    Mat3_t It;  // Inertia tensor of the polyhedron
    PolyhedraMP(V,F,vol,CM,It);
    Particles.Push (new Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));
    Particles[Particles.Size()-1]->x       = CM;
    Particles[Particles.Size()-1]->Props.V = vol;
    Particles[Particles.Size()-1]->Props.m = vol*rho;
    Vec3_t I;
    Quaternion_t Q;
    Vec3_t xp,yp,zp;
    Eig(It,I,xp,yp,zp);
    CheckDestroGiro(xp,yp,zp);
    I *= rho;
    Q(0) = 0.5*sqrt(1+xp(0)+yp(1)+zp(2));
    Q(1) = (yp(2)-zp(1))/(4*Q(0));
    Q(2) = (zp(0)-xp(2))/(4*Q(0));
    Q(3) = (xp(1)-yp(0))/(4*Q(0));
    Q = Q/norm(Q);
    Particles[Particles.Size()-1]->I     = I;
    Particles[Particles.Size()-1]->Q     = Q;
    double Dmax = Distance(CM,V[0])+R;
    for (size_t i=1; i<V.Size(); ++i)
    {
        if (Distance(CM,V[i])+R > Dmax) Dmax = Distance(CM,V[i])+R;
    }
    Particles[Particles.Size()-1]->Ekin = 0.0;
    Particles[Particles.Size()-1]->Erot = 0.0;
    Particles[Particles.Size()-1]->Dmax  = Dmax;
    Particles[Particles.Size()-1]->PropsReady = true;
    Particles[Particles.Size()-1]->Index = Particles.Size()-1;

    Vec3_t Y = X;
    Particles[Particles.Size()-1]->Translate(Y);

}

inline void Domain::AddRice (int Tag, const Vec3_t & X, double R, double L, double rho, double Angle, Vec3_t * Axis)
{
    // vertices
    Array<Vec3_t> V(2);
    V[0] = 0.0, 0.0,  L/2;
    V[1] = 0.0, 0.0, -L/2;

    // edges
    Array<Array <int> > E(1);
    E[0].Resize(2);
    E[0] = 0, 1;

    // faces
    Array<Array <int> > F(0); // no faces

    // calculate the rotation
    bool ThereisanAxis = true;
    if (Axis==NULL)
    {
        Angle   = (1.0*rand())/RAND_MAX*2*M_PI;
        Axis = new Vec3_t((1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX);
        ThereisanAxis = false;
    }
    Quaternion_t q;
    NormalizeRotation (Angle,(*Axis),q);
    for (size_t i=0; i<V.Size(); i++)
    {
        Vec3_t t;
        Rotation (V[i],q,t);
        V[i] = t+X;
    }

    // add particle
    Particles.Push (new Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));

    Particles[Particles.Size()-1]->Q          = q;
    Particles[Particles.Size()-1]->Props.V    = (4.0/3.0)*M_PI*R*R*R + M_PI*L*R*R;
    Particles[Particles.Size()-1]->Props.m    = rho*Particles[Particles.Size()-1]->Props.V;
    Particles[Particles.Size()-1]->I          = (1.0/3.0)*M_PI*rho*R*R*R*L*L + (1.0/12.0)*M_PI*rho*R*R*L*L*L + (3.0/4.0)*M_PI*rho*R*R*R*R*L + (8.0/15.0)*M_PI*rho*R*R*R*R*R,
                                                (1.0/3.0)*M_PI*rho*R*R*R*L*L + (1.0/12.0)*M_PI*rho*R*R*L*L*L + (3.0/4.0)*M_PI*rho*R*R*R*R*L + (8.0/15.0)*M_PI*rho*R*R*R*R*R,
                                                0.5*M_PI*rho*R*R*R*R*L + (8.0/15.0)*M_PI*rho*R*R*R*R*R;
    Particles[Particles.Size()-1]->x          = X;
    Particles[Particles.Size()-1]->Ekin       = 0.0;
    Particles[Particles.Size()-1]->Erot       = 0.0;
    Particles[Particles.Size()-1]->Dmax       = 0.5*L + R;
    Particles[Particles.Size()-1]->PropsReady = true;
    Particles[Particles.Size()-1]->Index      = Particles.Size()-1;

    // clean up
    if (!ThereisanAxis) delete Axis;
}

inline void Domain::AddPlane (int Tag, const Vec3_t & X, double R, double Lx, double Ly, double rho, double Angle, Vec3_t * Axis)
{
    // vertices
    Array<Vec3_t> V(4);
    double lx = Lx/2.0, ly = Ly/2.0;
    V[0] = -lx, -ly, 0.0;
    V[1] =  lx, -ly, 0.0;
    V[2] =  lx,  ly, 0.0;
    V[3] = -lx,  ly, 0.0;

    // edges
    Array<Array <int> > E(4);
    for (size_t i=0; i<4; ++i) E[i].Resize(2);
    E[ 0] = 0, 1;
    E[ 1] = 1, 2;
    E[ 2] = 2, 3;
    E[ 3] = 3, 0;

    // faces
    Array<Array <int> > F(1);
    F[0].Resize(4);
    F[0] = 0, 3, 2, 1;

    bool ThereisanAxis = true;
    if (Axis==NULL)
    {
        Angle   = 0.;
        Axis = new Vec3_t((1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX);
        ThereisanAxis = false;
    }
    Quaternion_t q;
    NormalizeRotation (Angle,(*Axis),q);
    for (size_t i=0; i<V.Size(); i++)
    {
        Vec3_t t;
        Rotation (V[i],q,t);
        V[i] = t+X;
    }

    // add particle
    Particles.Push (new Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));
    Particles[Particles.Size()-1]->Q          = q;
    Particles[Particles.Size()-1]->Props.V    = Lx*Ly*2*R;
    Particles[Particles.Size()-1]->Props.m    = rho*Lx*Ly*2*R;
    Particles[Particles.Size()-1]->I          = (1.0/12.0)*(Ly*Ly+4*R*R),(1.0/12.0)*(Lx*Lx+4*R*R),(1.0/12.0)*(Lx*Lx+Ly*Ly);
    Particles[Particles.Size()-1]->I         *= Particles[Particles.Size()-1]->Props.m;
    Particles[Particles.Size()-1]->x          = X;
    Particles[Particles.Size()-1]->Ekin       = 0.0;
    Particles[Particles.Size()-1]->Erot       = 0.0;
    Particles[Particles.Size()-1]->Dmax       = sqrt(Lx*Lx+Ly*Ly)+R;
    Particles[Particles.Size()-1]->PropsReady = true;
    Particles[Particles.Size()-1]->Index      = Particles.Size()-1;
    Particles[Particles.Size()-1]->Closed     = false;
    // clean up
    if (!ThereisanAxis) delete Axis;
}

inline void Domain::AddVoroCell (int Tag, voro::voronoicell & VC, double R, double rho, bool Erode, Vec3_t nv)
{
    Array<Vec3_t> V(VC.p);
    Array<Array <int> > E;
    Array<int> Eaux(2);
    for(int i=0;i<VC.p;i++) 
    {
        V[i] = Vec3_t(0.5*VC.pts[3*i]*nv(0),0.5*VC.pts[3*i+1]*nv(1),0.5*VC.pts[3*i+2]*nv(2));
        for(int j=0;j<VC.nu[i];j++) 
        {
            int k=VC.ed[i][j];
            if (VC.ed[i][j]<i) 
            {
                Eaux[0] = i;
                Eaux[1] = k;
                E.Push(Eaux);
            }
        }
    }
    Array<Array <int> > F;
    Array<int> Faux;
    for(int i=0;i<VC.p;i++) 
    {
        for(int j=0;j<VC.nu[i];j++) 
        {
            int k=VC.ed[i][j];
            if (k>=0) 
            {
                Faux.Push(i);
                VC.ed[i][j]=-1-k;
                int l=VC.cycle_up(VC.ed[i][VC.nu[i]+j],k);
                do 
                {
                    Faux.Push(k);
                    int m=VC.ed[k][l];
                    VC.ed[k][l]=-1-m;
                    l=VC.cycle_up(VC.ed[k][VC.nu[k]+l],m);
                    k=m;
                } while (k!=i);
                Array<int> Faux2(Faux.Size());
                for (size_t l = 0; l < Faux.Size();l++)
                {
                    Faux2[l] = Faux[Faux.Size()-1-l];
                }

                F.Push(Faux2);
                Faux.Clear();
                Faux2.Clear();
            }
        }
    }
    //VC.reset_edges();
    double vol; // volume of the polyhedron
    Vec3_t CM;  // Center of mass of the polyhedron
    Mat3_t It;  // Inertia tensor of the polyhedron
    PolyhedraMP(V,F,vol,CM,It);
    if (Erode) Erosion(V,E,F,R);
    // add particle
    Particles.Push (new Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));
    if (Erode) Particles[Particles.Size()-1]->Eroded = true;
    Particles[Particles.Size()-1]->x       = CM;
    Particles[Particles.Size()-1]->Props.V = vol;
    Particles[Particles.Size()-1]->Props.m = vol*rho;
    Vec3_t I;
    Quaternion_t Q;
    Vec3_t xp,yp,zp;
    Eig(It,I,xp,yp,zp);
    CheckDestroGiro(xp,yp,zp);
    I *= rho;
    Q(0) = 0.5*sqrt(1+xp(0)+yp(1)+zp(2));
    Q(1) = (yp(2)-zp(1))/(4*Q(0));
    Q(2) = (zp(0)-xp(2))/(4*Q(0));
    Q(3) = (xp(1)-yp(0))/(4*Q(0));
    Q = Q/norm(Q);
    Particles[Particles.Size()-1]->I     = I;
    Particles[Particles.Size()-1]->Q     = Q;
    double Dmax = Distance(CM,V[0])+R;
    for (size_t i=1; i<V.Size(); ++i)
    {
        if (Distance(CM,V[i])+R > Dmax) Dmax = Distance(CM,V[i])+R;
    }
    Particles[Particles.Size()-1]->Ekin = 0.0;
    Particles[Particles.Size()-1]->Erot = 0.0;
    Particles[Particles.Size()-1]->Dmax  = Dmax;
    Particles[Particles.Size()-1]->PropsReady = true;
    Particles[Particles.Size()-1]->Index = Particles.Size()-1;
}

inline void Domain::AddTorus (int Tag, Vec3_t const & X, Vec3_t const & N, double Rmax, double R, double rho)
{
    // Normalize normal vector
    Vec3_t n = N/norm(N);

    // Create the 2 vertices that define the torus
    Vec3_t P1 = OrthoSys::e0 - dot(OrthoSys::e0,n)*n;
    if (norm(P1)<1.0e-12) P1 = OrthoSys::e1 - dot(OrthoSys::e1,n)*n;
    P1       /= norm(P1);
    Vec3_t P2 = cross(n,P1);
    P2       /= norm(P2);
    Array<Vec3_t > V(2);
    V[0] = X + Rmax*P1;
    V[1] = X + Rmax*P2;

    // the torus has no edges or faces
    Array<Array <int> > E(0);
    Array<Array <int> > F(0);

    //Add the particle just with two vertices
    Particles.Push (new Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));

    //Input all the mass properties
    Vec3_t xp,yp,zp;
    xp = P1;
    zp = P2;
    yp = cross(zp,xp);
    CheckDestroGiro(xp,yp,zp);
    Quaternion_t q;
    q(0) = 0.5*sqrt(1+xp(0)+yp(1)+zp(2));
    q(1) = (yp(2)-zp(1))/(4*q(0));
    q(2) = (zp(0)-xp(2))/(4*q(0));
    q(3) = (xp(1)-yp(0))/(4*q(0));
    q = q/norm(q);

    if (std::isnan(norm(q))) q = 1.0,0.0,0.0,0.0;

    Particles[Particles.Size()-1]->Q          = q;
    Particles[Particles.Size()-1]->Props.V    = 2*M_PI*M_PI*R*R*Rmax;
    Particles[Particles.Size()-1]->Props.m    = rho*2*M_PI*M_PI*R*R*Rmax;
    Particles[Particles.Size()-1]->I          = 0.125*(4*Rmax*Rmax + 5*R*R), (Rmax*Rmax+0.75*R*R), 0.125*(4*Rmax*Rmax + 5*R*R);
    Particles[Particles.Size()-1]->I         *= Particles[Particles.Size()-1]->Props.m;
    Particles[Particles.Size()-1]->x          = X;
    Particles[Particles.Size()-1]->Ekin       = 0.0;
    Particles[Particles.Size()-1]->Erot       = 0.0;
    Particles[Particles.Size()-1]->Dmax       = Rmax + R;
    Particles[Particles.Size()-1]->PropsReady = true;
    Particles[Particles.Size()-1]->Index      = Particles.Size()-1;
    
    Particles[Particles.Size()-1]->Tori.Push(new Torus(&Particles[Particles.Size()-1]->x,Particles[Particles.Size()-1]->Verts[0],Particles[Particles.Size()-1]->Verts[1]));

}

inline void Domain::AddCylinder (int Tag, Vec3_t const & X0, double R0, Vec3_t const & X1, double R1, double R, double rho)
{
    // Define all key variables of the cylinder
    Vec3_t n = X1 - X0;
    n /= norm(n);
    Vec3_t P1 = OrthoSys::e0 - dot(OrthoSys::e0,n)*n;
    if (norm(P1)<1.0e-12) P1 = OrthoSys::e2 - dot(OrthoSys::e2,n)*n;
    P1       /= norm(P1);
    Vec3_t P2 = cross(n,P1);
    P2       /= norm(P2);

    //The cylinder is defined by 6 vertices
    Array<Vec3_t > V(6);
    V[0] = X0 + R0*P1;
    V[1] = X0 + R0*P2;
    V[2] = X0 - R0*P1;
    V[3] = X1 + R1*P1;
    V[4] = X1 + R1*P2;
    V[5] = X1 - R1*P1;

    // It has no edges or faces
    Array<Array <int> > E(0);
    Array<Array <int> > F(0);

    //Add the particle just with the vertices
    //
    Particles.Push (new Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));

    //Input all the mass properties
    Vec3_t xp,yp,zp;
    xp = P1;
    zp = P2;
    yp = cross(zp,xp);
    CheckDestroGiro(xp,yp,zp);
    Quaternion_t q;
    q(0) = 0.5*sqrt(1+xp(0)+yp(1)+zp(2));
    q(1) = (yp(2)-zp(1))/(4*q(0));
    q(2) = (zp(0)-xp(2))/(4*q(0));
    q(3) = (xp(1)-yp(0))/(4*q(0));
    q = q/norm(q);

    if (std::isnan(norm(q))) q = 1.0,0.0,0.0,0.0;

    Particles[Particles.Size()-1]->Q          = q;
    Particles[Particles.Size()-1]->Props.V    = 4.0*M_PI*R0*R*norm(X1-X0);
    Particles[Particles.Size()-1]->Props.m    = rho*4.0*M_PI*R0*R*norm(X1-X0);
    Particles[Particles.Size()-1]->I          = 1.0, 1.0, 1.0;
    Particles[Particles.Size()-1]->x          = 0.5*(X0 + X1);
    Particles[Particles.Size()-1]->Ekin       = 0.0;
    Particles[Particles.Size()-1]->Erot       = 0.0;
    Particles[Particles.Size()-1]->Dmax       = sqrt(0.25*dot(X1-X0,X1-X0)+std::max(R0,R1)*std::max(R0,R1)) + R;
    Particles[Particles.Size()-1]->PropsReady = true;
    Particles[Particles.Size()-1]->Index      = Particles.Size()-1;
    Particles[Particles.Size()-1]->Tori.Push     (new Torus(&X0,Particles[Particles.Size()-1]->Verts[0],Particles[Particles.Size()-1]->Verts[1]));
    Particles[Particles.Size()-1]->Tori.Push     (new Torus(&X1,Particles[Particles.Size()-1]->Verts[3],Particles[Particles.Size()-1]->Verts[4]));
    Particles[Particles.Size()-1]->Cylinders.Push(new Cylinder(Particles[Particles.Size()-1]->Tori[0],Particles[Particles.Size()-1]->Tori[1],Particles[Particles.Size()-1]->Verts[2],Particles[Particles.Size()-1]->Verts[5]));

    //std::cout << *Particles[Particles.Size()-1] << std::endl;
}

inline void Domain::AddFromJson (int Tag, char const * Filename, double R, double rho, double scale, bool Erode)
{

    Array<Vec3_t>       V;
    Array<Array <int> > E;
    Array<Array <int> > F;

    try {
        boost::property_tree::ptree pt;
        boost::property_tree::read_json(Filename, pt);
        BOOST_FOREACH(boost::property_tree::ptree::value_type & a, pt.get_child("verts")) {
            Vec3_t coords;
            int i = 0;
            BOOST_FOREACH(boost::property_tree::ptree::value_type & b, a.second.get_child("c")) {
                coords[i] = scale * boost::lexical_cast<double>(b.second.data());
                i++;
            }
            V.Push(coords);
        }
        BOOST_FOREACH(boost::property_tree::ptree::value_type & a, pt.get_child("edges")) {
            Array<int> vids(2);
            int i = 0;
            BOOST_FOREACH(boost::property_tree::ptree::value_type & b, a.second.get_child("verts")) {
                vids[i] = boost::lexical_cast<int>(b.second.data());
                i++;
            }
            E.Push(vids);
        }
        BOOST_FOREACH(boost::property_tree::ptree::value_type & a, pt.get_child("faces")) {
            Array<int>     vids;
            BOOST_FOREACH(boost::property_tree::ptree::value_type & b, a.second.get_child("verts")) {
                int vid = boost::lexical_cast<int>(b.second.data());
                vids .Push(vid);
            }
            F.Push(vids);
        }
        printf("[1;32mDEM::domain.h ConstructFromJson: finished[0m\n");
    } catch (std::exception & e) {
        throw new Fatal("DEM::domain.h: ConstructFromJson failed:\n\t%s", e.what());
    }
    double vol; // volume of the polyhedron
    Vec3_t CM;  // Center of mass of the polyhedron
    Mat3_t It;  // Inertia tensor of the polyhedron
    PolyhedraMP(V,F,vol,CM,It);
    if (Erode) Erosion(V,E,F,R);
    // add particle
    Particles.Push (new Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));
    if (Erode) Particles[Particles.Size()-1]->Eroded = true;
    Particles[Particles.Size()-1]->x       = CM;
    Particles[Particles.Size()-1]->Props.V = vol;
    Particles[Particles.Size()-1]->Props.m = vol*rho;
    Vec3_t I;
    Quaternion_t Q;
    Vec3_t xp,yp,zp;
    Eig(It,I,xp,yp,zp);
    CheckDestroGiro(xp,yp,zp);
    I *= rho;
    Q(0) = 0.5*sqrt(1+xp(0)+yp(1)+zp(2));
    Q(1) = (yp(2)-zp(1))/(4*Q(0));
    Q(2) = (zp(0)-xp(2))/(4*Q(0));
    Q(3) = (xp(1)-yp(0))/(4*Q(0));
    Q = Q/norm(Q);
    Particles[Particles.Size()-1]->I     = I;
    Particles[Particles.Size()-1]->Q     = Q;
    double Dmax = Distance(CM,V[0])+R;
    for (size_t i=1; i<V.Size(); ++i)
    {
        if (Distance(CM,V[i])+R > Dmax) Dmax = Distance(CM,V[i])+R;
    }
    Particles[Particles.Size()-1]->Ekin = 0.0;
    Particles[Particles.Size()-1]->Erot = 0.0;
    Particles[Particles.Size()-1]->Dmax  = Dmax;
    Particles[Particles.Size()-1]->PropsReady = true;
    Particles[Particles.Size()-1]->Index = Particles.Size()-1;
}

inline void Domain::AddFromOBJ  (int Tag, char const * Filename, double R, double rho, double scale, bool Erode)
{

    Array<Vec3_t>       V;
    Array<Array <int> > E;
    Array<Array <int> > F;

    String fn(Filename);
    if (!Util::FileExists(fn)) throw new Fatal("File <%s> not found",fn.CStr());
    ifstream fi(fn.CStr());
    
    std::string line;
    while(getline(fi,line))
    {
        std::istringstream ss( line );
        char c;
        double x,y,z;
	  	if (line[0]=='v' && line[1]==' ')
	  	{
	  		ss >> c >> x >> y >> z;
	        Vec3_t v(x,y,z);
            V.Push(v);
	  	}
	  	else if (line[0]=='f')
	  	{
            std::vector<std::string> ls;
            std::string lse;
			while (ss>>lse)	ls.push_back(lse);
            Array<int> face;
	        for (size_t m=1; m<ls.size(); ++m)
	        {
	        	 face.Push(stoi(ls[m])-1);
	        }
            F.Push(face);
        }
    }
    
    std::set<std::pair <int,int> > edges;
    for (size_t nf=0;nf<F.Size();nf++)
    {
        for (size_t nff=0;nff<F[nf].Size();nff++)
        {
            std::pair <int,int> p1,p2;
            p1 = std::make_pair(F[nf][(nff+1)%F[nf].Size()],F[nf][nff]);
            p2 = std::make_pair(F[nf][nff],F[nf][(nff+1)%F[nf].Size()]);
            if (!edges.count(p1)&&!edges.count(p2))
            {
                Array<int> edge;
                edge.Push(F[nf][nff]);
                edge.Push(F[nf][(nff+1)%F[nf].Size()]);
                E.Push(edge);
                edges.insert(p1);
            }
        }
    }
    

    double vol; // volume of the polyhedron
    Vec3_t CM;  // Center of mass of the polyhedron
    Mat3_t It;  // Inertia tensor of the polyhedron
    PolyhedraMP(V,F,vol,CM,It);
    if (Erode) Erosion(V,E,F,R);
    // add particle
    Particles.Push (new Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));
    if (Erode) Particles[Particles.Size()-1]->Eroded = true;
    Particles[Particles.Size()-1]->x       = CM;
    Particles[Particles.Size()-1]->Props.V = vol;
    Particles[Particles.Size()-1]->Props.m = vol*rho;
    Vec3_t I;
    Quaternion_t Q;
    Vec3_t xp,yp,zp;
    Eig(It,I,xp,yp,zp);
    CheckDestroGiro(xp,yp,zp);
    I *= rho;
    Q(0) = 0.5*sqrt(1+xp(0)+yp(1)+zp(2));
    Q(1) = (yp(2)-zp(1))/(4*Q(0));
    Q(2) = (zp(0)-xp(2))/(4*Q(0));
    Q(3) = (xp(1)-yp(0))/(4*Q(0));
    Q = Q/norm(Q);
    Particles[Particles.Size()-1]->I     = I;
    Particles[Particles.Size()-1]->Q     = Q;
    double Dmax = Distance(CM,V[0])+R;
    for (size_t i=1; i<V.Size(); ++i)
    {
        if (Distance(CM,V[i])+R > Dmax) Dmax = Distance(CM,V[i])+R;
    }
    Particles[Particles.Size()-1]->Ekin = 0.0;
    Particles[Particles.Size()-1]->Erot = 0.0;
    Particles[Particles.Size()-1]->Dmax  = Dmax;
    Particles[Particles.Size()-1]->PropsReady = true;
    Particles[Particles.Size()-1]->Index = Particles.Size()-1;
        
    printf("[1;32mDEM::domain.h ConstructFromOBJ: finished[0m\n");
}
#endif
