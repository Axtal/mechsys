/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raúl D. D. Farfan             *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

// Std Lib
#include <iostream>

// MechSys
#include <mechsys/util/maps.h>
#include <mechsys/vtk/win.h>
#include <mechsys/vtk/axes.h>
#include <mechsys/vtk/sgrid.h>
#include <mechsys/vtk/sphere.h>
#include <mechsys/vtk/spheres.h>
#include <mechsys/vtk/cube.h>
#include <mechsys/vtk/text2d.h>
#include <mechsys/util/colors.h>

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
    // input
    String filekey;
    int    proc_ini     = 0;
    int    proc_fin     = 0;
    bool   with_control = true;
    bool   show_ids     = true;
    bool   shadow       = true;
    if (argc>1) filekey      =      argv[1];
    if (argc>2) proc_ini     = atoi(argv[2]);
    if (argc>3) proc_fin     = atoi(argv[3]);
    if (argc>4) with_control = atoi(argv[4]);
    if (argc>5) show_ids     = atoi(argv[5]);
    if (argc>6) shadow       = atoi(argv[6]);
    if (argc<2) throw new Fatal("filename is needed as argument");

    // filekey or filename
    bool with_proc = true;
    bool with_stp  = true;
    String ext(".res");
    String buf;  buf.Printf("%s_proc_%d_%08d.res", filekey.CStr(), proc_ini, 0);
    if (!Util::FileExists(buf))
    {
        with_proc = false;
        buf.Printf("%s_%08d.res", filekey.CStr(), 0);
        if (!Util::FileExists(buf))
        {
            with_stp = false;
            if (Util::FileExists(filekey))
            {
                buf = filekey;
                filekey.resize(buf.size()-4);
                for (size_t i=0; i<filekey.size(); ++i) filekey[i] = buf[i];
                for (size_t i=0; i<4;              ++i) ext    [i] = buf[filekey.size()+i];
                cout << "Using file <" << filekey << ext << ">\n";
            }
            else throw new Fatal("Could not find any file named '%s_proc_%d_%08d.res', '%s_%08d.res', or '%s'",filekey.CStr(),proc_ini,0, filekey.CStr(),0, filekey.CStr());
        }
    }

    // parse control file
    int            nout = 1; // number of time output
    Array<int>     N(3);     // number of cells along each axis
    Array<double>  L(6);     // limits
    Array<int>     proc;     // processor ID of each cell
    Array<String>  clrs;     // colors of each processor
    buf.Printf("%s_control%s", filekey.CStr(), ext.CStr());
    if (!Util::FileExists(buf)) with_control = false;
    if (with_control)
    {
        std::fstream fi(buf.CStr(), std::ios::in);
        if (!fi.is_open()) throw new Fatal("Could not open file <%s>",buf.CStr());
        String line,key;
        double val;
        while (!fi.eof())
        {
            std::getline (fi,line);
            std::istringstream iss(line);
            if (iss >> key >> val)
            {
                if      (key=="fkey") {}
                else if (key=="dt")   {}
                else if (key=="nout") nout = val;
                else if (key=="nx")   N[0] = val;
                else if (key=="ny")   N[1] = val;
                else if (key=="nz")   N[2] = val;
                else if (key=="lxmi") L[0] = val;
                else if (key=="lxma") L[1] = val;
                else if (key=="lymi") L[2] = val;
                else if (key=="lyma") L[3] = val;
                else if (key=="lzmi") L[4] = val;
                else if (key=="lzma") L[5] = val;
                else if (key=="proc")
                {
                    if (proc.Size()>0) throw new Fatal("Error with input file");
                    int ncells = val;
                    proc.Resize (ncells);
                    int max_proc = 0;
                    for (int i=0; i<ncells; ++i)
                    {
                        iss >> proc[i];
                        if (proc[i]>max_proc) max_proc = proc[i];
                    }
                    int nprocs = max_proc+1;
                    Colors::GetRandom (nprocs, clrs);
                }
                else throw new Fatal("Key %s is wrong",key.CStr());
            }
        }
        fi.close();
    }

    // window
    VTK::Win win;
    win.GetCamera()->SetViewUp     (0,1,0);
    win.GetCamera()->SetPosition   (0.5,0.5,3);
    win.GetCamera()->SetFocalPoint (0.5,0.5,0);
    win.GetCamera()->ParallelProjectionOn();

    // axes
    VTK::Axes axe(1.1);
    axe.XLabel().SetOrientation (0,0,0);
    axe.YLabel().SetOrientation (0,0,0);
    axe.ZLabel().SetOrientation (0,0,0);
    axe.AddTo (win);

    // grid
    if (with_control)
    {
        Array<int> Nlines(3);
        for (size_t i=0; i<3; ++i) Nlines[i] = N[i]+1;
        if (with_proc)
        {
            double dx = (L[1]-L[0])/static_cast<double>(N[0]);
            double dy = (L[3]-L[2])/static_cast<double>(N[1]);
            double dz = (L[5]-L[4])/static_cast<double>(N[2]);
            int    nx = N[0];
            int    ny = N[1];
            int    nc = N[0]*N[1]*N[2];
            for (int n=0; n<nc; ++n)
            {
                int    K    =  n / (nx*ny);
                int    J    = (n % (nx*ny)) / nx;
                int    I    = (n % (nx*ny)) % nx;
                double xmin = L[0] +  I   *dx;
                double xmax = L[0] + (I+1)*dx;
                double ymin = L[2] +  J   *dy;
                double ymax = L[2] + (J+1)*dy;
                double zmin = L[4] +  K   *dz;
                double zmax = L[4] + (K+1)*dz;
                Vec3_t cen((xmin+xmax)/2.0, (ymin+ymax)/2.0, (zmin+zmax)/2.0);
                VTK::Cube cube(cen, xmax-xmin, ymax-ymin, zmax-zmin);
                //printf("n=%d, proc[n]=%d, clrs[proc[n]]=%s\n",n,proc[n],clrs[proc[n]].CStr());
                cube.SetColor (clrs[proc[n]].CStr(), 0.1);
                cube.AddTo    (win);
            }
        }
        else
        {
            VTK::SGrid grd(Nlines.GetPtr(), L.GetPtr());
            grd.SetColor ("black", 0.2);
            //grd.ShowIds  (90,90,45,0.003,8,false);
            grd.AddTo    (win);
        }
    }

    // text
    VTK::Text2D txt(0.,0.,"");
    txt.AddTo (win);

    // spheres
    int nproc = (with_proc ? proc_fin-proc_ini+1 : 1);
    Array<String> clr;
    Colors::GetRandom (nproc, clr);
    Array<Array<VTK::Spheres*> > sph(nproc);
    for (int stp=0; stp<nout; ++stp)
    {
        int idx_proc = 0;
        String buf2;
        if (with_control) buf2.Printf("%s_proc{",filekey.CStr());
        else              buf2 = filekey;
        for (int proc=proc_ini; proc<=proc_fin; ++proc)
        {
            // filename
            if (with_proc) buf.Printf("%s_proc_%d_%08d.res", filekey.CStr(), proc, stp);
            else
            {
                if (with_stp) buf.Printf("%s_%08d.res", filekey.CStr(), stp);
                else          buf.Printf("%s%s", filekey.CStr(), ext.CStr());
            }
            if (proc==proc_fin) buf2.Printf("%s%d",  buf2.CStr(), proc);
            else                buf2.Printf("%s%d,", buf2.CStr(), proc);
            if (!Util::FileExists(buf)) cout << "Could not find file <" << buf << ">\n";

            // read data
            Table tab;
            tab.Read (buf.CStr());
            Array<double> const & idd = tab("id");
            Array<double> const & xc  = tab("xc");
            Array<double> const & yc  = tab("yc");
            Array<double> const & zc  = tab("zc");
            Array<double> const & ra  = tab("ra");
            Array<int> id(idd.Size());
            for (size_t i=0; i<idd.Size(); ++i) id[i]=static_cast<int>(idd[i]);

            // spheres
            Array<Vec3_t> X(xc.Size());
            for (size_t i=0; i<xc.Size(); ++i) X[i] = xc[i], yc[i], zc[i];
            if (shadow)
            {
                if (sph[idx_proc].Size()>0) sph[idx_proc].Last()->SetColor ("black",0.03);
                sph[idx_proc].Push (new VTK::Spheres(X,ra));
            }
            else
            {
                if (sph[idx_proc].Size()>0) sph[idx_proc].Last()->DelFrom (win);
                sph[idx_proc].Push (new VTK::Spheres(X,ra));
            }
            sph[idx_proc].Last()->Ids = id.GetPtr();
            double opac = 1.0;
            if (show_ids) { sph[idx_proc].Last()->ShowIds(0,0,0,0.003,12,false);  opac=0.9; }
            sph[idx_proc].Last()->SetColor (clr[idx_proc].CStr(),opac);
            sph[idx_proc].Last()->AddTo    (win, /*rstcam*/(stp==0));
            idx_proc++;

            for (size_t i=0; i<id.Size(); ++i)
            {
                if (id[i]==316 || id[i]==333 || id[i]==125 || id[i]==200 || id[i]==273 || id[i]==154 || id[i]==13 || id[i]==84 || id[i]==152 || id[i]==186)
                {
                    Vec3_t x(xc[i],yc[i],zc[i]);
                    VTK::Sphere s(x,ra[i]*1.1);
                    s.SetColor ("yellow",0.4);
                    s.AddTo(win, /*rstcam*/false);
                }
            }

        }

        // text
        buf2.Printf ("%s}:stp=%d (stp_fin=%d)",buf2.CStr(),stp,nout-1);
        txt.SetText (buf2.CStr());

        // window
        win.Show();
    }

    // clean up
    for (int    i=0; i<nproc;         ++i)
    for (size_t j=0; j<sph[i].Size(); ++j) delete sph[i][j];

    // end
    return 0;
}
MECHSYS_CATCH
