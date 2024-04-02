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

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
    String fkey;
    int    nprocs  = 1;
    int    stp_ini = 0;
    int    stp_fin = 0;
    if (argc>1) fkey    =      argv[1];
    if (argc>2) nprocs  = atoi(argv[2]);
    if (argc>3) stp_ini = atoi(argv[3]);
    if (argc>4) stp_fin = atoi(argv[4]);
    if (argc<2) throw new Fatal("Filekey must be provided as argument");

    printf("  running from stp=%d to stp=%d\n",stp_ini,stp_fin);
    for (int stp=stp_ini; stp<=stp_fin; ++stp)
    {
        // results
        std::map<int,double> id2xc;
        std::map<int,double> id2yc;
        std::map<int,double> id2zc;
        std::map<int,double> id2ra;
        std::map<int,double> id2vx;
        std::map<int,double> id2vy;
        std::map<int,double> id2vz;

        // join
        int max_id = 0;
        printf("  ====> Step %d, Proc: ",stp);
        for (int proc=0; proc<nprocs; ++proc)
        {
            // read data
            String buf;  buf.Printf("%s_proc_%d_%08d.res", fkey.CStr(), proc, stp);
            Table tab;
            tab.Read (buf.CStr());

            // auxiliar variables
            Array<double> const & idd = tab("id");
            Array<double> const & ctd = tab("ct");
            Array<int> id(idd.Size());
            Array<int> ct(ctd.Size());
            for (size_t k=0; k<idd.Size(); ++k) { id[k]=static_cast<int>(idd[k]);  ct[k]=static_cast<int>(ctd[k]); }
            Array<double> const & xc = tab("xc");
            Array<double> const & yc = tab("yc");
            Array<double> const & zc = tab("zc");
            Array<double> const & ra = tab("ra");
            Array<double> const & vx = tab("vx");
            Array<double> const & vy = tab("vy");
            Array<double> const & vz = tab("vz");

            // results
            printf("%d ",proc);
            for (int i=0; i<static_cast<int>(id.Size()); ++i)
            {
                if (ct[i]==-4) continue; // ignore particles not in proc
                if (id2xc.count(id[i])>0) { if (fabs(id2xc[id[i]]-xc[i])>1.0e-15) throw new Fatal("\nId=%d, ct=%d: xc=%g in %s should be equal to %g",id[i],ct[i],xc[i],buf.CStr(),id2xc[id[i]]); }
                if (id2yc.count(id[i])>0) { if (fabs(id2yc[id[i]]-yc[i])>1.0e-15) throw new Fatal("\nId=%d, ct=%d: yc=%g in %s should be equal to %g",id[i],ct[i],yc[i],buf.CStr(),id2yc[id[i]]); }
                if (id2zc.count(id[i])>0) { if (fabs(id2zc[id[i]]-zc[i])>1.0e-15) throw new Fatal("\nId=%d, ct=%d: zc=%g in %s should be equal to %g",id[i],ct[i],zc[i],buf.CStr(),id2zc[id[i]]); }
                if (id2ra.count(id[i])>0) { if (fabs(id2ra[id[i]]-ra[i])>1.0e-15) throw new Fatal("\nId=%d, ct=%d: ra=%g in %s should be equal to %g",id[i],ct[i],ra[i],buf.CStr(),id2ra[id[i]]); }
                if (id2vx.count(id[i])>0) { if (fabs(id2vx[id[i]]-vx[i])>1.0e-15) throw new Fatal("\nId=%d, ct=%d: vx=%g in %s should be equal to %g",id[i],ct[i],vx[i],buf.CStr(),id2vx[id[i]]); }
                if (id2vy.count(id[i])>0) { if (fabs(id2vy[id[i]]-vy[i])>1.0e-15) throw new Fatal("\nId=%d, ct=%d: vy=%g in %s should be equal to %g",id[i],ct[i],vy[i],buf.CStr(),id2vy[id[i]]); }
                if (id2vz.count(id[i])>0) { if (fabs(id2vz[id[i]]-vz[i])>1.0e-15) throw new Fatal("\nId=%d, ct=%d: vz=%g in %s should be equal to %g",id[i],ct[i],vz[i],buf.CStr(),id2vz[id[i]]); }
                id2xc[id[i]] = xc[i];
                id2yc[id[i]] = yc[i];
                id2zc[id[i]] = zc[i];
                id2ra[id[i]] = ra[i];
                id2vx[id[i]] = vx[i];
                id2vy[id[i]] = vy[i];
                id2vz[id[i]] = vz[i];
                if (id[i]>max_id) max_id = id[i];
            }
        }
        printf("\n");

        // header
        Array<String> keys("id", "xc", "yc", "zc", "ra", "vx", "vy", "vz", "ct");
        std::ostringstream oss;
        oss << Util::_6 << keys[0];
        for (size_t i=1; i<keys.Size()-1; ++i) { oss << Util::_8s << keys[i]; }
        oss << Util::_6 << keys.Last() << "\n";

        // values
        for (int the_id=0; the_id<=max_id; ++the_id)
        {
            oss << Util::_6  << the_id;
            oss << Util::_8s << id2xc[the_id];
            oss << Util::_8s << id2yc[the_id];
            oss << Util::_8s << id2zc[the_id];
            oss << Util::_8s << id2ra[the_id];
            oss << Util::_8s << id2vx[the_id];
            oss << Util::_8s << id2vy[the_id];
            oss << Util::_8s << id2vz[the_id];
            oss << Util::_6  << -1;
            oss << "\n";
        }

        // open file and save data
        String buf;  buf.Printf("%s_join_%08d.res", fkey.CStr(), stp);
        std::ofstream of(buf.CStr(), std::ios::out);
        of << oss.str();
        of.close();
    }

    return 0;
}
MECHSYS_CATCH
