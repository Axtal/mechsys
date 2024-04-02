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

#ifndef MECHSYS_INPFILE_H
#define MECHSYS_INPFILE_H

// Std Lib
#include <fstream>

// wxWidgets
#ifdef USE_WXWIDGETS
  #include <mechsys/gui/common.h>
  #include <mechsys/gui/wxdict.h>
  #include <mechsys/gui/wxsipair.h>
  #include <mechsys/gui/wxarrayint.h>
#endif

// MechSys
#include <mechsys/util/maps.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>
#include <mechsys/linalg/matvec.h>
#include <mechsys/fem/element.h>   // for PROB
#include <mechsys/fem/geomelem.h>  // for GEOM
#include <mechsys/matfile.h>
#include <mechsys/models/model.h>

/*
struct PathIncs
{
    double dsx, dsy, dsz, dsxy, dsyz, dszx; // stress increments
    double dex, dey, dez, dexy, deyz, dezx; // strain increments (percentage)
    double lode, dp;                        // path given Lode angle (deg), dpoct and dez (percentage)
    bool   zPath;                           // use lode, dp and dez ?
    int    ninc;                            // number of increments for this path. -1 => use general
    double k;                               // path given Lode, k=dqoct/dpoct, and dez
    bool   kPath;                           // with k=Dq/Dp
    double dpw, dSw;                        // increment of pore-water pressure and water saturation
    bool   HasDpw, HasDSw;                  // has dpw or dSw ?
    PathIncs () : dsx(0.),dsy(0.),dsz(0.),dsxy(0.),dsyz(0.),dszx(0.), 
                  dex(0.),dey(0.),dez(0.),dexy(0.),deyz(0.),dezx(0.),
                  lode(0.),dp(0.),zPath(false),ninc(-1),k(0.),kPath(false),
                  dpw(0.),dSw(0.),HasDpw(false),HasDSw(false) {}
};
*/

#ifdef USE_WXWIDGETS
class InpFile : public wxWindow
#else
class InpFile
#endif
{
public:
    // Constructor & Destructor
#ifdef USE_WXWIDGETS
     InpFile (wxFrame * Parent);
#else
     InpFile ();
#endif
    ~InpFile ();

    // Methods
    void Defaults    ();
    void Read        (char const * FileName);
    void SetPrmsInis (MatFile const & Mat, bool ForceGTY=false);
    void GetIncs     (int PathKey, double Div, Vec_t & dsig, Vec_t & deps, Array<bool> & PrescDeps, double & dpw, double & dSw, bool & PrescDpw, bool & PrescDSw) const;
    void SetSUp      (Model const * Mdl, Model::StressUpdate::pDbgFun pFun=NULL, void * UserData=NULL) const;
    void SetSolFlags (SDPair & Flags) const;

    // Data
    int    ninc;        ///<  1 general number of increments (for all load-unload paths)
    bool   cdrift;      ///<  2 correct YS drift
    double stol;        ///<  3 local error tolerance
    bool   ssout;       ///<  4 output substeps ?
    bool   ctetg;       ///<  5 Constant stiffness (linear) ?
    bool   fem;         ///<  6 use one Hex8 FEM element instead of point integration
    bool   dyn;         ///<  7 dynamic analysis ?
    bool   hm;          ///<  8 HydroMech ?
    double tf;          ///<  9 final time
    double dt;          ///< 10 time step
    double dtout;       ///< 11 output time step
    double tsw;         ///< 12 switch time (for dynamic simulation)
    int    ndiv;        ///< 13 mesh number of divisions
    int    nip;         ///< 14 Number of integration points in element
    bool   o2;          ///< 15 quadratic elements ?
    bool   ray;         ///< 16 Rayleigh damping ?
    double am;          ///< 17 Damping Am
    double ak;          ///< 18 Damping Ak
    bool   rk;          ///< 19 Runge-Kutta instead of GN22 ?
    String rkscheme;    ///< 20 Runge-Kutta scheme 
    double rkstol;      ///< 21 Runge-Kutta tolerance
    String refdat;      ///< 22 reference data file
    String refsim;      ///< 23 reference simulation file
    String refana;      ///< 24 reference analytical solution file
    int    idxvert1;    ///< 25 index of vertex # 1 for output
    int    idxvert2;    ///< 26 index of vertex # 1 for output
    int    idxvert3;    ///< 27 index of vertex # 1 for output
    double optdbl1;     ///< 28 optional double 1
    double optdbl2;     ///< 29 optional double 2
    double optdbl3;     ///< 30 optional double 3
    bool   hasoptdbl1, hasoptdbl2, hasoptdbl3;
    int    nldt_nsml;   ///< 31 nonlinear timesteps Nsml
    int    nldt_nn;     ///< 32 nonlinear timesteps N
    int    nldt_n;      ///< 33 nonlinear timesteps n
    double nldt_ll;     ///< 34 nonlinear timesteps denominator
    int    nldt_sch;    ///< 35 nonlinear timesteps scheme (0 or 1)
    double nldt_m;      ///< 36 nonlinear timesteps multiplier for larger timesteps in sch==1
    int    maxit;       ///< 37 max num of iterations
    double tolr;        ///< 38 tolerance for residual
    String fnkey;       ///< 39 filename key
    double pcam0;       ///< 40 pcam0
    bool   haspcam0;    ///< has pcam0 ?
    String scheme;      ///< 41 solver scheme
    bool   vtufile;     ///< 42 write vtu file ?
    String suscheme;    ///< 43 stress-update scheme
    double sustol;      ///< 44 stress-update STOL
    String surkscheme;  ///< 45 stress-update RK scheme
    size_t dcmaxit;     ///< 46 drift correction max iterations
    double dcftol;      ///< 47 drift correction f tolerance
    double pw0;         ///< 48 pw0
    bool   haspw0;      ///< has pw0 ?
    bool   rkdyncte;    ///< 49 rk scheme dyn cte M and C
    bool   uwp;         ///< 50 u-w-p TPM formulation ?

    Array<String> AllUKeys, AllUKeysBCF;
    Array<String> AllFKeys, AllFKeysBCF;
    Array<String> AllEKeys, AllEKeysBCF;

    // Additional data
    Dict * Prms; ///< parameters (set by SetMat)
    Dict * Inis; ///< initial values (set by SetMat)

#ifdef USE_WXWIDGETS
    // Methods
    void Sync (bool Dat2Ctrl=false) { if (Dat2Ctrl) TransferDataToWindow(); else TransferDataFromWindow(); } ///< Synchronise (validate/transfer) data in controls

    // Data
    wxAuiManager       Aui;    ///< Aui manager
    wxString           LstDir; ///< Last accessed directory
    String             FName;  ///< Input file (.inp) filename
    GUI::WxDictTable * Path;   ///< Path increments
    GUI::WxDict      * GPath;  ///< Grid for editing the path

    // Additional data
    GUI::WxDictTable     * Prps;       ///< elements properties
    GUI::WxArrayIntTable * OutNods;    ///< output nodes
    GUI::WxDict          * GPrps;      ///< grid: elements properties
    GUI::WxArrayInt      * GOutNods;   ///< grid: output nodes
    GUI::WxSIPairTable   * Tag2MatID;  ///< element tag to material ID
    GUI::WxSIPairTable   * Tag2XMatID; ///< element tag to extra (flow) material ID
    GUI::WxSIPair        * GMatId2Tag; ///< material ID to element tag

    // Boundary conditions
    Array<GUI::WxDictTable*> Stages;  ///< boundary conditions
    //Array<GUI::WxDict*>      GStages; ///< grid: boundary conditions

    // Events
    void OnLoad (wxCommandEvent & Event);
    void OnSave (wxCommandEvent & Event);
    DECLARE_EVENT_TABLE();

#else
    // Data
    Dict * Path; ///< Path increments

    // Additional data
    Dict       * Prps;       ///< elements properties
    Array<int> * OutNods;    ///< output nodes
    SIPair     * Tag2MatID;  ///< element tag to material ID
    SIPair     * Tag2XMatID; ///< element tag to extra (flow) material ID

    // Boundary conditions
    Array<Dict*> Stages; ///< boundary conditions
#endif

#ifdef USE_BOOST_PYTHON
    BPy::tuple PyReadPrmIni (char const * MatFN, char const * InpFN, int Tag)
    {
        BPy::dict prm, ini;
        MatFile   mat;
        mat.Read    (MatFN);
        Read        (InpFN);
        SetPrmsInis (mat);
        Prms->PyGet (Tag, prm);
        Inis->PyGet (Tag, ini);
        return BPy::make_tuple (prm, ini);
    }
#endif
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


#ifdef USE_WXWIDGETS

inline InpFile::~InpFile ()
{
    Aui.UnInit ();
    delete Prms;
    delete Inis;
}

#else

inline InpFile::InpFile ()
{
    Defaults ();
    Path       = new Dict;
    Prps       = new Dict;
    OutNods    = new Array<int>;
    Tag2MatID  = new SIPair;
    Tag2XMatID = new SIPair;
    Prms       = new Dict;
    Inis       = new Dict;
}

inline InpFile::~InpFile ()
{
    delete Path;
    delete Prps;
    delete OutNods;
    delete Tag2MatID;
    delete Tag2XMatID;
    delete Prms;
    delete Inis;
    for (size_t i=0; i<Stages.Size(); ++i) delete Stages[i];
}

#endif

inline void InpFile::Defaults ()
{
    ninc       = -1;     //  1
    cdrift     = false;  //  2
    stol       = -1;     //  3
    ssout      = false;  //  4
    ctetg      = false;  //  5
    fem        = false;  //  6
    dyn        = false;  //  7
    hm         = false;  //  8
    tf         = -1;     //  9
    dt         = -1;     // 10
    dtout      = -1;     // 11
    tsw        = -1;     // 12
    ndiv       = -1;     // 13
    nip        = -1;     // 14
    o2         = false;  // 15
    ray        = false;  // 16
    am         = -1;     // 17
    ak         = -1;     // 18
    rk         = false;  // 19
    rkscheme   = "";     // 20
    rkstol     = -1;     // 21
    refdat     = "";     // 22
    refsim     = "";     // 23
    refana     = "";     // 24
    idxvert1   = -1;     // 25
    idxvert2   = -1;     // 26
    idxvert3   = -1;     // 27
    optdbl1    = 0;      // 28
    optdbl2    = 0;      // 29
    optdbl3    = 0;      // 30
    hasoptdbl1 = false; hasoptdbl2=false; hasoptdbl3=false;
    nldt_nsml  = -1;     // 31
    nldt_nn    = -1;     // 32
    nldt_n     = -1;     // 33
    nldt_ll    = -1;     // 34
    nldt_sch   = -1;     // 35
    nldt_m     = -1;     // 36
    maxit      = -1;     // 37
    tolr       = -1;     // 38
    fnkey      = "";     // 39
    pcam0      = 0;      // 40
    haspcam0   = false;
    scheme     = "";     // 41
    vtufile    = false;  // 42
    suscheme   = "";     // 43
    sustol     = -1;     // 44
    surkscheme = "";     // 45
    dcmaxit    = 0;      // 46
    dcftol     = -1;     // 47
    pw0        = 0;      // 48
    haspw0     = false;
    rkdyncte   = true;   // 49
    uwp        = false;  // 50

    // U and F keys
    for (FEM::ElementVarKeys_t::const_iterator it=FEM::ElementVarKeys.begin(); it!=FEM::ElementVarKeys.end(); ++it)
    {
        // U keys
        Array<String> keys;
        Util::Keys2Array (it->second.first, keys);
        for (size_t i=0; i<keys.Size(); ++i)
        {
            AllUKeys.XPush (keys[i]);
            keys[i].ToUpper();
            AllUKeysBCF.XPush (keys[i]);
        }

        // F keys
        keys.Clear();
        Util::Keys2Array (it->second.second, keys);
        for (size_t i=0; i<keys.Size(); ++i)
        {
            AllFKeys.XPush (keys[i]);
            keys[i].ToUpper();
            AllFKeysBCF.XPush (keys[i]);
        }
    }
    /*
    std::cout << "AllUKeys    = " << AllUKeys    << std::endl;
    std::cout << "AllUKeysBCF = " << AllUKeysBCF << std::endl;
    std::cout << "AllFKeys    = " << AllFKeys    << std::endl;
    std::cout << "AllFKeysBCF = " << AllFKeysBCF << std::endl;
    */

    // Extra keys
    for (FEM::ElementExtraKeys_t::const_iterator it=FEM::ElementExtraKeys.begin(); it!=FEM::ElementExtraKeys.end(); ++it)
    {
        for (size_t i=0; i<it->second.Size(); ++i)
        {
            AllEKeys   .XPush (it->second[i]);
            AllEKeysBCF.XPush (it->second[i].ToUpperCpy());
        }
    }
}

inline void InpFile::Read (char const * FileName)
{
    // parse input file
    std::fstream inp_file(FileName, std::ios::in);
    if (!inp_file.is_open()) throw new Fatal("InpFile::Read: Could not open file <%s>",FileName);
    bool   reading_path   = false;
    bool   reading_eprps  = false;
    bool   reading_stages = false;
    bool   reading_bcs    = false;
    int    npath          = 0;
    int    nelemprps      = 0;
    int    nstages        = 0;
    int    nbcs           = -1;
    int    ndat           = -1;
    size_t line_num       = 1;
    int    idxdat         = 0;
    int    idxbcs         = 0;
    int    idxpath        = 0;
    int    idxeprps       = 0;
    int    idxstage       = 0;
    int    elemtag        = -1;
    int    bcstag         = 0;
    SDPair elemprps;
    SDPair bcs;
    Path      -> clear();
    Prps      -> clear();
    OutNods   -> Resize(0);
    Tag2MatID -> clear();
    Tag2XMatID-> clear();
    Prms      -> clear();
    Inis      -> clear();
    for (size_t i=0; i<Stages.Size(); ++i) delete Stages[i];
    Stages.Resize(0);
    while (!inp_file.eof())
    {
        String line,key,equal,str_val;
        double val;
        std::getline (inp_file,line);
        std::istringstream iss(line);
        if (iss >> key >> equal >> str_val)
        {
            val = atof(str_val.CStr());
            if (key[0]=='#') { line_num++; continue; }
            if (reading_path)
            {
                if      (key=="ndat") ndat = atoi(str_val.CStr());
                else if (ndat<0) throw new Fatal("InpFile::Read: Reading path. Error in file <%s> at line # %d: key 'ndat' must come before data. '%s' is in the wrong place",FileName,line_num,key.CStr());
                else if (key=="kcam")  { Path->Set(idxpath, "kcam" , val     ); idxdat++; }
                else if (key=="dpcam") { Path->Set(idxpath, "dpcam", val     ); idxdat++; }
                else if (key=="lode")  { Path->Set(idxpath, "lode" , val     ); idxdat++; if (val<30. || val>90.) throw new Fatal("InpFile::Read: Error in file <%s> at line # %d: Lode angle alpha must be inside [30,90]. Alpha==%g is invalid",FileName,line_num,val); }
                else if (key=="dex")   { Path->Set(idxpath, "dex"  , val/100.); idxdat++; }
                else if (key=="dey")   { Path->Set(idxpath, "dey"  , val/100.); idxdat++; }
                else if (key=="dez")   { Path->Set(idxpath, "dez"  , val/100.); idxdat++; }
                else if (key=="dexy")  { Path->Set(idxpath, "dexy" , val/100.); idxdat++; }
                else if (key=="deyz")  { Path->Set(idxpath, "deyz" , val/100.); idxdat++; }
                else if (key=="dezx")  { Path->Set(idxpath, "dezx" , val/100.); idxdat++; }
                else if (key=="dsx")   { Path->Set(idxpath, "dsx"  , val     ); idxdat++; }
                else if (key=="dsy")   { Path->Set(idxpath, "dsy"  , val     ); idxdat++; }
                else if (key=="dsz")   { Path->Set(idxpath, "dsz"  , val     ); idxdat++; }
                else if (key=="dsxy")  { Path->Set(idxpath, "dsxy" , val     ); idxdat++; }
                else if (key=="dsyz")  { Path->Set(idxpath, "dsyz" , val     ); idxdat++; }
                else if (key=="dszx")  { Path->Set(idxpath, "dszx" , val     ); idxdat++; }
                else if (key=="dpw")   { Path->Set(idxpath, "dpw"  , val     ); idxdat++; }
                else if (key=="dSw")   { Path->Set(idxpath, "dSw"  , val     ); idxdat++; }
                else if (key=="ninc")  { Path->Set(idxpath, "ninc" , val     ); idxdat++; }
                else throw new Fatal("InpFile::Read: Reading path. Error in file <%s> at line # %d when reading data of Path # %d. Key==%s is invalid or in the wrong place",FileName,line_num,idxpath,key.CStr());
                if (idxdat==ndat)
                {
                    ndat   = -1;
                    idxdat = 0;
                    idxpath++;
                    if (idxpath==npath) reading_path = false;
                }
            }
            else if (reading_eprps)
            {
                if      (key=="ndat") ndat = atoi(str_val.CStr());
                else if (ndat<0) throw new Fatal("InpFile::Read: Reading elements properties. Error in file <%s> at line # %d: key 'ndat' must come before data. '%s' is in the wrong place",FileName,line_num,key.CStr());
                else if (key=="elemtag") { elemtag = atoi(str_val.CStr());                       idxdat++; }
                else if (key=="prob")    { elemprps.Set (key.CStr(), FEM::PROB(str_val.CStr())); idxdat++; }
                else if (key=="geom")    { elemprps.Set (key.CStr(), FEM::GEOM(str_val.CStr())); idxdat++; }
                else if (key=="psa")     { elemprps.Set (key.CStr(), atof(str_val.CStr()));      idxdat++; }
                else if (key=="pse")     { elemprps.Set (key.CStr(), atof(str_val.CStr()));      idxdat++; }
                else if (key=="fra")     { elemprps.Set (key.CStr(), atof(str_val.CStr()));      idxdat++; }
                else if (key=="d2d")     { elemprps.Set (key.CStr(), atof(str_val.CStr()));      idxdat++; }
                else if (key=="d3d")     { elemprps.Set (key.CStr(), atof(str_val.CStr()));      idxdat++; }
                else if (key=="rho")     { elemprps.Set (key.CStr(), atof(str_val.CStr()));      idxdat++; }
                else if (key=="geosta")  { elemprps.Set (key.CStr(), atof(str_val.CStr()));      idxdat++; }
                else if (key=="pospw")   { elemprps.Set (key.CStr(), atof(str_val.CStr()));      idxdat++; }
                else if (key=="K0")      { elemprps.Set (key.CStr(), atof(str_val.CStr()));      idxdat++; }
                else if (key=="surf")    { elemprps.Set (key.CStr(), atof(str_val.CStr()));      idxdat++; }
                else if (key=="water")   { elemprps.Set (key.CStr(), atof(str_val.CStr()));      idxdat++; }
                else if (key=="E")       { elemprps.Set (key.CStr(), atof(str_val.CStr()));      idxdat++; }
                else if (key=="A")       { elemprps.Set (key.CStr(), atof(str_val.CStr()));      idxdat++; }
                else if (key=="Izz")     { elemprps.Set (key.CStr(), atof(str_val.CStr()));      idxdat++; }
                else throw new Fatal("InpFile::Read: Reading elements properties. Error in file <%s> at line # %d when reading data of Properties # %d. Key==%s is invalid or in the wrong place",FileName,line_num,idxeprps,key.CStr());
                if (idxdat==ndat)
                {
                    Prps->Set (elemtag, elemprps);
                    ndat    = -1;
                    idxdat  = 0;
                    elemtag = -1;
                    elemprps.clear();
                    idxeprps++;
                    if (idxeprps==nelemprps) reading_eprps = false;
                }
            }
            else if (reading_stages && !reading_bcs)
            {
#ifdef USE_WXWIDGETS
                if (key=="nbcs") { nbcs = atoi(str_val.CStr());  Stages.Push(new GUI::WxDictTable);  reading_bcs=true; }
#else
                if (key=="nbcs") { nbcs = atoi(str_val.CStr());  Stages.Push(new Dict);  reading_bcs=true; }
#endif
                else throw new Fatal("InpFile::Read: Reading boundary conditions (stages). Error in file <%s> at line # %d: key '%s' is in the wrong place",FileName,line_num,key.CStr());
            }
            else if (reading_bcs)
            {
                if      (key=="ndat") ndat = atoi(str_val.CStr());
                else if (ndat<0) throw new Fatal("InpFile::Read: Reading boundary conditions (stages). Error in file <%s> at line # %d: key 'ndat' must come after 'nbcs' and before data. '%s' is in the wrong place",FileName,line_num,key.CStr());
                else if (key=="tag")           { bcstag = atoi(str_val.CStr());              idxdat++; }
                else if (AllUKeys   .Has(key)) { bcs.Set (key.CStr(), atof(str_val.CStr())); idxdat++; }
                else if (AllFKeys   .Has(key)) { bcs.Set (key.CStr(), atof(str_val.CStr())); idxdat++; }
                else if (AllEKeys   .Has(key)) { bcs.Set (key.CStr(), atof(str_val.CStr())); idxdat++; }
                else if (AllUKeysBCF.Has(key)) { bcs.Set (key.CStr(), atof(str_val.CStr())); idxdat++; }
                else if (AllFKeysBCF.Has(key)) { bcs.Set (key.CStr(), atof(str_val.CStr())); idxdat++; }
                else if (AllEKeysBCF.Has(key)) { bcs.Set (key.CStr(), atof(str_val.CStr())); idxdat++; }
                else if (key=="fgrav")         { bcs.Set (key.CStr(), atof(str_val.CStr())); idxdat++; }
                else throw new Fatal("InpFile::Read: Reading boundary conditions (stages). Error in file <%s> at line # %d when reading data of Stage # %d. Key==%s is invalid or in the wrong place",FileName,line_num,idxstage,key.CStr());
                if (idxdat==ndat)
                {
                    Stages[idxstage]->Set (bcstag, bcs);
                    ndat   = -1;
                    idxdat = 0;
                    bcstag = 0;
                    bcs.clear ();
                    idxbcs++;
                    if (idxbcs==nbcs)
                    {
                        reading_bcs = false;
                        nbcs        = -1;
                        idxbcs      = 0;
                        idxstage++;
                        if (idxstage==nstages) reading_stages = false;
                    }
                }
            }
            else
            {
                if      (key=="ninc")       ninc      = atoi(str_val.CStr());       //  1
                else if (key=="cdrift")     cdrift    = (bool)atoi(str_val.CStr()); //  2
                else if (key=="stol")       stol      = val;                        //  3
                else if (key=="ssout")      ssout     = (bool)atoi(str_val.CStr()); //  4
                else if (key=="ctetg")      ctetg     = (bool)atoi(str_val.CStr()); //  5
                else if (key=="fem")        fem       = val;                        //  6
                else if (key=="dyn")        dyn       = (bool)atoi(str_val.CStr()); //  7
                else if (key=="hm")         hm        = (bool)atoi(str_val.CStr()); //  8
                else if (key=="tf")         tf        = val;                        //  9
                else if (key=="dt")         dt        = val;                        // 10
                else if (key=="dtout")      dtout     = val;                        // 11
                else if (key=="tsw")        tsw       = val;                        // 12
                else if (key=="ndiv")       ndiv      = atoi(str_val.CStr());       // 13
                else if (key=="nip")        nip       = atoi(str_val.CStr());       // 14
                else if (key=="o2")         o2        = (bool)atoi(str_val.CStr()); // 15
                else if (key=="ray")        ray       = (bool)atoi(str_val.CStr()); // 16
                else if (key=="am")         am        = val;                        // 17
                else if (key=="ak")         ak        = val;                        // 18
                else if (key=="rk")         rk        = (bool)atoi(str_val.CStr()); // 19
                else if (key=="rkscheme")   rkscheme  = str_val;                    // 20
                else if (key=="rkstol")     rkstol    = val;                        // 21
                else if (key=="refdat")     refdat    = str_val;                    // 22
                else if (key=="refsim")     refsim    = str_val;                    // 23
                else if (key=="refana")     refana    = str_val;                    // 24
                else if (key=="idxvert1")   idxvert1  = atoi(str_val.CStr());       // 25
                else if (key=="idxvert2")   idxvert2  = atoi(str_val.CStr());       // 26
                else if (key=="idxvert3")   idxvert3  = atoi(str_val.CStr());       // 27
                else if (key=="optdbl1")  { optdbl1   = val;   hasoptdbl1=true; }   // 28
                else if (key=="optdbl2")  { optdbl2   = val;   hasoptdbl2=true; }   // 29
                else if (key=="optdbl3")  { optdbl3   = val;   hasoptdbl3=true; }   // 30
                else if (key=="nldt_nsml")  nldt_nsml = atoi(str_val.CStr());       // 31
                else if (key=="nldt_nn")    nldt_nn   = atoi(str_val.CStr());       // 32
                else if (key=="nldt_n")     nldt_n    = atoi(str_val.CStr());       // 33
                else if (key=="nldt_ll")    nldt_ll   = val;                        // 34
                else if (key=="nldt_sch")   nldt_sch  = atoi(str_val.CStr());       // 35
                else if (key=="nldt_m")     nldt_m    = val;                        // 36
                else if (key=="maxit")      maxit     = atoi(str_val.CStr());       // 37
                else if (key=="tolr")       tolr      = val;                        // 38
                else if (key=="fnkey")      fnkey     = str_val;                    // 39
                else if (key=="pcam0")    { pcam0     = val;     haspcam0 = true; } // 40
                else if (key=="scheme")     scheme    = str_val;                    // 41
                else if (key=="vtufile")    vtufile   = (int)val;                   // 42
                else if (key=="suscheme")   suscheme  = str_val;                    // 43
                else if (key=="sustol")     sustol    = val;                        // 44
                else if (key=="surkscheme") surkscheme= str_val;                    // 45
                else if (key=="dcmaxit")    dcmaxit   = (int)val;                   // 46
                else if (key=="dcftol")     dcftol    = val;                        // 47
                else if (key=="pw0")      { pw0       = val;     haspw0 = true; }   // 48
                else if (key=="rkdyncte")   rkdyncte  = (bool)atoi(str_val.CStr()); // 49
                else if (key=="uwp")        uwp       = (bool)atoi(str_val.CStr()); // 50
                else if (key=="npath")    { npath     = (int)val;  reading_path   = true; }
                else if (key=="nelemprps"){ nelemprps = (int)val;  reading_eprps  = true; }
                else if (key=="nstages")  { nstages   = (int)val;  reading_stages = true; }
                else if (key=="matids" || key=="xmatids")
                {
                    String left, right, str_id, str_tag;
                    line.Split (left, right, "=");
                    std::istringstream subiss(right);
                    bool pair_found = false;
                    while (subiss >> str_id >> str_tag)
                    {
                        pair_found = true;
                        int id  = atoi(str_id.CStr());
                        int tag = atoi(str_tag.CStr());
                        if (id<0 || tag>=0) throw new Fatal("InpFile::Read: Error in file <%s> @ line # %d with Key==%s. Material ids must be zero or positive and element tags must be negative. Ex.: matids = 0 -1  1 -2  2 -3",FileName,line_num,key.CStr());
                        if (key=="xmatids") Tag2XMatID->Set (str_tag.CStr(), id);
                        else                Tag2MatID ->Set (str_tag.CStr(), id);
                    }
                    if (!pair_found) throw new Fatal("InpFile::Read: Error in file <%s> @ line # %d with Key==%s. Material ids must be zero or positive and element tags must be negative. Ex.: matids = 0 -1  1 -2  2 -3",FileName,line_num,key.CStr());
                }
                else if (key=="outnods")
                {
                    String left, right, nod;
                    line.Split (left, right, "=");
                    std::istringstream subiss(right);
                    while (subiss >> nod) OutNods->Push (atoi(nod.CStr()));
                }
                else throw new Fatal("InpFile::Read: Error in file <%s> @ line # %d: Key==%s in invalid (as general data)",FileName,line_num,key.CStr());
            }
        }
        line_num++;
    }
    if ((size_t)idxpath!=Path->Keys.Size()) throw new Fatal("InpFile::Read: Error in file <%s>: Not all Path data could be read for npath==%zd",FileName,npath);

    // filename key
    if (fnkey=="")
    {
        String buf(FileName);
        buf.GetFNKey (fnkey);
    }

#ifdef USE_WXWIDGETS
    Sync (/*Dat2Ctrl*/true);
    GPath      -> ReBuild ();
    GPrps      -> ReBuild ();
    GOutNods   -> ReBuild ();
    GMatId2Tag -> ReBuild ();
#endif
}

inline void InpFile::SetPrmsInis (MatFile const & Mat, bool ForceGTY)
{
    for (size_t i=0; i<Tag2MatID->Keys.Size(); ++i)
    {
        int            tag  = atoi(Tag2MatID->Keys[i].CStr());
        int            id   = (*Tag2MatID)(Tag2MatID->Keys[i]);
        SDPair const & prms = (*Mat.ID2Prms)(id);
        SDPair const & inis = (*Mat.ID2Inis)(id);
        Prms->Set (tag, prms);
        Inis->Set (tag, inis);
        if (ForceGTY)
        {
            if (Prps->Keys.Size()==0) throw new Fatal("InpFile::SetPrmsInis: Dictionary of properties (Prps) is empty => ForceGTY (geometry type) cannot be applied. You may have to specify 'nelemprps = 1' in your Input file.");
            SDPair const & prps = (*Prps)(tag);
            Array<String> gtypes("d3d", "d2d", "psa", "pse", "fra");
            for (size_t i=0; i<gtypes.Size(); ++i)
            {
                if (prps.HasKey(gtypes[i])) (*Prms)(tag).Set (gtypes[i].CStr(), 1.0);
            }
        }
        if (haspcam0) (*Inis)(tag).Set ("sx sy sz", -pcam0, -pcam0, -pcam0);
    }
    for (size_t i=0; i<Tag2XMatID->Keys.Size(); ++i)
    {
        int            tag  = atoi(Tag2XMatID->Keys[i].CStr());
        int            id   = (*Tag2XMatID)(Tag2XMatID->Keys[i]);
        SDPair         prms = (*Mat.ID2Prms)(id);
        SDPair const & inis = (*Mat.ID2Inis)(id);
        double nam = prms("name");
        prms.Del("name");
        prms.Set("xname", nam);
        Prms->Set (tag, prms);
        Inis->Set (tag, inis);
        if (haspw0) (*Inis)(tag).Set ("pw", pw0);
    }
}

inline void InpFile::GetIncs (int PathKey, double Div, Vec_t & dsig, Vec_t & deps, Array<bool> & PrescDeps, double & dpw, double & dSw, bool & PrescDpw, bool & PrescDSw) const
{
    // number of increments
    SDPair const & path = (*Path)(Path->Keys[PathKey]);

    // set prescribed increments
    set_to_zero (dsig);
    set_to_zero (deps);
    PrescDeps = false,false,false, false,false,false;
    dpw       = 0.0;
    dSw       = 0.0;
    PrescDpw  = false;
    PrescDSw  = false;
    if (path.HasKey("dsx")) { dsig(0) = path("dsx")/Div;                       }
    if (path.HasKey("dsy")) { dsig(1) = path("dsy")/Div;                       }
    if (path.HasKey("dsz")) { dsig(2) = path("dsz")/Div;                       }
    if (path.HasKey("dex")) { deps(0) = path("dex")/Div;                       }
    if (path.HasKey("dey")) { deps(1) = path("dey")/Div;  PrescDeps[1] = true; }
    if (path.HasKey("dez")) { deps(2) = path("dez")/Div;  PrescDeps[2] = true; }
    if (path.HasKey("dpw")) { dpw     = path("dpw")/Div;  PrescDpw     = true; }
    if (path.HasKey("dSw")) { dSw     = path("dSw")/Div;  PrescDSw     = true; }
}

//inline void InpFile::SetSolver (FEM::Solver & Sol) const
//{
    //if (ssout    ) Sol.SSOut  = ssout;
    //if (ctetg    ) Sol.CteTg  = ctetg;
    //if (hm       ) Sol.DampTy = FEM::Solver::HMCoup_t;
    //if (ray      ) Sol.DampTy = FEM::Solver::Rayleigh_t;
    //if (am    >=0) Sol.DampAm = am;
    //if (ak    >=0) Sol.DampAk = ak;
    //if (rk       )
    //{
        //Sol.DScheme = FEM::Solver::RK_t;
        //if (rkscheme!="") Sol.RKScheme = rkscheme;
        //if (rkstol   >=0) Sol.RKSTOL   = rkstol;
    //}
    //if (maxit  >0 ) Sol.MaxIt = maxit;
    //if (tolr   >=0) Sol.TolR  = tolr;
    //if (scheme!="") Sol.SetScheme (scheme.CStr());
//
    //Sol.NLSteps.clear();
    //if (nldt_nsml>0) Sol.NLSteps.Set("nsml",(double)nldt_nsml);
    //if (nldt_nn  >0) Sol.NLSteps.Set("nn"  ,(double)nldt_nn  );
    //if (nldt_n   >0) Sol.NLSteps.Set("n"   ,(double)nldt_n   );
    //if (nldt_ll  >0) Sol.NLSteps.Set("ll"  ,(double)nldt_ll  );
    //if (nldt_sch >0) Sol.NLSteps.Set("sch" ,(double)nldt_sch );
    //if (nldt_m   >0) Sol.NLSteps.Set("m"   ,(double)nldt_m   );
//}
//
//inline void InpFile::SetSolver (FEM::RKSolver & Sol) const
//{
    //if (ctetg       ) Sol.LinProb = ctetg;
    //if (hm          ) Sol.DampTy  = FEM::RKSolver::HMCoup_t;
    //if (ray         ) Sol.DampTy  = FEM::RKSolver::Rayleigh_t;
    //if (am    >=0   ) Sol.DampAm  = am;
    //if (ak    >=0   ) Sol.DampAk  = ak;
    //if (rkscheme!="") Sol.Scheme  = rkscheme;
    //if (rkstol   >=0) Sol.STOL    = rkstol;
    //Sol.DynCteM = rkdyncte;
//}

inline void InpFile::SetSUp (Model const * Mdl, Model::StressUpdate::pDbgFun pFun, void * UserData) const
{
    Mdl->SUp.CDrift = cdrift;
    Mdl->SUp.DbgFun = pFun;
    Mdl->SUp.DbgDat = UserData;
    if (suscheme  !="") Mdl->SUp.SetScheme (suscheme);
    if (surkscheme!="") Mdl->SUp.RKScheme = surkscheme;
    if (sustol      >0) Mdl->SUp.STOL     = sustol;
}

inline void InpFile::SetSolFlags (SDPair & Flags) const
{
}

std::ostream & operator<< (std::ostream & os, InpFile const & IF)
{
    if (IF.ninc        >=0) os << "ninc      = " << IF.ninc      << "\n"; //  1
    if (IF.cdrift         ) os << "cdrift    = " << IF.cdrift    << "\n"; //  2
    if (IF.stol        >=0) os << "stol      = " << IF.stol      << "\n"; //  3
    if (IF.ssout          ) os << "ssout     = " << IF.ssout     << "\n"; //  4
    if (IF.ctetg          ) os << "ctetg     = " << IF.ctetg     << "\n"; //  5
    if (IF.fem            ) os << "fem       = " << IF.fem       << "\n"; //  6
    if (IF.dyn            ) os << "dyn       = " << IF.dyn       << "\n"; //  7
    if (IF.hm             ) os << "hm        = " << IF.hm        << "\n"; //  8
    if (IF.tf          >=0) os << "tf        = " << IF.tf        << "\n"; //  9
    if (IF.dt          >=0) os << "dt        = " << IF.dt        << "\n"; // 10
    if (IF.dtout       >=0) os << "dtout     = " << IF.dtout     << "\n"; // 11
    if (IF.tsw         >=0) os << "tsw       = " << IF.tsw       << "\n"; // 12
    if (IF.ndiv        >=0) os << "ndiv      = " << IF.ndiv      << "\n"; // 13
    if (IF.nip         >=0) os << "nip       = " << IF.nip       << "\n"; // 14
    if (IF.o2             ) os << "o2        = " << IF.o2        << "\n"; // 15
    if (IF.ray            ) os << "ray       = " << IF.ray       << "\n"; // 16
    if (IF.am          >=0) os << "am        = " << IF.am        << "\n"; // 17
    if (IF.ak          >=0) os << "ak        = " << IF.ak        << "\n"; // 18
    if (IF.rk             ) os << "rk        = " << IF.rk        << "\n"; // 19
    if (IF.rkscheme   !="") os << "rkscheme  = " << IF.rkscheme  << "\n"; // 20
    if (IF.rkstol      >=0) os << "rkstol    = " << IF.rkstol    << "\n"; // 21
    if (IF.refdat     !="") os << "refdat    = " << IF.refdat    << "\n"; // 22
    if (IF.refsim     !="") os << "refsim    = " << IF.refsim    << "\n"; // 23
    if (IF.refana     !="") os << "refana    = " << IF.refana    << "\n"; // 24
    if (IF.idxvert1    >=0) os << "idxvert1  = " << IF.idxvert1  << "\n"; // 25
    if (IF.idxvert2    >=0) os << "idxvert2  = " << IF.idxvert2  << "\n"; // 26
    if (IF.idxvert3    >=0) os << "idxvert3  = " << IF.idxvert3  << "\n"; // 27
    if (IF.hasoptdbl1     ) os << "optdbl1   = " << IF.optdbl1   << "\n"; // 28
    if (IF.hasoptdbl2     ) os << "optdbl2   = " << IF.optdbl2   << "\n"; // 29
    if (IF.hasoptdbl3     ) os << "optdbl3   = " << IF.optdbl3   << "\n"; // 30
    if (IF.nldt_nsml   >0 ) os << "nldt_nsml = " << IF.nldt_nsml << "\n"; // 31
    if (IF.nldt_nn     >0 ) os << "nldt_nn   = " << IF.nldt_nn   << "\n"; // 32
    if (IF.nldt_n      >0 ) os << "nldt_n    = " << IF.nldt_n    << "\n"; // 33
    if (IF.nldt_ll     >0 ) os << "nldt_ll   = " << IF.nldt_ll   << "\n"; // 34
    if (IF.nldt_sch    >0 ) os << "nldt_sch  = " << IF.nldt_sch  << "\n"; // 35
    if (IF.nldt_m      >0 ) os << "nldt_m    = " << IF.nldt_m    << "\n"; // 36
    if (IF.maxit       >0 ) os << "maxit     = " << IF.maxit     << "\n"; // 37
    if (IF.tolr        >=0) os << "tolr      = " << IF.tolr      << "\n"; // 38
    if (IF.fnkey.size()>0 ) os << "fnkey     = " << IF.fnkey     << "\n"; // 39
    if (IF.haspcam0       ) os << "pcam0     = " << IF.pcam0     << "\n"; // 40
    if (IF.scheme     !="") os << "scheme    = " << IF.scheme    << "\n"; // 41
    if (IF.vtufile        ) os << "vtufile   = " << IF.vtufile   << "\n"; // 42
    if (IF.suscheme   !="") os << "suscheme  = " << IF.suscheme  << "\n"; // 43
    if (IF.sustol      >=0) os << "sustol    = " << IF.sustol    << "\n"; // 44
    if (IF.surkscheme !="") os << "surkscheme= " << IF.surkscheme<< "\n"; // 45
    if (IF.dcmaxit      >0) os << "dcmaxit   = " << IF.dcmaxit   << "\n"; // 46
    if (IF.dcftol       >0) os << "dcftol    = " << IF.dcftol    << "\n"; // 47
    if (IF.haspw0         ) os << "pw0       = " << IF.pw0       << "\n"; // 48
    if (IF.rkdyncte       ) os << "rkdyncte  = " << IF.rkdyncte  << "\n"; // 49
    if (IF.uwp            ) os << "uwp       = " << IF.uwp       << "\n"; // 50
    os << "\nPath:\n"                        << (*IF.Path)      << "\n";
    os << "\nElements properties:\n"         << (*IF.Prps)      << "\n";
    os << "\nOutput nodes:\n"                << (*IF.OutNods)   << "\n";
    os << "\nTags => Material IDs:\n"        << (*IF.Tag2MatID) << "\n";
    os << "\nTags => E(x)tra Material IDs:\n"<< (*IF.Tag2XMatID)<< "\n";
    os << "\nParameters:\n"                  << (*IF.Prms)      << "\n";
    os << "\nInitial values:\n"              << (*IF.Inis)      << "\n";
    os << "\nBoundary conditions (stages):\n";
    for (size_t i=0; i<IF.Stages.Size(); ++i)
    {
        os << "Stage # " << i << ":\n";
        os << (*IF.Stages[i]) << "\n";
    }
    return os;
}

#ifdef USE_WXWIDGETS

enum
{
    ID_INPFILE_LOAD = wxID_HIGHEST+1000,
    ID_INPFILE_SAVE,
};

BEGIN_EVENT_TABLE(InpFile, wxWindow)
    EVT_BUTTON (ID_INPFILE_LOAD, InpFile::OnLoad)
    EVT_BUTTON (ID_INPFILE_SAVE, InpFile::OnSave)
END_EVENT_TABLE()

inline InpFile::InpFile (wxFrame * Parent)
    : wxWindow (Parent, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxFULL_REPAINT_ON_RESIZE)
{
    // default values
    Defaults();

    // force validation of child controls
    SetExtraStyle (wxWS_EX_VALIDATE_RECURSIVELY);

    // tell wxAuiManager to manage this window
    Aui.SetManagedWindow (this);

    // control panel
    ADD_WXPANEL  (pnl, szt, szr, 1, 2);
    ADD_WXBUTTON (pnl, szr, ID_INPFILE_LOAD, c0, "Load");
    ADD_WXBUTTON (pnl, szr, ID_INPFILE_SAVE, c1, "Save");

    // main
    ADD_WXPANEL     (p_mai, sz_mai_t, sz_mai, 7, 2);
    ADD_WXNUMINPUT2 (p_mai, sz_mai, wxID_ANY, c_ninc      , "ninc     ", ninc     ); //   1
    ADD_WXCHECKBOX2 (p_mai, sz_mai, wxID_ANY, c_fem       , "fem      ", fem      ); //   6
    ADD_WXCHECKBOX2 (p_mai, sz_mai, wxID_ANY, c_dyn       , "dyn      ", dyn      ); //   7
    ADD_WXCHECKBOX2 (p_mai, sz_mai, wxID_ANY, c_hm        , "hm       ", hm       ); //   8
    ADD_WXNUMINPUT2 (p_mai, sz_mai, wxID_ANY, c_pcam0     , "pcam0    ", pcam0    ); //  40
                                                                                            
    // local integration                                                                    
    ADD_WXPANEL     (p_loc, sz_loc_t, sz_loc, 2, 2);                                        
    ADD_WXCHECKBOX2 (p_loc, sz_loc, wxID_ANY, c_cdrift    , "cdrift   ", cdrift   ); //   2
    ADD_WXNUMINPUT2 (p_loc, sz_loc, wxID_ANY, c_stol      , "stol     ", stol     ); //   3
                                                                                            
    // fem solution                                                                         
    ADD_WXPANEL     (p_fem, sz_fem_t, sz_fem, 18, 2);                             
    ADD_WXCHECKBOX2 (p_fem, sz_fem, wxID_ANY, c_ssout     , "ssout    ", ssout    ); //   4
    ADD_WXCHECKBOX2 (p_fem, sz_fem, wxID_ANY, c_ctetg     , "ctetg    ", ctetg    ); //   5
    ADD_WXNUMINPUT2 (p_fem, sz_fem, wxID_ANY, c_tf        , "tf       ", tf       ); //   9
    ADD_WXNUMINPUT2 (p_fem, sz_fem, wxID_ANY, c_dt        , "dt       ", dt       ); //  10
    ADD_WXNUMINPUT2 (p_fem, sz_fem, wxID_ANY, c_dtout     , "dtout    ", dtout    ); //  11
    ADD_WXNUMINPUT2 (p_fem, sz_fem, wxID_ANY, c_tsw       , "tsw      ", tsw      ); //  12
    ADD_WXNUMINPUT2 (p_fem, sz_fem, wxID_ANY, c_ndiv      , "ndiv     ", ndiv     ); //  13
    ADD_WXNUMINPUT2 (p_fem, sz_fem, wxID_ANY, c_nip       , "nip      ", nip      ); //  14
    ADD_WXCHECKBOX2 (p_fem, sz_fem, wxID_ANY, c_o2        , "o2       ", o2       ); //  15
    ADD_WXCHECKBOX2 (p_fem, sz_fem, wxID_ANY, c_ray       , "ray      ", ray      ); //  16
    ADD_WXNUMINPUT2 (p_fem, sz_fem, wxID_ANY, c_am        , "am       ", am       ); //  17
    ADD_WXNUMINPUT2 (p_fem, sz_fem, wxID_ANY, c_ak        , "ak       ", ak       ); //  18
    ADD_WXCHECKBOX2 (p_fem, sz_fem, wxID_ANY, c_rk        , "rk       ", rk       ); //  19
    ADD_WXTEXTCTRL2 (p_fem, sz_fem, wxID_ANY, c_rkscheme  , "rkscheme ", rkscheme ); //  20
    ADD_WXNUMINPUT2 (p_fem, sz_fem, wxID_ANY, c_rkstol    , "rkstol   ", rkstol   ); //  21
    ADD_WXNUMINPUT2 (p_fem, sz_fem, wxID_ANY, c_maxit     , "maxit    ", maxit    ); //  37
    ADD_WXNUMINPUT2 (p_fem, sz_fem, wxID_ANY, c_tolr      , "tolr     ", tolr     ); //  38
    ADD_WXTEXTCTRL2 (p_fem, sz_fem, wxID_ANY, c_scheme    , "scheme   ", scheme   ); //  41
    ADD_WXCHECKBOX2 (p_fem, sz_fem, wxID_ANY, c_vtufile   , "vtufile  ", vtufile  ); //  42
                                                                                            
    // reference files                                                                      
    ADD_WXPANEL     (p_rfi, sz_rfi_t, sz_rfi, 3, 2);                                        
    ADD_WXTEXTCTRL2 (p_rfi, sz_rfi, wxID_ANY, c_refdat    , "refdat   ", refdat   ); //  22
    ADD_WXTEXTCTRL2 (p_rfi, sz_rfi, wxID_ANY, c_refsim    , "refsim   ", refsim   ); //  23
    ADD_WXTEXTCTRL2 (p_rfi, sz_rfi, wxID_ANY, c_refana    , "refana   ", refana   ); //  24
    c_refdat->SetMinSize (wxSize(200,20));
    c_refsim->SetMinSize (wxSize(200,20));
    c_refana->SetMinSize (wxSize(200,20));
                                                                                            
    // nonlinear steps                                                                      
    ADD_WXPANEL     (p_nls, sz_nls_t, sz_nls, 6, 2);                                        
    ADD_WXNUMINPUT2 (p_nls, sz_nls, wxID_ANY, c_nldt_nsml , "nldt_nsml", nldt_nsml); //  31
    ADD_WXNUMINPUT2 (p_nls, sz_nls, wxID_ANY, c_nldt_nn   , "nldt_nn  ", nldt_nn  ); //  32
    ADD_WXNUMINPUT2 (p_nls, sz_nls, wxID_ANY, c_nldt_n    , "nldt_n   ", nldt_n   ); //  33
    ADD_WXNUMINPUT2 (p_nls, sz_nls, wxID_ANY, c_nldt_ll   , "nldt_ll  ", nldt_ll  ); //  34
    ADD_WXNUMINPUT2 (p_nls, sz_nls, wxID_ANY, c_nldt_sch  , "nldt_sch ", nldt_sch ); //  35
    ADD_WXNUMINPUT2 (p_nls, sz_nls, wxID_ANY, c_nldt_m    , "nldt_m   ", nldt_m   ); //  36
                                                                                            
    // others                                                                               
    ADD_WXPANEL     (p_oth, sz_oth_t, sz_oth, 6, 2);                                        
    ADD_WXNUMINPUT2 (p_oth, sz_oth, wxID_ANY, c_idxvert1  , "idxvert1 ", idxvert1 ); //  25
    ADD_WXNUMINPUT2 (p_oth, sz_oth, wxID_ANY, c_idxvert2  , "idxvert2 ", idxvert2 ); //  26
    ADD_WXNUMINPUT2 (p_oth, sz_oth, wxID_ANY, c_idxvert3  , "idxvert3 ", idxvert3 ); //  27
    ADD_WXNUMINPUT2 (p_oth, sz_oth, wxID_ANY, c_optdbl1   , "optdbl1  ", optdbl1  ); //  28
    ADD_WXNUMINPUT2 (p_oth, sz_oth, wxID_ANY, c_optdbl2   , "optdbl2  ", optdbl2  ); //  29
    ADD_WXNUMINPUT2 (p_oth, sz_oth, wxID_ANY, c_optdbl3   , "optdbl3  ", optdbl3  ); //  30

    // path
    Path  = new GUI::WxDictTable;         Path->Transposed = false;
    GPath = new GUI::WxDict (this, Path); GPath->FitCol    = true;

    // additional data
    Prps       = new GUI::WxDictTable;
    OutNods    = new GUI::WxArrayIntTable;
    Tag2MatID  = new GUI::WxSIPairTable;
    GPrps      = new GUI::WxDict     (this, Prps);
    GOutNods   = new GUI::WxArrayInt (this, OutNods);
    GMatId2Tag = new GUI::WxSIPair   (this, Tag2MatID);
    Prms       = new Dict;
    Inis       = new Dict;

    // notebook
    ADD_WXNOTEBOOK (this, nbk0);
    ADD_WXNOTEBOOK (this, nbk1);
    ADD_WXNOTEBOOK (this, nbk2);
    nbk0->AddPage  (p_mai,      "Main",                 false);
    nbk0->AddPage  (p_loc,      "Local Integration",    false);
    nbk0->AddPage  (p_oth,      "Others",               false);
    nbk2->AddPage  (p_fem,      "FEM Solution",         false);
    nbk2->AddPage  (p_nls,      "Nonlinear Steps",      false);
    nbk0->AddPage  (p_rfi,      "Reference Files",      false);
    nbk1->AddPage  (GPath,      "Path",                 false);
    nbk2->AddPage  (GPrps,      "Elements Properties",  false);
    nbk2->AddPage  (GOutNods,   "Output Nodes",         false);
    nbk2->AddPage  (GMatId2Tag, "Material IDs => Tags", false);

    // commit all changes to wxAuiManager
    Aui.AddPane (pnl,  wxAuiPaneInfo().Name("cpnl").Caption("cpnl").Top().MinSize(wxSize(100,50)).DestroyOnClose(false).CaptionVisible(false) .CloseButton(false));
    Aui.AddPane (nbk0, wxAuiPaneInfo().Name("nbk0").Caption("General Input Data").Centre().Position(0).DestroyOnClose(false).CaptionVisible(true).CloseButton(false));
    Aui.AddPane (nbk1, wxAuiPaneInfo().Name("nbk1").Caption("Stress/Strain Path").Centre().Position(1).DestroyOnClose(false).CaptionVisible(true).CloseButton(false));
    Aui.AddPane (nbk2, wxAuiPaneInfo().Name("nbk2").Caption("FEM Input Data")    .Centre().Position(2).DestroyOnClose(false).CaptionVisible(true).CloseButton(false));
    Aui.Update  ();
}

inline void InpFile::OnLoad (wxCommandEvent & Event)
{
    wxFileDialog fd(this, "Load input (.inp) file", LstDir, "", "*.inp");
    if (fd.ShowModal()==wxID_OK)
    {
        LstDir = fd.GetDirectory ();
        try { Read (fd.GetPath().ToStdString().c_str()); }
        catch (Fatal * e) { WxError(e->Msg().CStr()); }
    }
}

inline void InpFile::OnSave (wxCommandEvent & Event)
{
    Sync ();
    wxFileDialog fd(this, "Save input (.inp) file", LstDir, "", "*.inp", wxFD_SAVE|wxFD_OVERWRITE_PROMPT);
    if (fd.ShowModal()==wxID_OK)
    {
        std::cout << "fem = " << fem << std::endl;
        std::fstream of(fd.GetPath().ToStdString().c_str(), std::ios::out);
        of << "ninc       = " << ninc       << std::endl; //  1
        of << "cdrift     = " << cdrift     << std::endl; //  2
        of << "stol       = " << stol       << std::endl; //  3
        of << "ssout      = " << ssout      << std::endl; //  4
        of << "ctetg      = " << ctetg      << std::endl; //  5
        of << "fem        = " << fem        << std::endl; //  6
        of << "dyn        = " << dyn        << std::endl; //  7
        of << "hm         = " << hm         << std::endl; //  8
        of << "tf         = " << tf         << std::endl; //  9
        of << "dt         = " << dt         << std::endl; // 10
        of << "dtout      = " << dtout      << std::endl; // 11
        of << "tsw        = " << tsw        << std::endl; // 12
        of << "ndiv       = " << ndiv       << std::endl; // 13
        of << "nip        = " << nip        << std::endl; // 14
        of << "o2         = " << o2         << std::endl; // 15
        of << "ray        = " << ray        << std::endl; // 16
        of << "am         = " << am         << std::endl; // 17
        of << "ak         = " << ak         << std::endl; // 18
        of << "rk         = " << rk         << std::endl; // 19
        of << "rkscheme   = " << rkscheme   << std::endl; // 20
        of << "rkstol     = " << rkstol     << std::endl; // 21
        of << "refdat     = " << refdat     << std::endl; // 22
        of << "refsim     = " << refsim     << std::endl; // 23
        of << "refana     = " << refana     << std::endl; // 24
        of << "idxvert1   = " << idxvert1   << std::endl; // 25
        of << "idxvert2   = " << idxvert2   << std::endl; // 26
        of << "idxvert3   = " << idxvert3   << std::endl; // 27
        of << "optdbl1    = " << optdbl1    << std::endl; // 28
        of << "optdbl2    = " << optdbl2    << std::endl; // 29
        of << "optdbl3    = " << optdbl3    << std::endl; // 30
        of << "nldt_nsml  = " << nldt_nsml  << std::endl; // 31
        of << "nldt_nn    = " << nldt_nn    << std::endl; // 32
        of << "nldt_n     = " << nldt_n     << std::endl; // 33
        of << "nldt_ll    = " << nldt_ll    << std::endl; // 34
        of << "nldt_sch   = " << nldt_sch   << std::endl; // 35
        of << "nldt_m     = " << nldt_m     << std::endl; // 36
        of << "maxit      = " << maxit      << std::endl; // 37
        of << "tolr       = " << tolr       << std::endl; // 39
        of << "pcam0      = " << pcam0      << std::endl; // 40
        of << "scheme     = " << scheme     << std::endl; // 41
        of << "vtufile    = " << vtufile    << std::endl; // 42
        of << "suscheme   = " << suscheme   << std::endl; // 43
        of << "sustol     = " << sustol     << std::endl; // 44
        of << "surkscheme = " << surkscheme << std::endl; // 45
        of << "dcmaxit    = " << dcmaxit    << std::endl; // 46
        of << "dcftol     = " << dcftol     << std::endl; // 47
        of << "pw0        = " << pw0        << std::endl; // 48
        of << "rkdyncte   = " << rkdyncte   << std::endl; // 49
        of << "uwp        = " << uwp        << std::endl; // 50
        of << "npath      = " << Path->Keys.Size() << std::endl;

        String buf;
        for (size_t i=0; i<Path->Keys.Size(); ++i)
        {
            SDPair const & path = (*Path)(Path->Keys[i]);
            of << std::endl;
            of << "ndat = " << path.Keys.Size() << std::endl;
            for (size_t j=0; j<path.Keys.Size(); ++j)
            {
                buf.Printf("  %-5s = %g\n", path.Keys[j].CStr(), path(path.Keys[j]));
                of << buf;
            }
        }
        of.close();
    }
}

#endif


/////////////////////////////////////////////////////////////////////////////////////////// Macros ///////


#define INIT_MAT_INP(argc, argv, inpfn, matfn, verbose, forcegty,   MAT, INP) \
    if (argc>1) inpfn   =      argv[1];                                       \
    if (argc>2) matfn   =      argv[2];                                       \
    if (argc>3) verbose = atoi(argv[3]);                                      \
    MatFile MAT;                                                              \
    InpFile INP;                                                              \
    MAT.Read        (matfn.CStr());                                           \
    INP.Read        (inpfn.CStr());                                           \
    INP.SetPrmsInis (MAT, forcegty);                                          \
    if (verbose)                                                              \
    {                                                                         \
        printf("\n%s--- Materials data <%s> ----------------------------------------------------%s\n",TERM_CLR1,matfn.CStr(),TERM_RST); \
        std::cout << MAT << std::endl;                                        \
        printf("\n%s--- Input data <%s> --------------------------------------------------------%s\n",TERM_CLR1,inpfn.CStr(),TERM_RST); \
        std::cout << INP << std::endl;                                        \
    }

#define INIT_MAT_INP_(argc, argv,   MAT, INP) \
    String inpfn("input.inp");                \
    String matfn("materials.mat");            \
    bool   verbose  = true;                   \
    bool   forcegty = false;                  \
    INIT_MAT_INP(argc, argv, inpfn, matfn, verbose, forcegty,   MAT, INP);


#endif // MECHSYS_INPFILE_H
