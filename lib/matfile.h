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

#ifndef MECHSYS_MATFILE_H
#define MECHSYS_MATFILE_H

// Std Lib
#include <sstream>
#include <fstream>

// wxWidgets
#ifdef USE_WXWIDGETS
  #include <mechsys/gui/wxdict.h>
  #include <mechsys/gui/common.h>
#endif

// MechSys
#include <mechsys/util/maps.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>
#define INCLUDE_MODELS_ONLY
#include <mechsys/fem/fem.h>
#undef  INCLUDE_MODELS_ONLY

#ifdef USE_WXWIDGETS
class MatFile : public wxWindow
#else
class MatFile
#endif
{
public:
    // Constructor
#ifdef USE_WXWIDGETS
     MatFile (wxFrame * Parent);
    ~MatFile () { Aui.UnInit(); }
#else
     MatFile () { ID2Prms = new Dict;  ID2Inis = new Dict; }
    ~MatFile () { delete ID2Prms;  delete ID2Inis; }
#endif

    // Methods
    void Read (char const * FileName);
    void Save (char const * FileName);

#ifdef USE_WXWIDGETS
    // Methods
    void AddMdl ();

    // Data
    GUI::WxDictTable * ID2Prms;  ///< maps ID to parameters
    GUI::WxDictTable * ID2Inis;  ///< maps ID to initial values
    GUI::WxDict      * DPrms;    ///< grid view for ID2Prms
    GUI::WxDict      * DInis;    ///< grid view for ID2Inis
    wxAuiManager       Aui;      ///< Aui manager
    wxString           LstDir;   ///< Last accessed directory
    wxComboBox       * CbxMdl;   ///< Model names
    String             FName;    ///< Material file (.mat) filename
    wxArrayString      MNames;   ///< Model names

    // Events
    void OnLoad (wxCommandEvent & Event);
    void OnSave (wxCommandEvent & Event);
    void OnAdd  (wxCommandEvent & Event) { AddMdl (); }
    void OnDel  (wxCommandEvent & Event);
    DECLARE_EVENT_TABLE();
#else
    // Data
    Dict * ID2Prms; ///< maps ID to parameters
    Dict * ID2Inis; ///< maps ID to initial values
#endif

#ifdef USE_BOOST_PYTHON
    void PyGetPrmIni (int ID, BPy::dict & Prm, BPy::dict & Ini)
    {
        ID2Prms->PyGet (ID, Prm);
        ID2Inis->PyGet (ID, Ini);
    }
#endif

};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline void MatFile::Read (char const * FileName)
{
    // open material file
    std::fstream mat_file(FileName, std::ios::in);
    if (!mat_file.is_open()) throw new Fatal("MatFile::Read: Could not open file <%s>",FileName);

    // parse
    int  idxcte       = 0;
    int  idxprm       = 0;
    int  idxini       = 0;
    int  nctes        = 0;
    int  nprms        = 0;
    int  ninis        = 0;
    int  line_num     = 1;
    bool reading_name = false;
    bool reading_ctes = false;
    bool reading_prms = false;
    bool reading_inis = false;
    int  ID           = 0;
    String model_name;
    Array<String> const * prm_names = NULL;
    Array<String> const * ivs_names = NULL;
    ID2Prms->clear ();
    ID2Inis->clear ();
    while (!mat_file.eof())
    {
        String line,key,equal,strval;
        std::getline (mat_file,line);
        std::istringstream iss(line);
        if (iss >> key >> equal >> strval)
        {
            if (key[0]=='#') { line_num++; continue; }
            else if (key=="ID")
            {
                ID = atoi(strval.CStr());
                if (ID<0) throw new Fatal("MatFile::Read: Error in <%s> file at line # %d: IDs must be greater than zero or equal to zero. %d is invalid",FileName,line_num,ID);
                reading_name = true;
            }
            else if (reading_name)
            {
                if (key=="name")
                {
                    // iterators to prm and ivs names
                    Str2ArrayStr_t::const_iterator ita = MODEL_PRM_NAMES.find(strval);
                    Str2ArrayStr_t::const_iterator itb = MODEL_IVS_NAMES.find(strval);

                    // check
                    if (ID2Prms->HasKey(ID))             throw new Fatal("MatFile::Read: Error in <%s> file at line # %d: IDs must be unique. %d is repeated",FileName,line_num,ID);
                    if (MODEL.find(strval)==MODEL.end()) throw new Fatal("MatFile::Read: Error in <%s> file at line # %d: Model 'name' = %s is not available in MODEL",FileName,line_num,strval.CStr());
                    if (ita==MODEL_PRM_NAMES.end())      throw new Fatal("MatFile::Read: Error in <%s> file at line # %d: Model 'name' = %s is not available in MODEL_PRM_NAMES",FileName,line_num,strval.CStr());
                    if (itb==MODEL_IVS_NAMES.end())      throw new Fatal("MatFile::Read: Error in <%s> file at line # %d: Model 'name' = %s is not available in MODEL_IVS_NAMES",FileName,line_num,strval.CStr());

                    // set
                    ID2Prms->Set (ID, "name", MODEL(strval));
                    model_name   = strval;
                    reading_name = false;
                    reading_ctes = true;
                    nctes        = 0;
                    idxcte       = -1;
                    prm_names    = &ita->second;
                    ivs_names    = &itb->second;
                }
                else throw new Fatal("MatFile::Read: Error in <%s> file at line # %d: 'name' must follow 'ID'. '%s' is invalid or in the wrong place",FileName,line_num,key.CStr());
            }
            else if (reading_ctes)
            {
                if (nctes==0)
                {
                    if (key=="nctes") nctes = atoi(strval.CStr());
                    else throw new Fatal("MatFile::Read: Error in <%s> file at line # %d: 'nctes' must follow 'name'. '%s' is invalid or in the wrong place",FileName,line_num,key.CStr());
                    if (nctes==0)
                    {
                        reading_ctes = false;
                        reading_prms = true;
                        nprms        = 0;
                        idxprm       = -1;
                    }
                }
                else if (key=="nprms")
                {
                    if (idxcte==nctes-1)
                    {
                        reading_ctes = false;
                        reading_prms = true;
                        nprms        = atoi(strval.CStr());
                        idxprm       = -1;
                    }
                    else throw new Fatal("MatFile::Read: Error in <%s> file at line # %d: 'nprms' must appear after all 'nctes' constants are read. nctes=%d and %d were read so far",FileName,line_num,nctes,idxcte+1);
                }
                else if (idxcte<nctes)
                {
                    if (MODEL_CTE_NAMES.Has(key))
                    {
                        ID2Prms->Set (ID, key.CStr(), atof(strval.CStr()));
                        idxcte++;
                    }
                    else throw new Fatal("MatFile::Read: Error in <%s> file at line # %d: constant named '%s' is not available for model '%s'",FileName,line_num,key.CStr(),model_name.CStr());
                }
                else throw new Fatal("MatFile::Read: Error in <%s> file at line # %d: there are more constants than what specified by nctes=%d. The reading of constants finishes when 'nprms' is found. '%s' is invalid or in the wrong place (idxcte=%d)",FileName,line_num,nctes,key.CStr(),idxcte);
            }
            else if (reading_prms)
            {
                if (nprms==0)
                {
                    if (key=="nprms") nprms = atoi(strval.CStr());
                    else throw new Fatal("MatFile::Read: Error in <%s> file at line # %d: 'nprms' must follow constants. '%s' is invalid or in the wrong place",FileName,line_num,key.CStr());
                }
                else if (key=="ninis")
                {
                    if (idxprm==nprms-1)
                    {
                        ninis        = atoi(strval.CStr());
                        idxini       = -1;
                        reading_prms = false;
                        if (ninis>0) reading_inis = true;
                        else
                        {
                            ID2Inis->Set (ID, SDPair());
                            reading_inis = false;
                        }
                    }
                    else throw new Fatal("MatFile::Read: Error in <%s> file at line # %d: 'ninis' must appear after all 'nprms' parameters are read. nprms=%d and %d were read so far",FileName,line_num,nprms,idxprm+1);
                }
                else if (idxprm<nprms)
                {
                    if (prm_names->Has(key))
                    {
                        ID2Prms->Set (ID, key.CStr(), atof(strval.CStr()));
                        idxprm++;
                    }
                    else throw new Fatal("MatFile::Read: Error in <%s> file at line # %d: parameter named '%s' is not available for model '%s'",FileName,line_num,key.CStr(),model_name.CStr());
                }
                else throw new Fatal("MatFile::Read: Error in <%s> file at line # %d: there are more parameters than what specified by nprms=%d. The reading of parameters finishes when 'ninis' is found. '%s' is invalid or in the wrong place (idxprm=%d)",FileName,line_num,nprms,key.CStr(),idxprm);
            }
            else if (reading_inis)
            {
                if (key=="ID" || key=="name" || key=="nctes" || key=="nprms" || key=="ninis") throw new Fatal("MatFile: Error in <%s> file at line # %d: There are not enough initial values corresponding to ninis=%d. '%s' is the wrong place",FileName,line_num,ninis,key.CStr());
                else if (idxini==ninis-1)
                {
                    reading_name = false;
                    reading_ctes = false;
                    reading_prms = false;
                    reading_inis = false;
                }
                else if (idxini<ninis)
                {
                    if (ivs_names->Has(key))
                    {
                        ID2Inis->Set (ID, key.CStr(), atof(strval.CStr()));
                        idxini++;
                    }
                    else throw new Fatal("MatFile::Read: Error in <%s> file at line # %d: initial value '%s' is not available for model '%s'",FileName,line_num,key.CStr(),model_name.CStr());
                }
                else throw new Fatal("MatFile::Read: Error in <%s> file at line # %d: there are more initial values than what specified by ninis=%d. The reading of initial values finishes when %d values are read. '%s' is invalid or in the wrong place (idxini)",FileName,line_num,ninis,ninis,key.CStr(),idxini);
            }
            else throw new Fatal("MatFile::Read: Problem with <%s> file at line # %d. '%s' is invalid or in the wrong place",FileName,line_num,key.CStr());
        }
        line_num++;
    }

#ifdef USE_WXWIDGETS
    DPrms->ReBuild  ();
    DInis->ReBuild  ();
    TransferDataToWindow ();
#endif
}

inline void MatFile::Save (char const * FileName)
{
    std::ostringstream oss;
    oss << "############ Materials ##############\n\n";
    String buf;
    for (size_t i=0; i<ID2Prms->Keys.Size(); ++i)
    {
        int            id   = ID2Prms->Keys[i];
        SDPair const & prms = (*ID2Prms)(id);
        SDPair const & inis = (*ID2Inis)(id);
        String name;
        MODEL.Val2Key (prms("name"), name);
        buf.Printf("%-8s = %d\n", "ID",    id);                  oss<<buf;
        buf.Printf("%-8s = %s\n", "name",  name.CStr());         oss<<buf;
        buf.Printf("%-8s = %d\n", "nprms", prms.Keys.Size()-1);  oss<<buf;
        for (size_t j=0; j<prms.Keys.Size(); ++j)
        {
            if (prms.Keys[j]!="name")
            {
                buf.Printf("%-8s = %g\n", prms.Keys[j].CStr(), prms(prms.Keys[j]));
                oss << buf;
            }
        }
        buf.Printf("%-8s = %d\n", "ninis", inis.Keys.Size());
        oss << buf;
        for (size_t j=0; j<inis.Keys.Size(); ++j)
        {
            buf.Printf("%-8s = %g\n", inis.Keys[j].CStr(), inis(inis.Keys[j]));
            oss << buf;
        }
        oss << std::endl;
    }
    std::fstream of(FileName, std::ios::out);
    of << oss.str();
    of.close();
}

std::ostream & operator<< (std::ostream & os, MatFile const & MF)
{
    os << "ID2Prms =\n";
    os << (*MF.ID2Prms) << std::endl;
    os << "ID2Inis =\n";
    os << (*MF.ID2Inis) << std::endl;
    return os;
}

#ifdef USE_WXWIDGETS

enum
{
    ID_MATFILE_LOAD = wxID_HIGHEST+2000,
    ID_MATFILE_SAVE ,
    ID_MATFILE_ADD  ,
    ID_MATFILE_DEL  ,
};

BEGIN_EVENT_TABLE(MatFile, wxWindow)
    EVT_BUTTON (ID_MATFILE_LOAD, MatFile::OnLoad)
    EVT_BUTTON (ID_MATFILE_SAVE, MatFile::OnSave)
    EVT_BUTTON (ID_MATFILE_ADD,  MatFile::OnAdd)
    EVT_BUTTON (ID_MATFILE_DEL,  MatFile::OnDel)
END_EVENT_TABLE()

inline MatFile::MatFile (wxFrame * Parent)
    : wxWindow (Parent, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxFULL_REPAINT_ON_RESIZE)
{
    // force validation of child controls
    SetExtraStyle (wxWS_EX_VALIDATE_RECURSIVELY);

    // tell wxAuiManager to manage this window
    Aui.SetManagedWindow (this);

    // control panel
    ADD_WXPANEL    (pnl, szt, szr, 1, 5);
    ADD_WXBUTTON   (pnl, szr, ID_MATFILE_LOAD, c0,     "Load");
    ADD_WXBUTTON   (pnl, szr, ID_MATFILE_SAVE, c1,     "Save");
    ADD_WXCOMBOBOX (pnl, szr, wxID_ANY,        CbxMdl, "Select Model");
    ADD_WXBUTTON   (pnl, szr, ID_MATFILE_ADD,  c2,     "Add" );
    ADD_WXBUTTON   (pnl, szr, ID_MATFILE_DEL,  c3,     "Del" );
    CbxMdl->SetMinSize (wxSize(200,12));

    // grids and models
    ID2Prms = new GUI::WxDictTable;
    ID2Inis = new GUI::WxDictTable;
    DPrms   = new GUI::WxDict (this, ID2Prms);
    DInis   = new GUI::WxDict (this, ID2Inis);
    for (ModelFactory_t::iterator it=ModelFactory.begin(); it!=ModelFactory.end(); ++it)
    {
        String const & model_name = it->first;
        Str2ArrayStr_t::const_iterator iprm = MODEL_PRM_NAMES.find(model_name);
        Str2ArrayStr_t::const_iterator iivs = MODEL_IVS_NAMES.find(model_name);
        if (iprm==MODEL_PRM_NAMES.end()) WxError("MatFile::MatFile: __internal_error__ Model named <%s> is not in map: MODEL_PRM_NAMES",model_name.CStr());
        if (iivs==MODEL_IVS_NAMES.end()) WxError("MatFile::MatFile: __internal_error__ Model named <%s> is not in map: MODEL_IVS_NAMES",model_name.CStr());
        MNames.Add (model_name);
    }
    ID2Prms->Val2Name = &MODEL;
    CbxMdl->Set (MNames);
    if (MNames.size()>0) CbxMdl->SetValue(MNames[0]);

    // commit all changes to wxAuiManager
    Aui.AddPane (pnl,   wxAuiPaneInfo().Name("cpnl") .Caption("Material file: Control Panel")   .Top   ().            DestroyOnClose(false).CaptionVisible(true).CloseButton(false).MinSize(wxSize(100,50)).Floatable(false));
    Aui.AddPane (DPrms, wxAuiPaneInfo().Name("dprms").Caption("Material file: Model parameters").Centre().Position(1).DestroyOnClose(false).CaptionVisible(true).CloseButton(false));
    Aui.AddPane (DInis, wxAuiPaneInfo().Name("dinis").Caption("Material file: Initial values")  .Centre().Position(2).DestroyOnClose(false).CaptionVisible(true).CloseButton(false));
    Aui.Update  ();
}

inline void MatFile::OnLoad (wxCommandEvent & Event)
{
    wxFileDialog fd(this, "Load material (.mat) file", LstDir, "", "*.mat");
    if (fd.ShowModal()==wxID_OK)
    {
        LstDir = fd.GetDirectory ();
        try { Read (fd.GetPath().ToStdString().c_str()); }
        catch (Fatal * e) { WxError(e->Msg().CStr()); }
    }
}

inline void MatFile::OnSave (wxCommandEvent & Event)
{
    TransferDataFromWindow(); // Synchronise (validate/transfer) data in controls
    wxFileDialog fd(this, "Save material (.mat) file", LstDir, "", "*.mat", wxFD_SAVE|wxFD_OVERWRITE_PROMPT);
    if (fd.ShowModal()==wxID_OK) Save (fd.GetPath().ToStdString().c_str());
}

inline void MatFile::AddMdl ()
{
    String model_name = CbxMdl->GetValue().ToStdString();
    int max_id = -1;
    for (size_t i=0; i<ID2Prms->Keys.Size(); ++i) if (ID2Prms->Keys[i]>max_id) max_id = ID2Prms->Keys[i];
    int id = max_id + 1;
    Str2ArrayStr_t::const_iterator iprm = MODEL_PRM_NAMES.find(model_name);
    Str2ArrayStr_t::const_iterator iivs = MODEL_IVS_NAMES.find(model_name);
    ID2Prms->Set     (id, "name", MODEL(model_name));
    ID2Prms->SetZero (id, iprm->second);
    ID2Inis->SetZero (id, iivs->second);
    DPrms->ReBuild   ();
    DInis->ReBuild   ();
}

inline void MatFile::OnDel (wxCommandEvent & Event)
{
    wxArrayInt sel = (DPrms->Tab->Transposed ? DPrms->Grd->GetSelectedRows() : DPrms->Grd->GetSelectedCols());
    if (DPrms->Tab->Transposed) DPrms->Tab->DeleteRows (sel);
    else                        DPrms->Tab->DeleteCols (sel);
    if (DInis->Tab->Transposed) DInis->Tab->DeleteRows (sel);
    else                        DInis->Tab->DeleteCols (sel);
    DPrms->ReBuild ();
    DInis->ReBuild ();
}

#endif

#endif // MECHSYS_MATFILE_H
