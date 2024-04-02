/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo                                    *
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

#ifndef MECHSYS_DATFILES_H
#define MECHSYS_DATFILES_H

// STL
#include <cmath>   // for ceil and floor
#include <cfloat>  // for DBL_EPSILON

// MechSys
#include <mechsys/linalg/matvec.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/array.h>
#include <mechsys/util/maps.h>
#ifdef USE_WXWIDGETS
  #include <mechsys/gui/plotxy.h>
  #include <mechsys/gui/common.h>
#endif

#ifdef USE_WXWIDGETS
class DatFiles : public wxWindow
#else
class DatFiles
#endif
{
public:
    // Constructor
#ifdef USE_WXWIDGETS
     DatFiles (wxFrame * Parent);
    ~DatFiles () { Aui.UnInit(); }
#else
    DatFiles () {}
#endif

    // Methods
#ifdef USE_WXWIDGETS
    void Read (wxArrayString const & FNames, wxArrayString const & ShortNames); ///< Read files
#else
    void Read (Array<String> const & FNames); ///< Read files
#endif

    // Data
    Array<Array<Vec_t> >  sig;                  ///< Stresses
    Array<Array<Vec_t> >  eps;                  ///< Strains [%]
    Array<Array<double> > ed, ev;               ///< Cambridge strain invariants [%]
    Array<Array<double> > p, q, t, lp, qp, mqp; ///< Octahedral invariants

    // Data corresponding to simulations
    Array<Array<double> > ED, EV;      ///< Cambridge strain invariants [%]
    Array<Array<double> > LP, QP, MQP; ///< Octahedral invariants

#ifdef USE_WXWIDGETS
    // Methods
    void Sync    () { TransferDataFromWindow(); } ///< Synchronise (validate/transfer) data in controls
    void ReBuild ();                              ///< Re-build plots

    // Data
    wxAuiManager   Aui;    ///< Aui Manager
    GUI::PlotXY  * qped;   ///< q/p versus Ed plot
    GUI::PlotXY  * qpev;   ///< q/p versus Ev plot
    GUI::PlotXY  * eved;   ///< Ev versus Ed plot
    GUI::PlotXY  * evlp;   ///< Ev versus log(p) plot
    bool           Multq;  ///< multiply q according to Lode angle ?
    bool           PltAll; ///< plot all data at the same time ?
    wxComboBox   * CbxFNs; ///< data filenames
    wxString       LstDir; ///< Last accessed directory

    // Events
    void OnLoad    (wxCommandEvent & Event);
    void OnReBuild (wxCommandEvent & Event) { ReBuild (); }
    DECLARE_EVENT_TABLE()
#endif
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


#ifdef USE_WXWIDGETS
inline void DatFiles::Read (wxArrayString const & FNames, wxArrayString const & ShortNames)
#else
inline void DatFiles::Read (Array<String> const & FNames)
#endif
{
    // resize data arrays
    size_t nfiles = FNames.size();
    sig.Resize (nfiles);
    eps.Resize (nfiles);
    ed .Resize (nfiles);
    ev .Resize (nfiles);
    p  .Resize (nfiles);
    q  .Resize (nfiles);
    t  .Resize (nfiles);
    lp .Resize (nfiles);
    qp .Resize (nfiles);
    mqp.Resize (nfiles);

    // load data
    for (size_t i=0; i<FNames.size(); ++i)
    {
        // read tabular data
        Table dat;
        dat.Read (FNames[i].ToStdString().c_str());
        size_t nrow    = 0;
        bool   has_sxy = false;
        bool   has_syz = false;
        bool   has_szx = false;
        bool   has_exy = false;
        bool   has_eyz = false;
        bool   has_ezx = false;
        Array<String> skeys(6);
        Array<String> ekeys(6);
        double m = 1.0; // multiplier for sig
        double n = 1.0; // multiplier for eps to convert to %
        for (size_t j=0; j<dat.Keys.Size(); ++j)
        {
            if      (dat.Keys[j]=="Sx")  { skeys="Sx","Sy","Sz","Sxy","Syz","Szx"; ekeys="Ex","Ey","Ez","Exy","Eyz","Ezx";  m=-1.0;  nrow=dat("Sx").Size(); } // old_data_file_xyz
            else if (dat.Keys[j]=="Sa")  { skeys="Sr","St","Sa","Srt","Sat","Sar"; ekeys="Er","Et","Ea","Ert","Eat","Ear";  m=-1.0;  nrow=dat("Sa").Size(); } // old_data_file_art
            else if (dat.Keys[j]=="sx")  { skeys="sx","sy","sz","sxy","syz","szx"; ekeys="ex","ey","ez","exy","eyz","ezx";  n=100.0; nrow=dat("sx").Size(); } // new_data_file_xyz
            else if (dat.Keys[j]=="Sxy" || dat.Keys[j]=="Srt" || dat.Keys[j]=="sxy") has_sxy = true;
            else if (dat.Keys[j]=="Syz" || dat.Keys[j]=="Sat" || dat.Keys[j]=="syz") has_syz = true;
            else if (dat.Keys[j]=="Szx" || dat.Keys[j]=="Sar" || dat.Keys[j]=="szx") has_szx = true;
            else if (dat.Keys[j]=="Exy" || dat.Keys[j]=="Ert" || dat.Keys[j]=="exy") has_exy = true;
            else if (dat.Keys[j]=="Eyz" || dat.Keys[j]=="Eat" || dat.Keys[j]=="eyz") has_eyz = true;
            else if (dat.Keys[j]=="Ezx" || dat.Keys[j]=="Ear" || dat.Keys[j]=="ezx") has_ezx = true;
        }
        if (nrow==0) throw new Fatal("DatFiles::OnLoad: Could not find (sx,sy,sz,sxy) columns in file %s",FNames[i].ToStdString().c_str());

        // resize sub arrays
        sig[i].Resize (nrow);
        eps[i].Resize (nrow);
        ed [i].Resize (nrow);
        ev [i].Resize (nrow);
        p  [i].Resize (nrow);
        q  [i].Resize (nrow);
        t  [i].Resize (nrow);
        lp [i].Resize (nrow);
        qp [i].Resize (nrow);
        mqp[i].Resize (nrow);

        // calculate invariants
        for (size_t j=0; j<nrow; ++j)
        {
            sig[i][j].change_dim (6);
            sig[i][j] = m*dat(skeys[0])[j], m*dat(skeys[1])[j], m*dat(skeys[2])[j],
                        (has_sxy ? m*dat(skeys[3])[j]*sqrt(2.0) : 0.0),
                        (has_syz ? m*dat(skeys[4])[j]*sqrt(2.0) : 0.0),
                        (has_szx ? m*dat(skeys[5])[j]*sqrt(2.0) : 0.0);

            eps[i][j].change_dim (6);
            eps[i][j] = m*n*dat(ekeys[0])[j], m*n*dat(ekeys[1])[j], m*n*dat(ekeys[2])[j],
                        (has_exy ? m*n*dat(ekeys[3])[j]*sqrt(2.0) : 0.0),
                        (has_eyz ? m*n*dat(ekeys[4])[j]*sqrt(2.0) : 0.0),
                        (has_ezx ? m*n*dat(ekeys[5])[j]*sqrt(2.0) : 0.0);

            OctInvs (sig[i][j], p[i][j], q[i][j], t[i][j]);
            ev [i][j] = Calc_evoct (eps[i][j]);
            ed [i][j] = Calc_edoct (eps[i][j]);
            lp [i][j] = log(p[i][j]);
            qp [i][j] = q[i][j]/p[i][j];
            mqp[i][j] = (t[i][j]<0.0 ? -qp[i][j] : qp[i][j]);
        }
    }

#ifdef USE_WXWIDGETS
    // refresh cbx => replot
    CbxFNs->Set          (ShortNames);
    CbxFNs->SetSelection (0);
#endif
}

#ifdef USE_WXWIDGETS

enum
{
    ID_LOAD = wxID_HIGHEST+1,
    ID_SELDATA,
    ID_MULTQ,
    ID_PLTALL,
};

BEGIN_EVENT_TABLE(DatFiles, wxWindow)
    EVT_BUTTON   (ID_LOAD,    DatFiles::OnLoad)
    EVT_COMBOBOX (ID_SELDATA, DatFiles::OnReBuild)
    EVT_CHECKBOX (ID_MULTQ,   DatFiles::OnReBuild)
    EVT_CHECKBOX (ID_PLTALL,  DatFiles::OnReBuild)
END_EVENT_TABLE()

inline DatFiles::DatFiles (wxFrame * Parent)
    : wxWindow (Parent, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxFULL_REPAINT_ON_RESIZE),
      Multq  (true),
      PltAll (false)
{
    // force validation of child controls
    SetExtraStyle (wxWS_EX_VALIDATE_RECURSIVELY);

    // tell wxAuiManager to manage this frame
    Aui.SetManagedWindow (this);

    // plots
    qped = new GUI::PlotXY (this, "Deviatoric stress/strains",             "ed [%]", "q/p");
    qpev = new GUI::PlotXY (this, "Deviatoric stress - volumetric strain", "ev [%]", "q/p");
    eved = new GUI::PlotXY (this, "Volumetric strain - deviatoric strain", "ed [%]", "ev [%]");
    evlp = new GUI::PlotXY (this, "Volumetric strain - natural log (p)",   "ln(p)",  "ev [%]");
    qped->ShowLastY = false;
    qpev->ShowLastY = false;
    eved->ShowLastY = false;
    evlp->ShowLastY = false;

    // control panel
    ADD_WXPANEL    (pnl, szt, sz, 1, 4);
    ADD_WXBUTTON   (pnl, sz, ID_LOAD,    c0,     "Load data");
    ADD_WXCOMBOBOX (pnl, sz, ID_SELDATA, CbxFNs, "Data files");
    ADD_WXCHECKBOX (pnl, sz, ID_MULTQ,   c2,     "Lode multiply q", Multq);
    ADD_WXCHECKBOX (pnl, sz, ID_PLTALL,  c3,     "Plot all",        PltAll);
    CbxFNs->SetMinSize (wxSize(200,20));

    // commit all changes to wxAuiManager
    Aui.AddPane (pnl,  wxAuiPaneInfo().Name("cpnl").Caption("Data").Top   ().MinSize(wxSize(100,50)) .DestroyOnClose(false).CaptionVisible(true) .CloseButton(false));
    Aui.AddPane (qped, wxAuiPaneInfo().Name("qped").Caption("qped").Centre().Position(0)             .DestroyOnClose(false).CaptionVisible(false).CloseButton(false));
    Aui.AddPane (eved, wxAuiPaneInfo().Name("eved").Caption("eved").Centre().Position(1)             .DestroyOnClose(false).CaptionVisible(false).CloseButton(false));
    Aui.AddPane (evlp, wxAuiPaneInfo().Name("evlp").Caption("evlp").Right ().MinSize(wxSize(200,100)).DestroyOnClose(false).CaptionVisible(false).CloseButton(false));
    Aui.AddPane (qpev, wxAuiPaneInfo().Name("qpev").Caption("qpev").Right ().MinSize(wxSize(200,100)).DestroyOnClose(false).CaptionVisible(false).CloseButton(false));
    Aui.Update  ();

}

inline void DatFiles::ReBuild ()
{
    // update control's data
    Sync ();

    // filenames
    wxArrayString fnames = CbxFNs->GetStrings ();

    // disconnect plots
    qped->DelCurves ();
    qpev->DelCurves ();
    eved->DelCurves ();
    evlp->DelCurves ();

    // reconnect plots
    bool   mul = Multq;
    bool   all = PltAll;
    size_t ini = (all ? 0         : CbxFNs->GetSelection());
    size_t num = (all ? ed.Size() : ini+1);
    for (size_t i=ini; i<num; ++i)
    {
        // data
        GUI::CurveProps & c0 = qped->AddCurve (&ed[i], (mul ? &mqp[i] : &qp[i]), fnames[i].ToStdString().c_str());
        GUI::CurveProps & c1 = qpev->AddCurve (&ev[i], (mul ? &mqp[i] : &qp[i]), fnames[i].ToStdString().c_str());
        GUI::CurveProps & c2 = eved->AddCurve (&ed[i],                  &ev[i],  fnames[i].ToStdString().c_str());
        GUI::CurveProps & c3 = evlp->AddCurve (&lp[i],                  &ev[i],  fnames[i].ToStdString().c_str());
        c0.Typ=GUI::CT_BOTH;  c0.Psz=4;  c0.Pen.Set((all ? GUI::LinClr(i) : "black"), (i>5 ? "dash" : "solid"), 1);
        c1.Typ=GUI::CT_BOTH;  c1.Psz=4;  c1.Pen.Set((all ? GUI::LinClr(i) : "black"), (i>5 ? "dash" : "solid"), 1);
        c2.Typ=GUI::CT_BOTH;  c2.Psz=4;  c2.Pen.Set((all ? GUI::LinClr(i) : "black"), (i>5 ? "dash" : "solid"), 1);
        c3.Typ=GUI::CT_BOTH;  c3.Psz=4;  c3.Pen.Set((all ? GUI::LinClr(i) : "black"), (i>5 ? "dash" : "solid"), 1);

        // simulations
        if (ED.Size()>i)
        {
            wxString buf;  buf.Printf("sim:%s",fnames[i].ToStdString().c_str());
            GUI::CurveProps & C0 = qped->AddCurve (&ED[i], (mul ? &MQP[i] : &QP[i]), buf.ToStdString().c_str());
            GUI::CurveProps & C1 = qpev->AddCurve (&EV[i], (mul ? &MQP[i] : &QP[i]), buf.ToStdString().c_str());
            GUI::CurveProps & C2 = eved->AddCurve (&ED[i],                  &EV[i],  buf.ToStdString().c_str());
            GUI::CurveProps & C3 = evlp->AddCurve (&LP[i],                  &EV[i],  buf.ToStdString().c_str());
            C0.Typ=GUI::CT_LINES;  C0.Pen.Set((all ? GUI::LinClr(i) : "red"), (i>5 ? "dash" : "solid"), 1);
            C1.Typ=GUI::CT_LINES;  C1.Pen.Set((all ? GUI::LinClr(i) : "red"), (i>5 ? "dash" : "solid"), 1);
            C2.Typ=GUI::CT_LINES;  C2.Pen.Set((all ? GUI::LinClr(i) : "red"), (i>5 ? "dash" : "solid"), 1);
            C3.Typ=GUI::CT_LINES;  C3.Pen.Set((all ? GUI::LinClr(i) : "red"), (i>5 ? "dash" : "solid"), 1);
        }
    }

    // redraw
    qped->Redraw ();
    qpev->Redraw ();
    eved->Redraw ();
    evlp->Redraw ();
}

inline void DatFiles::OnLoad (wxCommandEvent & Event)
{
    wxFileDialog fd(this, "Select data files (.dat)", LstDir, "", "*.dat", wxFD_MULTIPLE);
    if (fd.ShowModal()==wxID_OK)
    {
        // get filenames
        wxArrayString paths, fnames;
        fd.GetPaths     (paths);
        fd.GetFilenames (fnames);

        // disconnect plots
        qped->DelCurves ();
        qpev->DelCurves ();
        eved->DelCurves ();
        evlp->DelCurves ();
        CbxFNs->Clear   ();

        // load data
        LstDir = fd.GetDirectory ();
        try { Read (paths, fnames); }
        catch (Fatal * e) { WxError(e->Msg().CStr()); }
    }
}

#endif // USE_WXWIDGETS

#endif // MECHSYS_DATFILES_H
