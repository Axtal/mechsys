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

// Boost-Python
#include <boost/python.hpp> // this includes everything

#define USE_BOOST_PYTHON

namespace BPy = boost::python;

// MechSys
#include <mechsys/matfile.h>
#include <mechsys/inpfile.h>
#include <mechsys/linalg/jacobirot.h>
#include <mechsys/numerical/odesolver.h>
#include <mechsys/models/smpinvs.h>

// MechSys -- FEM
#include <mechsys/fem/fem.h>

// functions overloadings
BOOST_PYTHON_FUNCTION_OVERLOADS (FUN_PHI2M,        Phi2M,        1, 2)
BOOST_PYTHON_FUNCTION_OVERLOADS (FUN_M2PHI,        M2Phi,        1, 2)
BOOST_PYTHON_FUNCTION_OVERLOADS (FUN_JACOBIROT,    PyJacobiRot,  3, 4)

// member overloadings
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (OD_Init,         Init,         2, 6)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (SI_Calc,         PyCalc,       1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (IN_SetPrmsInis,  SetPrmsInis,  1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MG_ReadMesh,     ReadMesh,     1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MG_SetVert,      SetVert,      4, 5)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MG_WriteVTU,     WriteVTU,     1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MG_WriteMPY,     WriteMPY,     1, 5)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MS_Generate,     PyGenerate,   1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MS_GenBox,       GenBox,       0, 7)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MU_Generate,     Generate,     0, 4)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MU_GenBox,       GenBox,       0, 5)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (MU_WritePLY,     WritePLY,     1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (DO_PrintResults, PrintResults, 0, 4)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (DO_WriteVTU,     WriteVTU,     1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (SO_Solve,        Solve,        0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS (SO_DynSolve,     DynSolve,     3, 4)

// module
BOOST_PYTHON_MODULE (mechsys)
{

//////////////////////////////////////////////////////////////////////////////////// util /////

// String
BPy::class_<String>("String")
    .def("PyStr", &String::PyStr)
    .def(BPy::self_ns::str(BPy::self))
    ;

// SDPair
BPy::class_<SDPair>("SDPair")
    .def("Set", &SDPair::PySet)
    .def(BPy::self_ns::str(BPy::self))
    ;

// Dict
BPy::class_<Dict>("Dict")
    .def("Set", &Dict::PySet)
    .def("Get", &Dict::PyGet)
    .def(BPy::self_ns::str(BPy::self))
    ;

// Table
BPy::class_<Table>("Table")
    .def("Set", &Table::PySet)
    .def(BPy::self_ns::str(BPy::self))
    ;

// Fatal
BPy::register_exception_translator<Fatal *>(&PyExceptTranslator);

// Util
BPy::def("FindBestSquare", Util::PyFindBestSquare);

/////////////////////////////////////////////////////////////////////////////// Numerical /////

BPy::class_<Numerical::PyODESolver>("ODESolver", "ODESolver", BPy::init<BPy::object &, char const *, char const *>())
    .def("Init",        &Numerical::PyODESolver::Init, OD_Init())
    .def("Evolve",      &Numerical::PyODESolver::Evolve)
    .def("EvolveOut",   &Numerical::PyODESolver::EvolveOut)
    .def_readwrite("Y", &Numerical::PyODESolver::y)
    ;

/////////////////////////////////////////////////////////////////////////////////// linalg ////

BPy::def("pqth2L",    Pypqth2L);
BPy::def("Phi2M",     Phi2M,       FUN_PHI2M());
BPy::def("M2Phi",     M2Phi,       FUN_M2PHI());
BPy::def("JacobiRot", PyJacobiRot, FUN_JACOBIROT());
BPy::def("EigenProjAnalytic", PyEigenProjAnalytic);

/////////////////////////////////////////////////////////////////////////////// Invariants ////

BPy::class_<SMPInvs>("SMPInvs")
    .def("Calc", &SMPInvs::PyCalc, SI_Calc())
    .def_readwrite("b", &SMPInvs::b)
    .def(BPy::self_ns::str(BPy::self))
    ;

///////////////////////////////////////////////////////////////////////////////// MatFile /////

BPy::class_<MatFile>("MatFile")
    .def("Read",      &MatFile::Read)
    .def("GetPrmIni", &MatFile::PyGetPrmIni)
    .def(BPy::self_ns::str(BPy::self))
    ;

///////////////////////////////////////////////////////////////////////////////// InpFile /////

BPy::class_<InpFile>("InpFile")
    .def("Read",        &InpFile::Read)
    .def("SetPrmsInis", &InpFile::SetPrmsInis, IN_SetPrmsInis())
    .def("ReadPrmIni",  &InpFile::PyReadPrmIni)
    .def(BPy::self_ns::str(BPy::self))
    .def_readwrite("ninc"       , &InpFile::ninc      ) //   1 
    .def_readwrite("cdrift"     , &InpFile::cdrift    ) //   2
    .def_readwrite("stol"       , &InpFile::stol      ) //   3
    .def_readwrite("ssout"      , &InpFile::ssout     ) //   4
    .def_readwrite("ctetg"      , &InpFile::ctetg     ) //   5
    .def_readwrite("fem"        , &InpFile::fem       ) //   6
    .def_readwrite("dyn"        , &InpFile::dyn       ) //   7
    .def_readwrite("hm"         , &InpFile::hm        ) //   8
    .def_readwrite("tf"         , &InpFile::tf        ) //   9
    .def_readwrite("dt"         , &InpFile::dt        ) //  10
    .def_readwrite("dtout"      , &InpFile::dtout     ) //  11
    .def_readwrite("tsw"        , &InpFile::tsw       ) //  12
    .def_readwrite("ndiv"       , &InpFile::ndiv      ) //  13
    .def_readwrite("nip"        , &InpFile::nip       ) //  14
    .def_readwrite("o2"         , &InpFile::o2        ) //  15
    .def_readwrite("ray"        , &InpFile::ray       ) //  16
    .def_readwrite("am"         , &InpFile::am        ) //  17
    .def_readwrite("ak"         , &InpFile::ak        ) //  18
    .def_readwrite("rk"         , &InpFile::rk        ) //  19
    .def_readwrite("rkscheme"   , &InpFile::rkscheme  ) //  20
    .def_readwrite("rkstol"     , &InpFile::rkstol    ) //  21
    .def_readwrite("refdat"     , &InpFile::refdat    ) //  22
    .def_readwrite("refsim"     , &InpFile::refsim    ) //  23
    .def_readwrite("refana"     , &InpFile::refana    ) //  24
    .def_readwrite("idxvert1"   , &InpFile::idxvert1  ) //  25
    .def_readwrite("idxvert2"   , &InpFile::idxvert2  ) //  26
    .def_readwrite("idxvert3"   , &InpFile::idxvert3  ) //  27
    .def_readwrite("optdbl1"    , &InpFile::optdbl1   ) //  28
    .def_readwrite("optdbl2"    , &InpFile::optdbl2   ) //  29
    .def_readwrite("optdbl3"    , &InpFile::optdbl3   ) //  30
    .def_readwrite("hasoptdbl1" , &InpFile::hasoptdbl1) //  28b
    .def_readwrite("hasoptdbl2" , &InpFile::hasoptdbl2) //  29b
    .def_readwrite("hasoptdbl3" , &InpFile::hasoptdbl3) //  30b
    .def_readwrite("nldt_nsml"  , &InpFile::nldt_nsml ) //  31
    .def_readwrite("nldt_nn"    , &InpFile::nldt_nn   ) //  32
    .def_readwrite("nldt_n"     , &InpFile::nldt_n    ) //  33
    .def_readwrite("nldt_ll"    , &InpFile::nldt_ll   ) //  34
    .def_readwrite("nldt_sch"   , &InpFile::nldt_sch  ) //  35
    .def_readwrite("nldt_m"     , &InpFile::nldt_m    ) //  36
    .def_readwrite("maxit"      , &InpFile::maxit     ) //  37
    .def_readwrite("tolr"       , &InpFile::tolr      ) //  38
    .def_readwrite("fnkey"      , &InpFile::fnkey     ) //  39
    .def_readwrite("pcam0"      , &InpFile::pcam0     ) //  40
    .def_readwrite("haspcam0"   , &InpFile::haspcam0  ) //  40b
    .def_readwrite("scheme"     , &InpFile::scheme    ) //  41
    .def_readwrite("vtufile"    , &InpFile::vtufile   ) //  42
    .def_readwrite("suscheme"   , &InpFile::suscheme  ) //  43
    .def_readwrite("sustol"     , &InpFile::sustol    ) //  44
    .def_readwrite("surkscheme" , &InpFile::surkscheme) //  45
    .def_readwrite("dcmaxit"    , &InpFile::dcmaxit   ) //  46
    .def_readwrite("dcftol"     , &InpFile::dcftol    ) //  47
    .def_readwrite("pw0"        , &InpFile::pw0       ) //  48
    .def_readwrite("rkdyncte"   , &InpFile::rkdyncte  ) //  49
    .def_readwrite("uwp"        , &InpFile::uwp       ) //  50
    ;

//////////////////////////////////////////////////////////////////////////////////// mesh /////

// Block
BPy::class_<Mesh::Block>("Block")
    .def("Set", &Mesh::Block::PySet)
    .def(BPy::self_ns::str(BPy::self))
    ;

// Generic
BPy::class_<Mesh::Generic>("Generic","generic mesh", BPy::init<int>())
    .def("ReadMesh",      &Mesh::Generic::ReadMesh, MG_ReadMesh())
    .def("SetSize",       &Mesh::Generic::SetSize)
    .def("SetVert",       &Mesh::Generic::SetVert,  MG_SetVert())
    .def("SetCell",       &Mesh::Generic::PySetCell)
    .def("SetBryTag",     &Mesh::Generic::SetBryTag)
    .def("AddLinCells",   &Mesh::Generic::PyAddLinCells)
    .def("WriteVTU",      &Mesh::Generic::WriteVTU, MG_WriteVTU())
    .def("WriteMPY",      &Mesh::Generic::WriteMPY, MG_WriteMPY())
    .def("GetVertsEdges", &Mesh::Generic::PyGetVertsEdges)
    .def(BPy::self_ns::str(BPy::self))
    ;

// Structured
BPy::class_<Mesh::Structured, BPy::bases<Mesh::Generic> >("Structured","structured mesh", BPy::init<int>())
    .def("Generate", &Mesh::Structured::PyGenerate, MS_Generate())
    .def("GenBox",   &Mesh::Structured::GenBox,     MS_GenBox())
    .def(BPy::self_ns::str(BPy::self))
    ;

// Unstructured
BPy::class_<Mesh::Unstructured, BPy::bases<Mesh::Generic> >("Unstructured","Unstructured mesh", BPy::init<int>())
    .def("Set",      &Mesh::Unstructured::PySet)
    .def("Generate", &Mesh::Unstructured::Generate, MU_Generate())
    .def("GenBox",   &Mesh::Unstructured::GenBox,   MU_GenBox())
    .def("WritePLY", &Mesh::Unstructured::WritePLY, MU_WritePLY())
    .def(BPy::self_ns::str(BPy::self))
    ;

///////////////////////////////////////////////////////////////////////////////////// fem /////

// PROB, GEOM, and MODEL
BPy::def("PROB",  PyPROB);
BPy::def("GEOM",  PyGEOM);
BPy::def("MODEL", PyMODEL);

// Domain
BPy::class_<FEM::Domain>("FEM_Domain", "FEM domain", BPy::init<Mesh::Generic const &, Dict const &, Dict const &, Dict const &>())
    .def("SetBCs",       &FEM::Domain::SetBCs)
    .def("PrintResults", &FEM::Domain::PrintResults, DO_PrintResults())
    .def("WriteVTU",     &FEM::Domain::WriteVTU,     DO_WriteVTU())
    .def(BPy::self_ns::str(BPy::self))
    ;

// Solver
//BPy::class_<FEM::Solver>("Solver","generic solver", BPy::init<FEM::Domain &, SDPair>())
    //;

// STDSolver
//BPy::class_<FEM::STDSolver, BPy::bases<FEM::Solver> >("FEM_STDSolver", "FEM (standard) solver", BPy::init<FEM::Domain &, SDPair>())
    //.def("Solve",     &FEM::STDSolver::Solve,    SO_Solve())
    //.def("DynSolve",  &FEM::STDSolver::DynSolve, SO_DynSolve())
    //.def_readwrite("MaxIt", &FEM::STDSolver::MaxIt)
    //.def_readwrite("TolR",  &FEM::STDSolver::TolR)
    //;

} // BOOST_PYTHON_MODULE
