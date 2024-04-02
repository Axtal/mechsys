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

/**
  @mainpage MechSys: Reference Manual

  @section s_convlibs Convenience Libraries
  @subsection ss_util Utilities
  <ul>
  <li>namespaces: Util
  <li>util/array.h</li>
  <li>util/colors.h</li>
  <li>util/fatal.h</li>
  <li>util/maps.h</li>
  <li>util/numstreams.h</li>
  <li>util/stopwatch.h</li>
  <li>util/string.h</li>
  <li>util/tree.h</li>
  <li>util/util.h</li>
  </ul>

  @subsection ss_linalg Linear algebra and tensors
  <ul>
  <li>namespaces: LinAlg
  <li>linalg/jacobirot.h</li>
  <li>linalg/laexpr.h</li>
  <li>linalg/matrix.h</li>
  <li>linalg/matvec.h</li>
  <li>linalg/mumps.h</li>
  <li>linalg/sparse_crmatrix.h</li>
  <li>linalg/sparse_matrix.h</li>
  <li>linalg/sparse_triplet.h</li>
  <li>linalg/superlud.h</li>
  <li>linalg/superlu.h</li>
  <li>linalg/umfpack.h</li>
  <li>linalg/vector.h</li>
  </ul>

  @subsection ss_vtk Visualization toolkit
  <ul>
  <li>namespaces: VTK
  <li>vtk/arrow.h</li>
  <li>vtk/arrows.h</li>
  <li>vtk/axes.h</li>
  <li>vtk/cube.h</li>
  <li>vtk/cylinder.h</li>
  <li>vtk/disk.h</li>
  <li>vtk/isosurf.h</li>
  <li>vtk/line.h</li>
  <li>vtk/plane.h</li>
  <li>vtk/sgrid.h</li>
  <li>vtk/sphere.h</li>
  <li>vtk/spheres.h</li>
  <li>vtk/text2d.h</li>
  <li>vtk/text.h</li>
  <li>vtk/triangulate.h</li>
  <li>vtk/ugrid.h</li>
  <li>vtk/win.h</li>
  </ul>

  @subsection ss_py Python auxiliary scripts
  <ul>
  <li>python/msys_ana.py</li>
  <li>python/msys_drawmesh.py</li>
  <li>python/msys_fcrits.py</li>
  <li>python/msys_fig.py</li>
  <li>python/msys_invs.py</li>
  <li>python/msys_linfit.py</li>
  <li>python/msys_matvec.py</li>
  <li>python/msys_plt.py</li>
  </ul>

  @section s_mesh Mesh Generation
  <ul>
  <li>namespaces: Mesh</li>
  <li>mesh/alphashape.h</li>
  <li>mesh/paragrid3d.h</li>
  </ul>
  @subsection ss_generic Generic
  <ul>
  <li>mesh/mesh.h</li>
  </ul>
  @subsection ss_struct Structured
  <ul>
  <li>structured.h</li>
  </ul>
  @subsection ss_unstruct Unstructured
  <ul>
  <li>unstructured.h</li>
  </ul>

  @section s_blendergui Blender/CAD and GUI
  @subsection ss_blender Blender/CAD
  <ul>
  <li>blender/msys_blender_3dlink.py</li>
  <li>blender/msys_blender_cad.py</li>
  <li>blender/msys_blender_dict.py</li>
  <li>blender/msys_blender_fem.py</li>
  <li>blender/msys_blender_gui.py</li>
  <li>blender/msys_blender_main.py</li>
  <li>blender/msys_blender_mesh.py</li>
  <li>blender/msys_blender_mex3dlink.py</li>
  <li>blender/msys_blender_mexpt.py</li>
  <li>blender/msys_blender_shandler.py</li>
  </ul>

  @subsection ss_gui GUI: wxWidgets and FLTK
  <ul>
  <li>namespaces: GUI</li>
  <li>gui/common.h</li>
  <li>gui/plotxy.h</li>
  <li>gui/wxarrayint.h</li>
  <li>gui/wxdict.h</li>
  <li>gui/wxmyapp.h</li>
  <li>gui/wxnuminput.h</li>
  <li>gui/wxsipair.h</li>
  <li>gui/wxstringvalidator.h</li>
  </ul>

  @section s_models Material Models
  @section ss_general General interface
  <ul>
  <li>models/anisoinvs.h</li>
  <li>models/driver.h</li>
  <li>models/model.h</li>
  <li>models/smpinvs.h</li>
  </ul>
  @subsection ss_mech Mechanical models
  <ul>
  <li>models/bbmx.h</li>
  <li>models/camclay.h</li>
  <li>models/elastoplastic.h</li>
  <li>models/equilibstate.h</li>
  <li>models/linelastic.h</li>
  <li>models/neohookean.h</li>
  <li>models/nlelastic.h</li>
  <li>models/problemep.h</li>
  <li>models/strainupdate.h</li>
  <li>models/stressupdate.h</li>
  <li>models/unconv01.h</li>
  <li>models/unconv02.h</li>
  <li>models/unconv03.h</li>
  <li>models/unconv04.h</li>
  </ul>
  @subsection ss_flow Flow models
  <ul>
  <li>models/flowstate.h</li>
  <li>models/flowupdate.h</li>
  <li>models/linflow.h</li>
  </ul>
  @subsection ss_coupled Coupled models
  <ul>
  <li>models/hmstressupdate.h</li>
  <li>models/unsatflow.h</li>
  <li>models/unsatflowstate.h</li>
  </ul>

  @section s_num Numerical Methods
  <ul>
  <li>namespaces: Numerical</li>
  <li>numerical/min.h</li>
  <li>numerical/montecarlo.h</li>
  <li>numerical/multipleroots.h</li>
  <li>numerical/numdiff.h</li>
  <li>numerical/odesolver.h</li>
  <li>numerical/quadrature.h</li>
  <li>numerical/root.h</li>
  </ul>

  @section s_fem Finite Element Method (FEM)
  <ul>
  <li>namespaces: FEM</li>
  <li>fem/beam.h</li>
  <li>fem/domain.h</li>
  <li>fem/element.h</li>
  <li>fem/equilibelem.h</li>
  <li>fem/fem.h</li>
  <li>fem/flowelem.h</li>
  <li>fem/geomelem.h</li>
  <li>fem/hydromechelem.h</li>
  <li>fem/nlrod.h</li>
  <li>fem/node.h</li>
  <li>fem/quadrature.h</li>
  <li>fem/rod.h</li>
  <li>fem/solver.h</li>
  <li>fem/usigcondelem.h</li>
  <li>fem/usigelem.h</li>
  <li>fem/usigepselem.h</li>
  <li>fem/uwpelem.h</li>
  </ul>
  @subsection ss_elems Geometric elements
  <ul>
  <li>fem/elems/hex20.h</li>
  <li>fem/elems/hex8.h</li>
  <li>fem/elems/lin2.h</li>
  <li>fem/elems/quad4.h</li>
  <li>fem/elems/quad8.h</li>
  <li>fem/elems/tet10.h</li>
  <li>fem/elems/tri15.h</li>
  <li>fem/elems/tri3.h</li>
  <li>fem/elems/tri6.h</li>
  </ul>
  @subsection ss_solvers FE Solvers
  <ul>
  <li>fem/solvers/rksolver.h</li>
  <li>fem/solvers/stdsolver.h</li>
  <li>fem/solvers/uwpsolver.h</li>
  </ul>

  @section s_dem Discrete Element Method (DEM)
  <ul>
  <li>namespaces: DEM</li>
  <li>dem/basic_functions.h</li>
  <li>dem/cylinder.h</li>
  <li>dem/distance.h</li>
  <li>dem/domain.h</li>
  <li>dem/edge.h</li>
  <li>dem/face.h</li>
  <li>dem/graph.h</li>
  <li>dem/interacton.h</li>
  <li>dem/particle.h</li>
  <li>dem/quaternion.h</li>
  <li>dem/special_functions.h</li>
  <li>dem/torus.h</li>
  <li>dem/visualise.h</li>
  </ul>

  @section s_lbm Lattice Boltzmann Method (LBM)
  <ul>
  <li>namespaces: LBM</li>
  <li>lbm/Cell.h</li>
  <li>lbm/Dem.h</li>
  <li>lbm/Domain.h</li>
  <li>lbm/Interacton.h</li>
  <li>lbm/Lattice.h</li>
  </ul>

  @section s_sph Smoothed Particle Hydrodynamics (SPH)
  <ul>
  <li>namespaces: SPH</li>
  <li>sph/domain.h</li>
  <li>sph/interacton.h</li>
  <li>sph/particle.h</li>
  <li>sph/special_functions.h</li>
  </ul>

  @section s_mpm Material Point Method (MPM)
  <ul>
  <li>namespaces: MPM</li>
  <li>mpm/chull2d.h</li>
  <li>mpm/colormaps.h</li>
  <li>mpm/defs.h</li>
  <li>mpm/drawarea2d.h</li>
  <li>mpm/fmtnum.h</li>
  <li>mpm/grid2d.h</li>
  <li>mpm/infobox.h</li>
  <li>mpm/jacobirot.h</li>
  <li>mpm/mpoints2d.h</li>
  <li>mpm/output.h</li>
  <li>mpm/plotxy.h</li>
  <li>mpm/problems.h</li>
  <li>mpm/stress_update.h</li>
  <li>mpm/tensors.h</li>
  <li>mpm/tiled.h</li>
  </ul>

 */



/**
  @namespace Util
  @brief Convenience utilities
*/

/**
  @namespace Numerical
  @brief %Numerical methods
*/

/**
  @namespace FEM
  @brief Finite %Element Method
*/

/**
  @namespace DEM
  @brief Discrete %Element Method
*/

/**
  @namespace LBM
  @brief %Lattice Boltzmann Method
*/

/**
  @namespace SPH
  @brief Smoothed Particles Hydrodynamics
*/

/**
  @namespace MPM
  @brief Material Point Method
*/

/**
  @namespace LinAlg
  @brief Linear algebra routines
*/

/**
  @namespace ClrMap
  @brief %Colour map routines
*/

/**
  @namespace MPL
  @brief Matplotlib routines
*/

/**
  @namespace MPM::VTU
  @brief Paraview vtu files generation
*/

/**
  @namespace MPM::VTU::Out
  @brief Paraview output files
*/

/**
  @namespace Sparse
  @brief %Sparse matrices
*/

/**
  @namespace SuperLU
  @brief Linear solver for sparse systems
*/

/**
  @namespace SuperLUd
  @brief Parallel linear solver
*/

/**
  @namespace UMFPACK
  @brief Linear solver for sparse systems
*/

/**
  @namespace VTK
  @brief Visualization toolkit wrapper
*/

/**
  @namespace Mesh
  @brief %Mesh generation
*/

/**
  @namespace GUI
  @brief Graphical user interfaces
*/

/**
  @namespace Colors
  @brief Set of colors
*/

/**
  @namespace MUMPS
  @brief %MUMPS linear solver
*/

/**
  @namespace OrthoSys
  @brief Orthogonal Cartesian system
*/


/** @namespace msys_ana
 *  @brief Python auxiliary script */
/** @namespace msys_drawmesh
 *  @brief Python auxiliary script */
/** @namespace msys_fcrits
 *  @brief Python auxiliary script */
/** @namespace msys_fig
 *  @brief Python auxiliary script */
/** @namespace msys_invs
 *  @brief Python auxiliary script */
/** @namespace msys_linfit
 *  @brief Python auxiliary script */
/** @namespace msys_matvec
 *  @brief Python auxiliary script */
/** @namespace msys_plt
 *  @brief Python auxiliary script */

/** @namespace msys_blender_3dlink
 *  @brief Blender graphical interface module */
/** @namespace msys_blender_cad
 *  @brief Blender graphical interface module */
/** @namespace msys_blender_dict
 *  @brief Blender graphical interface module */
/** @namespace msys_blender_fem
 *  @brief Blender graphical interface module */
/** @namespace msys_blender_gui
 *  @brief Blender graphical interface module */
/** @namespace msys_blender_main
 *  @brief Blender graphical interface module */
/** @namespace msys_blender_mesh
 *  @brief Blender graphical interface module */
/** @namespace msys_blender_mex3dlink
 *  @brief Blender graphical interface module */
/** @namespace msys_blender_mexpt
 *  @brief Blender graphical interface module */
/** @namespace msys_blender_shandler
 *  @brief Blender graphical interface module */
