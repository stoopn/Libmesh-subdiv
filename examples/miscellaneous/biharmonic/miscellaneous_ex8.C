/* The Next Great Finite Element Library. */
/* Copyright (C) 2003  Benjamin S. Kirk */
// C++ include files that we need
#include <iostream>

// LibMesh include files.
#include "libmesh.h"
#include "mesh.h"
#include "mesh_refinement.h"
#include "mesh_modification.h"
#include "mesh_tools.h"
#include "linear_implicit_system.h"
#include "equation_systems.h"
#include "fe.h"
#include "quadrature.h"
#include "node.h"
#include "elem.h"
#include "vector_value.h"
#include "tensor_value.h"
#include "dense_matrix.h"
#include "dense_submatrix.h"
#include "dense_vector.h"
#include "dense_subvector.h"
#include "sparse_matrix.h"
#include "numeric_vector.h"
#include "vtk_io.h"
#include "exodusII_io.h"

// These are the include files typically needed for subdivision elements.
#include "dof_map_subdiv.h"
#include "face_tri3_sd.h"
#include "mesh_subdiv_support.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Function prototype.  This is the function that will assemble
// the stiffness matrix and the right-hand-side vector ready
// for solution.
void assemble_shell (EquationSystems& es, const std::string& system_name);

// Begin the main program.
int main (int argc, char** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);
  
  // Skip this 3D example if libMesh was compiled as 1D/2D-only.
  libmesh_example_assert(3 == LIBMESH_DIM, "3D support");
  
  // Tell the user what we are doing.
  std::cout << "Running example miscellaneous_ex8\n" << std::endl;        
  
  // Create the mesh for the 2D plate.
  Mesh mesh (2);

  // Read the coarse square mesh.
  //mesh.read ("square_mesh.off");
    mesh.read ("disc2.off");
    

  // Resize the square plate to edge length L.
  const Real L = 100.;
  MeshTools::Modification::scale(mesh, L, L, L);

  // Quadrisect the mesh triangles a few times to obtain a
  // finer mesh.  Subdivision surface elements require the
  // refinement data to be removed afterwards.
  MeshRefinement mesh_refinement (mesh);
  mesh_refinement.uniformly_refine (1);
  MeshTools::Modification::flatten (mesh);

  // Write the mesh before the ghost elements are added.
#if defined(LIBMESH_HAVE_VTK)
  VTKIO(mesh).write ("without_ghosts.pvtu");
#endif
#if defined(LIBMESH_HAVE_EXODUS_API)
  ExodusII_IO(mesh).write ("without_ghosts.e");
#endif

  // Print information about the triangulated mesh to the screen.
  mesh.print_info();

  // Turn the triangulated mesh into a subdivision mesh
  // and add an additional row of "ghost" elements around
  // it in order to complete the extended local support of
  // the triangles at the boundaries.  If the second
  // argument is set to true, the outhermost existing
  // elements are converted into ghost elements, and the
  // actual physical mesh is thus getting smaller.
  MeshTools::Subdiv::prepare_subdiv_mesh (mesh, false);
  
  // Print information about the subdivision mesh to the screen.
  mesh.print_info();

  // Write the mesh with the ghost elements added.
  // Compare this to the original mesh to see the difference.
#if defined(LIBMESH_HAVE_VTK)
  VTKIO(mesh).write ("with_ghosts.pvtu");
#endif
#if defined(LIBMESH_HAVE_EXODUS_API)
  ExodusII_IO(mesh).write ("with_ghosts.e");
#endif
  
  // Create an equation systems object.
  EquationSystems equation_systems (mesh);
  
  // Declare the system and its variables.
  // Create an explicit system named "Shell".
  LinearImplicitSystem & system = equation_systems.add_system<LinearImplicitSystem> ("Shell");

  // Subdivision elements have a larger local support than
  // conventional elements, which requires the use of a
  // special subdivision DofMap.  We thus need to replace
  // the default DofMap by DofMapSubdiv.
  system.replace_dof_map (AutoPtr<DofMap>(new DofMapSubdiv(system.number())));

  // Add the three translational deformation variables
  // "u", "v", "w" to "Shell".  Since subdivision shell
  // elements meet the C1-continuity requirement, no
  // rotational or other auxiliary variables are needed.
  // Loop Subdivision Elements are always interpolated
  // by quartic box splines, hence the order must always
  // be \p FOURTH.
  system.add_variable ("u", FOURTH, SUBDIV);
    
  // Give the system a pointer to the matrix and rhs assembly
  // function.
  system.attach_assemble_function (assemble_shell);

  // Use the parameters of the equation systems object to
  // tell the shell system about the material properties, the
  // shell thickness, and the external load.
  
  // Initialize the data structures for the equation system.
  equation_systems.init();
  
  // Print information about the system to the screen.
  equation_systems.print_info();

  // Solve the linear system.
  system.solve();

  // After solving the system, write the solution to a VTK
  // or ExodusII output file ready for import in, e.g.,
  // Paraview.
#if defined(LIBMESH_HAVE_VTK)
  VTKIO(mesh).write_equation_systems ("out.pvtu", equation_systems);
#endif
#if defined(LIBMESH_HAVE_EXODUS_API)
  ExodusII_IO(mesh).write_equation_systems ("out.e", equation_systems);
#endif

  
  // All done.
  return 0;
}

typedef struct {
    RealVectorValue a[2];  //Tangent vectors (2D surface)
    DenseMatrix<Real> gcov;  //Covariant metric tensor components
    DenseMatrix<Real> gcon;  //Contravariant metric tensor components
    RealVectorValue acovd[2][2];  //Covariant derivatives of tangent vectors a_(a,b)
    double g_abc[2][2][2];  //Covariant derivatives of metric tensor g_(ab,c)
    std::vector<Real> LaplaceSphi;
    RealVectorValue n;  //Surface normal (not normalized)
    double J1;    //Length of surface normal = Jacobian of mapping from parametric to physical space
    double invgdet;  // 1/det(gcov)   => needed for construction of contravariant metric
} SurfaceData;


double CHR(int m, int i, int j,  SurfaceData &SD)
{
    int k; double s = 0.0;
    for (k=0; k<2; k++)
    {
        s += 0.5*SD.gcon(k,m)*(SD.g_abc[i][k][j] + SD.g_abc[j][k][i] - SD.g_abc[i][j][k]);
    }
    return s;
}


//CHECKED AND CORRECT
void setup_surface_data(int qp, SurfaceData& SD, AutoPtr<FEBase>& element_fe)
{
	START_LOG("setup_surface_metrics", "ShellElement");
    
    SD.gcov.resize(2,2);
    SD.gcon.resize(2,2);
    int n_dofs = element_fe->get_phi().size();
    SD.LaplaceSphi.resize(n_dofs);
    
	const std::vector<std::vector<RealGradient> >& dphi = element_fe->get_dphi();
	const std::vector<std::vector<RealTensor> >& d2phi = element_fe->get_d2phi();
	const std::vector< RealGradient >& dxyzdxi = element_fe->get_dxyzdxi();
	const std::vector< RealGradient >& dxyzdeta = element_fe->get_dxyzdeta();
	const std::vector< RealGradient >& d2xyzdxi2 = element_fe->get_d2xyzdxi2();
	const std::vector< RealGradient >& d2xyzdeta2 = element_fe->get_d2xyzdeta2();
	const std::vector< RealGradient >& d2xyzdxideta = element_fe->get_d2xyzdxideta();
	RealGradient ddispdxi, ddispdeta, d2dispdxi2, d2dispdeta2, d2dispdxideta;
	Real det;
    
	SD.a[0] = dxyzdxi[qp];
	SD.a[1] = dxyzdeta[qp];
	SD.n = SD.a[0].cross(SD.a[1]);
	//SD.invgdet = sqrt(ar.acov[2]*ar.acov[2]);
	//SD.invgdet = 1/(SD.n.size());
	
	SD.acovd[0][0] = d2xyzdxi2[qp];
	SD.acovd[0][1] = SD.acovd[1][0] = d2xyzdxideta[qp];
	SD.acovd[1][1] = d2xyzdeta2[qp];
	
	// Covariant metric tensor:
    SD.gcov(0,0) = SD.a[0]*SD.a[0];
    SD.gcov(1,0) = SD.gcov(0,1) = SD.a[0]*SD.a[1];
    SD.gcov(1,1) = SD.a[1]*SD.a[1];
	//Contravariant metric tensor:
	det = (SD.gcov(0,0)*SD.gcov(1,1) - SD.gcov(1,0)*SD.gcov(1,0));
	SD.invgdet = 1./det;
    SD.gcon(0,0) = SD.invgdet*SD.gcov(1,1);
    SD.gcon(1,0) = SD.gcon(0,1) = -SD.invgdet*SD.gcov(0,1);
    SD.gcon(1,1) = SD.invgdet*SD.gcov(0,0);
	
    //covariant derivatives of metric tensor g_(ab,c)
    SD.g_abc[0][0][0] = 2. * SD.a[0]*SD.acovd[0][0];
    SD.g_abc[0][0][1] = 2. * SD.a[0]*SD.acovd[0][1];
    SD.g_abc[0][1][0] = SD.g_abc[1][0][0] = SD.a[1]*SD.acovd[0][0] + SD.acovd[1][0]*SD.a[0];
    SD.g_abc[0][1][1] = SD.g_abc[1][0][1] = SD.a[1]*SD.acovd[0][1] + SD.acovd[1][1]*SD.a[0];
    SD.g_abc[1][1][0] = 2. * SD.a[1]*SD.acovd[1][0];
    SD.g_abc[1][1][1] = 2. * SD.a[1]*SD.acovd[1][1];
    
    for (int a=0; a<n_dofs; a++) {
        double laplace = 0.0;
        //Contract with contravariant metric:
        int j,k;
        for (j=0; j<2; j++)
        {
            for (k=0; k<2; k++) {
                // laplace += SD.gcon(j,k)*N2[a][j][k] - SD.gcon(j,k)*(CHR(0,j,k,SD.g_abc,SD.gcon)*N1[a][0]
                //                                                       + CHR(1,j,k,SD.g_abc,SD.gcon)*N1[a][1]);
                laplace += SD.gcon(j,k)*d2phi[a][qp](j,k) -SD.gcon(j,k)*(CHR(0,j,k,SD)*dphi[a][qp](0)
                                                                         + CHR(1,j,k,SD)*dphi[a][qp](1));
            }
        }
        SD.LaplaceSphi[a]=laplace;
        //  PetscPrintf(PETSC_COMM_SELF,"LaplaceSN=%g\n",laplace);
    }
    
}

// We now define the matrix and rhs vector assembly function
// for the shell system.  This function implements the
// linear Kirchhoff-Love theory for thin shells.  At the
// end we also take into account the boundary conditions
// here, using the penalty method.
void assemble_shell (EquationSystems& es, const std::string& system_name)
{
  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert(system_name == "Shell");

  // Get a constant reference to the mesh object.
  const MeshBase& mesh = es.get_mesh();

  // Get a reference to the shell system object.
  LinearImplicitSystem & system = es.get_system<LinearImplicitSystem> ("Shell");


  // Numeric ids corresponding to each variable in the system.
  const unsigned int u_var = system.variable_number ("u");

    // Get the Finite Element type for "u".  Note this will be
  // the same as the type for "v" and "w".
  FEType fe_type = system.variable_type (u_var);

  // Build a Finite Element object of the specified type.
  AutoPtr<FEBase> fe (FEBase::build(2, fe_type));
  
  // A Gauss quadrature rule for numerical integration.
  // For subdivision shell elements, a single Gauss point per
  // element is sufficient, hence we use extraorder = 0.
  const int extraorder = 0;
  AutoPtr<QBase> qrule (fe_type.default_quadrature_rule (2, extraorder));

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule (qrule.get());

  // The element Jacobian * quadrature weight at each integration point.   
  const std::vector<Real>& JxW = fe->get_JxW();

  // The surface tangents in both directions at the quadrature points. 
  const std::vector<RealGradient>& dxyzdxi  = fe->get_dxyzdxi();
  const std::vector<RealGradient>& dxyzdeta = fe->get_dxyzdeta();

  // The second partial derivatives at the quadrature points.
  const std::vector<RealGradient>& d2xyzdxi2    = fe->get_d2xyzdxi2();
  const std::vector<RealGradient>& d2xyzdeta2   = fe->get_d2xyzdeta2();
  const std::vector<RealGradient>& d2xyzdxideta = fe->get_d2xyzdxideta();

  // The element shape function and its derivatives evaluated at the
  // quadrature points.
  const std::vector<std::vector<Real> >&          phi = fe->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
  const std::vector<std::vector<RealTensor> >&  d2phi = fe->get_d2phi();

  // A reference to the \p DofMap object for this system.  The \p DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.
  const DofMap & dof_map = system.get_dof_map();

  // Define data structures to contain the element stiffness matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe".
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

   // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<unsigned int> dof_indices;
  std::vector<unsigned int> dof_indices_u;
  
  // Now we will loop over all the elements in the mesh.  We will
  // compute the element matrix and right-hand-side contribution.
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for (; el != end_el; ++el)
  {
    // Store a pointer to the element we are currently
    // working on.  This allows for nicer syntax later.
    const Elem* elem = *el;

    // The ghost elements at the boundaries need to be excluded
    // here, as they don't belong to the physical shell,
    // but serve for a proper boundary treatment only.
    libmesh_assert(elem->type() == TRI3SD);
    const Tri3SD* sd_elem = static_cast<const Tri3SD*> (elem);
    if (sd_elem->is_ghost())
      continue;

    // Get the degree of freedom indices for the
    // current element.  These define where in the global
    // matrix and right-hand-side this element will
    // contribute to.
    dof_map.dof_indices (elem, dof_indices);
    dof_map.dof_indices (elem, dof_indices_u, u_var);
    
    const unsigned int n_dofs   = dof_indices.size();
    const unsigned int n_u_dofs = dof_indices_u.size();

    // Compute the element-specific data for the current
    // element.  This involves computing the location of the
    // quadrature points and the shape functions
    // (phi, dphi, d2phi) for the current element.
    fe->reinit (elem);

    // Zero the element matrix and right-hand side before
    // summing them.  We use the resize member here because
    // the number of degrees of freedom might have changed from
    // the last element.
    Ke.resize (n_dofs, n_dofs);
    Fe.resize (n_dofs);
  
    // Now we will build the element matrix and right-hand-side.
    for (unsigned int qp=0; qp<qrule->n_points(); ++qp)
    {
      // First, we compute the external force resulting
      // from a load q distributed uniformly across the plate.
      // Since the load is supposed to be transverse to the plate,
      // it affects the z-direction, i.e. the "w" variable.
      for (unsigned int i=0; i<n_u_dofs; ++i)
        Fe(i) += JxW[qp] * phi[i][qp] * 0.0000001;
      
      // Next, we assemble the stiffness matrix.  This is only valid
      // for the linear theory, i.e., for small deformations, where
      // reference and deformed surface metrics are indistinguishable.

        SurfaceData SD;
        setup_surface_data(qp,SD,fe);
        
      // Loop over all pairs of nodes I,J.
      for (unsigned int a=0; a<n_u_dofs; ++a)
      {
        for (unsigned int b=0; b<n_u_dofs; ++b)
        {    //Laplace term
           /* Ke(a,b) = JxW[qp]*(  (SD.gcon(0,0)*dphi[a][qp](0)*dphi[b][qp](0)
                                + SD.gcon(0,1)*dphi[a][qp](0)*dphi[b][qp](1)
                                + SD.gcon(1,0)*dphi[a][qp](1)*dphi[b][qp](0)
                                + SD.gcon(1,1)*dphi[a][qp](1)*dphi[b][qp](1))
                               );*/
            //Biharmonic term
            Ke(a,b) = JxW[qp]*1.0*(SD.LaplaceSphi[a]*SD.LaplaceSphi[b]);
        }
      }

    } // end of the quadrature point qp-loop

    // The element matrix and right-hand-side are now built
    // for this element.  Add them to the global matrix and
    // right-hand-side vector.  The \p NumericMatrix::add_matrix()
    // and \p NumericVector::add_vector() members do this for us.
    system.matrix->add_matrix (Ke, dof_indices_u);
    system.rhs->add_vector    (Fe, dof_indices_u);
  } // end of non-ghost element loop
  
  // Next, we apply the boundary conditions.  In this case,
  // all boundaries are clamped by the penalty method, using
  // the special "ghost" nodes along the boundaries.  Note
  // that there are better ways to implement boundary conditions
  // for subdivision shells.  We use the simplest way here,
  // which is known to be overly restrictive and will lead to
  // a slightly too small deformation of the plate.
  el = mesh.active_local_elements_begin();

  for (; el != end_el; ++el)
  {
    // Store a pointer to the element we are currently
    // working on.  This allows for nicer syntax later.
    const Elem* elem = *el;

    // For the boundary conditions, we only need to loop over
    // the ghost elements.
    libmesh_assert(elem->type() == TRI3SD);
    const Tri3SD* gh_elem = static_cast<const Tri3SD*> (elem);
    if (!gh_elem->is_ghost())
      continue;

    // Find the side which is part of the physical plate boundary,
    // that is, the boundary of the original mesh without ghosts.
    for (unsigned int s=0; s<elem->n_sides(); ++s)
    {
      const Tri3SD* nb_elem = static_cast<const Tri3SD*> (elem->neighbor(s));
      if (nb_elem == NULL || nb_elem->is_ghost())
        continue;

      // Determine the four nodes involved in the boundary
      // condition treatment of this side.  The \p MeshTools::Subdiv
      // namespace provides lookup tables \p next and \p prev
      // for an efficient determination of the next and previous
      // nodes of an element, respectively.
      //
      //      n4
      //     /  \
      //    / gh \
      //  n2 ---- n3
      //    \ nb /
      //     \  /
      //      n1
      Node* nodes [4]; // n1, n2, n3, n4
      nodes[1] = gh_elem->get_node(s); // n2
      nodes[2] = gh_elem->get_node(MeshTools::Subdiv::next[s]); // n3
      nodes[3] = gh_elem->get_node(MeshTools::Subdiv::prev[s]); // n4

      // The node in the interior of the domain, \p n1, is the
      // hardest to find.  Walk along the edges of element \p nb until
      // we have identified it.
      unsigned int n = 0;
      nodes[0] = nb_elem->get_node(0);
      while (nodes[0]->id() == nodes[1]->id() || nodes[0]->id() == nodes[2]->id())
        nodes[0] = nb_elem->get_node(++n);

      // The penalty value.  \f$ \frac{1}{\epsilon} \f$
      const Real penalty = 1.e10;

      // With this simple method, clamped boundary conditions are
      // obtained by penalizing the displacements of all four nodes.
      // This ensures that the displacement field vanishes on the
      // boundary side \p s.
      for (unsigned int n=0; n<4; ++n)
      {
        const unsigned int u_dof = nodes[n]->dof_number (system.number(), u_var, 0);
        system.matrix->add (u_dof, u_dof, penalty);
      }
    }
  } // end of ghost element loop
}
