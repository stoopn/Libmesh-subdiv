// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



// C++ Includes -------------------------------------
#include <set>
#include <algorithm>

// Local Includes -----------------------------------
#include "coupling_matrix.h"
#include "dense_matrix.h"
#include "dense_vector_base.h"
#include "dof_map_subdiv.h"
#include "elem.h"
#include "face_tri3_sd.h"
#include "fe_interface.h"
#include "fe_type.h"
#include "libmesh_logging.h"
#include "mesh_base.h"
#include "mesh_tools.h"
#include "mesh_subdiv_support.h"
#include "numeric_vector.h"
#include "parallel.h"
#include "sparse_matrix.h"
#include "sparsity_pattern.h"
#include "string_to_enum.h"
#include "threads_allocators.h"


namespace libMesh
{

// ------------------------------------------------------------
// DofMap member functions


void DofMapSubdiv::reinit(MeshBase& mesh)
{
  libmesh_assert (mesh.is_prepared());

  START_LOG("reinit()", "DofMapSubdiv");

  //this->clear();

  const unsigned int n_var = this->n_variables();

  //------------------------------------------------------------
  // Then set the number of variables for each \p DofObject
  // equal to n_variables() for this system.  This will
  // handle new \p DofObjects that may have just been created
  {
    // All the nodes
    MeshBase::node_iterator       node_it  = mesh.nodes_begin();
    const MeshBase::node_iterator node_end = mesh.nodes_end();

    for ( ; node_it != node_end; ++node_it)
      (*node_it)->set_n_vars(this->sys_number(),n_var);

    // All the elements
    MeshBase::element_iterator       elem_it  = mesh.elements_begin();
    const MeshBase::element_iterator elem_end = mesh.elements_end();

    for ( ; elem_it != elem_end; ++elem_it)
      (*elem_it)->set_n_vars(this->sys_number(),n_var);
  }

  //------------------------------------------------------------
  // Next allocate space for the DOF indices
  for (unsigned int var=0; var<this->n_variables(); var++)
    {
      const Variable &var_description =	this->variable(var);
      const FEType& base_fe_type              = this->variable_type(var);

      // For all the active elements
      MeshBase::element_iterator       elem_it  = mesh.active_elements_begin();
      const MeshBase::element_iterator elem_end = mesh.active_elements_end();

      // Count vertex degrees of freedom first
      for ( ; elem_it != elem_end; ++elem_it)
	{
	  Elem* elem  = *elem_it;
	  libmesh_assert (elem != NULL);

	  // Skip the numbering if this variable is
	  // not active on this element's subdoman
	  if (!var_description.active_on_subdomain(elem->subdomain_id()))
	    continue;

	  const ElemType type = elem->type();
          const unsigned int dim = elem->dim();

          FEType fe_type = base_fe_type;

	  fe_type.order = static_cast<Order>(fe_type.order +
                                             elem->p_level());

	  // Allocate the vertex DOFs
	  for (unsigned int n=0; n<elem->n_nodes(); n++)
	    {
	      Node* node = elem->get_node(n);

	      if (elem->is_vertex(n))
	        {
	          const unsigned int old_node_dofs =
	            node->n_comp(this->sys_number(), var);

		  const unsigned int vertex_dofs =
		    std::max(FEInterface::n_dofs_at_node(dim, fe_type,
                                                         type, n),
                             old_node_dofs);

		  // Some discontinuous FEs have no vertex dofs
		  if (vertex_dofs > old_node_dofs)
		    {
		      node->set_n_comp(this->sys_number(), var,
				       vertex_dofs);
		      // Abusing dof_number to set a "this is a
		      // vertex" flag
		      node->set_dof_number(this->sys_number(),
					   var, 0, vertex_dofs);
		    }
	        }
	    }
	} // done counting vertex dofs

      // count edge & face dofs next
      elem_it = mesh.active_elements_begin();
    }

  //------------------------------------------------------------
  // Finally, clear all the current DOF indices
  // (distribute_dofs expects them cleared!)
  this->invalidate_dofs(mesh);

  STOP_LOG("reinit()", "DofMapSubdiv");
}



void DofMapSubdiv::distribute_local_dofs_var_major(unsigned int &next_free_dof,
                                             MeshBase& mesh)
{
  const unsigned int sys_num = this->sys_number();
  const unsigned int n_vars  = this->n_variables();

  // We now only add remote dofs to the _send_list
  // unsigned int send_list_size = 0;

  // We will cache the first local index for each variable
  this->first_local_dof().clear();

  //-------------------------------------------------------------------------
  // First count and assign temporary numbers to local dofs
  for (unsigned var=0; var<n_vars; var++)
    {
      this->first_local_dof().push_back(next_free_dof);

      const Variable var_description = this->variable(var);

      MeshBase::element_iterator       elem_it  = mesh.active_local_elements_begin();
      const MeshBase::element_iterator elem_end = mesh.active_local_elements_end();

      for ( ; elem_it != elem_end; ++elem_it)
        {
          // Only number dofs connected to active
          // elements on this processor.
          Elem* elem  = *elem_it;

	  // ... and only variables which are active on
	  // on this element's subdomain
	  if (!var_description.active_on_subdomain(elem->subdomain_id()))
	    continue;

          const unsigned int n_nodes = elem->n_nodes();

          // First number the nodal DOFS
          for (unsigned int n=0; n<n_nodes; n++)
            {
              Node* node = elem->get_node(n);

              // assign dof numbers (all at once) if this is
              // our node and if they aren't already there
              if ((node->n_comp(sys_num,var) > 0) &&
                  (node->processor_id() == libMesh::processor_id()) &&
                  (node->dof_number(sys_num,var,0) ==
                   DofObject::invalid_id))
                {
                  node->set_dof_number(sys_num,
                                       var,
                                       0,
                                       next_free_dof);
                  next_free_dof += node->n_comp(sys_num,var);
                }
            }
        } // end loop on elements

      // we may have missed assigning DOFs to nodes that we own
      // but to which we have no connected elements matching our
      // variable restriction criterion.  this will happen, for example,
      // if variable V is restricted to subdomain S.  We may not own
      // any elements which live in S, but we may own nodes which are
      // *connected* to elements which do.  in this scenario these nodes
      // will presently have unnumbered DOFs. we need to take care of
      // them here since we own them and no other processor will touch them.
      {
	MeshBase::node_iterator       node_it  = mesh.local_nodes_begin();
	const MeshBase::node_iterator node_end = mesh.local_nodes_end();

	for (; node_it != node_end; ++node_it)
	  {
	    Node *node = *node_it;
	    libmesh_assert(node);

	    if (node->n_comp(sys_num,var))
	      if (node->dof_number(sys_num,var,0) == DofObject::invalid_id)
		{
		  node->set_dof_number (sys_num,
					var,
					0,
					next_free_dof);

		  next_free_dof += node->n_comp(sys_num,var);
		}
	  }
      }
    } // end loop on variables

  // Cache the last local dof number too
  this->first_local_dof().push_back(next_free_dof);
}



void DofMapSubdiv::dof_indices (const Elem* const elem,
			  std::vector<unsigned int>& di,
			  const unsigned int vn) const
{
  START_LOG("dof_indices()", "DofMapSubdiv");

  libmesh_assert (elem != NULL);
  libmesh_assert (elem->type() == TRI3SD);
  const Tri3SD* sd_elem = static_cast<const Tri3SD*>(elem);
  
  if (sd_elem->is_ghost())
    return;
  
  // Create the node patch for this element
  std::vector<Node *> nodes;
  libMesh::MeshTools::Subdiv::find_one_ring(sd_elem, nodes);

  const ElemType type        = elem->type();
  const unsigned int sys_num = this->sys_number();
  const unsigned int n_vars  = this->n_variables();
  const unsigned int dim     = elem->dim();

  // Clear the DOF indices vector
  di.clear();

  // Get the dof numbers
  for (unsigned int v=0; v<n_vars; v++)
    if ((v == vn) || (vn == libMesh::invalid_uint))
    {
      if (this->variable(v).active_on_subdomain(elem->subdomain_id()))
	{ // Do this for all the variables if one was not specified
	  // or just for the specified variable

	  // Increase the polynomial order on p refined elements
	  FEType fe_type = this->variable_type(v);

	  // Get the node-based DOF numbers
	  for (unsigned int n=0; n<nodes.size(); n++)
	    {
	      const Node* node      = nodes[n];

	      const unsigned int nc = FEInterface::n_dofs_at_node (dim,
								   fe_type,
								   type,
								   n);

		for (unsigned int i=0; i<nc; i++)
		  {
		    libmesh_assert (node->dof_number(sys_num,v,i) !=
				    DofObject::invalid_id);
		    di.push_back(node->dof_number(sys_num,v,i));
		  }
	    }
	}
      } // end loop over variables

  STOP_LOG("dof_indices()", "DofMapSubdiv");
}

} // namespace libMesh
