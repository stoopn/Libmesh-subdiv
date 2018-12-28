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



#ifndef __dof_map_subdiv_h__
#define __dof_map_subdiv_h__

// Local Includes -----------------------------------
#include "dof_map.h"

// C++ Includes   -----------------------------------

namespace libMesh
{

// ------------------------------------------------------------
// DofMap class definition for subdivision surfaces

/**
 * This class handles the numbering of degrees of freedom on a
 * subdivision surface mesh.
 */
class DofMapSubdiv : public DofMap
{
public:

  /**
   * Constructor.  Requires the number of the system for which we
   * will be numbering degrees of freedom.
   */
  DofMapSubdiv(const unsigned int sys_number) : DofMap(sys_number) {};

  /**
   * Destructor.
   */
  virtual ~DofMapSubdiv() {};

  /**
   * Distributes the global degrees of freedom, for dofs on
   * this processor.  In this format the local
   * degrees of freedom are in a contiguous block for each
   * variable in the system.
   * Starts at index next_free_dof, and increments it to
   * the post-final index.
   */
  virtual void distribute_local_dofs_var_major (unsigned int& next_free_dof,
						MeshBase& mesh);

  /**
   * Fills the vector \p di with the global degree of freedom indices
   * for the element. If no variable number is specified then all
   * variables are returned.
   */
  virtual void dof_indices (const Elem* const elem,
			    std::vector<unsigned int>& di,
			    const unsigned int vn = libMesh::invalid_uint) const;

  /**
   * Reinitialize the underlying data strucures conformal to the current mesh.
   */
  virtual void reinit (MeshBase& mesh);
};

} // namespace libMesh

#endif // __dof_map_subdiv_h__
