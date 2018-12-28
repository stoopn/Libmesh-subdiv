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

#include "fe_subdiv_map.h"
#include "mesh_subdiv_support.h"

void FESubdivMap::resize_quadrature_map_vectors(const unsigned int dim, unsigned int n_qp)
{
    FEMap::resize_quadrature_map_vectors(dim, n_qp);
    d3xyzdxi3_map.resize(n_qp);
    d3xyzdeta3_map.resize(n_qp);
    d3xyzdxi2eta_map.resize(n_qp);
    d3xyzdeta2xi_map.resize(n_qp);
}

void FESubdivMap::compute_map( const unsigned int dim,
				const std::vector<Real>& qw,
				const Elem* elem )
{
   // Start logging the map computation.
  START_LOG("compute_map()", "FESubdivMap");

	libmesh_assert (dim == 2);
	libmesh_assert (elem != NULL);
	libmesh_assert (elem->type() == TRI3SD);
	const Tri3SD* sd_elem = static_cast<const Tri3SD*>(elem);

  // The number of quadrature points.
  const unsigned int n_qp = qw.size();

 	// create the 1-ring around element elem
	std::vector<Node *> nodes;
	MeshTools::Subdiv::find_one_ring(sd_elem, nodes);

	this->resize_quadrature_map_vectors (dim, n_qp);

	for (unsigned int p = 0; p < n_qp; ++p)
	{
		// Finally, build the local manifold metric and other things
		xyz[p].zero();
		dxyzdxi_map[p].zero();
		dxyzdeta_map[p].zero();
		d2xyzdxi2_map[p].zero();
		d2xyzdeta2_map[p].zero();
		d2xyzdxideta_map[p].zero();
        
        d3xyzdxi3_map[p].zero();
        d3xyzdeta3_map[p].zero();
        d3xyzdxi2eta_map[p].zero();
        d3xyzdeta2xi_map[p].zero();
        
		for (unsigned int j = 0; j < phi_map.size(); ++j)
		{
			xyz[p].add_scaled             (*nodes[j], phi_map[j][p]);
			dxyzdxi_map[p].add_scaled     (*nodes[j], dphidxi_map[j][p]);
			dxyzdeta_map[p].add_scaled    (*nodes[j], dphideta_map[j][p]);
			d2xyzdxi2_map[p].add_scaled   (*nodes[j], d2phidxi2_map[j][p]);
			d2xyzdeta2_map[p].add_scaled  (*nodes[j], d2phideta2_map[j][p]);
			d2xyzdxideta_map[p].add_scaled(*nodes[j], d2phidxideta_map[j][p]);
            
            d3xyzdxi3_map[p].add_scaled   (*nodes[j], d3phidxi3_map[j][p]);
            d3xyzdeta3_map[p].add_scaled  (*nodes[j], d3phideta3_map[j][p]);
            d3xyzdxi2eta_map[p].add_scaled(*nodes[j], d3phidxi2eta_map[j][p]);
            d3xyzdeta2xi_map[p].add_scaled(*nodes[j], d3phideta2xi_map[j][p]);
		}

		JxW[p] = dxyzdxi_map[p].cross(dxyzdeta_map[p]).size() * qw[p];
	}

  // Stop logging the map computation.
  STOP_LOG("compute_map()", "FESubdivMap");

  return;
}
