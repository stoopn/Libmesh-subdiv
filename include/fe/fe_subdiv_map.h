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



#ifndef __fe_subdiv_map_h__
#define __fe_subdiv_map_h__

#include "fe_map.h"
#include "libmesh_logging.h"
#include "elem.h"

namespace libMesh
{

  class FESubdivMap : public FEMap
  {
  public:
    
    FESubdivMap()
      : FEMap(){};

    virtual ~FESubdivMap(){};

    /**
     * Special implementation for subdivision surface finite elements
     */
    virtual void compute_map( const unsigned int dim,
			      const std::vector<Real>& qw,
			      const Elem* elem );
      
      /**
       * @returns the reference to physical map 3rd derivative
       */
      std::vector<std::vector<Real> >& get_d3phidxi3_map()
      { return d3phidxi3_map; }
      std::vector<std::vector<Real> >& get_d3phideta3_map()
      { return d3phideta3_map; }
      std::vector<std::vector<Real> >& get_d3phidxi2eta_map()
      { return d3phidxi2eta_map; }
      std::vector<std::vector<Real> >& get_d3phideta2xi_map()
      { return d3phideta2xi_map; }
      /**
       * @returns the third partial derivatives.
       */
      const std::vector<RealGradient>& get_d3xyzdxi3() const
      { return d3xyzdxi3_map; }
      const std::vector<RealGradient>& get_d3xyzdeta3() const
      { return d3xyzdeta3_map; }
      const std::vector<RealGradient>& get_d3xyzdxi2eta() const
      { return d3xyzdxi2eta_map; }
      const std::vector<RealGradient>& get_d3xyzdeta2xi() const
      { return d3xyzdeta2xi_map; }

  protected:
      void resize_quadrature_map_vectors(const unsigned int dim, unsigned int n_qp);
      /**
       * Map for the third derivatives.
       */
      std::vector<std::vector<Real> >   d3phidxi3_map;
      std::vector<std::vector<Real> >   d3phideta3_map;
      std::vector<std::vector<Real> >   d3phidxi2eta_map;
      std::vector<std::vector<Real> >   d3phideta2xi_map;
      /**
       * Vector of second partial derivatives in xi:
       * d^2(x)/d(xi)^2, d^2(y)/d(xi)^2, d^2(z)/d(xi)^2
       */
      std::vector<RealGradient> d3xyzdxi3_map;
      std::vector<RealGradient> d3xyzdeta3_map;
      std::vector<RealGradient> d3xyzdxi2eta_map;
      std::vector<RealGradient> d3xyzdeta2xi_map;
      
      
      
      
  }; // class FESubdivMap
} // namespace libMesh
#endif //__fe_subdiv_map_h__
