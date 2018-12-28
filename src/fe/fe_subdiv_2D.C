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



// Local includes
#include "fe.h"
#include "libmesh_logging.h"
#include "fe_type.h"
#include "quadrature.h"
#include "face_tri3_sd.h"
#include "fe_macro.h"
#include "dense_matrix.h"
#include "utility.h"
#include "fe_subdiv_map.h"


namespace libMesh
{

FESubdiv::FESubdiv(const FEType& fet) :
  FE<2,SUBDIV>(fet)
{
	// Only 2D meshes in 3D space are supported
	libmesh_assert(LIBMESH_DIM == 3);
}



void FESubdiv::init_subdiv_matrix(DenseMatrix<Real> &A,
			  unsigned int valence)
{
	A.resize(valence + 12, valence + 12);

	// A = (S11 0; S21 S22), see Cirak et al.,
	// Int. J. Numer. Meth. Engng. 2000; 47:2039-2072, Appendix A.2.

	// First, set the static S21 part
	A(valence+ 1,0        ) = 0.125;
	A(valence+ 1,1        ) = 0.375;
	A(valence+ 1,valence  ) = 0.375;
	A(valence+ 2,0        ) = 0.0625;
	A(valence+ 2,1        ) = 0.625;
	A(valence+ 2,2        ) = 0.0625;
	A(valence+ 2,valence  ) = 0.0625;
	A(valence+ 3,0        ) = 0.125;
	A(valence+ 3,1        ) = 0.375;
	A(valence+ 3,2        ) = 0.375;
	A(valence+ 4,0        ) = 0.0625;
	A(valence+ 4,1        ) = 0.0625;
	A(valence+ 4,valence-1) = 0.0625;
	A(valence+ 4,valence  ) = 0.625;
	A(valence+ 5,0        ) = 0.125;
	A(valence+ 5,valence-1) = 0.375;
	A(valence+ 5,valence  ) = 0.375;
	A(valence+ 6,1        ) = 0.375;
	A(valence+ 6,valence  ) = 0.125;
	A(valence+ 7,1        ) = 0.375;
	A(valence+ 8,1        ) = 0.375;
	A(valence+ 8,2        ) = 0.125;
	A(valence+ 9,1        ) = 0.125;
	A(valence+ 9,valence  ) = 0.375;
	A(valence+10,valence  ) = 0.375;
	A(valence+11,valence-1) = 0.125;
	A(valence+11,valence  ) = 0.375;

	// Next, set the static S22 part
	A(valence+ 1,valence+1) = 0.125;
	A(valence+ 2,valence+1) = 0.0625;
	A(valence+ 2,valence+2) = 0.0625;
	A(valence+ 2,valence+3) = 0.0625;
	A(valence+ 3,valence+3) = 0.125;
	A(valence+ 4,valence+1) = 0.0625;
	A(valence+ 4,valence+4) = 0.0625;
	A(valence+ 4,valence+5) = 0.0625;
	A(valence+ 5,valence+5) = 0.125;
	A(valence+ 6,valence+1) = 0.375;
	A(valence+ 6,valence+2) = 0.125;
	A(valence+ 7,valence+1) = 0.125;
	A(valence+ 7,valence+2) = 0.375;
	A(valence+ 7,valence+3) = 0.125;
	A(valence+ 8,valence+2) = 0.125;
	A(valence+ 8,valence+3) = 0.375;
	A(valence+ 9,valence+1) = 0.375;
	A(valence+ 9,valence+4) = 0.125;
	A(valence+10,valence+1) = 0.125;
	A(valence+10,valence+4) = 0.375;
	A(valence+10,valence+5) = 0.125;
	A(valence+11,valence+4) = 0.125;
	A(valence+11,valence+5) = 0.375;

	// Last, set the S11 part: first row
	std::vector<Real> weights;
	get_limit_mask(weights, valence);
	for (unsigned int i = 0; i <= valence; ++i)
		A(0,i) = weights[i];

	// second row
	A(1,0) = 0.375;
	A(1,1) = 0.375;
	A(1,2) = 0.125;
	A(1,valence) = 0.125;

	// third to second-to-last rows
	for (unsigned int i = 2; i < valence; ++i)
	{
		A(i,0  ) = 0.375;
		A(i,i-1) = 0.125;
		A(i,i  ) = 0.375;
		A(i,i+1) = 0.125;
	}

	// last row
	A(valence,0) = 0.375;
	A(valence,1) = 0.125;
	A(valence,valence-1) = 0.125;
	A(valence,valence  ) = 0.375;
}



Real FESubdiv::regular_shape(const unsigned int i,
			  const Real v,
			  const Real w)
{
	// These are the 12 quartic box splines, see Cirak et al.,
	// Int. J. Numer. Meth. Engng. 2000; 47:2039-2072, Appendix A.1.

	const Real u = 1 - v - w;
	libmesh_assert(0 <= v);
	libmesh_assert(0 <= w);
	libmesh_assert(0 <= u);

	using Utility::pow;
	const Real factor = 1. / 12;

	switch (i)
	{
		case 0:
			return factor*(pow<4>(u) + 2*u*u*u*v);
		case 1:
			return factor*(pow<4>(u) + 2*u*u*u*w);
		case 2:
			return factor*(pow<4>(u) + 2*u*u*u*w + 6*u*u*u*v + 6*u*u*v*w + 12*u*u*v*v + 6*u*v*v*w + 6*u*v*v*v +
				2*v*v*v*w + pow<4>(v));
		case 3:
			return factor*(6*pow<4>(u) + 24*u*u*u*w + 24*u*u*w*w + 8*u*w*w*w + pow<4>(w) + 24*u*u*u*v + 
				60*u*u*v*w + 36*u*v*w*w + 6*v*w*w*w + 24*u*u*v*v + 36*u*v*v*w + 12*v*v*w*w + 8*u*v*v*v +
				6*v*v*v*w + pow<4>(v));
		case 4:
			return factor*(pow<4>(u) + 6*u*u*u*w + 12*u*u*w*w + 6*u*w*w*w + pow<4>(w) + 2*u*u*u*v + 6*u*u*v*w +
				6*u*v*w*w + 2*v*w*w*w);
		case 5:
			return factor*(2*u*v*v*v + pow<4>(v));
		case 6:
			return factor*(pow<4>(u) + 6*u*u*u*w + 12*u*u*w*w + 6*u*w*w*w + pow<4>(w) + 8*u*u*u*v + 36*u*u*v*w +
				36*u*v*w*w + 8*v*w*w*w + 24*u*u*v*v + 60*u*v*v*w + 24*v*v*w*w + 24*u*v*v*v + 24*v*v*v*w + 6*pow<4>(v));
		case 7:
			return factor*(pow<4>(u) + 8*u*u*u*w + 24*u*u*w*w + 24*u*w*w*w + 6*pow<4>(w) + 6*u*u*u*v + 36*u*u*v*w +
				60*u*v*w*w + 24*v*w*w*w + 12*u*u*v*v + 36*u*v*v*w + 24*v*v*w*w + 6*u*v*v*v + 8*v*v*v*w + pow<4>(v));
		case 8:
			return factor*(2*u*w*w*w + pow<4>(w));
		case 9:
			return factor*(2*v*v*v*w + pow<4>(v));
		case 10:
			return factor*(2*u*w*w*w + pow<4>(w) + 6*u*v*w*w + 6*v*w*w*w + 6*u*v*v*w + 12*v*v*w*w + 2*u*v*v*v +
				6*v*v*v*w + pow<4>(v));
		case 11:
			return factor*(pow<4>(w) + 2*v*w*w*w);

		default:
			libmesh_error();
	}

	libmesh_error();
	return 0.; 
}



Real FESubdiv::regular_shape_deriv(const unsigned int i,
			  const unsigned int j,
			  const Real v,
			  const Real w)
{
	const Real u = 1 - v - w;
	const Real factor = 1. / 12;

	switch (j) // j=0: xi-directional derivative, j=1: eta-directional derivative
	{
		case 0: // xi derivatives
		{	      
			switch (i) // shape function number
			{
				case 0:
					return factor*(-6*v*u*u - 2*u*u*u);
				case 1:
					return factor*(-4*u*u*u - 6*u*u*w);
				case 2:
					return factor*(-2*v*v*v - 6*v*v*u + 6*v*u*u + 2*u*u*u);
				case 3:
					return factor*(-4*v*v*v - 24*v*v*u - 24*v*u*u - 18*v*v*w - 48*v*u*w - 12*u*u*w -
						12*v*w*w - 12*u*w*w - 2*w*w*w);
				case 4:
					return factor*(-6*v*u*u - 2*u*u*u - 12*v*u*w-12*u*u*w - 6*v*w*w - 18*u*w*w - 4*w*w*w);
				case 5:
					return factor*(2*v*v*v + 6*v*v*u);
				case 6:
					return factor*(24*v*v*u + 24*v*u*u + 4*u*u*u + 12*v*v*w + 48*v*u*w + 18*u*u*w +
						12*v*w*w + 12*u*w*w + 2*w*w*w);
				case 7:
					return factor*(-2*v*v*v - 6*v*v*u + 6*v*u*u + 2*u*u*u - 12*v*v*w + 12*u*u*w -
						12*v*w*w + 12*u*w*w);
				case 8:
					return -w*w*w/6;
				case 9:
					return factor*(4*v*v*v + 6*v*v*w);
				case 10:
					return factor*(2*v*v*v + 6*v*v*u + 12*v*v*w + 12*v*u*w + 18*v*w*w + 6*u*w*w + 4*w*w*w);
				case 11:
					return w*w*w/6;
				default:
					libmesh_error();
			}
		}
		case 1: // eta derivatives
		{
			switch (i) // shape function number
			{
				case 0:
					return factor*(-6*v*u*u - 4*u*u*u);
				case 1:
					return factor*(-2*u*u*u - 6*u*u*w);
				case 2:
					return factor*(-4*v*v*v - 18*v*v*u - 12*v*u*u - 2*u*u*u - 6*v*v*w - 12*v*u*w -
						6*u*u*w);
				case 3:
					return factor*(-2*v*v*v-12*v*v*u - 12*v*u*u - 12*v*v*w - 48*v*u*w - 24*u*u*w -
						18*v*w*w - 24*u*w*w - 4*w*w*w);
				case 4:
					return factor*(2*u*u*u + 6*u*u*w - 6*u*w*w - 2*w*w*w);
				case 5:
					return -v*v*v/6;
				case 6:
					return factor*(12*v*v*u + 12*v*u*u + 2*u*u*u - 12*v*v*w + 6*u*u*w - 12*v*w*w -
						6*u*w*w - 2*w*w*w);
				case 7:
					return factor*(2*v*v*v + 12*v*v*u + 18*v*u*u + 4*u*u*u + 12*v*v*w + 48*v*u*w +
						24*u*u*w + 12*v*w*w + 24*u*w*w);
				case 8:
					return factor*(6*u*w*w + 2*w*w*w);
				case 9:
					return v*v*v/6;
				case 10:
					return factor*(4*v*v*v + 6*v*v*u + 18*v*v*w + 12*v*u*w + 12*v*w*w + 6*u*w*w +
						2*w*w*w);
				case 11:
					return factor*(6*v*w*w + 4*w*w*w);
				default: 
					libmesh_error();
			}
		}
		default:
			libmesh_error();
	}

	libmesh_error();
	return 0.;
}



Real FESubdiv::regular_shape_second_deriv(const unsigned int i,
			  const unsigned int j,
			  const Real v,
			  const Real w)
{
	const Real u = 1 - v - w;
	const Real factor = 1. / 12;

	switch (j)
	{
		case 0: // xi-xi derivative
		{
			switch (i) // shape function number
			{
				case 0:
					return v*u;
				case 1:
					return u*u + u*w;
				case 2:
					return -2*v*u;
				case 3:
					return v*v - 2*u*u + v*w - 2*u*w;
				case 4:
					return v*u + v*w + u*w + w*w;
				case 5:
					return v*u;
				case 6:
					return factor*(-24*v*v + 12*u*u - 24*v*w + 12*u*w);
				case 7:
					return -2*v*u - 2*v*w - 2*u*w - 2*w*w;
				case 8:
					return 0.;
				case 9:
					return v*v + v*w;
				case 10:
					return v*u + v*w + u*w + w*w;
				case 11:
					return 0.;
				default:
					libmesh_error();
			}
		}
		case 1: //eta-xi derivative
		{
			switch (i)
			{
				case 0:
					return factor*(12*v*u + 6*u*u);
				case 1:
					return factor*(6*u*u + 12*u*w);
				case 2:
					return factor*(6*v*v - 12*v*u - 6*u*u);
				case 3:
					return factor*(6*v*v - 12*u*u + 24*v*w + 6*w*w);
				case 4:
					return factor*(-6*u*u - 12*u*w + 6*w*w);
				case 5:
					return -v*v/2.;
				case 6:
					return factor*(-12*v*v + 6*u*u - 24*v*w - 12*u*w - 6*w*w);
				case 7:
					return factor*(-6*v*v - 12*v*u + 6*u*u - 24*v*w - 12*w*w);
				case 8:
					return -w*w/2.;
				case 9:
					return v*v/2.;
				case 10:
					return factor*(6*v*v + 12*v*u + 24*v*w + 12*u*w + 6*w*w);			
				case 11:
					return w*w/2.;
				default:
					libmesh_error();
			}
		}
		case 2: // eta-eta derivative
		{
		    switch (i)
			{
				case 0:
					return v*u + u*u;
				case 1:
					return u*w;
				case 2:
					return v*v + v*u + v*w + u*w;
				case 3:
					return -2*v*u - 2*u*u + v*w + w*w;
				case 4:
					return -2*u*w;
				case 5:
					return 0.;
				case 6:
					return -2*v*v - 2*v*u - 2*v*w - 2*u*w;
				case 7:
					return v*u + u*u - 2*v*w - 2*w*w;
				case 8:
					return u*w;
				case 9:
					return 0.;
				case 10:
					return v*v + v*u + v*w + u*w;
				case 11:
					return v*w + w*w;
				default:
					libmesh_error();
			}
		} 
		default:
			libmesh_error();
	}

	libmesh_error();
	return 0.;
}


    
    
Real FESubdiv::regular_shape_third_deriv(const unsigned int i,
                                              const unsigned int j,
                                              const Real v,
                                              const Real w)
{
        const Real u = 1 - v - w;
        const Real factor = 1. / 12;
        
        switch (j)
        {
            case 0: // xi-xi-xi derivative   xi=v, eta=w
            {
                switch (i) // shape function number
                {
                    case 0:
                        return 1.-2*v-w;
                    case 1:
                        return -2+ 2*v + w;
                    case 2:
                        return 2*v -2*u;
                    case 3:
                        return 4-2*v-w;
                    case 4:
                        return 1 - 2*v - w;
                    case 5:
                        return 1 - 2*v - w;
                    case 6:
                        return -2*v - 2 - w;
                    case 7:
                        return -2 + 4*v + 2*w;
                    case 8:
                        return 0.;
                    case 9:
                        return 2*v + w;
                    case 10:
                        return 1-2*v - w;
                    case 11:
                        return 0.;
                    default:
                        libmesh_error();
                }
            }
        case 1: //eta-xi-xi derivative (eta=w, xi=v)
            {
                switch (i)
                {
                    case 0:
                        return -v;
                    case 1:
                        return v-1;
                    case 2:
                        return 2*v;
                    case 3:
                        return 2-v;
                    case 4:
                        return 1 - v;
                    case 5:
                        return -v;
                    case 6:
                        return -v - 1;
                    case 7:
                        return 2*v - 2;
                    case 8:
                        return 0;
                    case 9:
                        return v;
                    case 10:
                        return -v + 1;
                    case 11:
                        return 0;
                    default:
                        libmesh_error();
                }
            }
                
        case 2: //eta-xi-eta derivative (eta=w, xi=v)
            {
                switch (i)
                {
                    case 0:
                        return -1 + w;
                    case 1:
                        return -w;
                    case 2:
                        return 1 - w;
                    case 3:
                        return 2 - w;
                    case 4:
                        return 2*w;
                    case 5:
                        return 0;
                    case 6:
                        return -2 + 2*w;
                    case 7:
                        return -1 - w;
                    case 8:
                        return -w;
                    case 9:
                        return 0;
                    case 10:
                        return -w + 1;
                    case 11:
                        return w;
                    default:
                        libmesh_error();
                }
            }
        case 3: // eta-eta-eta derivative (eta=w)
            {
                switch (i)
                {
                    case 0:
                        return v - 2 + 2*w;
                    case 1:
                        return -2*w + 1 - v;
                    case 2:
                        return -2*w + 1 - v;
                    case 3:
                        return -v + 4 - 2*w;
                    case 4:
                        return 4*w - 2 + 2*v;
                    case 5:
                        return 0.;
                    case 6:
                        return 4*w - 2 + 2*v;
                    case 7:
                        return -v - 2 - 2*w;
                    case 8:
                        return -2*w + 1 - v;
                    case 9:
                        return 0.;
                    case 10:
                        return -2*w + 1 - v;
                    case 11:
                        return v+ 2*w;
                    default:
                        libmesh_error();
                }
            } 
            default:
                libmesh_error();
        }
        
        libmesh_error();
        return 0.;
}

    

void FESubdiv::get_limit_mask(std::vector<Real> & weights,
			  const unsigned int valence)
{
	libmesh_assert(valence > 0);
	const Real cs = std::cos(2 * libMesh::pi / valence);
	const Real nb_weight = (0.625 - Utility::pow<2>(0.375 + 0.25 * cs)) / valence;
	weights.resize(1 + valence, nb_weight);
	weights[0] = 1.0 - valence * nb_weight;
}



void FESubdiv::init_shape_functions(const std::vector<Point> &qp,
			  const Elem *e)
{
	libmesh_assert(e != NULL);
	libmesh_assert(e->type() == TRI3SD);
	const Tri3SD* sd_elem = static_cast<const Tri3SD*>(e);

	START_LOG("init_shape_functions()", "FESubdiv");

	calculations_started = true;

	// If the user forgot to request anything, we'll be safe and calculate everything:
	if (!calculate_phi && !calculate_dphi && !calculate_d2phi)
		calculate_phi = calculate_dphi = calculate_d2phi = true;

	const unsigned int valence = sd_elem->get_ordered_valence(0);
	const unsigned int n_qp = qp.size();
	const unsigned int n_approx_shape_functions = valence + 6;

	// resize the vectors to hold current data
	phi.resize         (n_approx_shape_functions);
	dphi.resize        (n_approx_shape_functions);
	dphidxi.resize     (n_approx_shape_functions);
	dphideta.resize    (n_approx_shape_functions);
	d2phi.resize       (n_approx_shape_functions);
	d2phidxi2.resize   (n_approx_shape_functions);
	d2phidxideta.resize(n_approx_shape_functions);
	d2phideta2.resize  (n_approx_shape_functions);
    d3phidxi3.resize   (n_approx_shape_functions);
    d3phideta3.resize  (n_approx_shape_functions);
    d3phidxi2eta.resize(n_approx_shape_functions);
    d3phideta2xi.resize(n_approx_shape_functions);

	for (unsigned int i = 0; i < n_approx_shape_functions; ++i)
	{
		phi[i].resize         (n_qp);
		dphi[i].resize        (n_qp);
		dphidxi[i].resize     (n_qp);
		dphideta[i].resize    (n_qp);
		d2phi[i].resize       (n_qp);
		d2phidxi2[i].resize   (n_qp);
		d2phidxideta[i].resize(n_qp);
		d2phideta2[i].resize  (n_qp);
        
        d3phidxi3[i].resize   (n_qp);
        d3phideta3[i].resize  (n_qp);
        d3phidxi2eta[i].resize(n_qp);
        d3phideta2xi[i].resize(n_qp);
	}

	// Renumbering of the shape functions
	static const unsigned int cvi[12] = {3,6,2,0,1,4,7,10,9,5,11,8};

	if (valence == 6) // This means that all vertices are regular, i.e. we have 12 shape functions
	{
		for (unsigned int i = 0; i < n_approx_shape_functions; ++i)
		{
			for (unsigned int p = 0; p < n_qp; ++p)
			{
				phi[i][p]          = FE<2,SUBDIV>::shape             (e, fe_type.order, cvi[i],    qp[p]);
				dphidxi[i][p]      = FE<2,SUBDIV>::shape_deriv       (e, fe_type.order, cvi[i], 0, qp[p]);
				dphideta[i][p]     = FE<2,SUBDIV>::shape_deriv       (e, fe_type.order, cvi[i], 1, qp[p]);
				dphi[i][p](0)      = dphidxi[i][p];
				dphi[i][p](1)      = dphideta[i][p];
				d2phidxi2[i][p]    = FE<2,SUBDIV>::shape_second_deriv(e, fe_type.order, cvi[i], 0, qp[p]);
				d2phidxideta[i][p] = FE<2,SUBDIV>::shape_second_deriv(e, fe_type.order, cvi[i], 1, qp[p]);
				d2phideta2[i][p]   = FE<2,SUBDIV>::shape_second_deriv(e, fe_type.order, cvi[i], 2, qp[p]);
				d2phi[i][p](0,0)   = d2phidxi2[i][p];
				d2phi[i][p](0,1)   = d2phi[i][p](1,0) = d2phidxideta[i][p];
				d2phi[i][p](1,1)   = d2phideta2[i][p];
                
                d3phidxi3[i][p]    = shape_third_deriv(e, fe_type.order, cvi[i], 0, qp[p]);
                d3phideta3[i][p]   = shape_third_deriv(e, fe_type.order, cvi[i], 3, qp[p]);
                d3phidxi2eta[i][p] = shape_third_deriv(e, fe_type.order, cvi[i], 1, qp[p]);
                d3phideta2xi[i][p] = shape_third_deriv(e, fe_type.order, cvi[i], 2, qp[p]);
			}
		}
	}
	else // vertex 0 is irregular by construction of the mesh
	{
		static const Real eps = 1e-10;

		// temporary values
		std::vector<Real> tphi(12);
		std::vector<Real> tdphidxi(12);
		std::vector<Real> tdphideta(12);   
		std::vector<Real> td2phidxi2(12);
		std::vector<Real> td2phidxideta(12);
		std::vector<Real> td2phideta2(12);
        
        std::vector<Real> td3phidxi3(12);
        std::vector<Real> td3phidxi2eta(12);
        std::vector<Real> td3phideta2xi(12);
        std::vector<Real> td3phideta3(12);
        
		for (unsigned int p = 0; p < n_qp; ++p)
		{
			// evaluate the number of the required subdivisions
			Real v = qp[p](0);
			Real w = qp[p](1);
			Real u = 1 - v - w;
			Real min = 0, max = 0.5;
			int n = 0;
			while (!(u > min-eps && u < max+eps))
			{
				++n;
				min = max;
				max += std::pow((Real)(2), -n-1);
			}

			// transform u, v and w according to the number of subdivisions required.
			const Real pow2 = std::pow((Real)(2), n);
			v *= pow2;
			w *= pow2;
			u = 1 - v - w;
			libmesh_assert(u < 0.5 + eps && u > -eps);

			// find out in which subdivided patch we are and setup the "selection matrix" P and the transformation Jacobian
			// (see Int. J. Numer. Meth. Engng. 2000; 47:2039-2072, Appendix A.2.)
			const int k = n+1;
			Real jfac; // the additional factor per derivative order
			DenseMatrix<Real> P(12, valence+12);
			if (v > 0.5 - eps)
			{
				v = 2*v - 1;
				w = 2*w;
				jfac = std::pow((Real)(2), k);
				P( 0,2        ) = 1;
				P( 1,0        ) = 1;
				P( 2,valence+3) = 1;
				P( 3,1        ) = 1;
				P( 4,valence  ) = 1;
				P( 5,valence+8) = 1;
				P( 6,valence+2) = 1;
				P( 7,valence+1) = 1;
				P( 8,valence+4) = 1;
				P( 9,valence+7) = 1;
				P(10,valence+6) = 1;
				P(11,valence+9) = 1;
			}
			else if (w > 0.5 - eps)
			{
				v = 2*v;
				w = 2*w - 1;
				jfac = std::pow((Real)(2), k);
				P( 0,0         ) = 1;
				P( 1,valence- 1) = 1;
				P( 2,1         ) = 1;
				P( 3,valence   ) = 1;
				P( 4,valence+ 5) = 1;
				P( 5,valence+ 2) = 1;
				P( 6,valence+ 1) = 1;
				P( 7,valence+ 4) = 1;
				P( 8,valence+11) = 1;
				P( 9,valence+ 6) = 1;
				P(10,valence+ 9) = 1;
				P(11,valence+10) = 1;
			}
			else
			{
				v = 1 - 2*v;
				w = 1 - 2*w;
				jfac = std::pow((Real)(-2), k);
				P( 0,valence+9) = 1;
				P( 1,valence+6) = 1;
				P( 2,valence+4) = 1;
				P( 3,valence+1) = 1;
				P( 4,valence+2) = 1;
				P( 5,valence+5) = 1;
				P( 6,valence  ) = 1;
				P( 7,1        ) = 1;
				P( 8,valence+3) = 1;
				P( 9,valence-1) = 1;
				P(10,0        ) = 1;
				P(11,2        ) = 1;
			} 

			u = 1 - v - w;
			if ((u > 1 + eps) || (u < -eps))
			{
				std::cout << "SUBDIV irregular patch: u is outside valid range!\n";
				libmesh_error();
			}

			DenseMatrix<Real> A;
			init_subdiv_matrix(A, valence);

			// compute P*A^k
			if (k > 1)
			{
				DenseMatrix<Real> Acopy(A);
				for (int e = 1; e < k; ++e)
				A.right_multiply(Acopy);
			}    
			P.right_multiply(A);

			const Point transformed_p(v,w);

			for (unsigned int i = 0; i < 12; ++i)
			{
				tphi[i]          = FE<2,SUBDIV>::shape             (e, fe_type.order, i,    transformed_p);
				tdphidxi[i]      = FE<2,SUBDIV>::shape_deriv       (e, fe_type.order, i, 0, transformed_p);
				tdphideta[i]     = FE<2,SUBDIV>::shape_deriv       (e, fe_type.order, i, 1, transformed_p);
				td2phidxi2[i]    = FE<2,SUBDIV>::shape_second_deriv(e, fe_type.order, i, 0, transformed_p);
				td2phidxideta[i] = FE<2,SUBDIV>::shape_second_deriv(e, fe_type.order, i, 1, transformed_p);
				td2phideta2[i]   = FE<2,SUBDIV>::shape_second_deriv(e, fe_type.order, i, 2, transformed_p);
                
                td3phidxi3[i]    = shape_third_deriv(e, fe_type.order, i, 0, transformed_p);
                td3phideta3[i]   = shape_third_deriv(e, fe_type.order, i, 3, transformed_p);
                td3phidxi2eta[i] = shape_third_deriv(e, fe_type.order, i, 1, transformed_p);
                td3phideta2xi[i] = shape_third_deriv(e, fe_type.order, i, 2, transformed_p);
			}

			// Finally, we can compute the irregular shape functions as the product of P
			// and the regular shape functions:
			Real sum1, sum2, sum3, sum4, sum5, sum6,  sum7, sum8, sum9, sum10;
			for (unsigned int j = 0; j < n_approx_shape_functions; ++j)
			{
				sum1 = sum2 = sum3 = sum4 = sum5 = sum6 = sum7 = sum8 = sum9 = sum10 = 0;
				for (unsigned int i = 0; i < 12; ++i)
				{
					sum1 += P(i,j) * tphi[i];
					sum2 += P(i,j) * tdphidxi[i];
					sum3 += P(i,j) * tdphideta[i];
					sum4 += P(i,j) * td2phidxi2[i];
					sum5 += P(i,j) * td2phidxideta[i];
					sum6 += P(i,j) * td2phideta2[i];
                    
                    sum7  += P(i,j) * td3phidxi3[i];
                    sum8  += P(i,j) * td3phideta3[i];
                    sum9  += P(i,j) * td3phidxi2eta[i];
                    sum10 += P(i,j) * td3phideta2xi[i];
				}
				phi[j][p]          = sum1;
				dphidxi[j][p]      = sum2 * jfac;
				dphideta[j][p]     = sum3 * jfac;
				dphi[j][p](0)      = dphidxi[j][p];
				dphi[j][p](1)      = dphideta[j][p];
				d2phidxi2[j][p]    = sum4 * jfac * jfac;
				d2phidxideta[j][p] = sum5 * jfac * jfac;
				d2phideta2[j][p]   = sum6 * jfac * jfac;
				d2phi[j][p](0,0)   = d2phidxi2[j][p];
				d2phi[j][p](0,1)   = d2phi[j][p](1,0) = d2phidxideta[j][p];
				d2phi[j][p](1,1)   = d2phideta2[j][p];
                
                d3phidxi3[j][p]    = sum7 * jfac * jfac * jfac;
                d3phideta3[j][p]   = sum8 * jfac * jfac * jfac;
                d3phidxi2eta[j][p] = sum9 * jfac * jfac * jfac;
                d3phideta2xi[j][p] = sum10 * jfac * jfac * jfac;

                
			}
		} // end quadrature loop
	} // end irregular vertex

  // Let the FESubdivMap use the same initialized shape functions
	this->_fe_map->get_phi_map()          = phi;
	this->_fe_map->get_dphidxi_map()      = dphidxi;
	this->_fe_map->get_dphideta_map()     = dphideta;
	this->_fe_map->get_d2phidxi2_map()    = d2phidxi2;
	this->_fe_map->get_d2phideta2_map()   = d2phideta2;
	this->_fe_map->get_d2phidxideta_map() = d2phidxideta;
    
    FESubdivMap* fesd_map = dynamic_cast<FESubdivMap*>(this->_fe_map.get());
  
//    AutoPtr<FESubdivMap> fesd_map = dynamic_cast< AutoPtr<FESubdivMap> >(this->_fe_map);
    fesd_map->get_d3phidxi3_map()    = d3phidxi3;
    fesd_map->get_d3phideta3_map()   = d3phideta3;
    fesd_map->get_d3phidxi2eta_map() = d3phidxi2eta;
    fesd_map->get_d3phideta2xi_map() = d3phideta2xi;
    
	STOP_LOG("init_shape_functions()", "FESubdiv");
}



void FESubdiv::attach_quadrature_rule(QBase *q)
{
	libmesh_assert(q != NULL);

	// currently, only Gauss quadrature is supported
	libmesh_assert(q->type() == QGAUSS);

	qrule = q;
	// make sure we don't cache results from a previous quadrature rule
	elem_type = INVALID_ELEM;
	return;
}



void FESubdiv::reinit(const Elem* elem,
			  const std::vector<Point>* const pts,
			  const std::vector<Real>* const)
{
	libmesh_assert(elem != NULL);
	libmesh_assert(elem->type() == TRI3SD);
	const Tri3SD* sd_elem = static_cast<const Tri3SD*>(elem);

	START_LOG("reinit()", "FESubdiv");

	libmesh_assert(!sd_elem->is_ghost());
	libmesh_assert(sd_elem->is_subdiv_updated());

	// check if vertices 1 and 2 are regular
	libmesh_assert(sd_elem->get_ordered_valence(1) == 6);
	libmesh_assert(sd_elem->get_ordered_valence(2) == 6);

	// no custom quadrature support
	libmesh_assert(pts == NULL);
	libmesh_assert(qrule != NULL);
	qrule->init(elem->type());

	// Initialize the shape functions
	this->init_shape_functions(this->qrule->get_points(), elem);

	// The shape functions correspond to the qrule
	shapes_on_quadrature = true;

  // Compute the map for this element.
  this->_fe_map->compute_map (this->dim, this->qrule->get_weights(), elem);

	STOP_LOG("reinit()", "FESubdiv");
}

    
Real FESubdiv::shape_third_deriv(const ElemType type,
                                          const Order order,
                                          const unsigned int i,
                                          const unsigned int j,
                                          const Point& p)
    {
        switch (order)
        {
            case FOURTH:
            {
                switch (type)
                {
                    case TRI3SD:
                        libmesh_assert(i < 12);
                        return FESubdiv::regular_shape_third_deriv(i,j,p(0),p(1));
                    default:
                        std::cerr << "ERROR: Unsupported element type!" << std::endl;
                        libmesh_error();
                }
            }
            default:
                std::cerr << "ERROR: Unsupported polynomial order!" << std::endl;
                libmesh_error();
        }
        
        libmesh_error();
        return 0.;
    }
    
    
Real FESubdiv::shape_third_deriv(const Elem* elem,
                                          const Order order,
                                          const unsigned int i,
                                          const unsigned int j,
                                          const Point& p)
    {
        libmesh_assert(elem != NULL);
        return FESubdiv::shape_third_deriv(elem->type(), order, i, j, p);
    }
    
    
template <>
Real FE<2,SUBDIV>::shape(const ElemType type,
			  const Order order,
			  const unsigned int i,
			  const Point& p)
{
	switch (order)
	{      
		case FOURTH:
		{
			switch (type)
			{
				case TRI3SD:
					libmesh_assert(i < 12);
					return FESubdiv::regular_shape(i,p(0),p(1));
				default:
					std::cerr << "ERROR: Unsupported element type!" << std::endl;
					libmesh_error();
			}
		}
		default:
			std::cerr << "ERROR: Unsupported polynomial order!" << std::endl;
			libmesh_error();
	}

	libmesh_error();
	return 0.;
}


template <>
Real FE<2,SUBDIV>::shape(const Elem* elem,
			  const Order order,
			  const unsigned int i,
			  const Point& p)
{
	libmesh_assert(elem != NULL);
	return FE<2,SUBDIV>::shape(elem->type(), order, i, p);
}


template <>
Real FE<2,SUBDIV>::shape_deriv(const ElemType type,
			  const Order order,
			  const unsigned int i,
			  const unsigned int j,
			  const Point& p)
{
	switch (order)
	{      
		case FOURTH:
		{
			switch (type)
			{
				case TRI3SD:
					libmesh_assert(i < 12);
					return FESubdiv::regular_shape_deriv(i,j,p(0),p(1));
				default:
					std::cerr << "ERROR: Unsupported element type!" << std::endl;
					libmesh_error();
			}
		}
		default:
			std::cerr << "ERROR: Unsupported polynomial order!" << std::endl;
			libmesh_error();
	}

	libmesh_error();
	return 0.;
}



template <>
Real FE<2,SUBDIV>::shape_deriv(const Elem* elem,
			  const Order order,
			  const unsigned int i,
			  const unsigned int j,
			  const Point& p)
{
	libmesh_assert(elem != NULL);
	return FE<2,SUBDIV>::shape_deriv(elem->type(), order, i, j, p);
}


template <>
Real FE<2,SUBDIV>::shape_second_deriv(const ElemType type,
			  const Order order,
			  const unsigned int i,
			  const unsigned int j,
			  const Point& p)
{
	switch (order)
	{      
		case FOURTH:
		{
			switch (type)
			{
				case TRI3SD:
					libmesh_assert(i < 12);
					return FESubdiv::regular_shape_second_deriv(i,j,p(0),p(1));
				default:
					std::cerr << "ERROR: Unsupported element type!" << std::endl;
					libmesh_error();
			}
		}
		default:
			std::cerr << "ERROR: Unsupported polynomial order!" << std::endl;
			libmesh_error();
	}
	
	libmesh_error();
	return 0.;
}


template <>
Real FE<2,SUBDIV>::shape_second_deriv(const Elem* elem,
			  const Order order,
			  const unsigned int i,
			  const unsigned int j,
			  const Point& p)
{
	libmesh_assert(elem != NULL);
	return FE<2,SUBDIV>::shape_second_deriv(elem->type(), order, i, j, p);
}

    



template <>
void FE<2,SUBDIV>::nodal_soln(const Elem* elem,
			  const Order,
			  const std::vector<Number>& elem_soln,
			  std::vector<Number>& nodal_soln)
{
	libmesh_assert(elem != NULL);
	libmesh_assert(elem->type() == TRI3SD);
	const Tri3SD* sd_elem = static_cast<const Tri3SD*>(elem);
	
	nodal_soln.resize(3); // three nodes per element

	// Ghost nodes are auxiliary.
	if (sd_elem->is_ghost())
	{
		nodal_soln[0] = 0;
		nodal_soln[1] = 0;
		nodal_soln[2] = 0;
		return;
	}

	// First node (node 0 in the element patch):
	unsigned int j = sd_elem->local_node_number(sd_elem->get_ordered_node(0)->id());
	nodal_soln[j] = elem_soln[0];

	// Second node (node 1 in the element patch):
	j = sd_elem->local_node_number(sd_elem->get_ordered_node(1)->id());
	nodal_soln[j] = elem_soln[1];

	// Third node (node 'valence' in the element patch):
	j = sd_elem->local_node_number(sd_elem->get_ordered_node(2)->id());
	nodal_soln[j] = elem_soln[sd_elem->get_ordered_node(0)->valence()];
}



// the empty template specializations below are needed to avoid
// linker reference errors, but should never get called
template <>
void FE<2,SUBDIV>::side_map(const Elem*,
			  const Elem*,
			  const unsigned int,
			  const std::vector<Point>&,
			  std::vector<Point>&)
{
	libmesh_error();
}

template <>
void FE<2,SUBDIV>::edge_reinit(Elem const*,
			  unsigned int,
			  Real,
			  const std::vector<Point>* const,
			  const std::vector<Real>* const)
{
	libmesh_error();
}

template <>
Point FE<2,SUBDIV>::inverse_map(const Elem*,
			  const Point&,
			  const Real,
			  const bool)
{
	libmesh_error();
}

template <>
void FE<2,SUBDIV>::inverse_map(const Elem*,
			  const std::vector<Point>&,
			  std::vector<Point>&,
			  Real,
			  bool)
{
	libmesh_error();
}



// Loop subdivision elements are triangles
template <> unsigned int FE<2,SUBDIV>::n_dofs(const ElemType, const Order) { return 3; }

// Loop subdivision elements have only a single dof per node
template <> unsigned int FE<2,SUBDIV>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 1; }

// Subdivision FEMs have dofs only at the nodes
template <> unsigned int FE<2,SUBDIV>::n_dofs_per_elem(const ElemType, const Order) { return 0; }
  
// Subdivision FEMs have dofs only at the nodes
template <> void FE<2,SUBDIV>::dofs_on_side(const Elem *const, const Order, unsigned int, std::vector<unsigned int> &di) { di.resize(0); }
template <> void FE<2,SUBDIV>::dofs_on_edge(const Elem *const, const Order, unsigned int, std::vector<unsigned int> &di) { di.resize(0); }

// Subdivision FEMs are C^1 continuous
template <> FEContinuity FE<2,SUBDIV>::get_continuity() const { return C_ONE; }

// Subdivision FEMs are not hierarchic
template <> bool FE<2,SUBDIV>::is_hierarchic() const { return false; }

// Subdivision FEM shapes need reinit
template <> bool FE<2,SUBDIV>::shapes_need_reinit() const { return true; }

} // namespace libMesh
