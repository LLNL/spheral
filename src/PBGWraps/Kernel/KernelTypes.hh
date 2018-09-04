#ifndef __PBGWRAPS_KERNELTYPES__
#define __PBGWRAPS_KERNELTYPES__

#include "Geometry/Dimension.hh"
#include "Kernel/BSplineKernel.hh"
#include "Kernel/W4SplineKernel.hh"
#include "Kernel/GaussianKernel.hh"
#include "Kernel/SuperGaussianKernel.hh"
#include "Kernel/PiGaussianKernel.hh"
#include "Kernel/HatKernel.hh"
#include "Kernel/SincKernel.hh"
#include "Kernel/NSincPolynomialKernel.hh"
#include "Kernel/NBSplineKernel.hh"
#include "Kernel/QuarticSplineKernel.hh"
#include "Kernel/QuinticSplineKernel.hh"
#include "Kernel/TableKernel.hh"
#include "Kernel/WendlandC2Kernel.hh"
#include "Kernel/WendlandC4Kernel.hh"
#include "Kernel/WendlandC6Kernel.hh"
#include "Kernel/ExpInvKernel.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Names!
//------------------------------------------------------------------------------
typedef BSplineKernel<Dim<1> > BSplineKernel1d;
typedef BSplineKernel<Dim<2> > BSplineKernel2d;
typedef BSplineKernel<Dim<3> > BSplineKernel3d;
        
typedef W4SplineKernel<Dim<1> > W4SplineKernel1d;
typedef W4SplineKernel<Dim<2> > W4SplineKernel2d;
typedef W4SplineKernel<Dim<3> > W4SplineKernel3d;

typedef GaussianKernel<Dim<1> > GaussianKernel1d;
typedef GaussianKernel<Dim<2> > GaussianKernel2d;
typedef GaussianKernel<Dim<3> > GaussianKernel3d;

typedef SuperGaussianKernel<Dim<1> > SuperGaussianKernel1d;
typedef SuperGaussianKernel<Dim<2> > SuperGaussianKernel2d;
typedef SuperGaussianKernel<Dim<3> > SuperGaussianKernel3d;

typedef PiGaussianKernel<Dim<1> > PiGaussianKernel1d;
typedef PiGaussianKernel<Dim<2> > PiGaussianKernel2d;
typedef PiGaussianKernel<Dim<3> > PiGaussianKernel3d;

typedef HatKernel<Dim<1> > HatKernel1d;
typedef HatKernel<Dim<2> > HatKernel2d;
typedef HatKernel<Dim<3> > HatKernel3d;

typedef SincKernel<Dim<1> > SincKernel1d;
typedef SincKernel<Dim<2> > SincKernel2d;
typedef SincKernel<Dim<3> > SincKernel3d;

typedef NSincPolynomialKernel<Dim<1> > NSincPolynomialKernel1d;
typedef NSincPolynomialKernel<Dim<2> > NSincPolynomialKernel2d;
typedef NSincPolynomialKernel<Dim<3> > NSincPolynomialKernel3d;

typedef NBSplineKernel<Dim<1> > NBSplineKernel1d;
typedef NBSplineKernel<Dim<2> > NBSplineKernel2d;
typedef NBSplineKernel<Dim<3> > NBSplineKernel3d;

typedef QuarticSplineKernel<Dim<1> > QuarticSplineKernel1d;
typedef QuarticSplineKernel<Dim<2> > QuarticSplineKernel2d;
typedef QuarticSplineKernel<Dim<3> > QuarticSplineKernel3d;

typedef QuinticSplineKernel<Dim<1> > QuinticSplineKernel1d;
typedef QuinticSplineKernel<Dim<2> > QuinticSplineKernel2d;
typedef QuinticSplineKernel<Dim<3> > QuinticSplineKernel3d;

typedef TableKernel<Dim<1> > TableKernel1d;
typedef TableKernel<Dim<2> > TableKernel2d;
typedef TableKernel<Dim<3> > TableKernel3d;

typedef WendlandC2Kernel<Dim<1> > WendlandC2Kernel1d;
typedef WendlandC2Kernel<Dim<2> > WendlandC2Kernel2d;
typedef WendlandC2Kernel<Dim<3> > WendlandC2Kernel3d;
    
typedef WendlandC4Kernel<Dim<1> > WendlandC4Kernel1d;
typedef WendlandC4Kernel<Dim<2> > WendlandC4Kernel2d;
typedef WendlandC4Kernel<Dim<3> > WendlandC4Kernel3d;
    
typedef WendlandC6Kernel<Dim<1> > WendlandC6Kernel1d;
typedef WendlandC6Kernel<Dim<2> > WendlandC6Kernel2d;
typedef WendlandC6Kernel<Dim<3> > WendlandC6Kernel3d;

typedef ExpInvKernel<Dim<1> > ExpInvKernel1d;
typedef ExpInvKernel<Dim<2> > ExpInvKernel2d;
typedef ExpInvKernel<Dim<3> > ExpInvKernel3d;

}

#endif
