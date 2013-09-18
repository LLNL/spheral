//---------------------------------Spheral++----------------------------------//
// Kernel -- The interpolation kernel for use in smoothed field estimates.
//
// Created by JMO, Thu Apr 22 20:48:23 PDT 1999
//----------------------------------------------------------------------------//

// #include "Kernel.hh"
// #include "BSplineKernel.hh"
// #include "W4SplineKernel.hh"
// #include "GaussianKernel.hh"
// #include "SuperGaussianKernel.hh"
// #include "PiGaussianKernel.hh"

// // BSplineKernel.
// const double Kernel<Dim<1>, BSplineKernel<Dim<1> > >::mVolumeNormalization = 2.0/3.0;
// const double Kernel<Dim<2>, BSplineKernel<Dim<2> > >::mVolumeNormalization = 10.0/(7.0*M_PI);
// const double Kernel<Dim<3>, BSplineKernel<Dim<3> > >::mVolumeNormalization = 1.0/M_PI;

// const double Kernel<Dim<1>, BSplineKernel<Dim<1> > >::mKernelExtent = 2.0;
// const double Kernel<Dim<2>, BSplineKernel<Dim<2> > >::mKernelExtent = 2.0;
// const double Kernel<Dim<3>, BSplineKernel<Dim<3> > >::mKernelExtent = 2.0;

// // W4SplineKernel
// const double Kernel<Dim<1>, W4SplineKernel<Dim<1> > >::mVolumeNormalization = 1.0;
// const double Kernel<Dim<2>, W4SplineKernel<Dim<2> > >::mVolumeNormalization = 30.0/(7.0*M_PI);
// const double Kernel<Dim<3>, W4SplineKernel<Dim<3> > >::mVolumeNormalization = 5.0/(6.0*M_PI);

// const double Kernel<Dim<1>, W4SplineKernel<Dim<1> > >::mKernelExtent = 2.0;
// const double Kernel<Dim<2>, W4SplineKernel<Dim<2> > >::mKernelExtent = 2.0;
// const double Kernel<Dim<3>, W4SplineKernel<Dim<3> > >::mKernelExtent = 2.0;

// // Gaussian Kernel
// const double Kernel<Dim<1>, GaussianKernel<Dim<1> > >::mVolumeNormalization = 1.0/sqrt(M_PI);
// const double Kernel<Dim<2>, GaussianKernel<Dim<2> > >::mVolumeNormalization = 1.0/M_PI;
// const double Kernel<Dim<3>, GaussianKernel<Dim<3> > >::mVolumeNormalization = 1.0/pow(M_PI, 1.5);

// const double Kernel<Dim<1>, GaussianKernel<Dim<1> > >::mKernelExtent = 3.0;
// const double Kernel<Dim<2>, GaussianKernel<Dim<2> > >::mKernelExtent = 3.0;
// const double Kernel<Dim<3>, GaussianKernel<Dim<3> > >::mKernelExtent = 3.0;

// // SuperGaussian Kernel
// const double Kernel<Dim<1>, SuperGaussianKernel<Dim<1> > >::mVolumeNormalization = 1.0/sqrt(M_PI);
// const double Kernel<Dim<2>, SuperGaussianKernel<Dim<2> > >::mVolumeNormalization = 1.0/(M_PI);
// const double Kernel<Dim<3>, SuperGaussianKernel<Dim<3> > >::mVolumeNormalization = 1.0/pow(M_PI, 1.5);

// const double Kernel<Dim<1>, SuperGaussianKernel<Dim<1> > >::mKernelExtent = 3.0;
// const double Kernel<Dim<2>, SuperGaussianKernel<Dim<2> > >::mKernelExtent = 3.0;
// const double Kernel<Dim<3>, SuperGaussianKernel<Dim<3> > >::mKernelExtent = 3.0;

// // PiGaussian Kernel
// const double Kernel<Dim<1>, PiGaussianKernel<Dim<1> > >::mVolumeNormalization = 2.0/3.6256099082;
// const double Kernel<Dim<2>, PiGaussianKernel<Dim<2> > >::mVolumeNormalization = 2.0/pow(M_PI, 1.5);
// const double Kernel<Dim<3>, PiGaussianKernel<Dim<3> > >::mVolumeNormalization = 1.0/(M_PI*1.2254167024);

// const double Kernel<Dim<1>, PiGaussianKernel<Dim<1> > >::mKernelExtent = 2.0;
// const double Kernel<Dim<2>, PiGaussianKernel<Dim<2> > >::mKernelExtent = 2.0;
// const double Kernel<Dim<3>, PiGaussianKernel<Dim<3> > >::mKernelExtent = 2.0;

