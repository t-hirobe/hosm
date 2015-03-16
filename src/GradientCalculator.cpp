#include "GradientCalculator.h"
#include "FFT2D.h"
#include <fstream>

using namespace std;



FFTGradientCalculator::FFTGradientCalculator(int nx, int ny, double dx, double dy)
    : GradientCalculator(nx, ny, dx, dy),
      invVal_(nx*(ny/2+1)),
      invValX_(nx*(ny/2+1)),
      invValY_(nx*(ny/2+1))
{
}

void FFTGradientCalculator::execute(const ubvector_d& val, const SpectrumProperty& specProp)
{
    const int nx = nx_;
    const int ny = ny_;
    FFT2D::forward(val, invVal_, nx, ny);
    invVal_[0] = complex_d(0,0);

// NG
//    ubvector_c middat = element_prod(invVal_, specProp.low_pass_filter()) * scale_ * complex_d(0, 1);
//cout << &(invValX_[0]) << endl;
//    invValX_ = element_prod(middat, specProp.norm_kxlist()) * specProp.dkx();
//cout << &(invValX_[0]) << endl;
//    invValY_ = element_prod(specProp.norm_kylist(), middat) * specProp.dky();
    for(int ix=0; ix<nx; ix++) {
        for(int iy=0; iy<ny/2+1; iy++) {
            int id = ix*(ny/2+1) + iy;
            if(specProp.low_pass_filter2()[id] == 0) {
                invValX_[id] = complex_d(0, 0);
                invValY_[id] = complex_d(0, 0);
            } else {
                double kx = specProp.norm_kxlist2()[id] * specProp.dkx();
                double ky = specProp.norm_kylist2()[id] * specProp.dky();
                invValX_[id] = kx * invVal_[id] * complex_d(0, 1);
                invValY_[id] = ky * invVal_[id] * complex_d(0, 1);
            }

        }
    }
    FFT2D::backward(&invValX_[0], &derivX_[0], nx, ny);
    FFT2D::backward(&invValY_[0], &derivY_[0], nx, ny);
}


FFTGradientCalculator::~FFTGradientCalculator()
{
}


void Central2ndGradientCalculator::execute(const ubvector_d& val, const SpectrumProperty& specProp) {
    for(int ix=0; ix<nx_; ix++) {
        for(int iy=0; iy<ny_; iy++) {
            int ixm1 = ix-1;
            int ixp1 = ix+1;
            int iym1 = iy-1;
            int iyp1 = iy+1;
            if(ix==    0) ixm1 = nx_-1;
            if(ix==nx_-1) ixp1 = 0;
            if(iy==    0) iym1 = ny_-1;
            if(iy==ny_-1) iyp1 = 0;
            int w = ixp1*ny_+iy;
            int e = ixm1*ny_+iy;
            int n = ix*ny_+iyp1;
            int s = ix*ny_+iym1;
            int c = ix*ny_+iy;
            derivX_[c] = (val[w]-val[e])/2./dx_;
            derivY_[c] = (val[n]-val[s])/2./dy_;
        }
    }
}
/*
void Upwind2ndGradientCalculator::execute(const ubvector_d& val, const SpectrumProperty& specProp) {
    for(int ix=0; ix<nx_; ix++) {
        for(int iy=0; iy<ny_; iy++) {
            int ixm1 = ix-1;
            int ixp1 = ix+1;
            int iym1 = iy-1;
            int iyp1 = iy+1;
            if(ix==    0) ixm1 = nx_-1;
            if(ix==nx_-1) ixp1 = 0;
            if(iy==    0) iym1 = ny_-1;
            if(iy==ny_-1) iyp1 = 0;
            int w = ixp1*ny_+iy;
            int e = ixm1*ny_+iy;
            int n = ix*ny_+iyp1;
            int s = ix*ny_+iym1;
            int c = ix*ny_+iy;
            derivX_[c] = (val[w]-val[e])/2./dx_;
            derivY_[c] = (val[n]-val[s])/2./dy_;
        }
    }
}
*/
