#include "WCreater.h"
#include "FFT2D.h"
#include <cmath>

using namespace std;


template<int M>
int WCreaterImpl<M>::key_(int i, int j)
{
    return i + (i+j-2)*(i+j+1)/2;
}

template<int M>
void WCreaterImpl<M>::setup(int nx, int ny)
{
    nx_ = nx;
    ny_ = ny;
    // memory allocation
    phi_.assign(M, ubvector_d(nx*ny) );
    w_.assign(M, ubvector_d(nx*ny) );
    c_.assign(M, ubvector_c(nx*(ny/2+1), complex_d(0, 0)) );
    int num = M*(M+1)/2;
    phidz_.assign(num, ubvector_d(nx*ny));
    phidz_k_.assign(num, ubvector_c(nx*(ny/2+1), complex_d(0, 0)) );
}

template<int M>
ubvector_d WCreaterImpl<M>::w(int level)
{
    ubvector_d tmp(nx_*ny_, 0);
    for(int l=0; l<level; l++) tmp += w_[l];
    return tmp;
}

template<int M>
ubvector_d WCreaterImpl<M>::w2(int level)
{
    ubvector_d tmp(nx_*ny_, 0);
    for(int l=0; l<level; l++) {
        for(int k=0; k<=l; k++) {
            tmp += element_prod(w_[k], w_[l-k]);
        }
    }

    return tmp;
}

template<int M>
void WCreaterImpl<M>::calcDerivPhi_(const ubvector_d& eta, const ubvector_d& phi, const SpectrumProperty& specProp)
{
    const int nx = nx_;
    const int ny = ny_;
    // copy
    for(int i=0; i<nx*ny; i++) phi_[0][i] = phi[i];

    for(int m=0; m<M; m++) {

        {
            // derive  (d/dz)...((d/dz)^(M-m))phi_[m] from phi_[m]
            FFT2D::forward(phi_[m], c_[m], nx, ny);

            ubvector_d powk(nx*ny, 1.0);
          // (d/dz)^k phi_[m]
            for(int k=1; k<=M-m; k++) {
                const int id = key_(k, m);
                for(int ix=0; ix<nx; ix++) {
                    for(int iy=0; iy<ny/2+1; iy++) {
                        int i = ix*(ny/2+1) + iy;
                        if(specProp.low_pass_filter2()[i] == 0) {
                            phidz_k_[id][i] = complex_d(0, 0);
                        } else {
                            powk[i] *= specProp.klist2()[i];
                            phidz_k_[id][i] = powk[i] * c_[m][i];
                        }
                    }
                }
                phidz_k_[id][specProp.key2(0, 0)] = complex_d(0, 0);

                FFT2D::backward(phidz_k_[id], phidz_[id], nx, ny);
            }
        }

        if(m+1 == M) break;

        {
            // derive phi_[m+1]
            ubvector_d poweta(nx*ny, 1.0);
            double fact = 1;
            for(int id=0; id<nx*ny; id++) phi_[m+1][id] = 0;
            for(int k=1; k<=m+1; k++) {
                fact *= k;
                for(int id=0; id<nx*ny; id++) {
                    poweta[id] *= eta[id];
                    phi_[m+1][id] += -poweta[id] * phidz_[key_(k, m-k)][id] / fact;
                }
            }
        }
    }

}

template<int M>
void WCreaterImpl<M>::execute(const ubvector_d& eta, const ubvector_d& phi, const SpectrumProperty& specProp)
{
    const int nx = nx_;
    const int ny = ny_;

    calcDerivPhi_(eta, phi, specProp);

    for(int m=0; m<M; m++) {
        ubvector_d poweta(nx*ny, 1.0);
        ubvector_d tmp(nx*ny, 0.0);
        double fact = 1;
        for(int k=0; k<=m; k++) {
            int id = key_(k+1, m-k);
            for(int i=0; i<nx*ny; i++) {
                tmp[i] += poweta[i] * phidz_[id][i] / fact;
                poweta[i] *= eta[i];
            }
            fact *= k+1;
        }
        for(int i=0; i<nx*ny; i++) w_[m][i] = tmp[i];
    }

}
