#include "SpectrumProperty.h"
#include <iostream>
#include <cstdio>

using namespace std;


string SpectrumProperty::keygen_(int kx, int ky) const
{
    char str[32];
    sprintf(str, "%d_%d", kx, ky);
    return string(str);
}

SpectrumProperty::SpectrumProperty(int nx, int ny, double lx, double ly, int M):
    ny_(ny), klist_(nx*ny), norm_kxlist_(nx*ny), norm_kylist_(nx*ny), low_pass_filter_(nx*ny), id_sort_kxky_(nx*ny),
    klist2_(nx*ny), norm_kxlist2_(nx*ny), norm_kylist2_(nx*ny), low_pass_filter2_(nx*ny), id_sort_kxky2_(nx*ny)
{
    dkx_ = 2.0*M_PI/lx;
    dky_ = 2.0*M_PI/ly;
    for(int ix=0; ix<nx; ix++) {
        for(int iy=0; iy<ny; iy++) {
            int id = ix*ny + iy;
            // important here
            int norm_kx = nx/2<=ix ? -nx+ix : ix;
            int norm_ky = ny/2<=iy ? -ny+iy : iy;
            double kx = norm_kx * dkx_;
            double ky = norm_ky * dky_;
            double k = sqrt(kx*kx+ky*ky);
            keymap_[keygen_(norm_kx, norm_ky)] = id;
            klist_[id] = k;
            norm_kxlist_[id] = norm_kx;
            norm_kylist_[id] = norm_ky;
            if(-nx/double(M+1) <= norm_kx && norm_kx <= nx/double(M+1) &&
               -ny/double(M+1) <= norm_ky && norm_ky <= ny/double(M+1)) {
                low_pass_filter_[id] =  1;
            } else {
                low_pass_filter_[id] =  0;
            }
        }
    }

    int id = 0;
    for(int norm_kx=-nx/2; norm_kx<nx/2; norm_kx++) {
        for(int norm_ky=-ny/2; norm_ky<ny/2; norm_ky++) {
            id_sort_kxky_[id] = key(norm_kx, norm_ky);
            id++;
        }
    }




    for(int ix=0; ix<nx; ix++) {
        for(int iy=0; iy<ny/2+1; iy++) {
            int id = ix*(ny/2+1) + iy;
            int norm_kx = nx/2<=ix ? -nx+ix : ix;
            int norm_ky = iy;
            double kx = norm_kx * dkx_;
            double ky = norm_ky * dky_;
            double k = sqrt(kx*kx+ky*ky);
            keymap2_[keygen_(norm_kx, norm_ky)] = id;
            klist2_[id] = k;
            norm_kxlist2_[id] = norm_kx;
            norm_kylist2_[id] = norm_ky;
            if(-nx/double(M+1) <= norm_kx && norm_kx <= nx/double(M+1) &&
               -ny/double(M+1) <= norm_ky && norm_ky <= ny/double(M+1)) {
                low_pass_filter2_[id] =  1;
            } else {
                low_pass_filter2_[id] =  0;
            }
        }
    }
    id = 0;
    for(int norm_kx=-nx/2; norm_kx<nx/2; norm_kx++) {
        for(int norm_ky=-ny/2; norm_ky<ny/2; norm_ky++) {
            id_sort_kxky2_[id] = key(norm_kx, norm_ky);
            id++;
        }
    }


}
