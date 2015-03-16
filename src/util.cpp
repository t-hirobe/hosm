#include "util.h"
#include <iostream>
#include <fstream>

using namespace std;

void writeGrid(const double* vals, string fileName, int nx, int ny) {
    ofstream out(fileName.c_str());
    for(int ix=0; ix<nx; ix++) {
        for(int iy=0; iy<ny; iy++) {
            int id = ix*ny+iy;
            out << ix << " "
                << iy << " "
                << vals[id] << endl;
        }
        out << endl;
    }
    out.close();
}

void writeSpectrum(const double* vals, string fileName, int nx, int ny, const SpectrumProperty& specProp) {
    vector<complex_d> output_c(nx*ny);
    vector<complex_d> input_c(nx*ny);

    fft(&vals[0], &output_c[0], nx, ny);
    ofstream out(fileName.c_str());
    for(int i=0; i<nx*ny; i++) {
        int id = specProp.id_sort_kxky()[i];
        int kx = specProp.norm_kxlist()[id];
        int ky = specProp.norm_kylist()[id];
        out << kx << " "
            << ky << " "
            << abs(output_c[id]) << " "
            << arg(output_c[id]) << endl;
        if(i % ny == ny-1) out << endl;
    }
    out.close();
}



void fft(const double* input,
         complex_d* output,
         int nx, int ny)
{
    double scale = 1./((double)nx*ny);
    vector<complex_d> tmp(nx*ny);
    fftw_plan plan =
        fftw_plan_dft_2d(
            nx, ny,
            reinterpret_cast<fftw_complex*>(&(tmp[0])),
            reinterpret_cast<fftw_complex*>(output),
            FFTW_FORWARD,
            FFTW_ESTIMATE
            );
    for(int i=0; i<nx*ny; i++) tmp[i] = complex_d(input[i], 0);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    for(int i=0; i<nx*ny; i++) output[i] *= scale;
}
