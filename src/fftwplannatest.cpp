#include <iostream>
#include <cstdlib>
#include <complex>
#include <fftw3.h>
#include <vector>
#include <cmath>

using namespace std;

typedef complex<double> complex_d;

void filter(vector<double>& input, int nx, int ny) {
    vector<complex_d> input_c(nx*ny);
    vector<complex_d> output_c(nx*ny);
//    fft2d_.forward(input_c_, output_c_, nx, ny);
    {
        double scale = 1./((double)nx*ny);
        fftw_plan plan =
            fftw_plan_dft_2d(
                nx, ny,
                reinterpret_cast<fftw_complex*>(&(input_c[0])),
                reinterpret_cast<fftw_complex*>(&(output_c[0])),
                FFTW_FORWARD,
                FFTW_MEASURE);
//                FFTW_ESTIMATE);
//                FFTW_PATIENT);
        for(int i=0; i<nx*ny; i++) input_c[i] = complex_d(input[i], 0);
        fftw_execute(plan);
        fftw_destroy_plan(plan);
        for(int i=0; i<nx*ny; i++) output_c[i] *= scale;
    }

//    fft2d_.backward(output_c_, input_c_, nx, ny);
    vector<complex_d> input_x(nx*ny);
    vector<complex_d> output_x(nx*ny);
    {
        fftw_plan plan =
            fftw_plan_dft_2d(
                nx, ny,
                reinterpret_cast<fftw_complex*>(&(output_x[0])),
                reinterpret_cast<fftw_complex*>(&(input_x[0])),
                FFTW_BACKWARD,
                FFTW_MEASURE);
//                FFTW_ESTIMATE);
//                FFTW_PATIENT);
        for(int i=0; i<nx*ny; i++) output_x[i] = output_c[i];
        fftw_execute(plan);
        fftw_destroy_plan(plan);
    }

    for(int i=0; i<nx*ny; i++) input[i] = input_x[i].real();

}


int main() {
    const int nx = 100;
    const int ny = 200;
    const int modes = 6;

    const int kxs[modes]      = {3,   7, -11, -3, 10, 0};
    const int kys[modes]      = {8, -12,  13, -9,  0, 5};
    const double ampls[modes] = {1,   2,   3,  4,  5, 6};
    double phases[modes];
    for(int mode=0; mode<modes; mode++) {
        const double phase = 2*M_PI * ((double)((rand() / ((double)RAND_MAX+1.0f))) - 0.5);
        phases[mode] = phase;
    }

    const double lx = 1.0;
    const double ly = 1.0;

    vector<double> h(nx*ny);
    vector<double> h0(nx*ny);
    for(int ix=0; ix<nx; ix++) {
        for(int iy=0; iy<ny; iy++) {
            int id = ix*ny + iy;
            h[id] = 0;
            const double x = ix*lx/nx;
            const double y = iy*ly/ny;
            for(int mode=0; mode<modes; mode++) {
                const double kx = 2*M_PI/lx* kxs[mode];
                const double ky = 2*M_PI/ly* kys[mode];
                h[id]      += ampls[mode] * cos(kx * x + ky * y + phases[mode]);
                h0[id] = h[id];
            }
        }
    }

    filter(h, nx, ny);
    double sum = 0;
    for(int i=0; i<nx*ny; i++) sum += fabs(h0[i] - h[i])/nx/ny;
    cout << sum << endl;


    return 0;
}
