#include "SpectrumCreater.h"
#include <fstream>
#include <string>
#include <cmath>
#include "Define.h"
#include "FFT2D.h"

using namespace std;



void createSpectrum(string inputFile, double* eta, double* phi, int nx, int ny, const SpectrumProperty& specProp, double gravity)
{
    ifstream in(inputFile.c_str());
    ubvector_c eta_k(nx*ny, complex_d(0,0));
    ubvector_c phi_k(nx*ny, complex_d(0,0));
    if(!in.fail()) {
        string line = "";
        while(getline(in, line)) {
            if(line == "") break;
            istringstream is(line);
            int ix, iy;
            double value;
            is >> ix >> iy >> value;
            int id = specProp.key(ix, iy);
            double kx = specProp.norm_kxlist()[id] * specProp.dkx();
            double ky = specProp.norm_kylist()[id] * specProp.dky();
            double k = sqrt(kx*kx + ky*ky);
            double w = sqrt(k*gravity);
            double ampl = value;
            double phase = 2*M_PI * (double)((rand() / ((double)RAND_MAX+1.0f)));
            cout << ix << "," <<
                iy << "," <<
                kx << "," <<
                ky << "," <<
                2*M_PI/kx << "," <<
                2*M_PI/ky << "," <<
                w << "," <<
                ampl << "," <<
                phase << endl;
            eta_k[id] = ampl * exp(complex_d(0, phase));
            phi_k[id] = -gravity / w * eta_k[id] * complex_d(0, 1);
        }
        ubvector_c eta_cmpl(nx*ny);
        ubvector_c phi_cmpl(nx*ny);
        FFT2D::backward(eta_k, eta_cmpl, nx, ny);
        FFT2D::backward(phi_k, phi_cmpl, nx, ny);
        for(int i=0; i<nx*ny; i++) eta[i] = eta_cmpl[i].real();
        for(int i=0; i<nx*ny; i++) phi[i] = phi_cmpl[i].real();
    }
    in.close();
}


double jonswap_psi(double w, double wp, double E, double gravity) {
    double alpha = 3.279 * E;
    double gamma = 3.3;
    double sigma = w < wp ? 0.07 : 0.09;

    return alpha * gravity * gravity * pow(w, -5) * exp(-5/4.0 * pow(w/wp, -4))
        * pow(gamma, exp(-pow(w-wp,2) / (2*sigma*sigma*wp*wp)) );
}

double jonswap_G(double theta) {
    double val = 0;
    if(fabs(theta) <= M_PI/2.0) val = (2/M_PI)*pow(cos(theta), 2);
    return val;
}

void createJONSWAPSpectrum(double* eta, double* phi, int nx, int ny, const SpectrumProperty& specProp, double E, double wp, double gravity)
{
    ubvector_c eta_k(nx*ny, complex_d(0,0));
    ubvector_c phi_k(nx*ny, complex_d(0,0));

    double dS = specProp.dkx() * specProp.dky();
    for(int i=0; i<nx*ny; i++) {
        int id = specProp.id_sort_kxky()[i];
        double kx = specProp.norm_kxlist()[id] * specProp.dkx();
        double ky = specProp.norm_kylist()[id] * specProp.dky();
        double k = specProp.klist()[id];
        if(k == 0) continue;
        double w = sqrt(gravity*k);
        double theta = atan2(ky, kx);
        double psi = jonswap_psi(w, wp, E, gravity);
        double G = jonswap_G(theta);
        double b2 = gravity*gravity/2.0/pow(w, 4) * psi * G * dS;
//        double ampl = 1.0/M_PI*sqrt(k/2.0/w) * sqrt(b2);
        double ampl = 1.0/M_PI*sqrt(k/2.0/w) * sqrt(b2) * 2;
        double phase = 2*M_PI * (double)((rand() / ((double)RAND_MAX+1.0f)));
        eta_k[id] = ampl * exp(complex_d(0, phase));
        phi_k[id] = -gravity / w * eta_k[id] * complex_d(0, 1);
    }
    ubvector_c eta_cmpl(nx*ny);
    ubvector_c phi_cmpl(nx*ny);
    FFT2D::backward(eta_k, eta_cmpl, nx, ny);
    FFT2D::backward(phi_k, phi_cmpl, nx, ny);
    for(int i=0; i<nx*ny; i++) eta[i] = eta_cmpl[i].real();
    for(int i=0; i<nx*ny; i++) phi[i] = phi_cmpl[i].real();

}
