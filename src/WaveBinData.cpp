#include <iostream>
#include <fstream>
#include "WaveBinData.h"
#include "FFT2D.h"

using namespace std;


void WaveBinData::write_int_(int dat) {
    file_.write((char *) &dat, sizeof(int));
}

void WaveBinData::write_float_(float dat) {
    file_.write((char *) &dat, sizeof(float));
}

void WaveBinData::write_grid_int_(const int *dat)
{
    for(int i=0; i<nx_*ny_; i++) {
        int id = specProp_.id_sort_kxky()[i];
        if(specProp_.low_pass_filter()[id] == 0) continue;
        write_int_(dat[id]);
    }
}

void WaveBinData::write_grid_complex_(const complex_d *dat)
{
    for(int i=0; i<nx_*ny_; i++) {
        int id = specProp_.id_sort_kxky()[i];
        if(specProp_.low_pass_filter()[id] == 0) continue;
        write_float_(abs(dat[id]));
    }
    for(int i=0; i<nx_*ny_; i++) {
        int id = specProp_.id_sort_kxky()[i];
        if(specProp_.low_pass_filter()[id] == 0) continue;
        write_float_(arg(dat[id]));
    }
}

WaveBinData::WaveBinData(string fileName, int nx, int ny, const SpectrumProperty& specProp)
    : fileName_(fileName), nx_(nx), ny_(ny), specProp_(specProp),
      eta_k_(nx_*ny_, complex_d(0, 0)),
      phi_k_(nx_*ny_, complex_d(0, 0)),
      eta_r_(nx_*ny_, complex_d(0, 0)),
      phi_r_(nx_*ny_, complex_d(0, 0))
{
    file_.open(fileName_.c_str(), ios::out | ios::binary | ios::trunc);
    isOpened_ = true;
}

WaveBinData::~WaveBinData() {
    if(isOpened_) close();
}

void WaveBinData::addData(const double *eta, const double *phi, double time) {
    for(int id=0; id<nx_*ny_; id++) eta_r_[id] = complex_d(eta[id], 0);
    for(int id=0; id<nx_*ny_; id++) phi_r_[id] = complex_d(phi[id], 0);
    FFT2D::forward(eta_r_, eta_k_, nx_, ny_);
    FFT2D::forward(phi_r_, phi_k_, nx_, ny_);
    write_float_(time);
    write_grid_complex_(&eta_k_[0]);
    write_grid_complex_(&phi_k_[0]);
}

void WaveBinData::addHeader(int level, double lx, double ly, double dt, double endTime, int outputStep, double gravity) {
    write_int_( level );
    write_int_( nx_ );
    write_int_( ny_ );
    write_float_( lx );
    write_float_( ly );
    write_float_( specProp_.dkx() );
    write_float_( specProp_.dky() );
    write_float_( dt );
    write_float_( endTime );
    write_int_( outputStep );
    write_float_( gravity );
    write_grid_int_( &(specProp_.norm_kxlist()[0]) );
    write_grid_int_( &(specProp_.norm_kylist()[0]) );
    file_.flush();
}

void WaveBinData::close() {
    file_.close();
    isOpened_ = false;
}
