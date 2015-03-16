#include "FFT2D.h"

using namespace std;

FFT2D* FFT2D::instance_ = new FFT2D();

bool operator==(const FFT2D::plan_key_& a, const FFT2D::plan_key_& b);
std::size_t hash_value(const FFT2D::plan_key_& p);


FFT2D::~FFT2D() {
    //destroy all registered plan
    hash_map::const_iterator it = plan_map_.begin();
    for(hash_map::const_iterator it = plan_map_.begin();
        it != plan_map_.end();
        it++) {
        fftw_destroy_plan( it->second );
    }
}

fftw_plan FFT2D::fft2d_plan_create_(complex_d* input,
                                    complex_d* output,
                                    int nx, int ny,
                                    int direction,
                                    unsigned flags)
{
    return fftw_plan_dft_2d(
        nx, ny,
        reinterpret_cast<fftw_complex*>(input),
        reinterpret_cast<fftw_complex*>(output),
        direction,
        flags
        );
}

fftw_plan FFT2D::fft2d_r2c_plan_create_(double* input,
                                        complex_d* output,
                                        int nx, int ny,
                                        unsigned flags)
{
    return fftw_plan_dft_r2c_2d(
        nx, ny,
        input,
        reinterpret_cast<fftw_complex*>(output),
        flags
        );
}

fftw_plan FFT2D::fft2d_c2r_plan_create_(complex_d* input,
                                        double* output,
                                        int nx, int ny,
                                        unsigned flags)
{
    return fftw_plan_dft_c2r_2d(
        nx, ny,
        reinterpret_cast<fftw_complex*>(input),
        output,
        flags
        );
}

void FFT2D::forward(const double *input, complex_d *output, int nx, int ny) {
    instance_->execute_plan_(input, output, nx, ny);
    double scale = 1./((double)nx*ny);
    for(int i=0; i<nx*(ny/2+1); i++) output[i] *= scale;
}

void FFT2D::backward(const complex_d *input, double *output, int nx, int ny) {
    instance_->execute_plan_(input, output, nx, ny);
}

void FFT2D::execute_plan_(const complex_d *input, double *output, int nx, int ny) {
    // check whather the plan is already created or not
    // key is input, output, nx and ny
    plan_key_ key(input, output, nx, ny, FFTW_BACKWARD);

    hash_map::const_iterator it = plan_map_.find(key);
    fftw_plan plan;
    if(it != plan_map_.end()) {
        plan = it->second;
    } else {
        vector<complex_d> tmp(nx*(ny/2+1));
        for(int i=0; i<nx*(ny/2+1); i++) tmp[i] = input[i];
        plan = fft2d_c2r_plan_create_(const_cast<complex_d*>(input), output, nx, ny, FFTW_MEASURE);
        for(int i=0; i<nx*(ny/2+1); i++) (const_cast<complex_d*>(input))[i] = tmp[i];
        plan_map_[key] = plan;
    }

    fftw_execute(plan);
}

void FFT2D::execute_plan_(const double *input, complex_d *output, int nx, int ny) {
    // check whather the plan is already created or not
    // key is input, output, nx and ny
    plan_key_ key(input, output, nx, ny, FFTW_FORWARD);

    hash_map::const_iterator it = plan_map_.find(key);
    fftw_plan plan;
    if(it != plan_map_.end()) {
        plan = it->second;
    } else {
        vector<double> tmp(nx*ny);
        for(int i=0; i<nx*ny; i++) tmp[i] = input[i];
        plan = fft2d_r2c_plan_create_(const_cast<double*>(input), output, nx, ny, FFTW_MEASURE);
        for(int i=0; i<nx*ny; i++) (const_cast<double*>(input))[i] = tmp[i];
        plan_map_[key] = plan;
    }

    fftw_execute(plan);
}


void FFT2D::forward(const complex_d *input, complex_d *output, int nx, int ny) {
    instance_->execute_plan_(input, output, nx, ny, FFTW_FORWARD);
    double scale = 1./((double)nx*ny);
    for(int i=0; i<nx*ny; i++) output[i] *= scale;
}

void FFT2D::backward(const complex_d *input, complex_d *output, int nx, int ny) {
    instance_->execute_plan_(input, output, nx, ny, FFTW_BACKWARD);
}


void FFT2D::execute_plan_(const complex_d *input, complex_d *output, int nx, int ny, int direction) {
    // check whather the plan is already created or not
    // key is input, output, nx and ny
    plan_key_ key(input, output, nx, ny, direction);

    hash_map::const_iterator it = plan_map_.find(key);
    fftw_plan plan;
    if(it != plan_map_.end()) {
        plan = it->second;
    } else {
        vector<complex_d> tmp(nx*ny);
        for(int i=0; i<nx*ny; i++) tmp[i] = input[i];
        plan = fft2d_plan_create_(const_cast<complex_d*>(input), output, nx, ny, direction, FFTW_MEASURE);
        for(int i=0; i<nx*ny; i++) (const_cast<complex_d*>(input))[i] = tmp[i];
        plan_map_[key] = plan;
    }

    fftw_execute(plan);
}


bool operator==(const FFT2D::plan_key_& a, const FFT2D::plan_key_& b) {
    return a.input_     == b.input_ &&
           a.output_    == b.output_ &&
           a.direction_ == b.direction_ &&
           a.nx_ == b.nx_ && a.ny_ == b.ny_;
}

std::size_t hash_value(const FFT2D::plan_key_& p) {
    std::size_t seed = 0;
    boost::hash_combine(seed, p.input_);
    boost::hash_combine(seed, p.output_);
    boost::hash_combine(seed, p.direction_);
    boost::hash_combine(seed, p.nx_);
    boost::hash_combine(seed, p.ny_);
    return seed;
}
