#ifndef FFT_2D__
#define FFT_2D__

#include <complex>
#include <fftw3.h>
#include <boost/unordered_map.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include "Define.h"


class FFT2D
{
public:
    virtual ~FFT2D();
    static void forward( const complex_d *input, complex_d *output, int nx, int ny);
    static void backward(const complex_d *input, complex_d *output, int nx, int ny);
    // rft
//    static void forward( double *input, complex_d *output, int nx, int ny);
    static void forward(const double *input, complex_d *output, int nx, int ny);
    static void backward(const complex_d *input, double *output, int nx, int ny);
    static void forward( const std::vector<double> &input, std::vector<complex_d> &output, int nx, int ny) {
        forward(&input[0], &output[0], nx, ny);
    }
    static void backward(const std::vector<complex_d> &input, std::vector<double> &output, int nx, int ny) {
        backward(&input[0], &output[0], nx, ny);
    }
    static void forward( const ubvector_d &input, ubvector_c &output, int nx, int ny) {
        forward(&input[0], &output[0], nx, ny);
    }
    static void backward(const ubvector_c &input, ubvector_d &output, int nx, int ny) {
        backward(&input[0], &output[0], nx, ny);
    }

    static void forward( const std::vector<complex_d> &input, std::vector<complex_d> &output, int nx, int ny) {
        forward(&input[0], &output[0], nx, ny);
    }
    static void backward(const std::vector<complex_d> &input, std::vector<complex_d> &output, int nx, int ny) {
        backward(&input[0], &output[0], nx, ny);
    }
    static void forward( const ubvector_c &input, ubvector_c &output, int nx, int ny) {
        forward(&input[0], &output[0], nx, ny);
    }
    static void backward(const ubvector_c &input, ubvector_c &output, int nx, int ny) {
        backward(&input[0], &output[0], nx, ny);
    }

    struct plan_key_
    {
        explicit plan_key_(const void* input, const void* output, int nx, int ny, int direction):
            input_(input), output_(output), nx_(nx), ny_(ny), direction_(direction) {};
        const void* input_;
        const void* output_;
        int nx_;
        int ny_;
        int direction_;
    };

private:
    FFT2D() {};
    FFT2D(const FFT2D& obj) {};
    static FFT2D* instance_;
    typedef boost::unordered_map<plan_key_, fftw_plan> hash_map;
    hash_map plan_map_;

    fftw_plan fft2d_plan_create_(complex_d* input,
                                 complex_d* output,
                                 int nx, int ny,
                                 int direction,
                                 unsigned flags = FFTW_MEASURE);
    fftw_plan fft2d_r2c_plan_create_(double* input,
                                     complex_d* output,
                                     int nx, int ny,
                                     unsigned flags = FFTW_MEASURE);
    fftw_plan fft2d_c2r_plan_create_(complex_d* input,
                                     double* output,
                                     int nx, int ny,
                                     unsigned flags = FFTW_MEASURE);
    void execute_plan_(const complex_d *input, complex_d *output, int nx, int ny, int direction);
    void execute_plan_(const double *input, complex_d *output, int nx, int ny);
    void execute_plan_(const complex_d *input, double *output, int nx, int ny);
};


#endif /* FFT_2D__ */
