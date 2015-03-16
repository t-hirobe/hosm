#ifndef UTIL_H__
#define UTIL_H__

#include <complex>
#include <string>
#include <fftw3.h>
#include "SpectrumProperty.h"
#include "Define.h"

void writeGrid(const double* vals, std::string fileName, int nx, int ny);
void writeSpectrum(const double* vals, std::string fileName, int nx, int ny, const SpectrumProperty& specProp);

void  fft(const double* input, complex_d* output, int nx, int ny);

#endif /* UTIL_H__ **/
