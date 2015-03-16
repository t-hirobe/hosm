#ifndef SPECTRUM_CREATER__
#define SPECTRUM_CREATER__

#include <string>
#include "SpectrumProperty.h"

void createSpectrum(std::string inputFile, double* eta, double* phi, int nx, int ny, const SpectrumProperty& specProp, double gravity);
void createJONSWAPSpectrum(double* eta, double* phi, int nx, int ny, const SpectrumProperty& specProp, double jonswap_E, double jonswap_wp, double gravity);

#endif /* SPECTRUM_CREATER__ */
