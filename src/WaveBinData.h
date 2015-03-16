#ifndef WAVE_BIN_DATA__
#define WAVE_BIN_DATA__

#include <string>
#include <fstream>
#include "SpectrumProperty.h"
#include "Define.h"

class WaveBinData {
public:
    explicit WaveBinData(std::string fileName, int nx, int ny, const SpectrumProperty& specProp);
    virtual ~WaveBinData();
    void addHeader(int level, double lx, double ly, double dt, double endTime, int outputStep, double gravity);
    void addData(const double *eta, const double *phi, double time);
    void close();
private:
    std::string fileName_;
    int nx_;
    int ny_;
    const SpectrumProperty& specProp_;
    std::vector<complex_d> eta_k_;
    std::vector<complex_d> phi_k_;
    std::vector<complex_d> eta_r_;
    std::vector<complex_d> phi_r_;
    std::ofstream file_;
    bool isOpened_;
    void write_int_(int dat);
    void write_float_(float dat);
    void write_grid_int_(const int *dat);
    void write_grid_complex_(const complex_d *dat);
};

#endif /* WAVE_BIN_DATA__ */
