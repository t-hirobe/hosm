#ifndef LOGGER_H__
#define LOGGER_H__

#include <vector>
#include <fstream>
#include <complex>
#include <boost/numeric/ublas/vector.hpp>
#include "SpectrumProperty.h"
#include "Define.h"


class Logger
{
public:
    explicit Logger() {}
    virtual ~Logger() { if(isOpened_) close(); }
    void open(std::string filename);
    void set(int nx, int ny);
    void addLog(double time, const ubvector_d& eta, const ubvector_d& phi, const SpectrumProperty& specProp);
    void close();

private:
    struct Dat
    {
        explicit Dat(int ix, int iy) : x(ix), y(iy) {}
        int x;
        int y;
    };

    int nx_, ny_;
    std::vector<complex_d> input_;
    std::vector<complex_d> output_;
    bool isOpened_;

    std::vector<Dat> targets_;
    std::vector<std::ofstream*> outTargets_;
    std::ofstream whole_;
};


#endif /* LOGGER_H__ **/
