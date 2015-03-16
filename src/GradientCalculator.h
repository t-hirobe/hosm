#ifndef GRADIEND_CALCULATOR__
#define GRADIEND_CALCULATOR__

#include <boost/numeric/ublas/vector.hpp>
#include <vector>
#include <complex>
#include "SpectrumProperty.h"
#include "Define.h"

class GradientCalculator
{
public:
    explicit GradientCalculator(int nx, int ny, double dx, double dy)
        : nx_(nx), ny_(ny), dx_(dx), dy_(dy),
        derivX_(nx*ny), derivY_(nx*ny) {}
    virtual ~GradientCalculator() {}
    virtual void execute(const ubvector_d& val, const SpectrumProperty& specProp) = 0;

    const ubvector_d& derivX() const { return derivX_; }
    const ubvector_d& derivY() const { return derivY_; }

protected:
    int nx_;
    int ny_;
    double dx_;
    double dy_;
    ubvector_d derivX_;
    ubvector_d derivY_;
};


class FFTGradientCalculator : public GradientCalculator
{
public:
    explicit FFTGradientCalculator(int nx, int ny, double dx, double dy);
    virtual ~FFTGradientCalculator();
    void execute(const ubvector_d& val, const SpectrumProperty& specProp);

private:
    ubvector_c invVal_;
    ubvector_c invValX_;
    ubvector_c invValY_;
};

class Central2ndGradientCalculator : public GradientCalculator
{
public:
    explicit Central2ndGradientCalculator(int nx, int ny, double dx, double dy)
        : GradientCalculator(nx, ny, dx, dy) {};
    virtual ~Central2ndGradientCalculator(){}
    void execute(const ubvector_d& val, const SpectrumProperty& specProp);

};
/*
class Upwind2ndGradientCalculator : public GradientCalculator
{
public:
    explicit Upwind2ndGradientCalculator(int nx, int ny, double dx, double dy)
        : GradientCalculator(nx, ny, dx, dy) {};
    virtual ~Upwind2ndGradientCalculator(){}
    void execute(const ubvector_d& val, const SpectrumProperty& specProp);

};
*/
#endif /* GRADIEND_CALCULATOR__ **/
