#ifndef W_CREATER_H__
#define W_CREATER_H__

#include <complex>
#include <vector>
#include <boost/numeric/ublas/vector.hpp>
#include "SpectrumProperty.h"
#include "Define.h"

class WCreater
{
public:
    explicit WCreater() {}
    virtual ~WCreater() {}
    virtual void setup(int nx, int ny) = 0;
    virtual void execute(const ubvector_d& eta, const ubvector_d& phi, const SpectrumProperty& specProp) = 0;

    virtual const ubvector_d& wComponent(int level) = 0;
    virtual ubvector_d w(int level) = 0;
    virtual ubvector_d w2(int level) = 0;
};

template<int Lv>
class WCreaterImpl : public WCreater
{
public:
    explicit WCreaterImpl() : WCreater() {}
    explicit WCreaterImpl(int nx, int ny) : WCreater() { setup(nx, ny);}
    virtual ~WCreaterImpl() {}
    void setup(int nx, int ny);
    void execute(const ubvector_d& eta, const ubvector_d& phi, const SpectrumProperty& specProp);

    const ubvector_d& wComponent(int level) { return w_[level]; }
    ubvector_d w(int level);
    ubvector_d w2(int level);

private:
    int nx_, ny_;
    typedef std::vector< ubvector_c > map_imag_c;
    typedef std::vector< ubvector_d > map_imag_d;
    map_imag_d phi_;
    map_imag_d phidz_;
    map_imag_c c_;
    map_imag_c phidz_k_;
    std::vector<ubvector_d> w_;
    int key_(int i, int j);
    void calcDerivPhi_(const ubvector_d& eta, const ubvector_d& phi, const SpectrumProperty& specProp);

};

template class WCreaterImpl<1>;
template class WCreaterImpl<2>;
template class WCreaterImpl<3>;
template class WCreaterImpl<4>;
template class WCreaterImpl<5>;
template class WCreaterImpl<6>;

#endif /** W_CREATER_H__ */
