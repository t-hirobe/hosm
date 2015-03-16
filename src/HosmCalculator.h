#ifndef HOSM_CALCULATOR__
#define HOSM_CALCULATOR__

#include "WCreater.h"
#include "GradientCalculator.h"
#include "SpectrumProperty.h"
#include "Define.h"

class HosmCalculator
{
public:
    enum TimeIntegration {
        Euler = 0,
        RK4 = 1,
        RK3 = 2,
        RK2 = 3
    };

    explicit HosmCalculator(int nx, int ny, double dx, double dy, int level,
                            TimeIntegration ti, double gravity);
    virtual ~HosmCalculator(){ delete gradPhi_; delete gradEta_; delete wCreater_; }
    void setDx(double dx) { dx_ = dx; }
    void setDy(double dy) { dy_ = dy; }
    void next(ubvector_d& eta, ubvector_d& phi, double dt);
    const SpectrumProperty& specProp() const { return specProp_; }
    double getKineticEnergy() const { return kinetic_energy_; }
    double getPotentialEnergy() const { return potential_energy_; }

private:
    int nx_;
    int ny_;
    double dx_;
    double dy_;
    double grav_;
    int level_;
    SpectrumProperty specProp_;
    GradientCalculator *gradPhi_;
    GradientCalculator *gradEta_;
    std::vector<ubvector_d> k_eta_;
    std::vector<ubvector_d> k_phi_;
    WCreater *wCreater_;
    double kinetic_energy_;
    double potential_energy_;
    ubvector_c output_c_;

    void f_(ubvector_d& k_eta, ubvector_d& k_phi,
            const ubvector_d& eta, const ubvector_d& phi);
    void (HosmCalculator::*next_)(ubvector_d& eta, ubvector_d& phi, double dt);
    void nextEuler_(ubvector_d& eta, ubvector_d& phi, double dt);
    void nextRK4_(ubvector_d& eta, ubvector_d& phi, double dt);
    void nextRK3_(ubvector_d& eta, ubvector_d& phi, double dt); // RK3 Kutta
    void nextRK2_(ubvector_d& eta, ubvector_d& phi, double dt);
    void filter(ubvector_d& val, double dt);
};


#endif /** HOSM_CALCULATOR__ */
