#include "HosmCalculator.h"
#include "FFT2D.h"

using namespace std;
using namespace boost::numeric::ublas;


HosmCalculator::HosmCalculator(
    int nx, int ny, double dx, double dy, int level,
    TimeIntegration ti, double gravity)
    : nx_(nx), ny_(ny), dx_(dx), dy_(dy), grav_(gravity), level_(level),
      specProp_(nx, ny, nx*dx, ny*dy, level),
      output_c_(nx*(ny/2+1))
{
    gradPhi_ = new FFTGradientCalculator(nx, ny, dx, dy);
    gradEta_ = new FFTGradientCalculator(nx, ny, dx, dy);
//    gradPhi_ = new Upwind2ndGradientCalculator(nx, ny, dx, dy);
//    gradEta_ = new Upwind2ndGradientCalculator(nx, ny, dx, dy);
//    gradPhi_ = new Central2ndGradientCalculator(nx, ny, dx, dy);
//    gradEta_ = new Central2ndGradientCalculator(nx, ny, dx, dy);
    if(ti == RK4) {
        next_ = &HosmCalculator::nextRK4_;
        k_eta_.assign(4, ubvector_d(nx*ny));
        k_phi_.assign(4, ubvector_d(nx*ny));
    } else if(ti == RK3) {
        next_ = &HosmCalculator::nextRK3_;
        k_eta_.assign(3, ubvector_d(nx*ny));
        k_phi_.assign(3, ubvector_d(nx*ny));
    } else if(ti == RK2) {
        next_ = &HosmCalculator::nextRK2_;
        k_eta_.assign(2, ubvector_d(nx*ny));
        k_phi_.assign(2, ubvector_d(nx*ny));
    } else {
        next_ = &HosmCalculator::nextEuler_;
        k_eta_.assign(1, ubvector_d(nx*ny));
        k_phi_.assign(1, ubvector_d(nx*ny));
    }

    switch(level) {
    case 1:
        wCreater_ = new WCreaterImpl<1>(nx, ny);
        break;
    case 2:
        wCreater_ = new WCreaterImpl<2>(nx, ny);
        break;
    case 3:
        wCreater_ = new WCreaterImpl<3>(nx, ny);
        break;
    case 4:
        wCreater_ = new WCreaterImpl<4>(nx, ny);
        break;
    case 5:
        wCreater_ = new WCreaterImpl<5>(nx, ny);
        break;
    case 6:
        wCreater_ = new WCreaterImpl<6>(nx, ny);
        break;
    default:
        wCreater_ = new WCreaterImpl<3>(nx, ny);
        // TODO: throw error
        break;
    }
}

void HosmCalculator::f_(ubvector_d& k_eta, ubvector_d& k_phi,
                        const ubvector_d& eta, const ubvector_d& phi)
{
    gradEta_->execute(eta, specProp_);
    gradPhi_->execute(phi, specProp_);
    const ubvector_d& dhx = gradEta_->derivX();
    const ubvector_d& dhy = gradEta_->derivY();
    const ubvector_d& dpx = gradPhi_->derivX();
    const ubvector_d& dpy = gradPhi_->derivY();
    const ubvector_d& gradH = element_prod(dhx, dhx) + element_prod(dhy, dhy);
    const ubvector_d& gradP = element_prod(dpx, dpx) + element_prod(dpy, dpy);
    const ubvector_d& dhdp  = element_prod(dhx, dpx) + element_prod(dhy, dpy);
    wCreater_->execute(eta, phi, specProp_);

    k_eta =
        - dhdp + wCreater_->w(level_) + element_prod(wCreater_->w(level_-2), gradH);
    k_phi =
        - grav_*eta - 0.5*gradP +
                0.5*(
                    wCreater_->w2(level_) +
                    element_prod(wCreater_->w2(level_-2), gradH)
                    );

}


void HosmCalculator::next(ubvector_d& eta, ubvector_d& phi, double dt)
{
    // alias phi and eta
    filter(eta, dt);
    filter(phi, dt);
    return (this->*next_)(eta, phi, dt);
}

void HosmCalculator::filter(ubvector_d& input, double dt) {
    const int nx = nx_;
    const int ny = ny_;

    FFT2D::forward(input, output_c_, nx, ny);
    double nusqrt = sqrt(1.8e-6) * 0.1;
    for(int id=0; id<nx*(ny/2+1); id++) {
        if(specProp_.low_pass_filter2()[id] == 0) {
            output_c_[id] = complex_d(0, 0);
        } else {
//            double ampl = abs(output_c_[id]);
//            double phase = arg(output_c_[id]);
//            output_c_[id] = ampl*(1-2*specProp_.klist2()[id]*nusqrt*dt) * exp(complex_d(0, phase));
        }
    }
    output_c_[0] = complex_d(0, 0);  // remove singular point
    FFT2D::backward(output_c_, input, nx, ny);

}

void HosmCalculator::nextRK4_(ubvector_d& eta, ubvector_d& phi, double dt)
{
    f_(k_eta_[0], k_phi_[0], eta, phi);
    f_(k_eta_[1], k_phi_[1], eta + 0.5*dt*k_eta_[0], phi + 0.5*dt*k_phi_[0]);
    f_(k_eta_[2], k_phi_[2], eta + 0.5*dt*k_eta_[1], phi + 0.5*dt*k_phi_[1]);
    f_(k_eta_[3], k_phi_[3], eta +     dt*k_eta_[2], phi +     dt*k_phi_[2]);

    eta += dt/6.0*(k_eta_[0] + 2*k_eta_[1] + 2*k_eta_[2] + k_eta_[3]);
    phi += dt/6.0*(k_phi_[0] + 2*k_phi_[1] + 2*k_phi_[2] + k_phi_[3]);

    kinetic_energy_ = 0.5*inner_prod(phi, k_eta_[3]);
    potential_energy_ = 0.5 * grav_ * inner_prod(eta, eta);

}

void HosmCalculator::nextRK3_(ubvector_d& eta, ubvector_d& phi, double dt)
{
    f_(k_eta_[0], k_phi_[0], eta, phi);
    f_(k_eta_[1], k_phi_[1], eta + 0.5*dt*k_eta_[0], phi + 0.5*dt*k_phi_[0]);
    f_(k_eta_[2], k_phi_[2], eta - dt*k_eta_[0] + 2.0*dt*k_eta_[1], phi - dt*k_phi_[0] + 2.0*dt*k_phi_[1]);

    eta += dt/6.0*(k_eta_[0] + 4*k_eta_[1] + k_eta_[2]);
    phi += dt/6.0*(k_phi_[0] + 4*k_phi_[1] + k_phi_[2]);

    kinetic_energy_ = 0.5*inner_prod(phi, k_eta_[2]);
    potential_energy_ = 0.5 * grav_ * inner_prod(eta, eta);

}

void HosmCalculator::nextRK2_(ubvector_d& eta, ubvector_d& phi, double dt)
{
    f_(k_eta_[0], k_phi_[0], eta, phi);
    f_(k_eta_[1], k_phi_[1], eta + 0.5*dt*k_eta_[0], phi + 0.5*dt*k_phi_[0]);

    eta += dt*k_eta_[1];
    phi += dt*k_phi_[1];

    kinetic_energy_ = 0.5*inner_prod(phi, k_eta_[1]);
    potential_energy_ = 0.5 * grav_ * inner_prod(eta, eta);

}

void HosmCalculator::nextEuler_(ubvector_d& eta, ubvector_d& phi, double dt)
{
    f_(k_eta_[0], k_phi_[0], eta, phi);

    eta += dt*k_eta_[0];
    phi += dt*k_phi_[0];

    kinetic_energy_ = 0.5*inner_prod(phi, k_eta_[0]);
    potential_energy_ = 0.5 * grav_ * inner_prod(eta, eta);
}
