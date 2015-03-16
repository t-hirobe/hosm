#include "Logger.h"
#include "Define.h"
#include "FFT2D.h"
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <algorithm>

using namespace boost;
using namespace std;

void Logger::open(std::string filename)
{
    ifstream in(filename.c_str());
    string line;
    while(getline(in, line)) {
        istringstream is(line);
        double ix,iy;
        is >> ix >> iy;
        targets_.push_back(Dat(ix, iy));
    }
    in.close();

    BOOST_FOREACH(const Dat& dat, targets_) {
        string name = (format("kx_%d_ky_%d.out") % dat.x % dat.y).str();
        outTargets_.push_back(new ofstream(name.c_str()));
    }
    whole_.open("whole.out");
    isOpened_ = true;
}

void Logger::set(int nx, int ny)
{
    input_.assign(nx*ny, complex_d(0, 0));
    output_.assign(nx*ny, complex_d(0, 0));
    nx_ = nx;
    ny_ = ny;
}

void Logger::addLog(double time, const ubvector_d& eta, const ubvector_d& phi, const SpectrumProperty& specProp)
{
    for(int i=0; i<nx_*ny_; i++) input_[i] = eta[i];
    FFT2D::forward(input_, output_, nx_, ny_);
    for(unsigned int id=0; id<targets_.size();id++) {
        const Dat& dat = targets_[id];
        int kx = dat.x;
        int ky = dat.y;
        int i = specProp.key(kx, ky);
        complex_d tmp = output_[i];

        (*outTargets_[id]) <<
            format("%f %e %f")
            % time
            % abs(tmp) % arg(tmp)
                          << endl;
    }

    // for test
    double k1 = arg(output_[specProp.key(10, 2)]);
    double k2 = arg(output_[specProp.key(3, -2)]);
    double k3 = arg(output_[specProp.key( 8, 1)]);
    double k4 = arg(output_[specProp.key(5, -1)]);
    double fmismach = fmod(k1+k2-k3-k4, 2*M_PI);
    if(fmismach >  M_PI) fmismach -= 2*M_PI;
    if(fmismach < -M_PI) fmismach += 2*M_PI;
    whole_ <<
        format("%f %f")
        % time
        % fmismach
           << endl;
}

void Logger::close()
{
    for(unsigned int id=0; id<targets_.size(); id++) {
        outTargets_[id]->close();
        delete outTargets_[id];
    }
    whole_.close();
    isOpened_ = false;
}
