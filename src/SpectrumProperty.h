#ifndef SPECTRUM_PROPERTY__
#define SPECTRUM_PROPERTY__

#include <boost/numeric/ublas/vector.hpp>
#include <map>
#include <string>
#include "Define.h"

class SpectrumProperty
{
public:
    explicit SpectrumProperty(int nx, int ny, double lx, double ly, int M);
    const ubvector_d& klist() const { return klist_; }
//    const ubvector_d& kxlist() const { return norm_kxlist_ * dkx_; }
//    const ubvector_d& kylist() const { return norm_kylist_ * dky_; }
    const ubvector_i& norm_kxlist() const { return norm_kxlist_; }
    const ubvector_i& norm_kylist() const { return norm_kylist_; }
    const ubvector_d& low_pass_filter() const { return low_pass_filter_; }
    const ubvector_i& id_sort_kxky() const { return id_sort_kxky_; }
    const double dkx() const { return dkx_; }
    const double dky() const { return dky_; }
    int key(int kx, int ky) const { return keymap_.find(keygen_(kx, ky))->second; }

    // for real to complex fft
    const ubvector_d& klist2() const { return klist2_; }
    const ubvector_i& norm_kxlist2() const { return norm_kxlist2_; }
    const ubvector_i& norm_kylist2() const { return norm_kylist2_; }
    const ubvector_d& low_pass_filter2() const { return low_pass_filter2_; }
    const ubvector_i& id_sort_kxky2() const { return id_sort_kxky2_; }
    int key2(int kx, int ky) const { return keymap2_.find(keygen_(kx, ky))->second; }

private:
    int ny_;
    ubvector_d klist_;
    ubvector_i norm_kxlist_;
    ubvector_i norm_kylist_;
    ubvector_d low_pass_filter_;
    ubvector_i id_sort_kxky_;

    ubvector_d klist2_;
    ubvector_i norm_kxlist2_;
    ubvector_i norm_kylist2_;
    ubvector_d low_pass_filter2_;
    ubvector_i id_sort_kxky2_;
    double dkx_, dky_;
    std::map<std::string, int> keymap_;
    std::map<std::string, int> keymap2_;
    std::string keygen_(int kx, int ky) const;

};

#endif /* SPECTRUM_PROPERTY__ */
