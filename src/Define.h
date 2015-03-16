#ifndef DEFINE_H__
#define DEFINE_H__

#include <complex>
#include <boost/numeric/ublas/vector.hpp>

typedef std::complex<double> complex_d;
typedef boost::numeric::ublas::vector<int> ubvector_i;
typedef boost::numeric::ublas::vector<double, std::vector<double> > ubvector_d;
typedef boost::numeric::ublas::vector<complex_d, std::vector<complex_d> > ubvector_c;


#endif /* DEFINE_H__ */
