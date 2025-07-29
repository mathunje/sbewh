#ifndef SBE_WH_OMP_REDUCTIONS_H
#define SBE_WH_OMP_REDUCTIONS_H

#include<complex.h>


inline void sfp_omp_complex_add(std::complex<double> &lhs, const std::complex<double> &rhs){ lhs += rhs; }

#pragma omp declare reduction ( + :  std::complex<double> : sfp_omp_complex_add(omp_out, omp_in) ) initializer(omp_priv={0, 0})

#endif
