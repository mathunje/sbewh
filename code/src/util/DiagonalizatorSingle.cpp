#include "DiagonalizatorSingle.h"

DiagonalizatorSingle::DiagonalizatorSingle(DiagMode dm, unsigned matSize) :
    diagMode(dm), matSize(matSize)
{
    switch(diagMode){
        case jacobi: initJacobi(); break;
#ifdef SBE_WH_LAPACK
        case lzheev: initZheev(); break;
        case lzheevd: initZheevd(); break;
#endif
    }
}

int DiagonalizatorSingle::diag(std::complex<double> *A, double * eigenValues)
{
    switch(diagMode){
        case jacobi: return diagJacobi(A, eigenValues);
#ifdef SBE_WH_LAPACK
        case lzheev: return diagZheev(A, eigenValues);
        case lzheevd: return diagZheevd(A, eigenValues);
#endif
        default: return -1;
    }
}

double DiagonalizatorSingle::getError(const std::complex<double> *U, const std::complex<double> * origA) const
{
    double diagSum = 0;
    double offDiagSum = 0;
    std::vector<std::complex<double>> temp(matSize * matSize, 0);
    for(unsigned a=0; a<matSize; a++)
        for(unsigned b=0; b<matSize; b++)
            for(unsigned c=0; c<matSize; c++)
                temp[a*matSize+b] += U[a*matSize+c] * origA[c*matSize+b];
    for(unsigned a=0; a<matSize; a++)
        for(unsigned b=0; b<matSize; b++){
            std::complex<double> res_ab = 0;
            for(unsigned c=0; c<matSize; c++)
                res_ab += temp[a*matSize+c] * std::conj(U[b*matSize+c]);
            if (a == b){
                diagSum += std::abs(res_ab);
            } else {
                offDiagSum += std::abs(res_ab);
            }
        }
    return offDiagSum / diagSum;
}

/************************* JACOBI *************************/

void DiagonalizatorSingle::initJacobi()
{
    works.resize(matSize * matSize);
}

inline std::pair<double, double> DiagonalizatorSingle::csHalf(double y, double x)
{
    double t = y / ( fabs(x) + sqrt(x*x + y*y) );
    double c = 1 / sqrt(1 + t * t);
    double s = t * c;
    if (x >= 0)
        return {c, s};
    return { fabs(s), copysign(c, y) };
}

inline void DiagonalizatorSingle::jacobiRotate(unsigned p, unsigned q, std::complex<double> *A, std::complex<double> *U)
{
    std::complex<double> Apq = A[p * matSize + q];
    auto [c1, s1] = csHalf( -std::real(Apq), std::imag(Apq) );
    auto [c2, s2] = csHalf( 2 * abs(Apq), std::real(A[p*(matSize+1)])-std::real(A[q*(matSize+1)]) );
    std::complex<double> Rpp( -s1*s2-1, -c1*s2 );
    std::complex<double> Rpq( -c1*c2, -s1*c2 );
    std::complex<double> Rqp = -std::conj(Rpq);
    std::complex<double> Rqq = std::conj(Rpp);
    for(unsigned i=0; i<matSize; i++){
        std::complex<double> dUp = Rpp * U[p*matSize+i] + Rpq * U[q*matSize+i];
        std::complex<double> dUq = Rqp * U[p*matSize+i] + Rqq * U[q*matSize+i];
        U[p * matSize + i] += dUp;
        U[q * matSize + i] += dUq;
    }
    // accesses full matrix -- can be improved to access only upper triangular part
    for(unsigned i=0; i<matSize; i++){
        std::complex<double> dApq = Rpp * A[p*matSize+i] + Rpq * A[q*matSize+i];
        std::complex<double> dAqp = Rqp * A[p*matSize+i] + Rqq * A[q*matSize+i];
        A[p * matSize + i] += dApq;
        A[q * matSize + i] += dAqp;
    }
    for(unsigned i=0; i<matSize; i++){
        std::complex<double> dApq = Rqq * A[i*matSize+p] - Rqp * A[i*matSize+q];
        std::complex<double> dAqp = -Rpq * A[i*matSize+p] + Rpp * A[i*matSize+q];
        A[i*matSize+p] += dApq;
        A[i*matSize+q] += dAqp;
    }
    A[p * matSize + q] = 0;
    A[q * matSize + p] = 0;
}

int DiagonalizatorSingle::diagJacobi(std::complex<double> * A, double * eigVals)
{
    std::complex<double> * U = works.data();
    for(unsigned a=0; a<matSize; a++)
        for(unsigned b=0; b<matSize; b++)
            U[a*matSize+b] = (a==b);
    for(unsigned sweep=0; sweep<sweepCount; sweep++){
        double sum = 0;
        for(unsigned i=0; i<matSize; i++)
            for(unsigned j=i; j<matSize; j++)
                sum += fabs(std::real(A[i*matSize+j])) + fabs(std::imag(A[i*matSize+j]));
        sum /= 2;
        for(unsigned a=0; a<matSize; a++){
            for(unsigned b=a+1; b<matSize; b++)
                if (std::abs(A[a*matSize+b]) > eps * sum)
                    jacobiRotate(a, b, A, U);
        }
    }
    for(unsigned i=0; i<matSize; i++)
        eigVals[i] = std::real(A[i*(matSize+1)]);
    for(unsigned i=0; i<matSize*matSize; i++)
        A[i] = U[i];
    return 0;
}

#ifdef SBE_WH_LAPACK
/************************* ZHEEV **************************/

void DiagonalizatorSingle::initZheev()
{
    lapack_complex_double work_query;
    int queryInfo = LAPACKE_zheev_work(LAPACK_COL_MAJOR, 'V', 'U', matSize, NULL,
                  matSize, NULL, &work_query, -1, NULL);
    assert(queryInfo == 0);
    int lwork = work_query.real;
    int lrwork = 3 * matSize - 2;
    works.resize(lwork);
    rworks.resize(lrwork);
}


int DiagonalizatorSingle::diagZheev(std::complex<double> * A, double * eigVals)
{
    int lwork = works.size();
    int info = LAPACKE_zheev_work(LAPACK_COL_MAJOR, 'V', 'U', matSize, reinterpret_cast<lapack_complex_double*>(A), matSize,
                       eigVals,
                       reinterpret_cast<lapack_complex_double*>(works.data()), lwork, rworks.data() );
    return info;
}


/************************* ZHEEVD *************************/

void DiagonalizatorSingle::initZheevd()
{
    lapack_complex_double work_query({-1.0, 0.0});
    double work_query_r = -1;
    lapack_int work_query_i = -1;
    int queryInfo = LAPACKE_zheevd_work(LAPACK_COL_MAJOR, 'V', 'U', matSize, NULL, matSize, NULL,
                                        &work_query, -1, &work_query_r, -1, &work_query_i, -1);
    assert(queryInfo == 0);
    int lwork = work_query.real;
    int lrwork = work_query_r;
    int liwork = (int)work_query_i;
    works.resize(lwork);
    rworks.resize(lrwork);
    iworks.resize(liwork);
}

int DiagonalizatorSingle::diagZheevd(std::complex<double> *A, double * eigVals)
{
    int lwork = works.size();
    int lrwork = rworks.size();
    int liwork = iworks.size();
    int info = LAPACKE_zheevd_work(LAPACK_COL_MAJOR, 'V', 'U', matSize,
                       reinterpret_cast<lapack_complex_double*>(A), matSize,
                       eigVals,
                       reinterpret_cast<lapack_complex_double*>(works.data()), lwork,
                       rworks.data(), lrwork,
                       iworks.data(), liwork);
    return info;
}

#endif // SBE_WH_LAPACK
