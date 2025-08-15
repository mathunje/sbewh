#include "gtest/gtest.h"
#include "DiagonalizatorMultiple.h"

#include<random>

namespace{

double orthogonalUU(const std::complex<double> *U, unsigned matSize){
    double maxErr = 0;
    for(unsigned a=0; a<matSize; a++)
        for(unsigned b=0; b<matSize; b++){
            std::complex<double> p = 0;
            for(unsigned c=0; c<matSize; c++)
                p += std::conj(U[c*matSize+a]) * U[c*matSize+b];
            double cErr = std::abs(p - (double)(a==b) );
            if ( cErr > maxErr )
                maxErr = cErr;
        }
    return maxErr;
}

std::vector<std::complex<double>> transformA(const std::complex<double> *A, const std::complex<double> *U, unsigned matSize){
    std::vector<std::complex<double>> temp(matSize * matSize, 0);
    std::vector<std::complex<double>> res(matSize * matSize, 0);
    for(unsigned a=0; a<matSize; a++)
        for(unsigned b=0; b<matSize; b++)
            for(unsigned c=0; c<matSize; c++)
                temp[a*matSize+b] += U[a*matSize+c] * A[c*matSize+b];
    for(unsigned a=0; a<matSize; a++)
        for(unsigned b=0; b<matSize; b++)
            for(unsigned c=0; c<matSize; c++)
                res[a*matSize+b] += temp[a*matSize+c] * std::conj(U[b*matSize+c]);
    return res;
}

std::vector<std::complex<double>> createRandomHermitianMatrices(unsigned N, unsigned matCount, unsigned stride, unsigned offset=0)
{
    std::random_device rd;
    std::minstd_rand0 mt(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    std::vector<std::complex<double>> A(matCount*stride, -42);
    for(unsigned mc=0; mc<matCount; mc++)
        for(unsigned a=0; a<N; a++){
            A[mc*stride+a*N+a + offset] = dist(mt);
            for(unsigned b=a+1; b<N; b++){
                A[mc*stride+a*N+b + offset] = std::complex<double>(dist(mt), dist(mt));
                A[mc*stride+b*N+a + offset] = std::conj( A[mc*stride+a*N+b + offset] );
            }
        }
    return A;
}

TEST(Diagonalizator, JacobiStride){
    unsigned matCount = 1000;
    unsigned stride = 11;
    unsigned matSize = 2;
    unsigned offset = 1;
    std::vector<std::complex<double>> A(matCount*stride);
    for(unsigned i=0; i<matCount; i++){
        A[stride*i + offset] = 0;
        A[stride*i+1 + offset] = i;
        A[stride*i+2 + offset] = i;
        A[stride*i+3 + offset] = 0;
    }
    DiagonalizatorMultiple dm(DiagMode::jacobi, matCount, matSize, stride, offset);
    dm.diag(A.data());
    for(unsigned i=0; i<matCount; i++){
        auto ev = dm.getEigenValues(i);
        for(unsigned m=0; m<matSize; m++)
            EXPECT_NEAR( std::abs(ev[m]), i, 1e-10 );
    }
}

TEST(Diagonalizator, JacobiRandom){
    unsigned N = 8;
    unsigned stride = N * N;
    unsigned matCount = 22;
    unsigned sweepCount = 10;
    std::vector<std::complex<double>> A = createRandomHermitianMatrices(N, matCount, stride);
    std::vector<std::complex<double>> Aorig = A;
    DiagonalizatorMultiple dm(DiagMode::jacobi, matCount, N, stride);
    for(unsigned i=1; i<=sweepCount; i++){
        std::copy(Aorig.begin(), Aorig.end(), A.begin());
        dm.setSweepCount(i);
        dm.diag(A.data());
        printf("Max error after %u sweeps: %lg\n", i, dm.getMaxError(&Aorig[0]));
    }
    EXPECT_NEAR(dm.getMaxError(&Aorig[0]), 0, 1e-14);
    for(unsigned mc=0; mc<matCount; mc++){
        std::complex<double> * dA = A.data() + mc * stride;
        EXPECT_NEAR(orthogonalUU(dA, N), 0, 1e-14);
        auto uAu = transformA(&Aorig[mc*stride], dA, N);
        const double * ev = dm.getEigenValues(mc);
        for(unsigned a=0; a<N; a++)
            for(unsigned b=0; b<N; b++){
                if ( a==b ){
                    EXPECT_NEAR( std::real(uAu[a*N+b]), ev[a], 1e-13 );
                } else {
                    EXPECT_NEAR( std::real(uAu[a*N+b]), 0, 1e-13 );
                }
                EXPECT_NEAR( std::imag(uAu[a*N+b]), 0, 1e-13 );
            }
    }
}

#ifdef SFP_LAPACK

TEST(Diagonalizator, zheev){
    unsigned N = 8;
    unsigned matCount = 2;
    unsigned stride = 128;
    unsigned offset = 12;
    std::vector<std::complex<double>> A = createRandomHermitianMatrices(N, matCount, stride, offset);
    std::vector<std::complex<double>> Aorig = A;
    DiagonalizatorMultiple dm(DiagMode::zheev, matCount, N, stride, offset);
    ASSERT_EQ(dm.diag(A.data()), 0);
    for(unsigned mc=0; mc<matCount; mc++){
        std::complex<double> * dA = A.data() + mc * stride + offset;
        EXPECT_NEAR(orthogonalUU(dA, N), 0, 1e-14);
        auto uAu = transformA(&Aorig[mc*stride+offset], dA, N);
        const double * ev = dm.getEigenValues(mc);
        for(unsigned a=0; a<N; a++)
            for(unsigned b=0; b<N; b++){
                if ( a==b ){
                    EXPECT_NEAR( std::real(uAu[a*N+b]), ev[a], 1e-13 );
                } else {
                    EXPECT_NEAR( std::real(uAu[a*N+b]), 0, 1e-13 );
                }
                EXPECT_NEAR( std::imag(uAu[a*N+b]), 0, 1e-13 );
            }
    }
}

TEST(Diagonalizator, zheevd){
    unsigned N = 8;
    unsigned matCount = 2;
    unsigned stride = 128;
    unsigned offset = 12;
    std::vector<std::complex<double>> A = createRandomHermitianMatrices(N, matCount, stride, offset);
    std::vector<std::complex<double>> Aorig = A;
    DiagonalizatorMultiple dm(DiagMode::zheevd, matCount, N, stride, offset);
    ASSERT_EQ(dm.diag(A.data()), 0);
    for(unsigned mc=0; mc<matCount; mc++){
        std::complex<double> * dA = A.data() + mc * stride + offset;
        EXPECT_NEAR(orthogonalUU(dA, N), 0, 1e-14);
        auto uAu = transformA(&Aorig[mc*stride+offset], dA, N);
        const double * ev = dm.getEigenValues(mc);
        for(unsigned a=0; a<N; a++)
            for(unsigned b=0; b<N; b++){
                if ( a==b ){
                    EXPECT_NEAR( std::real(uAu[a*N+b]), ev[a], 1e-13 );
                } else {
                    EXPECT_NEAR( std::real(uAu[a*N+b]), 0, 1e-13 );
                }
                EXPECT_NEAR( std::imag(uAu[a*N+b]), 0, 1e-13 );
            }
    }
}

#endif


}
