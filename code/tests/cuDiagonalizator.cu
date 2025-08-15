#include "gtest/gtest.h"
#include "CuDiagonalizator.h"

#include<random>
#include<vector>
#include<complex>

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

std::vector<std::complex<double>> createRandomHermitianMatrices(unsigned N, unsigned matCount, unsigned offset=0)
{
    std::random_device rd;
    std::minstd_rand0 mt(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    std::vector<std::complex<double>> A(offset + matCount*N*N, -42);
    for(unsigned mc=0; mc<matCount; mc++)
        for(unsigned a=0; a<N; a++){
            A[mc*N*N+a*N+a + offset] = dist(mt);
            for(unsigned b=a+1; b<N; b++){
                A[mc*N*N+a*N+b + offset] = std::complex<double>(dist(mt), dist(mt));
                A[mc*N*N+b*N+a + offset] = std::conj( A[mc*N*N+a*N+b + offset] );
            }
        }
    return A;
}


TEST(CuDiagonalizator, JacobiStride){
    const int m = 25;
    const int batchSize = 5;
    const int offset = 10;
    const double eps = 1e-13;
    const int sweepCount = 10;
    std::vector< std::complex<double> > A = createRandomHermitianMatrices(m, batchSize, offset);
    std::vector< std::complex<double> > U(A.size(), 0);

    cudaStream_t stream = nullptr;
    cuDoubleComplex * d_A = nullptr;
    CUDA_CHECK(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));

    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_A), sizeof(std::complex<double>) * A.size()));
    CUDA_CHECK(
        cudaMemcpyAsync(d_A, A.data(), sizeof(std::complex<double>) * A.size(), cudaMemcpyHostToDevice, stream));
    {
        CuDiagonalizator cud(stream, batchSize, m, offset);
        cud.setTolerance(eps);
        cud.setSweepCount(sweepCount);
        cud.diag(d_A);
        std::vector<int> info = cud.getDiagInfo();
        for(int i : info)
            EXPECT_EQ(i, 0);
        CUDA_CHECK(
            cudaMemcpyAsync(U.data(), d_A, sizeof(std::complex<double>) * U.size(), cudaMemcpyDeviceToHost, stream));
        std::vector<double> evs = cud.getEigenValues();
        for(unsigned i=0; i<batchSize; i++){
            std::complex<double> *Ui = U.data() + offset + i * m*m;
            std::complex<double> *Ai = A.data() + offset + i * m*m;
            double * ei = evs.data() + i * m;
            EXPECT_NEAR( orthogonalUU(Ui, m), 0, 1e-14 );
            std::vector<std::complex<double>> R = transformA(Ai, Ui, m);
            for(unsigned a=0; a<m; a++)
                for(unsigned b=0; b<m; b++){
                    if ( a == b ){
                        EXPECT_NEAR(ei[a], std::real(R[a*m+b]), m * eps);
                        EXPECT_NEAR(0, std::imag(R[a*m+b]),  1e-14);
                    } else {
                        EXPECT_NEAR(0, std::real(R[a*m+b]), m * eps);
                        EXPECT_NEAR(0, std::imag(R[a*m+b]), m * eps);
                    }
                }
        }
    }

    CUDA_CHECK(cudaStreamDestroy(stream));
    CUDA_CHECK(cudaDeviceReset());
}


}
