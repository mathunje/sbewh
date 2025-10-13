#ifndef SBE_WH_MPI_UTIL_H
#define SBE_WH_MPI_UTIL_H

#include<string>
#include<array>
#include<vector>
#include<complex>

#include<mpi.h>

#include "sbewhConfig.h"
#include "MpiParameter.h"


void mpi_init(MpiParameter_t &p, int agrc, char **argv);
void mpi_finalize();




template <typename T> void mpi_bcast(T&, const MpiParameter_t&) = delete;

void mpi_bcast(std::string &s, const MpiParameter_t &p);
void mpi_bcast(std::vector<std::string> &v, const MpiParameter_t &p);


template<typename T> void mpi_bast_trivial(T&, const MpiParameter_t&) = delete;

template<typename T>
void mpi_bcast_trivial(std::vector<T> &v, const MpiParameter_t &p)
{
#ifdef SBE_WH_MPI
    int size = int(v.size());
    MPI_Bcast(&size, 1, MPI_INT, p.root, MPI_COMM_WORLD);
    if ( p.rank != p.root )
        v.resize(size);
    MPI_Bcast(reinterpret_cast<char*>(v.data()), size*sizeof(T), MPI_CHAR, p.root, MPI_COMM_WORLD);
#endif
}

template<typename T, long unsigned int N>
void mpi_bcast_trivial(std::array<T, N> &v, const MpiParameter_t &p)
{
#ifdef SBE_WH_MPI
    MPI_Bcast(reinterpret_cast<char*>(v.data()), N*sizeof(T), MPI_CHAR, p.root, MPI_COMM_WORLD);
#endif
}

template<typename T>
void mpi_bcast_trivial(T &v, const MpiParameter_t &p)
{
#ifdef SBE_WH_MPI
    T s = v;
    MPI_Bcast(reinterpret_cast<char*>(&s), sizeof(T), MPI_CHAR, p.root, MPI_COMM_WORLD);
    v = s;
#endif
}


template<typename T> void mpi_reduce(T*, size_t, int, MPI_Op) = delete;
inline void mpi_reduce(unsigned * data, size_t N, const MpiParameter_t &p, MPI_Op op){
#ifdef SBE_WH_MPI
    if ( p.rank == p.root ){
        MPI_Reduce(MPI_IN_PLACE, data, N, MPI_UNSIGNED, op, p.root, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(data, data, N, MPI_UNSIGNED, op, p.root, MPI_COMM_WORLD);
    }
#endif
}
inline void mpi_reduce(int *data, size_t N, const MpiParameter_t &p, MPI_Op op){
#ifdef SBE_WH_MPI
    if ( p.rank == p.root ){
        MPI_Reduce(MPI_IN_PLACE, data, N, MPI_INT, op, p.root, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(data, data, N, MPI_INT, op, p.root, MPI_COMM_WORLD);
    }
#endif
}
inline void mpi_reduce(float *data, size_t N, const MpiParameter_t &p, MPI_Op op){
#ifdef SBE_WH_MPI
    if ( p.rank == p.root ){
        MPI_Reduce(MPI_IN_PLACE, data, N, MPI_FLOAT, op, p.root, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(data, data, N, MPI_FLOAT, op, p.root, MPI_COMM_WORLD);
    }
#endif
}
inline void mpi_reduce(double *data, size_t N, const MpiParameter_t &p, MPI_Op op){
#ifdef SBE_WH_MPI
    if ( p.rank == p.root ){
        MPI_Reduce(MPI_IN_PLACE, data, N, MPI_DOUBLE, op, p.root, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(data, data, N, MPI_DOUBLE, op, p.root, MPI_COMM_WORLD);
    }
#endif
}

template<typename T>
inline void mpi_reduce(T &v, const MpiParameter_t &p, MPI_Op op){
#ifdef SBE_WH_MPI
    T s = v;
    mpi_reduce(&s, 1, p, op);
    v = s;
#endif
}

template<typename T>
inline void mpi_reduce(std::complex<T> &v, const MpiParameter_t &p, MPI_Op op){
#ifdef SBE_WH_MPI
    T s = v;
    mpi_reduce(reinterpret_cast<T*>(&s), 2, p, op);
    v = s;
#endif
}


template<typename T>
inline void mpi_reduce(std::vector<T> &v, const MpiParameter_t &p, MPI_Op op){
#ifdef SBE_WH_MPI
    mpi_reduce(v.data(), v.size(), p, op);
#endif
}

template<typename T>
inline void mpi_reduce(std::vector<std::complex<T> > &v, const MpiParameter_t &p, MPI_Op op){
#ifdef SBE_WH_MPI
    mpi_reduce(reinterpret_cast<T*>(v.data()), 2*v.size(), p, op);
#endif
}

template<typename T, size_t N>
inline void mpi_reduce(std::array<T, N> &v, const MpiParameter_t &p, MPI_Op op){
#ifdef SBE_WH_MPI
    mpi_reduce(v.data(), v.size(), p, op);
#endif
}
template<typename T, size_t N>
inline void mpi_reduce(std::array<std::complex<T>, N> &v, const MpiParameter_t &p, MPI_Op op){
#ifdef SBE_WH_MPI
    mpi_reduce(reinterpret_cast<T*>(v.data()), 2*v.size(), p, op);
#endif
}

template<typename T>
inline void mpi_sum(T &v, const MpiParameter_t &p){ mpi_reduce(v, p, MPI_SUM); }
template<typename T>
inline void mpi_min(T &v, const MpiParameter_t &p){ mpi_reduce(v, p, MPI_MIN); }
template<typename T>
inline void mpi_max(T &v, const MpiParameter_t &p){ mpi_reduce(v, p, MPI_MAX); }

#endif
