#include "mpiUtil.h"

void mpi_init(MpiParameter_t &p, int argc, char ** argv)
{
#ifdef SBE_WH_MPI
    MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &p.threadModeProvided);
    MPI_Comm_size(MPI_COMM_WORLD, &p.size);
    MPI_Comm_rank(MPI_COMM_WORLD, &p.rank);
    p.root = 0;
#endif
}

void mpi_finalize()
{
#ifdef SBE_WH_MPI
    MPI_Finalize();
#endif
}


void mpi_bcast(std::string &s, const MpiParameter_t &p)
{
#ifdef SBE_WH_MPI
    int size = s.size();
    MPI_Bcast(&size, 1, MPI_INT, p.root, MPI_COMM_WORLD);
    if ( p.rank != p.root )
        s.resize(size);
    if ( size > 0 )
        MPI_Bcast(reinterpret_cast<char*>(s.data()), size, MPI_CHAR, p.root, MPI_COMM_WORLD);
#endif
}


void mpi_bcast(std::vector<std::string> &v, const MpiParameter_t &p)
{
#ifdef SBE_WH_MPI
    int size = (int)v.size();
    MPI_Bcast(&size, 1, MPI_INT, p.root, MPI_COMM_WORLD);
    if ( p.rank != p.root )
        v.resize(size);
    for(unsigned i=0; i<v.size(); i++)
        mpi_bcast(v[i], p);
#endif
}
