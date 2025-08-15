#include "GeomVector.hpp"

void mpi_bcast(GeomVector3d &v, const MpiParameter_t &p)
{
#ifdef SBE_WH_MPI
    std::array<double, 3> data = v.getArray();
    mpi_bcast_trivial(data, p);
    if ( p.rank != p.root )
        v = GeomVector3d(data);
#endif
}

void mpi_bcast(CellIndex &v, const MpiParameter_t &p)
{
#ifdef SBE_WH_MPI
    std::array<int, 3> data = v.getArray();
    mpi_bcast_trivial(data, p);
    if ( p.rank != p.root )
        v = CellIndex(data);
#endif
}
