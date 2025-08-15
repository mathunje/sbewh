#ifndef SBE_WH_MPI_PARAMETER_H
#define SBE_WH_MPI_PARAMETER_H

struct MpiParameter_t {
    int rank;
    int size;
    int root;
    int threadModeProvided;
    MpiParameter_t() { rank = 0; size = 1; root = 0; threadModeProvided = 0; }
};

#endif
