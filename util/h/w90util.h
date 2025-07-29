#ifndef SBE_WH_W90_UTIL_H
#define SBE_WH_W90_UTIL_H

#include<stdio.h>
#include<cassert>
#include<string>
#include<complex>
#include<map>
#include<unordered_map>
#include<memory>

#include "GeomVector.hpp"
#include "unitConversion.h"


class W90_mat2D{
    unsigned colSize, rowSize;
    std::vector< std::complex<double> > mat;
public:
    W90_mat2D(unsigned size) : colSize(size), rowSize(size), mat(size*size, 0) { }
    W90_mat2D(unsigned rowSize, unsigned colSize) : colSize(colSize), rowSize(rowSize), mat(colSize*rowSize, 0) { }
    inline unsigned getColSize() const { return colSize; }
    inline unsigned getRowSize() const { return rowSize; }
    bool symmetrize();
    bool symmetrize(W90_mat2D &m); // both matrices are affected
    void toReal();
    inline std::complex<double>& operator()(unsigned a, unsigned b) { return mat[a*colSize+ b]; }
    inline std::complex<double> operator()(unsigned a, unsigned b) const { return mat.at(a*colSize + b); }
    inline std::complex<double> * data(){ return mat.data(); }
};

struct W90_wignerSeitzCell{
    W90_wignerSeitzCell(unsigned numWann) : H(numWann), D{numWann, numWann, numWann}, degeneracy(0)  {}
    W90_mat2D H;
    W90_mat2D D[3];
    CellIndex r;
    unsigned degeneracy;
};

struct W90_realSpaceOperators {
    W90_realSpaceOperators(unsigned numWann) : H(numWann), D{numWann, numWann, numWann} {}
    W90_mat2D H;
    W90_mat2D D[3];
};

// uses always atomic units
struct W90_tb {
    bool valid;
    std::string header;
    GeomVector3d a[3];
    unsigned numWann;
    std::vector<W90_wignerSeitzCell> wsCells;
};


struct W90_wsvecCell {
    CellIndex r;
    unsigned n; // wannier index
    unsigned m; // wannier index
    std::vector<CellIndex> T;
};

struct W90_wsvec {
    bool valid = false;
    std::string header;
    std::vector<W90_wsvecCell> cells;
};

struct W90_Umat {
    GeomVector3d k;
    W90_mat2D U;
};


W90_tb w90_read_tb(const char *fname, bool verbose=false);
W90_wsvec w90_read_wsvec(const char *fname, bool verbose=false);
bool w90_symmetrize(W90_tb &wtb, bool verbose=false);

// may be used to read U_k and U^disentangle_k.
// returns read data, success
std::pair< std::vector<W90_Umat>, bool> w90_read_Umat(const char *fname, bool verbose=false);

std::unordered_map<CellIndex, W90_realSpaceOperators> w90_calcRealSpaceOperators(const W90_tb &tb);
std::unordered_map<CellIndex, W90_realSpaceOperators> w90_calcRealSpaceOperators(const W90_tb &tb, const W90_wsvec &wsvec);

bool w90_symmetrizeRealSpaceOperators(std::unordered_map<CellIndex, W90_realSpaceOperators> &rso, bool discardImag=false);
void w90_discardImagOfRealSpaceOperators(std::unordered_map<CellIndex, W90_realSpaceOperators> &rso);
bool w90_projectRealSpaceOperators(std::unordered_map<CellIndex, W90_realSpaceOperators> &rso, unsigned dim);


#endif
