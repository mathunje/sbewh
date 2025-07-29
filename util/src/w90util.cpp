#include "w90util.h"

bool W90_mat2D::symmetrize(){
    if (colSize != rowSize)
        return false;
    for(unsigned i=0; i<colSize; i++){
        (*this)(i, i) = std::real((*this)(i, i));
        for(unsigned j=i+1; j<colSize; j++){
            auto nval = 0.5 * ((*this)(i,j) + std::conj((*this)(j,i)));
            (*this)(i, j) = nval;
            (*this)(j, i) = std::conj(nval);
        }
    }
    return true;
}

bool W90_mat2D::symmetrize(W90_mat2D &m){
    if (colSize != m.getRowSize() || rowSize != m.getColSize())
        return false;
    for(unsigned r=0; r<rowSize; r++)
        for(unsigned c=0; c<colSize; c++){
            auto nval = 0.5 * ((*this)(r,c) + std::conj(m(c,r)));
            (*this)(r, c) = nval;
            m(c, r) = std::conj(nval);
        }
    return true;
}


void W90_mat2D::toReal(){
    for(unsigned r=0; r<rowSize; r++)
        for(unsigned c=0; c<colSize; c++)
            (*this)(r, c) = std::real( (*this)(r, c) );
}


/********************************************************/


W90_tb w90_read_tb(const char *fname, bool verbose)
{
    W90_tb res;
    res.valid = false;
    FILE * f = fopen(fname, "r");
    if ( !f ){
        if ( verbose )
            printf("Could not open file %s\n", fname);
        return res;
    }
    const int bufSize = 1024;
    char buf[bufSize];
    // read  header line
    if( ! fgets(buf, sizeof(buf), f) ){
        if ( verbose )
            printf("Could not read header\n");
        fclose(f);
        return res;
    }
    res.header = std::string(buf);
    for(unsigned i=0; i<3; i++){
        double x[3];
        if ( fscanf(f, "%lf %lf %lf", x, x+1, x+2) != 3){
            if ( verbose )
                printf("Could not read lattice vector %u\n", i);
            fclose(f);
            return res;
        }
        for(int k=0; k<3; k++)
            x[k] = atomicUnits::from_A(x[k]);
        res.a[i] = GeomVector3d(x);
    }
    unsigned nrpts;
    if ( fscanf(f, "%u %u", &res.numWann, &nrpts) != 2 ){
        if ( verbose )
            printf("Could not read 'num_wann nrpts'\n");
        fclose(f);
        return res;
    }
    res.wsCells.resize(nrpts, res.numWann);
    for(unsigned i=0; i<nrpts; i++){
        if (fscanf(f, "%u", &res.wsCells[i].degeneracy) != 1){
            if ( verbose )
                printf("Could not read degeneracy for %u. wigner seitz cell point\n", i);
            fclose(f);
            return res;
        }
    }
    // read <0m|H|nR>
    for(unsigned i=0; i<nrpts; i++){
        if ( fscanf(f, "%d %d %d", &res.wsCells[i].r[0], &res.wsCells[i].r[1], &res.wsCells[i].r[2]) != 3){
            if ( verbose )
                printf("Could not read R vector of %u. WS cell during read of H\n", i);
            fclose(f);
            return res;
        }
        for(int n=0; n<res.numWann; n++)
            for(int m=0; m<res.numWann; m++){
                unsigned mread, nread;
                double hrReal_eV, hrImag_eV;
                if ( fscanf(f, "%d %d %lf %lf", &mread, &nread, &hrReal_eV, &hrImag_eV) != 4){
                    if ( verbose )
                        printf("Coud not read  <0n| H | Rm> for R=(%d, %d, %d) at %u. line\n",
                                res.wsCells[i].r[0], res.wsCells[i].r[1], res.wsCells[i].r[2],
                                n * res.numWann + m);
                    fclose(f);
                    return res;
                }
                if ( mread == 0 || mread > res.numWann || nread == 0 || nread > res.numWann){
                    if (verbose)
                        printf("Invalid wannierIndices %u %u\n", mread, nread);
                    fclose(f);
                    return res;
                }
                res.wsCells[i].H(mread-1, nread-1) = {atomicUnits::from_eV(hrReal_eV), atomicUnits::from_eV(hrImag_eV)};
            }
    }
    // read dipoles
    for(int i=0; i<nrpts; i++){
        int ciRaw[3];
        if ( fscanf(f, "%d %d %d", ciRaw, ciRaw+1, ciRaw+2) != 3){
            if ( verbose )
                printf("Could not read R vector of %u. WS cell during read of D\n", i);
            fclose(f);
            return res;
        }
        for(unsigned k=0; k<3; k++)
            if ( ciRaw[k] != res.wsCells[i].r[k] ){
                if (verbose )
                    printf("Order of R if different for D and H for %i. wigner seitz cell point\n", i);
                fclose(f);
                return res;
            }
        for(int n=0; n<res.numWann; n++)
            for(int m=0; m<res.numWann; m++){
                int mread, nread;
                double rxReal, rxImag, ryReal, ryImag, rzReal, rzImag;
                if ( fscanf(f, "%d %d %lf %lf %lf %lf %lf %lf", &mread, &nread,
                                            &rxReal, &rxImag, &ryReal, &ryImag, &rzReal, &rzImag) != 8){
                    if ( verbose )
                        printf("Coud not read  <0n| D | Rm> for R=(%d, %d, %d) at %u. line\n",
                                res.wsCells[i].r[0], res.wsCells[i].r[1], res.wsCells[i].r[2],
                                n * res.numWann + m);
                    fclose(f);
                    return res;
                }
                if ( mread == 0 || mread > res.numWann || nread == 0 || nread > res.numWann){
                    if (verbose)
                        printf("Invalid wannierIndices %u %u\n", mread, nread);
                    fclose(f);
                    return res;
                }
                res.wsCells[i].D[0](mread-1, nread-1) = {atomicUnits::from_A(rxReal), atomicUnits::from_A(rxImag)};
                res.wsCells[i].D[1](mread-1, nread-1) = {atomicUnits::from_A(ryReal), atomicUnits::from_A(ryImag)};
                res.wsCells[i].D[2](mread-1, nread-1) = {atomicUnits::from_A(rzReal), atomicUnits::from_A(rzImag)};
            }
    }
    fclose(f);
    res.valid = true;
    return res;
}

// returns read data, success
W90_wsvec w90_read_wsvec(const char *fname, bool verbose)
{
    W90_wsvec res;
    res.valid = false;
    FILE * f = fopen(fname, "r");
    if ( !f ){
        if ( verbose )
            printf("Could not open file %s\n", fname);
        return res;
    }
    const int bufSize = 1024;
    char buf[bufSize];
    // read  header line
    if( ! fgets(buf, sizeof(buf), f) ){
        if ( verbose )
            printf("Could not read header\n");
        fclose(f);
        return res;
    }
    res.header = std::string(buf);
    int r[3];
    unsigned m, n;
    while( fscanf(f, "%d %d %d %u %u", &r[0], &r[1], &r[2], &m, &n) == 5){
        W90_wsvecCell cc;
        cc.r = CellIndex(r);
        cc.m = m-1;
        cc.n = n-1;
        unsigned cnt;
        if ( fscanf(f, "%u", &cnt) != 1){
           if ( verbose )
               printf("Could not read supercell offset multiciplicity\n");
           return res;
        }
        cc.T.resize(cnt);
        for(unsigned i=0; i<cnt; i++){
            if ( fscanf(f, "%d %d %d", &r[0], &r[1], &r[2]) != 3){
                if ( verbose )
                    printf("Could not read superlattice offset vector\n");
                return res;
            }
            cc.T[i] = CellIndex(r);
        }
        res.cells.push_back(cc);
    }
    res.valid = true;
    return res;
}


bool w90_symmetrize(W90_tb &wtb, bool verbose)
{
    std::map<GeomVector<int, 3>, std::pair<unsigned, bool> > adaptedCells;
    for(unsigned i=0; i< wtb.wsCells.size(); i++)
        adaptedCells[wtb.wsCells[i].r] = {i, false};
    for(auto &[ci, ac] : adaptedCells){
        auto [index, adapted] = ac;
        if ( adapted )
            continue;
        if ( ci.at(0) == 0 && ci.at(1) == 0 && ci.at(2) == 0) {
            adaptedCells[ci].second = true;
            if ( ! wtb.wsCells[index].H.symmetrize()){
                if ( verbose )
                    printf("Coud not symmetrize H(0, 0, 0)\n");
                return false;
            }
            for(unsigned dir=0; dir<3; dir++)
                wtb.wsCells[index].D[dir].symmetrize();
            continue;
        }
        adaptedCells[ci].second = true;
        auto it = adaptedCells.find(-ci);
        if (it == adaptedCells.end()){
            if (verbose )
                printf("Can not symmetrize hamiltonian as for R=(%d, %d, %d) no counterpart -R exists\n",
                        ci.at(0), ci.at(1), ci.at(2));
            return false;
        }
        it->second.second = true;
        unsigned jndex = it->second.first;
        if ( ! wtb.wsCells[index].H.symmetrize(wtb.wsCells[jndex].H)) {
            if ( verbose )
                printf("Can not symmetrize hamiltonian for R=(%d, %d, %d) and -R\n",
                        ci.at(0), ci.at(1), ci.at(2));
            return false;
        }
        for(unsigned dir=0; dir<3; dir++)
            wtb.wsCells[index].D[dir].symmetrize(wtb.wsCells[jndex].D[dir]);
    }
    return true;
}


std::pair< std::vector<W90_Umat>, bool> w90_read_Umat(const char *fname, bool verbose)
{
    std::vector<W90_Umat> res;
    FILE * f = fopen(fname, "r");
    if ( ! f )
        return { res, false};
    const int bufSize = 1024;
    char buf[bufSize];
    // read  header line
    if( ! fgets(buf, sizeof(buf), f) ){
        if ( verbose )
            printf("Could not read header of '%s'\n", fname);
        fclose(f);
        return {res, false};
    }
    unsigned num_kpts, num_wann, num_bands;
    if ( fscanf(f, "%u %u %u", &num_kpts, &num_wann, &num_bands) != 3){
        if ( verbose )
            printf("Could not read dimension definitions of '%s'\n", fname);
        fclose(f);
        return {res, false};
    }
    for(unsigned ki=0; ki<num_kpts; ki++){
        double rk[3];
        if ( fscanf(f, "%lf %lf %lf", rk, rk+1, rk+2) != 3){
            if ( verbose )
                printf("Could not read k vector of %u. unitary matrix of '%s'\n", ki, fname);
            fclose(f);
            return {res, false};
        }
        W90_mat2D u(num_bands, num_wann);
        for(unsigned n=0; n<num_wann; n++)
            for(unsigned b=0; b<num_bands; b++){
                double re, im;
                if ( fscanf(f, "%lf %lf", &re, &im) != 2){
                    if ( verbose )
                        printf("Could not read (%u, %u) complex number in %u. unitary matrix of '%s'\n", n, b, ki, fname);
                    fclose(f);
                    return {res, false};
                }
                u(b, n) = re + std::complex<double>(0, 1) * im;
            }
        res.push_back( {GeomVector3d(rk), u});
    }
    fclose(f);
    return {res, true};
}


std::unordered_map<CellIndex, W90_realSpaceOperators> w90_calcRealSpaceOperators(const W90_tb &tb)
{
    std::unordered_map<CellIndex, W90_realSpaceOperators> res;
    assert ( tb.valid );
    unsigned numWann = tb.numWann;
    for(const W90_wignerSeitzCell &wsc : tb.wsCells){
        W90_realSpaceOperators rso(numWann);
        for(unsigned i=0; i<numWann; i++)
            for(unsigned j=0; j<numWann; j++){
                rso.H(i, j) = wsc.H(i, j) / (double)wsc.degeneracy;
                for(unsigned dir=0; dir<3; dir++)
                    rso.D[dir](i, j) = wsc.D[dir](i, j) / (double)wsc.degeneracy;
            }
        res.insert({wsc.r, rso});
    }
    return res;
}


std::unordered_map<CellIndex, W90_realSpaceOperators> w90_calcRealSpaceOperators(const W90_tb &tb, const W90_wsvec &wsvec)
{
    assert ( tb.valid );
    assert ( wsvec.valid );
    unsigned numWann = tb.numWann;
    std::unordered_map<CellIndex, std::map<std::tuple<unsigned, unsigned>, std::vector<CellIndex> >  > translations;
    for(const W90_wsvecCell c : wsvec.cells){
        for(CellIndex t : c.T)
            translations[c.r][{c.m, c.n}].push_back(t);
    }
    std::unordered_map<CellIndex, W90_realSpaceOperators> res;
    for(const W90_wignerSeitzCell &wsc : tb.wsCells){
        unsigned numWann = wsc.H.getColSize();
        const auto & subTrans = translations.at(wsc.r);
        for(unsigned i=0; i<numWann; i++)
            for(unsigned j=0; j<numWann; j++){
                const auto &tVec = subTrans.at({i, j});
                for(const CellIndex &t : tVec){
                    auto it = res.find(wsc.r + t);
                    if ( it == res.end())
                        it = res.insert(it, {wsc.r + t, W90_realSpaceOperators(numWann)});
                    it->second.H(i, j) += wsc.H(i, j) / (double)wsc.degeneracy / (double)tVec.size();
                    for(unsigned dir=0; dir<3; dir++)
                        it->second.D[dir](i, j) += wsc.D[dir](i, j) / (double)wsc.degeneracy / (double)tVec.size();
                }
            }
    }
    return res;
}

bool w90_symmetrizeRealSpaceOperators(std::unordered_map<CellIndex, W90_realSpaceOperators> &rso, bool discardImag)
{
    for(auto it=rso.begin(); it!= rso.end(); it++){
        CellIndex ci = it->first;
        W90_realSpaceOperators &r = it->second;
        if (ci.at(0) == 0 && ci.at(1) == 0 && ci.at(2) == 0){
            r.H.symmetrize();
            for(unsigned dir=0; dir<3; dir++)
                r.D[dir].symmetrize();
            if ( discardImag ){
                r.H.toReal();
                for(unsigned dir=0; dir<3; dir++)
                    r.D[dir].toReal();
            }
        }
        if ( (ci.at(0) < 0) ||
             (ci.at(0) == 0 && ci.at(1) < 0) ||
             (ci.at(0) == 0 && ci.at(1) == 0 && ci.at(2) < 0) ){
            auto itReflected = rso.find(-ci);
            if ( itReflected == rso.end() )
                return false;
            W90_realSpaceOperators &r2 = itReflected->second;
            r.H.symmetrize(r2.H);
            for(unsigned dir=0; dir<3; dir++)
                r.D[dir].symmetrize(r2.D[dir]);
            if ( discardImag ){
                r.H.toReal();
                r2.H.toReal();
                for(unsigned dir=0; dir<3; dir++){
                    r.D[dir].toReal();
                    r2.D[dir].toReal();
                }
            }
        }
    }
    return true;
}

void w90_discardImagOfRealSpaceOperators(std::unordered_map<CellIndex, W90_realSpaceOperators> &rso)
{
    for(auto it=rso.begin(); it!=rso.end(); it++){
        W90_realSpaceOperators & ops = it->second;
        ops.H.toReal();
        for(unsigned dir=0; dir<3; dir++)
            ops.D[dir].toReal();
    }
}


inline void add_w90_mat2D(W90_mat2D &lhs, const W90_mat2D &rhs)
{
    assert(lhs.getColSize() == rhs.getColSize() && lhs.getRowSize() == rhs.getRowSize());
    for(unsigned a=0; a<lhs.getRowSize(); a++)
        for(unsigned b=0; b<lhs.getColSize(); b++)
            lhs(a, b) += rhs(a, b);
}

inline void add_rso(W90_realSpaceOperators &lhs, const W90_realSpaceOperators &rhs)
{
     add_w90_mat2D(lhs.H, rhs.H);
     for(unsigned i=0; i<3; i++)
         add_w90_mat2D(lhs.D[i], rhs.D[i]);
}

bool w90_projectRealSpaceOperators(std::unordered_map<CellIndex, W90_realSpaceOperators> &rso, unsigned dim)
{
    std::unordered_map<CellIndex, W90_realSpaceOperators> newOps;
    for(auto it=rso.begin(); it!=rso.end();){
        const CellIndex &ci = it->first;
        CellIndex ciProj = ci;
        for(unsigned i=dim; i<3; i++)
            ciProj[i] = 0;
        if ( ci != ciProj ){
            auto pit = rso.find(ciProj);
            if ( pit == rso.end()){
                auto nit = newOps.find(ciProj);
                if ( nit == newOps.end() ){
                    newOps.insert({ciProj, it->second});
                } else {
                    add_rso(nit->second, it->second);
                }
            } else {
                add_rso(pit->second, it->second);
            }
            it = rso.erase(it);
        } else {
            it++;
        }
    }
    rso.merge(newOps);
    return true;
}
