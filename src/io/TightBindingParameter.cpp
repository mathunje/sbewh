#include "io/TightBindingParameter.h"

void TightBindingParameter_t::print(FILE *f, unsigned dim) const
{
    fprintf(f, "Tightbinding parameter\n");

    fprintf(f, "\tWannier tight binding file: %s\n", wannierTbFname.c_str());
    if ( wannierWsvecFname.compare("") )
        fprintf(f, "\tWannier wsvec file: %s\n", wannierWsvecFname.c_str());
    for(unsigned i=0; i<dim; i++)
        fprintf(f, "\tNw%u: %u\n", i+1, Nw[i]);

    fprintf(f, "\tLattice vectors [angstrom]:\n");
    for(unsigned i=0; i<dim; i++){
        fprintf(f, "\t");
        for(unsigned j=0; j<dim; j++)
            fprintf(f, "    %+8.5lf", atomicUnits::to_A(latticeVectors.at(i).at(j)));
        fprintf(f, "\n");
    }
    fprintf(f, "\tGround state temperature: %lf K\n", atomicUnits::to_K(gsTemp));
}


TightBindingParameterValidator::TightBindingParameterValidator(ParameterIO &paramIO, const std::string &fixedInputDir)
    : fixedInputDir(fixedInputDir)
{
    paramIO.startGroup("TightBinding", "Defines the (static) model parameters\n");
    ParameterString * pWannierSeed = new ParameterString("wannierSeed", "looks for seed_tb.dat and seed_wsvec.dat in format as created by wannier90.\n"
                                                                            "\tExpected units: energy -> eV, length -> A",
                                                           wannierSeed);
    paramIO.registerParam(pWannierSeed);
    paramIO.registerParam(new ParameterBool("symmetrizeHamiltonian", "if enabled: symmetry of wannier elements is enforced",
                                            symmetrizeHamiltonian, true));
    paramIO.registerParam(new ParameterBool("onlyRealMatrixElements", "if enabled: imaginary part of Hamiltonian and dipole matrix elements "
                                                                  "are ignored", onlyRealMatrixElements, false));

    paramIO.registerParam(new ParameterNumeric<unsigned>("occupiedBands",
                "number of bands occupied to define Fermi energy (if occupyBelowBandGap is not set)", occupiedBands, 0));
    paramIO.registerParam(new ParameterBool("occupiedBelowBandGap",
                "defines the Fermi energy in the middle of the lowest band gap", occupiedBelowBandGap, false));

    auto *pfl = new ParameterNumeric<double>("fermiLevel", "fixes Fermi level regardless of other options [eV]", customFermiLevel, 0);
    customFermiLevelSetPtr = pfl->getSetPtr();
    paramIO.registerParam(pfl);

    ParameterNumeric<double> * pT = new ParameterNumeric<double>("temperature", "electronic temperature [K] to define groundstate", gsTemp, 300);
    pT->setRange(0, 5000);
    paramIO.registerParam(pT);
    paramIO.finishGroup();
}



bool TightBindingParameterValidator::update(TightBindingParameter_t *p, unsigned dim)
{
    if ( fixedInputDir.size() > 0 )
        wannierSeed = std::filesystem::path(fixedInputDir) / std::filesystem::path(wannierSeed).filename();
    p->wannierTbFname = wannierSeed + "_tb.dat";
    p->wannierWsvecFname = wannierSeed + "_wsvec.dat";
    p->gsTemp = atomicUnits::from_K(gsTemp);
    p->occupiedBands = occupiedBands;
    p->occupiedBelowBandGap = occupiedBelowBandGap;
    p->customFermiLevel = atomicUnits::from_eV(customFermiLevel);
    p->useCustomFermiLevel = *customFermiLevelSetPtr;

    if ( wannierSeed.compare("") ){
        if ( !evalWannierFiles(p, dim) )
            return false;
    } else {
        Logger::error("No wannier seed set\n");
        return false;
    }
    p->minCellIndices = CellIndex({0, 0, 0});
    int maxCi[3] = {0, 0, 0};
    for(const auto &[ci, _] : p->realSpaceOps)
        for(unsigned i=0; i<3; i++){
            if(ci.at(i) < p->minCellIndices[i])
                p->minCellIndices[i] = ci.at(i);
            if(ci.at(i) > maxCi[i])
                maxCi[i] = ci.at(i);
        }
    unsigned tbDim = 1;
    for(unsigned i=0; i<3; i++){
        p->Nw[i] = maxCi[i] - p->minCellIndices[i] + 1;
        if (p->Nw[i] > 1)
            tbDim = 1 + i;
    }
    if ( tbDim < dim ){
        Logger::info("Dimensionality of Wannier file is %u and not %u\n", tbDim, dim);
    }
    if ( ! p->useCustomFermiLevel ){
        if ( ! occupiedBelowBandGap ){
            if ( occupiedBands == 0){
                Logger::error("Either occupiedBands > 0 or occupiedBelowBandGap must be set.\n");
                return false;
            }
            if ( occupiedBands >= p->numWann ){
                Logger::error("occupiedBands (%u) >= numWann (%u)\n", occupiedBands, p->numWann);
                return false;
            }
        }
    }
    return true;
}

bool TightBindingParameterValidator::evalWannierFiles(TightBindingParameter_t *p, unsigned dim)
{
    if ( ! std::filesystem::exists(p->wannierTbFname) ){
        Logger::error("The wannier file '%s' does not exist\n", p->wannierTbFname.c_str());
        return false;
    }
    W90_tb wtb = w90_read_tb(p->wannierTbFname.c_str());
    if ( ! wtb.valid ){
        Logger::error("Could not read '%s' properly\n", p->wannierTbFname.c_str());
        return false;
    }

    if ( std::filesystem::exists(p->wannierWsvecFname) ){
        Logger::info("Reading wsvec and using new interpolation scheme\n");
        W90_wsvec wsv = w90_read_wsvec(p->wannierWsvecFname.c_str());
        if ( ! wsv.valid ){
            Logger::error("Could not read '%s' properly\n", p->wannierWsvecFname.c_str());
            return false;
        }
        p->realSpaceOps = w90_calcRealSpaceOperators(wtb, wsv);
    } else {
        Logger::info("No wsvec found -- using old interpolation scheme\n");
        p->realSpaceOps = w90_calcRealSpaceOperators(wtb);
    }
    if ( symmetrizeHamiltonian ){
        w90_symmetrizeRealSpaceOperators(p->realSpaceOps, onlyRealMatrixElements);
    } else {
        if ( onlyRealMatrixElements )
            w90_discardImagOfRealSpaceOperators(p->realSpaceOps);
    }
    for(unsigned i=0; i<3; i++)
        p->latticeVectors[i] = wtb.a[i];
    if ( ! w90_projectRealSpaceOperators(p->realSpaceOps, dim) )
        return false;
    auto itOrig = p->realSpaceOps.find({0, 0, 0});
    if ( itOrig == p->realSpaceOps.end()){
        Logger::error("Projected hamiltonian does not contain R=(0, 0, 0)\n");
        return false;
    } else {
        p->numWann = itOrig->second.H.getColSize();
    }
    return true;
}

void mpi_bcast_TightBindingParameter(TightBindingParameter_t &p, const MpiParameter_t &mp)
{
    mpi_bcast_trivial(p.numWann, mp);
    mpi_bcast_trivial(p.gsTemp, mp);
    mpi_bcast_trivial(p.occupiedBands, mp);
    mpi_bcast_trivial(p.occupiedBelowBandGap, mp);
    for(unsigned i=0; i<p.latticeVectors.size(); i++)
        mpi_bcast(p.latticeVectors[i], mp);
    int size = (int)p.realSpaceOps.size();
    mpi_bcast_trivial(size, mp);
    auto it = p.realSpaceOps.begin();
    for(int i=0; i<size; i++){
        W90_realSpaceOperators rso(p.numWann);
        CellIndex ci;
        if( mp.rank == mp.root ){
            ci = it->first;
            rso = it->second;
            it++;
        }
        mpi_bcast(rso, mp);
        mpi_bcast(ci, mp);
        if ( mp.rank != mp.root )
           p.realSpaceOps.insert({ci, rso});
    }
    mpi_bcast(p.wannierTbFname, mp);
    mpi_bcast(p.wannierWsvecFname, mp);
    mpi_bcast_trivial(p.Nw, mp);
    mpi_bcast(p.minCellIndices, mp);
}
