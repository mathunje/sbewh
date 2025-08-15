#include "io/KspaceRegionParameter.h"

void KspaceRegionParameter_t::print(FILE *f) const
{
    if ( std::all_of(regions.begin(), regions.end(), [] (const KspaceRegion_t &r) { return std::isinf(r.fwhm); }) ){
         fprintf(f, "No custom KspaceRegionParameter set\n");
         return;
    }
    fprintf(f, "\tstore density matrix: %s\n", (storeDensityMatrix ? "yes" : "no"));
    if ( storeDensityMatrix ) {
        fprintf(f, "\tdensity matrix output stride: %u\n", densityOutputStride);
    }
    fprintf(f, "Custom KspaceRegionParameter:\n");
    fprintf(f, "\tname                            |  kfrac or dKFrac                  | dim or FWHM\n");
    fprintf(f, "\t---------------------------------------------------------------------------------\n");
    bool separated = true;
    for(unsigned ri=0; ri<regions.size(); ri++){
        const KspaceRegion_t &r = regions[ri];
        if ( std::isinf(r.fwhm) )
            continue;
        if ( ! separated && (r.sampleDim[0] != 1 || r.sampleDim[1] != 1 || r.sampleDim[2] != 1)){
            fprintf(f, "\t -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - \n");
        }
        fprintf(f, "\t%-20s  %-9s | %10lg %10lg %10lg  |%10lg\n", r.name.c_str(), (r.movingFrame ? "(moving)" :" (fixed)"),
                    r.kFrac.at(0), r.kFrac.at(1), r.kFrac.at(2), r.fwhm);
        if ( r.sampleDim[2] > 1 )
            fprintf(f, "\t%-20s  %-9s |\n","", "-> volume");
        if ( r.sampleDim[2] == 1 && r.sampleDim[1] > 1 )
            fprintf(f, "\t%-20s  %-9s |\n","", "-> area");
        if ( r.sampleDim[2] == 1 && r.sampleDim[1] == 1 && r.sampleDim[0] > 1 )
            fprintf(f, "\t%-20s  %-9s |\n","", "-> line");
        for(unsigned i=0; i<3; i++){
            if ( r.sampleDim[i] == 1 )
                continue;
            fprintf(f, "\t%-20s  %-9s | %10lg %10lg %10lg  | %u\n", "", "",
                    r.dkFrac[i].at(0), r.dkFrac[i].at(1), r.dkFrac[i].at(2), r.sampleDim[i]);
        }
        if ( ri+1 != regions.size() && (r.sampleDim[0] != 1 || r.sampleDim[1] != 1 || r.sampleDim[2] != 1) ){
            fprintf(f, "\t -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - \n");
            separated = true;
        } else {
            separated = false;
        }
    }
    fprintf(f, "\t---------------------------------------------------------------------------------\n\n");
}


KspaceRegionParameterValidator::KspaceRegionParameterValidator(ParameterIO &paramIO)
{

    paramIO.startGroup("KspaceRegions", "Defines subregions to track and store. Within a region the expectatation values "
                                        "are calculated via an gaussian convolution (of given fwhm) in fractional koorindates "
                                        "for a fixed (F) or comoving (M) central point.\n"
                                        "Possible definitions (format: name = ...):\n"
                                        "  Single Point:\n"
                                        "    {M|F} k1 k2 k3 width\n"
                                        "  Regions along a line:"
                                        "    {M|F} k_start_1 k_start_2 k_start_3 k_end_1, k_end_2 k_end_3 N {I|E} fwhm\n"
                                        "      with N: number of points, I: include endpoint, E: exclude endpoint\n"
                                        "  Regions on an area:\n"
                                        "    {M|F} ks_1 ks_2 ks_3 ka_1 ka_2 ka_3 kb_1 kb_2 kb_3 N_a N_b {I|E} fwhm\n"
                                        "      the region centers are calculated as k = ks + alpha * (ka-ks) + beta * (kb-ks)\n"
                                        "      where alpha, beta are calculated from number of regions N_a and N_b\n"
                                        "  Regions sampling a volume:\n"
                                        "    {M|F} ks_1 ks_2 ks_3 ka_1 ka_2 ka_3 kb_1 kb_2 kb_3 kc_1 kc_2 kc_3 N_a N_b N_c {I|E} fwhm\n"
                                        "      same as for area, but k = ks + alpha * (ka-ks) + beta * (kb-ks) + gamma * (kc-ks)");
    paramIO.allowCustomParameter(&params);
    paramIO.registerParam(new ParameterBool("storeDensityMatrix",
                "if set, the (interpolated) density matrix is stored for each kSpaceRegion.", storeDensityMatrix, false));
    paramIO.registerParam(new ParameterNumeric<unsigned>("densityOutputStride",
                "(interpolated) density matrix is stored only every the specified step.", densityOutputStride, 1));
    paramIO.finishGroup();
}



bool KspaceRegionParameterValidator::update(KspaceRegionParameter_t *p)
{
    if ( storeDensityMatrix && densityOutputStride == 0){
        Logger::error("densityOutputStride must be bigger than zero if storeDensityMatrix is enabled\n");
        return false;
    }
    p->storeDensityMatrix = storeDensityMatrix;
    p->densityOutputStride = densityOutputStride;
    // add default region, filling the complete k-space
    p->regions.push_back(KspaceRegion_t {.name = std::string("expectationValues"),
                                         .movingFrame =true,
                                         .fwhm = std::numeric_limits<double>::infinity(),
                                         .sampleDim = {1, 1, 1},
                                         .kFrac={0, 0, 0},
                                         .dkFrac={ {0, 0, 0}, {0, 0, 0}, {0, 0, 0} }
                                         });
    for(auto &[key, values] : params){
        char mode;
        double k[12];
        char rangeType = 'E';
        unsigned count[3] = {1, 1, 1};
        double fwhm;
        int dim = -1;
        // regions sampling volume
        if ( sscanf(values.c_str(), "%c %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %u %u %u %c %lg",
                    &mode, k, k+1, k+2, k+3, k+4, k+5, k+6, k+7, k+8, k+9, k+10, k+11,
                    count, count+1, count+2, &rangeType, &fwhm) == 18 )
            dim = 3;
        // regions sampling area
        if ( dim < 0 && sscanf(values.c_str(), "%c %lg %lg %lg %lg %lg %lg %lg %lg %lg %u %u %c %lg",
                    &mode, k, k+1, k+2, k+3, k+4, k+5, k+6, k+7, k+8,
                    count, count+1, &rangeType, &fwhm) == 14 )
            dim = 2;
        // regions sampling line
        if ( dim < 0 && sscanf(values.c_str(), "%c %lg %lg %lg %lg %lg %lg %u %c %lg",
                    &mode, k, k+1, k+2, k+3, k+4, k+5, count, &rangeType, &fwhm) == 10 )
            dim = 1;
        // single regions
        if ( dim < 0 && sscanf(values.c_str(), "%c %lg %lg %lg %lg", &mode, k, k+1, k+2, &fwhm) == 5 )
            dim = 0;

        if ( dim < 0 ){
            Logger::error("%s: could not properly scan KspaceRegionParamters '%s'\n", key.c_str(), values.c_str());
            Logger::error("See help for kSpaceRegions for definition of input format\n");
            return false;
        }
        if ( mode != 'F' &&  mode != 'M'){
            Logger::error("%s: mode must be either 'F' (fixed frame) or 'M' (moving frame)\n", key.c_str() );
            return false;
        }
        if ( dim > 0  && (rangeType != 'I' && rangeType != 'E') ){
            Logger::error("%s: range type must be either 'I' (inclusive) or 'E' (exclusive)\n", key.c_str() );
            return false;
        }
        for(unsigned i=0; i<(dim+1)*3; i++)
            if ( k[i] < -1e-10 || k[i] > 1 + 1e-10){
                Logger::error("%s: k[%d] must be in range [0, 1]\n", key.c_str(), (i%3)+1);
                return false;
            }
        if ( fwhm <= 0 || fwhm >= 1 ){
            Logger::error("%s: full width half maximum must be bigger 0 and smaller 1\n", key.c_str());
            return false;
        }
        GeomVector3d kStart({k[0], k[1], k[2]});
        GeomVector3d dk[3];
        for(unsigned i=0; i<dim; i++){
            if ( count[i] < 2 ){
                Logger::error("%s: specified sampling dimensions must all be bigger than 1\n", key.c_str());
                return false;
            }
            unsigned Neff = (rangeType == 'E') ? count[i] : count[i]-1;
            dk[i] = (1.0d/Neff) * ( GeomVector3d({k[3*i+3], k[3*i+4], k[3*i+5]}) - kStart );
        }
        for(unsigned i=dim; i<3; i++){
            count[dim] = 1;
            dk[i] = GeomVector3d({0, 0, 0});
        }
        p->regions.push_back(KspaceRegion_t {.name = key,
                                             .movingFrame = (mode == 'M'),
                                             .fwhm = fwhm,
                                             .sampleDim = {count[0], count[1], count[2]},
                                             .kFrac=kStart,
                                             .dkFrac={dk[0], dk[1], dk[2]},
                                            });
    }
    return true;
}

void mpi_bcast_KspaceRegionParameter(KspaceRegionParameter_t &p, const MpiParameter_t &mp)
{
    unsigned size = (unsigned)p.regions.size();
    mpi_bcast_trivial(size, mp);
    mpi_bcast_trivial(p.storeDensityMatrix, mp);
    mpi_bcast_trivial(p.densityOutputStride, mp);
    if ( mp.rank != mp.root )
        p.regions.resize(size);
    for(unsigned i=0; i<size; i++){
        mpi_bcast(p.regions[i].name, mp);
        mpi_bcast_trivial(p.regions[i].movingFrame, mp);
        mpi_bcast_trivial(p.regions[i].fwhm, mp);
        mpi_bcast_trivial(p.regions[i].sampleDim, mp);
        mpi_bcast(p.regions[i].kFrac, mp);
        mpi_bcast(p.regions[i].dkFrac[0], mp);
        mpi_bcast(p.regions[i].dkFrac[1], mp);
        mpi_bcast(p.regions[i].dkFrac[2], mp);
    }
}
