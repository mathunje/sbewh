#ifndef SBE_WH_REQUIRED_PARAMETER_H
#define SBE_WH_REQUIRED_PARAMETER_H

/* - compare to Parameter
 * - the here not listed members must be parsed */
struct RequiredParameter_t {
    bool par = false;
    bool tb = false;
    bool ft = false;
    bool diag = false;
    bool pulse = false;
    bool prop = false;
    bool ksr = false;
    bool out = false;
    bool multiConfigs = false;
};

#endif
