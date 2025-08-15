#ifndef SBE_WH_PARAMETER_BOOL_H
#define SBE_WH_PARAMETER_BOOL_H

#include "ParameterBase.h"
#include <cstring>

class ParameterBool : public ParameterBase{
    bool & data;
    bool buData;
    void parseData(const char * str) override;
    void backupData() override { buData = data; }
    void resetData() override { data = buData; }
public:
    ParameterBool(const std::string name, const std::string description, bool &ref);
    ParameterBool(const std::string name, const std::string description, bool &ref, bool defaultValue);
    virtual ~ParameterBool();
    void printValue(FILE *f) const override;
};

#endif
