#ifndef PARAMETER_STRING_H
#define PARAMETER_STRING_H

#include "ParameterBase.h"
#include <cstring>

class ParameterString : public ParameterBase{
    std::string & data;
    std::string buData;
    void parseData(const char * str) override;
    void backupData() override { buData = data; }
    void resetData() override { data = buData; }
public:
    ParameterString(const std::string name, const std::string description, std::string &ref);
    ParameterString(const std::string name, const std::string description, std::string &ref, std::string defaultValue);
    virtual ~ParameterString();
    void printValue(FILE *f) const override;

};

#endif
