#ifndef SBE_WH_PARAMETER_MULTICONFIG_H
#define SBE_WH_PARAMETER_MULTICONFIG_H

#include<filesystem>
#include<functional>

#include "ParameterBase.h"



class ParameterMultiConfig : public ParameterBase{
    std::string fname;
    std::function<bool(std::string)> callback;
    std::string fixedInputDir;
    void parseData(const char * str) override;
    void backupData() override {}
    void resetData() override {}
public:
    ParameterMultiConfig(const std::string name, const std::string description, std::function<bool(std::string)> callback, const std::string fixedInputDir="");
    virtual ~ParameterMultiConfig();
    bool isMultiConfig() const override{ return true; }
    void printValue(FILE *f) const override;
    const std::string &getFileName() const;
    std::function<bool(std::string)> getCallback() { return callback; }

};


#endif
