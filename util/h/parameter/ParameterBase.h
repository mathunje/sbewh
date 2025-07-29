#ifndef SBE_WH_PARAMETER_BASE_H
#define SBE_WH_PARAMETER_BASE_H

#include<string>
#include<stdio.h>

class ParameterBase{
protected:
    bool error;
    bool set;
    bool backupSet;
    const std::string name;
    const std::string description;
    std::string diagnose;
    void printDiagnoseStart(FILE *f) const;

    virtual void parseData(const char * str) = 0;
    virtual void backupData() = 0;
    virtual void resetData() = 0;
public:
    ParameterBase(const std::string name, const std::string description);
    virtual ~ParameterBase();
    void parse(const char * str);
    void backup();
    void reset();
    virtual void printValue(FILE *f) const = 0;
    virtual bool isMultiConfig() const { return false; }
    bool hasError() const;
    std::string getDiagnose() const;

    const std::string &getDescription() const;
    const std::string &getName() const;
    const bool * getSetPtr() const { return &set; }

    void print(FILE *f, bool withDescription, bool endLine=true) const;

};


#endif
