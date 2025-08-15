#include "ParameterBase.h"

ParameterBase::ParameterBase(const std::string name, const std::string description) :
    name(name), description(description)
{
    error = false;
    set = false;
    backupSet = false;
}

ParameterBase::~ParameterBase()
{

}

void ParameterBase::parse(const char * str)
{
    parseData(str);
    set = true;
}

void ParameterBase::backup()
{
    backupData();
    backupSet = set;
}

void ParameterBase::reset()
{
    resetData();
    set = backupSet;
}

const std::string& ParameterBase::getName() const
{
    return name;
}

const std::string& ParameterBase::getDescription() const
{
    return description;
}

bool ParameterBase::hasError() const
{
    return error;
}

std::string ParameterBase::getDiagnose() const{
    return diagnose;
}

void ParameterBase::print(FILE *f, bool withDescription, bool endLine) const
{
    fprintf(f, "%s = ", name.c_str());
    printValue(f);
    if ( withDescription){
        fprintf(f, "\n\t%s", description.c_str());
    }
    if (endLine)
        fprintf(f, "\n");
}
