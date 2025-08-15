#include "ParameterString.h"

ParameterString::ParameterString(const std::string name, const std::string description, std::string &ref) :
    ParameterBase(name, description), data(ref)
{

}

ParameterString::ParameterString(const std::string name, const std::string description, std::string &ref, const std::string defaultValue) :
    ParameterBase(name, description), data(ref)
{
    data = defaultValue;
}

ParameterString::~ParameterString()
{

}

void ParameterString::parseData(const char * str)
{
    data = std::string(str);
    error = false;
    diagnose = "";
}

void ParameterString::printValue(FILE *f) const
{
    fprintf(f, "\"%s\"", data.c_str());
}
