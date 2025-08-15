#include "ParameterBool.h"

ParameterBool::ParameterBool(const std::string name, const std::string description, bool &ref) :
    ParameterBase(name, description), data(ref)
{

}

ParameterBool::ParameterBool(const std::string name, const std::string description, bool &ref, bool defaultValue) :
    ParameterBase(name, description), data(ref)
{
    data = defaultValue;
}

ParameterBool::~ParameterBool()
{

}

void ParameterBool::parseData(const char * str)
{
    error = false;
    diagnose = "";
    if ( !strcmp(str, "false") ){
        data = false;
        return;
    }
    if ( !strcmp(str, "true") ){
        data = true;
        return;
    }
    int v = 0;
    if ( !sscanf(str, "%d", &v) ){
        error = true;
        diagnose = "Could not interpret '";
        diagnose.append(str);
        diagnose.append("' as boolean");
        return;
    }
    data = (v != 0);
}

void ParameterBool::printValue(FILE *f) const
{
    if ( data ){
        fprintf(f, "true");
    } else {
        fprintf(f, "false");
    }
}
