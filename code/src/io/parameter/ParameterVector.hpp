#ifndef SBE_WH_PARAMETER_VECTOR_HPP
#define SBE_WH_PARAMETER_VECTOR_HPP

#include "ParameterBase.h"
#include "ParameterNumeric.hpp"
#include <vector>
#include <cstring>

template<typename T>
class ParameterVector : public ParameterBase{
    std::vector<T> &vec;
    std::vector<T> backupVec;
    const char * delim;
    void parseData(const char * str) override;
    void backupData() override { backupVec = vec; }
    void resetData() override { vec = backupVec; }
public:
    ParameterVector(const std::string name, const std::string description, std::vector<T> &ref, const char *delimiter=",");
    ParameterVector(const std::string name, const std::string description, std::vector<T> &ref, const std::vector<T> defaultVec, const char * delimiter=",");
    virtual ~ParameterVector();
    void printValue(FILE *f) const override;

};


template<typename T>
ParameterVector<T>::ParameterVector(const std::string name, const std::string description, std::vector<T> &ref,
                                    const char * delimiter) :
    ParameterBase(name, description),  vec(ref), delim(delimiter)
{

}

template<typename T>
ParameterVector<T>::ParameterVector(const std::string name, const std::string description,
                                    std::vector<T> &ref, const std::vector<T> defaultVec,
                                    const char * delimiter) :
    ParameterBase(name, description),  vec(ref), delim(delimiter)
{
    vec = defaultVec;
}

template<typename T>
ParameterVector<T>::~ParameterVector()
{

}

template<typename T>
void ParameterVector<T>::parseData(const char * str)
{
    vec.clear();
    error = false;
    diagnose = "";
    char *strCopy = (char*) malloc(sizeof(char) * (1 + strlen(str)));
    strcpy(strCopy, str);
    char * token = strtok(strCopy, delim);
    while ( token ){
        T cval;
        if ( ! ParameterNumeric<T>::fromString(token, &cval) ){
            error = true;
            diagnose = "Could not read value of '";
            diagnose.append(token);
            diagnose.append("'");
            vec.clear();
            free(strCopy);
            return;
        }
        vec.push_back(cval);
        token = strtok(NULL, delim);
    }
    free(strCopy);
}

template<typename T>
void ParameterVector<T>::printValue(FILE *f) const
{
    for(unsigned i=0; i<vec.size(); i++){
        fprintf(f, "%s", ParameterNumeric<T>::toString(vec[i]).c_str());
        if ( i+1 < vec.size() )
            fprintf(f, ", ");
    }
}


#endif
