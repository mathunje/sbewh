#ifndef SBE_WH_PARAMETER_NUMERIC_HPP
#define SBE_WH_PARAMETER_NUMERIC_HPP

#include "ParameterBase.h"
#include <cstring>

template<typename T>
class ParameterNumeric : public ParameterBase{
    T & value;
    T backupValue;
    bool minSet;
    T minValue;
    bool maxSet;
    T maxValue;
    void updateValue(const T readValue);

    void parseData(const char * str) override;
    void backupData() override { backupValue = value; }
    void resetData() override { value = backupValue; }
public:
    ParameterNumeric(const std::string name, const std::string description, T &ref);
    ParameterNumeric(const std::string name, const std::string description, T &ref, const T defaultValue);
    void setMinValue(const T minValue);
    void setMaxValue(const T maxValue);
    void setRange(const T minValue, const T maxValue);
    virtual ~ParameterNumeric();
    void printValue(FILE *f) const override;

    static std::string toString(const T v);
    static bool fromString(const char * str, T *v);

};


template<typename T>
ParameterNumeric<T>::ParameterNumeric(const std::string name, const std::string description, T &ref) :
    ParameterBase(name, description),  value(ref)
{
    minSet = false;
    maxSet = false;
}

template<typename T>
ParameterNumeric<T>::ParameterNumeric(const std::string name, const std::string description, T &ref, const T defaultValue) :
    ParameterBase(name, description),  value(ref)
{
    value = defaultValue;
    minSet = false;
    maxSet = false;
}

template<typename T>
ParameterNumeric<T>::~ParameterNumeric()
{

}

template<typename T>
void ParameterNumeric<T>::setMinValue(const T minValue)
{
    minSet = true;
    this->minValue = minValue;
}

template<typename T>
void ParameterNumeric<T>::setMaxValue(const T maxValue)
{
    maxSet = true;
    this->maxValue = maxValue;
}

template<typename T>
void ParameterNumeric<T>::setRange(const T minValue, const T maxValue)
{
    minSet = true;
    this->minValue = minValue;
    maxSet = true;
    this->maxValue = maxValue;
}

template<typename T>
void ParameterNumeric<T>::parseData(const char * str)
{
    error = false;
    diagnose = "";
    T readValue;
    if ( ! fromString(str, &readValue) ){
        error = true;
        diagnose = "Could read value of '";
        diagnose.append(str);
        diagnose.append("'");
        return;
    }
    if ( minSet && ( readValue < minValue) ){
        error = true;
        diagnose = "Read value ";
        diagnose.append(toString(readValue));
        diagnose.append(" is smaller than minimum allowed value of ");
        diagnose.append(toString(minValue));
        return;
    }
    if ( maxSet && ( readValue > maxValue) ){
        error = true;
        diagnose = "Read value ";
        diagnose.append(toString(readValue));
        diagnose.append(" is greater than maximum allowed value of ");
        diagnose.append(toString(maxValue));
        return;
    }
    value = readValue;
}

template<typename T>
void ParameterNumeric<T>::printValue(FILE *f) const
{
    fprintf(f, "%s", toString(value).c_str());
}

/*********************** specializations ******************/
template<> bool ParameterNumeric<int>::fromString(const char * str, int *v);
template<> bool ParameterNumeric<unsigned int>::fromString(const char * str, unsigned int *v);
template<> bool ParameterNumeric<long long int>::fromString(const char * str, long long int *v);
template<> bool ParameterNumeric<unsigned long long int>::fromString(const char * str, unsigned long long int *v);
template<> bool ParameterNumeric<float>::fromString(const char * str, float *v);
template<> bool ParameterNumeric<double>::fromString(const char * str, double *v);

/// TODO: add writing + parsing for complex numbers (may contain empty spaces...)

#endif
