#include "ParameterNumeric.hpp"

template<> bool ParameterNumeric<int>::fromString(const char * str, int * v) { return sscanf(str, "%d", v) == 1; }
template<> std::string ParameterNumeric<int>::toString(const int v) { char buf[30]; sprintf(buf, "%d", v); return std::string(buf); }

template<> bool ParameterNumeric<unsigned int>::fromString(const char * str, unsigned int * v) { return sscanf(str, "%u", v) == 1; }
template<> std::string ParameterNumeric<unsigned int>::toString(const unsigned int v) { char buf[30]; sprintf(buf, "%u", v); return std::string(buf); }

template<> bool ParameterNumeric<long long int>::fromString(const char * str, long long int *v) { return sscanf(str, "%lld", v) == 1; }
template<> std::string ParameterNumeric<long long int>::toString(const long long int v) { char buf[30]; sprintf(buf, "%lld", v); return std::string(buf); }

template<> bool ParameterNumeric<unsigned long long int>::fromString(const char * str, unsigned long long int *v) { return sscanf(str, "%llu", v) == 1; }
template<> std::string ParameterNumeric<unsigned long long int>::toString(const unsigned long long int v) { char buf[30]; sprintf(buf, "%llu", v); return std::string(buf); }

template<> bool ParameterNumeric<float>::fromString(const char * str, float * v) { return sscanf(str, "%f", v) == 1; }
template<> std::string ParameterNumeric<float>::toString(const float v) { char buf[30]; sprintf(buf, "%g", v); return std::string(buf); }

template<> bool ParameterNumeric<double>::fromString(const char * str, double * v) { return sscanf(str, "%lg", v) == 1; }
template<> std::string ParameterNumeric<double>::toString(const double v) { char buf[30]; sprintf(buf, "%lg", v); return std::string(buf); }
