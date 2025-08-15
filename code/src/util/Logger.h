#ifndef SBE_WH_LOGGER_H
#define SBE_WH_LOGGER_H

#include<stdio.h>
#include<cstring>
#include<stdarg.h>

#include<chrono>
#include<ctime>



class Logger{
public:
    enum Severity {ERROR, WARNING, INFO, VERBOSE, VERYVERBOSE};
private:
    static Severity minSeverity;
    static std::chrono::time_point<std::chrono::steady_clock> startTime;

    static void printTime();
public:
    static void helloWorld(const char *name, const char *gitHash=nullptr, int boxWidth=50);
    static void setVerbosity(Severity s);
    static int print(const char *format, ...); // always printed
    static int error(const char *format, ...);
    static int warn(const char *format, ...);
    static int info(const char *format, ...);
    static int verbose(const char *format, ...);
    static int veryVerbose(const char *format, ...);
};

inline bool operator<(Logger::Severity a, Logger::Severity b){ return static_cast<int>(a) < static_cast<int>(b); }

#endif
