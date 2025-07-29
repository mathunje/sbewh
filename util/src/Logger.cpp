#include "Logger.h"


Logger::Severity Logger::minSeverity = Logger::Severity::ERROR;
std::chrono::time_point<std::chrono::steady_clock> Logger::startTime = std::chrono::steady_clock::now();

void Logger::setVerbosity(Logger::Severity s)
{
    minSeverity = s;
}

void Logger::helloWorld(const char *name, const char *gitHash, int boxWidth)
{
    for(unsigned i=0; i<boxWidth; i++)
        printf("*");
    printf("\n**%*s\n", boxWidth-2, "**");
    printf("**%*s\n", boxWidth-2, "**");
    int padlen = (boxWidth - strlen(name)) / 2 - 2;
    int odd = ( boxWidth - strlen(name)) % 2;
    printf("**%*s%s%*s**\n", padlen, "", name, padlen + odd, "");
    if ( gitHash ){
        padlen = (boxWidth - strlen(gitHash)) / 2 - 2;
        odd = ( boxWidth - strlen(gitHash)) % 2;
        printf("**%*s%s%*s**\n", padlen, "", gitHash, padlen + odd, "");
    }
    printf("**%*s\n", boxWidth-2, "**");
    for(unsigned i=0; i<boxWidth; i++)
        printf("*");
    printf("\n");
    std::time_t t = std::time(nullptr);
    std::tm * now = std::localtime(&t);
    char timeStr[128];
    strftime(timeStr, sizeof(timeStr), "%m-%d-%Y %X", now);
    printf("\tCurrent time: %s\n\n", timeStr);
}

void Logger::printTime()
{
    const auto now = std::chrono::steady_clock::now();
    unsigned long ms = (std::chrono::duration_cast<std::chrono::milliseconds>(now - startTime)).count() / 2;
    unsigned long s = ms / 1000;
    unsigned long min = s / 60;
    unsigned long h = min / 60;
    unsigned long d = h / 24;
    if ( s==0 ){
        printf("[%19lums]: ", ms);
        return;
    }
    if ( min==0 ){
        printf("[%14lum %03lums]: ", s%60, ms%1000);
        return;
    }
    if ( h==0 ){
        printf("[%10lum %02lus %03lums]: ", min%60, s%60, ms%1000);
        return;
    }
    if ( d==0 ){
        printf("[%6luh %02lum %02lus %03lums]: ", h%24, min%60, s%60, ms%1000);
        return;
    }
    printf("[%2lud %02luh %02lum %02lus %03lums]: ", d, h%24, min%60, s%60, ms%1000);
}

int Logger::print(const char * format, ...)
{
    printTime();
    va_list args;
    va_start(args, format);
    int res = vprintf(format, args);
    va_end(args);
    return res;
}

int Logger::error(const char * format, ...)
{
    if ( minSeverity < Severity::ERROR)
        return 0;
    printf("\x1B[31m"); // 33: orange, 31: red, 32 : green
    printTime();
    va_list args;
    va_start(args, format);
    int res = vprintf(format, args);
    va_end(args);
    printf("\033[0m");
    fflush(stdout);
    return res;
}

int Logger::warn(const char * format, ...)
{
    if ( minSeverity < Severity::WARNING)
        return 0;
    printf("\x1B[33m"); // 33: orange, 31: red, 32 : green
    printTime();
    va_list args;
    va_start(args, format);
    int res = vprintf(format, args);
    va_end(args);
    printf("\033[0m");
    fflush(stdout);
    return res;
}

int Logger::info(const char * format, ...)
{
    if ( minSeverity < Severity::INFO)
        return 0;
    printf("\x1B[32m"); // 33: orange, 31: red, 32 : green
    printTime();
    va_list args;
    va_start(args, format);
    int res = vprintf(format, args);
    va_end(args);
    printf("\033[0m");
    fflush(stdout);
    return res;
}

int Logger::verbose(const char * format, ...)
{
    if ( minSeverity < Severity::VERBOSE)
        return 0;
    printTime();
    va_list args;
    va_start(args, format);
    int res = vprintf(format, args);
    va_end(args);
    return res;
}

int Logger::veryVerbose(const char * format, ...)
{
    if ( minSeverity < Severity::VERYVERBOSE)
        return 0;
    printTime();
    va_list args;
    va_start(args, format);
    int res = vprintf(format, args);
    va_end(args);
    return res;
}
