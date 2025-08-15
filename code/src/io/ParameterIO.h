#ifndef SBE_WH_PARAMETER_IO_H
#define SBE_WH_PARAMETER_IO_H

#include <string>
#include <stdio.h>

#include<map>
#include<vector>
#include<list>
#include<algorithm>
#include<regex>

#include<cstring>
#include<cassert>

#include "parameter/ParameterBase.h"
#include "parameter/ParameterMultiConfig.h"
#include "detail/parse.h"


struct __paramCompare {
    inline bool operator() (const ParameterBase * a, const ParameterBase * b) const
                           { return a->getName() < b->getName(); }
};

struct __paramNode {
    std::string description;
    std::map<ParameterBase*,int, __paramCompare> params; // int captures the lineNumber
    std::map<std::string, __paramNode*> groups;
    std::map<std::string, std::string> *customParams;
    __paramNode * prev;
};


class ParameterIO {
    __paramNode * root;
    __paramNode * currentNode;
    __paramNode * referenceNode;
    std::vector<std::string> paramFnames;
    bool parseParameter(int lineNumber, std::list<std::string> &l);
    void parseGroups(std::list<std::string> &l, const std::string &str) const;
    bool updateCurrentGroup(char updateMode, std::list<std::string> &l, bool verbose=true);
    std::pair<__paramNode*, ParameterBase*> parseGroupString(const std::string s) const;

    inline bool matchParameter(const std::string &line, int lineNumber);
    inline bool matchSection(const std::string &line, bool verbose=true);

    bool parseMultiConfigFile(__paramNode *p, const char *fname, std::function<bool(std::string)> callback, bool reset);

    void printDefinedParams() const;
    void printHelp(const std::string name) const;
    void printHelpInfo() const;
public:
    enum ParseResult { Ok, Fail, HelpRequested };
    ParameterIO();
    ~ParameterIO();
    void startGroup(const std::string name, const std::string description);
    void finishGroup();
    void allowCustomParameter(std::map<std::string, std::string> *paramMap);
    void registerParam(ParameterBase * p);

    ParseResult parse(const std::string filename);
    ParseResult parse(int argc, const char * argv[], const char * defaultInputFile=NULL);

    bool parseMultiConfigs(bool reset=true);

    std::vector<std::string> getParamFileNames() const;
    void printDiagnose(FILE *f) const;
    void printParam(FILE *f, const std::string name = "", bool withDescription=false) const;
};

#endif
