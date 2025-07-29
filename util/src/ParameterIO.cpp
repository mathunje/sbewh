#include "ParameterIO.h"

ParameterIO::ParameterIO()
{
    root = new __paramNode;
    root->customParams = nullptr;
    root->prev = nullptr;
    currentNode = root;
    referenceNode = nullptr;
}

ParameterIO::~ParameterIO()
{
    std::list<__paramNode*> del;
    del.push_front(root);
    while(! del.empty() ){
        __paramNode * p = del.front();
        for(auto &git : p->groups)
            del.push_back(git.second);
        for(auto &pit : p->params)
            delete pit.first;
        del.pop_front();
        delete p;
    }
}

void ParameterIO::startGroup(const std::string name, const std::string description)
{
    auto it = currentNode->groups.find(name);
    if ( it != currentNode->groups.end()){
        currentNode = it->second;
        currentNode->customParams = nullptr;
        return;
    }
    __paramNode * newNode = new __paramNode;
    newNode->customParams = nullptr;
    currentNode->groups.insert({name, newNode});
    newNode->prev = currentNode;
    newNode->description = description;
    currentNode = newNode;
}

void ParameterIO::finishGroup()
{
    if ( currentNode->prev )
        currentNode = currentNode->prev;
}

void ParameterIO::registerParam(ParameterBase *p)
{
    if ( currentNode->params.find(p) != currentNode->params.end() ){
        printf("WARNING: Parameter '%s' already registered -- skipping\n", p->getName().c_str());
        return;
    }
    currentNode->params.insert({p, -1});
}

void ParameterIO::allowCustomParameter(std::map<std::string, std::string> *paramMap)
{
    currentNode->customParams = paramMap;
}


bool ParameterIO::parseParameter(int lineNumber, std::list<std::string> &l)
{
    std::string paramValue = l.back();
    l.pop_back();
    std::string paramName = l.back();
    l.pop_back();
    __paramNode * parseNode = currentNode;
    for(const auto &e : l){
        auto git = parseNode->groups.find(e);
        if ( git == parseNode->groups.end()){
            if ( lineNumber >= 0)
                printf("Line %d: ", lineNumber);
            printf("Subgroup '%s' not known\n", e.c_str());
            return false;
        }
        parseNode = git->second;
    }
    // hacky replacement of .find()
    auto pit = parseNode->params.end();
    for(auto it = parseNode->params.begin(); it != pit; it++)
        if ( it->first->getName() == paramName ){
            pit = it;
            break;
        }
    if ( pit == parseNode->params.end() ){
        if ( parseNode->customParams ){
            parseNode->customParams->insert({paramName, paramValue});
            return true;
        } else {
            if ( lineNumber >= 0)
                printf("Line %d: ", lineNumber);
            printf("Name '%s' not known\n", paramName.c_str());
        }
        return false;
    }
    ParameterBase * p = pit->first;
    pit->second = lineNumber;
    p->parse(paramValue.c_str());
    if ( p->hasError() )
        return false;
    return true;
}

bool ParameterIO::updateCurrentGroup(char updateMode, std::list<std::string> &l, bool verbose)
{
    switch(updateMode){
        case 0: currentNode = root; break;
        case '.': case ':':
            if (referenceNode){
                currentNode = referenceNode;
            } else {
                if (verbose)
                    printf("no reference node for [%c%s] defined\n", updateMode, l.begin()->c_str());
                return false;
            }
            break;
        default:
            break;
    }
    for(const auto &e : l){
        auto git = currentNode->groups.find(e);
        if ( git == currentNode->groups.end()){
            if (verbose)
                printf("subgroup '%s' not known\n", e.c_str());
            return false;
        }
        currentNode = git->second;
    }
    if ( updateMode != '.' )
        referenceNode = currentNode;
    return true;
}

void ParameterIO::parseGroups(std::list<std::string> &l, const std::string &str) const
{
    std::smatch sm;
    const std::regex regLastGroup("((?:.*\\.)*)(\\w+)\\.");
    std::string remainingStr(str);
    while ( std::regex_match(remainingStr, sm, regLastGroup) ){
        l.push_front(sm[2].str());
        remainingStr = sm[1].str();
    }
}

inline bool trimLine(char * line){
    char * pc = line;
    while( *pc && *pc != '#' && *pc != '\\' )
        pc++;
    char last = *pc;
    if ( last == '\\'){
        *(pc+1) = 0;
    } else {
        *pc = 0;
    }
    strtrim(line);
    pc = line;
    while(*pc && *pc != '\\')
        pc++;
    *pc = 0;
    return last == '\\';
}

inline bool readNextLine(FILE *f, std::string &line, int &lineNumber)
{
    const int maxLineSize = 8192;
    char lBuf[maxLineSize];
    if ( ! fgets(lBuf, maxLineSize, f) )
        return false;
    lineNumber++;
    if ( trimLine(lBuf) ){
        char * nextLine = lBuf + strlen(lBuf);
        while( fgets(nextLine, maxLineSize - (int)(lBuf-nextLine), f) && trimLine(nextLine)){
            nextLine += strlen(nextLine);
            lineNumber++;
        }
    }
    line = std::string(lBuf);
    return true;
}

inline bool ParameterIO::matchParameter(const std::string &line, int lineNumber)
{
    const std::regex regParameter ("((?:\\w+\\.)*)(\\w+)\\s*=\\s*(.*)");
    std::smatch sm;
    if ( std::regex_match(line, sm, regParameter) ){
        std::list<std::string> subStr;
        subStr.push_back(sm[2].str());
        subStr.push_back(sm[3].str());
        parseGroups(subStr, sm[1].str());
        if ( !parseParameter(lineNumber, subStr) )
            return false;
    }
    return true;
}

inline bool ParameterIO::matchSection(const std::string &line, bool verbose)
{
    const std::regex regSection ("\\[(\\.|:)?((?:\\w+\\.)*)(\\w+)\\]");
    std::smatch sm;
    if ( std::regex_match(line, sm, regSection) ){
        char updateMode = sm[1].str().c_str()[0];
        std::list<std::string> subStr;
        subStr.push_back(sm[3].str());
        parseGroups(subStr, sm[2].str());
        if ( !updateCurrentGroup(updateMode, subStr, verbose) )
            return false;
    }
    return true;
}

ParameterIO::ParseResult ParameterIO::parse(const std::string filename)
{
    paramFnames.push_back(filename);
    currentNode = root;
    FILE * f = fopen(filename.c_str(), "r");
    if ( !f ){
        printf("Could not open '%s'\n", filename.c_str());
        return ParseResult::Fail;
    }
    int lineNumber = 0;
    std::string line;
    while (readNextLine(f, line, lineNumber)){
        if ( ! matchParameter(line, lineNumber) ){
            fclose(f);
            return ParseResult::Fail;
        }
        if ( ! matchSection(line)){
            fclose(f);
            return ParseResult::Fail;
        }
    }
    fclose(f);
    return ParseResult::Ok;
}

ParameterIO::ParseResult ParameterIO::parse(int argc, const char *argv[], const char * defaultInputFile)
{
    bool noInput = true;
    for(int i=1; i<argc; i++){
        noInput &= (argv[i][0] == '-');
        if ( !strcmp(argv[i], "--help") || !strcmp(argv[i], "-help") ||
             !strcmp(argv[i], "--h") || !strcmp(argv[i], "-h") ){
            if ( i+1 < argc ){
                if ( ! strcmp(argv[i+1], "-h") || !strcmp(argv[i+1], "-help") || !strcmp(argv[i+1], "-par") ){
                    printHelpInfo();
                } else {
                    printHelp(argv[i+1]);
                }
            } else {
                printHelpInfo();
            }
            return ParseResult::HelpRequested;
        }
        if ( !strcmp(argv[i], "--par") || !strcmp(argv[i], "-par") ){
            printDefinedParams();
            return ParseResult::HelpRequested;
        }
    }
    if ( noInput ){
        if ( ! defaultInputFile ){
            printf("No [default] input file given\n");
            return ParseResult::Fail;
        } else {
            ParseResult pr = parse(defaultInputFile);
            if ( pr != ParseResult::Ok)
                return pr;
        }
    }
    for(unsigned i=1; i<argc; i++){
        currentNode = root;
        if (argv[i][0] != '-'){
            // parse all input files
            ParseResult pr = parse(argv[i]);
            if ( pr != ParseResult::Ok )
                return pr;
        } else {
            // parse option
            std::string inp(argv[i]);
            std::smatch sm;
            const std::regex regParameter ("-((?:\\w+\\.)*)(\\w+)\\s*=\\s*(.*)");
            if ( std::regex_match(inp, sm, regParameter) ){
                // found in definition
                std::list<std::string> subStr;
                subStr.push_back(sm[2].str());
                subStr.push_back(sm[3].str());
                parseGroups(subStr, sm[1].str());
                if ( ! parseParameter(-1, subStr) )
                    return ParseResult::Fail;
            } else {
                printf("Could not parse argument '%s'\n", argv[i]);
                printf("Use -help to get more information\n");
                printf("Use -par prints all defined parameters\n");
                return ParseResult::Fail;
            }
        }
    }
    return ParseResult::Ok;
}

bool ParameterIO::parseMultiConfigFile(__paramNode *p, const char *fname, std::function<bool(std::string)> callback, bool reset){
    FILE * f = fopen(fname, "r");
    if ( !f )
        return false;
    currentNode = p;
    std::string line;
    int lineNumber = 0;
    std::string configName = "";
    while (readNextLine(f, line, lineNumber)){
        if ( ! matchParameter(line, lineNumber) ){
            fclose(f);
            return false;
        }
        if ( ! matchSection(line, false) ){
            if ( configName.size() )
                callback(configName);
            if ( reset ){
                std::list<__paramNode*> pl;
                pl.push_front(p);
                while(! pl.empty() ){
                    __paramNode * pn = pl.front();
                    pl.pop_front();
                    for(auto &pit : pn->params)
                        pit.first->reset();
                    for(auto &git : pn->groups)
                        pl.push_back(git.second);
                }
            }
            // remove leading '[' and ending ']':
            configName = line.substr(1, line.size()-2);
            currentNode = p;
        }
    }
    fclose(f);
    if ( configName.size() )
        callback(configName);
    return true;
}

bool ParameterIO::parseMultiConfigs(bool reset){
    if ( reset ) {
        std::list<__paramNode*> pl;
        pl.push_front(root);
        while(! pl.empty() ){
            __paramNode * p = pl.front();
            pl.pop_front();
            for(auto &pit : p->params)
                pit.first->backup();
            for(auto &git : p->groups)
                pl.push_back(git.second);
        }
    }
    std::list<__paramNode*> pl;
    pl.push_front(root);
    while(! pl.empty() ){
        __paramNode * p = pl.front();
        pl.pop_front();
        for(auto &pit : p->params){
            auto pNode = pit.first;
            if ( pNode->isMultiConfig() ){
                auto pi = static_cast<ParameterMultiConfig*>(pNode);
                auto fname = pi->getFileName();
                auto callback = pi->getCallback();
                if ( fname.size() && ! parseMultiConfigFile(p, fname.c_str(), callback, reset) ){
                    printf("Could not parse multiconfiguration file '%s'\n", fname.c_str());
                    return false;
                }
            }
        }
        for(auto &git : p->groups)
            pl.push_back(git.second);
    }
    return true;
}

std::vector<std::string> ParameterIO::getParamFileNames() const
{
    return paramFnames;
}

std::pair<__paramNode*, ParameterBase*> ParameterIO::parseGroupString(const std::string s) const
{
    if ( s == "")
        return {root, nullptr};

    __paramNode * startNode = root;
    // parse name and adapt start ...
    std::smatch sm;
    const std::regex regLastGroup("(?:-)?((?:\\w+\\.)*)(\\w+)");
    if ( ! std::regex_match(s, sm, regLastGroup) )
        return {nullptr, nullptr};

    std::string lastPart = sm[2].str();
    // traverse tree until only last part remains
    std::list<std::string> groups;
    parseGroups(groups, sm[1].str());
    for(auto &g : groups){
        auto git = startNode->groups.find(g);
        if ( git == startNode->groups.end())
            return {nullptr, nullptr};
        startNode = git->second;
    }
    // last part -- as group
    auto git = startNode->groups.find(lastPart);
    if ( git != startNode->groups.end())
        return {git->second, nullptr};
    // last part as Parameter
    for(auto &p : startNode->params)
        if ( p.first->getName() == lastPart )
            return {nullptr, p.first};

    return {nullptr, nullptr};
}

void ParameterIO::printDiagnose(FILE *f) const
{
    std::list<__paramNode*> pl;
    std::list<std::tuple<__paramNode*, ParameterBase*, int> > errorList;
    pl.push_front(root);
    while(! pl.empty() ){
        __paramNode * p = pl.front();
        pl.pop_front();
        for(auto &pit : p->params)
            if ( pit.first->hasError() ){
                errorList.push_back(std::make_tuple(p, pit.first, pit.second));
            }
        for(auto &git : p->groups)
            pl.push_back(git.second);
    }
    if ( errorList.empty() ){
        fprintf(f, "Parameter diagnose indicated no error\n");
        return;
    }
    for(auto [pn, p, lineNumber] : errorList){
        std::list<std::string> groups;
        while(pn->prev){
            // hacky way of finding forward reference
            for(auto &git : pn->prev->groups)
                if ( git.second == pn )
                    groups.push_front(git.first);
            pn = pn->prev;
        }
        if ( lineNumber >= 0)
            fprintf(f, "[Line %d] ", lineNumber);
        for(auto & g : groups)
            fprintf(f, "%s.", g.c_str());
        fprintf(f, "%s: %s\n", p->getName().c_str(), p->getDiagnose().c_str());
    }
}

void ParameterIO::printDefinedParams() const
{
    if ( ! root )
        printf("No parameter defined\n");
    std::list< std::pair<__paramNode*, const std::string> > pl; // node, full-name
    pl.push_front({root, ""});
    while(! pl.empty() ){
        auto [p, cname]  = pl.front();
        pl.pop_front();
        for(auto &pit : p->params){
            std::string name = cname + pit.first->getName();
            printf("%s\n", name.c_str() );
        }
        for(auto iter = p->groups.rbegin(); iter != p->groups.rend(); ++iter)
            pl.push_front( {iter->second, cname + iter->first + "." } );
    }
}

void ParameterIO::printParam(FILE *f, const std::string name, bool withGroupDescription) const
{
    auto [startNode, param] = parseGroupString(name);
    if ( param ) {
        param->print(f, withGroupDescription);
        return;
    }
    if ( ! startNode ){
        fprintf(f, "Unrecognized parameter / group name '%s'", name.c_str());
        return;
    }
    // traverse / print in dfs order
    std::list< std::tuple<__paramNode*, const std::string, int> > pl; // node, name, level
    pl.push_front(std::make_tuple(startNode, name, 0));
    while(! pl.empty() ){
        auto [p, cname, level]  = pl.front();
        pl.pop_front();
        fprintf(f, "%*s%s\n", 2 * level, "", cname.c_str());
        for(auto &pit : p->params){
            fprintf(f, "%*s", 2*(level+1), "");
            pit.first->print(f, withGroupDescription);
        }
        for(auto iter = p->groups.rbegin(); iter != p->groups.rend(); ++iter)
            pl.push_front(std::make_tuple(iter->second, iter->first, level+1));
        if ( p->customParams ){
            fprintf(f, "%*sCustom parameters\n", 2*(level+1), "");
            for(const auto &[name, val] : *p->customParams )
                fprintf(f, "%*s%s: %s\n", 2*(level+2), "", name.c_str(), val.c_str());
        }
    }
}

void ParameterIO::printHelpInfo() const
{
    printf("Help\n");
    printf("\t-par          bare list of parameters\n");
    printf("\t-help         print this help function\n");
    printf("\t-help group   list defined subgroups (subgroup can be accessed via group.subgroup)\n");
    printf("\t-help name    show description of parameter\n");
    printf("List of main Groups:\n");
    for( auto [name, _] : root->groups)
        printf("\t%s\n", name.c_str());
}

void ParameterIO::printHelp(const std::string name) const
{
    auto [startNode, param] = parseGroupString(name);
    if ( !startNode ){
        if ( ! param ){
            printf("Cannot help: Unrecognized parameter / group name '%s'\n", name.c_str());
            printHelpInfo();
        } else {
            param->print(stdout, true);
        }
        return;
    }
    printf("Help for '%s'\n", name.c_str());
    printf("----\n");
    printf("%s\n", startNode->description.c_str());
    printf("----\n");
    if ( param ) {
        param->print(stdout, true);
        return;
    }
    if ( startNode->groups.size()){
        printf("List of subgroups:\n");
        for( auto [name, _] : startNode->groups)
            printf("\t%s\n", name.c_str());
    }
    printf("Parameter:\n");
    for(auto &pit : startNode->params){
        pit.first->print(stdout, true);
        printf("\n");
    }
}
