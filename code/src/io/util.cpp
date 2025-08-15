#include "io/util.h"


std::string getOutputDirectoryPath(const OutputParameter_t &out, const std::string fname)
{
    return std::filesystem::path(out.saveDir) / std::filesystem::path(fname);
}

bool createOutputDirectory(const OutputParameter_t &out)
{
    namespace fs = std::filesystem;
    fs::path sd(out.saveDir);
    if ( fs::exists(fs::status(sd)) ){
        Logger::warn("Save directory '%s' already exists\n", out.saveDir.c_str());
        return false;
    }
    if ( ! fs::create_directories(sd) ){
        Logger::error("Could not create save directory '%s' directory\n", out.saveDir.c_str());
        return false;
    }
    Logger::verbose("save to run '%s'\n", out.run.c_str());
    return true;
}

std::vector<std::string> loadArgsFromInputDirectory(const std::string &dir, bool *ok)
{
    *ok = true;
    std::filesystem::path p(dir);
    std::string argFname = p/"args.txt";
    FILE * f = fopen( argFname.c_str(), "r");
    if ( ! f ){
        Logger::error("Could not load '%s'\n", argFname.c_str());
        *ok = false;
        return {};
    }
    char *line = NULL;
    size_t len;
    for(unsigned i=0; i<2; i++)
        if (getline(&line, &len, f) < 0){
            Logger::error("Unexpected format of '%s'\n", argFname.c_str());
            free(line);
            fclose(f);
            *ok = false;
            return {};
        }
    std::vector<std::string> args;
    while (getline(&line, &len, f) >= 0){
        char *l = line;
        while( *l && *l!= '\r' && *l!= '\n')
            l++;
        *l = 0;
        if ( l == line )
            continue;
        args.push_back(line);
    }
    free(line);
    fclose(f);
    bool hasInput = false;
    unsigned i=0;
    while(i<args.size() && !hasInput)
        hasInput = args[i++][0] != '-';
    if ( ! hasInput )
        args.insert(args.begin(), p/"input.txt" );
    return args;
}
