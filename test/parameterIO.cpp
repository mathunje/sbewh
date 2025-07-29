#include "gtest/gtest.h"

#include "ParameterIO.h"
#include "ParameterNumeric.hpp"

ParameterIO paramIO(int &a, int &b, int &c, int &csg, int &d)
{
    ParameterIO p;
    p.registerParam(new ParameterNumeric<int>("c", "c_description", c));
    p.startGroup("group", "fancy");
        p.registerParam(new ParameterNumeric<int>("a", "a_description", a));
        p.startGroup("subgroup", "super fancy");
            p.registerParam(new ParameterNumeric<int>("b", "b_description", b));
        p.finishGroup();
        p.startGroup("sg", "ok");
            p.registerParam(new ParameterNumeric<int>("c", "c_sg_description", csg));
            p.startGroup("bottle", "for drinking");
                p.registerParam(new ParameterNumeric<int>("d", "d_description", d));
            p.finishGroup();
        p.finishGroup();
    p.finishGroup();
    return p;
}

TEST(ParameterIO, validFile){
    int a = 0, b = 0 , c = 0, csg = 0, d = 0;
    ParameterIO p = paramIO(a, b, c, csg, d);
    p.printParam(stdout, "", true);
    EXPECT_EQ(p.parse("../data/inputValid.txt"), ParameterIO::ParseResult::Ok);
    // p.printParam(stdout);
    EXPECT_EQ(a, 123248);
    EXPECT_EQ(b, 12);
    EXPECT_EQ(c, 42);
    EXPECT_EQ(csg, -12);
    EXPECT_EQ(d, 23);
}

TEST(ParameterIO, invalidFile){
    ParameterIO p;
    int a = 0;
    p.startGroup("g", "description");
    p.registerParam(new ParameterNumeric<int>("a", "", a));
    p.finishGroup();
    EXPECT_EQ(p.parse("../data/inputInvalid.txt"), ParameterIO::ParseResult::Fail);

    // p.printDiagnose(stdout);
}

TEST(ParameterIO, validArgv){
    ParameterIO p;
    int a = 0;
    p.startGroup("g", "description");
    p.registerParam(new ParameterNumeric<int>("a", "", a));
    p.finishGroup();
    const char * args[] = { "dummyname", "-g.a = 12", NULL };
    int argc = sizeof(args)/sizeof(args[0]) - 1;
    bool help = false;
    EXPECT_EQ(p.parse(argc, args, "../data/emptyInput.txt"), ParameterIO::ParseResult::Ok);
    EXPECT_EQ(a, 12);
}

TEST(ParameterIO, loadFileArgv){
    int a = 0, b = 0 , c = 0, csg = 0, d = 0;
    ParameterIO p = paramIO(a, b, c, csg, d);
    const char * args[] = { "dummyname", "../data/inputValid.txt", "-c=-1", NULL };
    int argc = sizeof(args)/sizeof(args[0]) - 1;
    EXPECT_EQ(p.parse(argc, args), ParameterIO::ParseResult::Ok);
    EXPECT_EQ(a, 123248);
    EXPECT_EQ(c, -1);
}

TEST(ParameterIO, help){
    int a = 0, b = 0 , c = 0, csg = 0, d = 0;
    ParameterIO p = paramIO(a, b, c, csg, d);
    const char * args[] = { "dummyname", "-help", "group.subgroup"};
    int argc = sizeof(args)/sizeof(args[0]);
    EXPECT_EQ(p.parse(argc, args), ParameterIO::ParseResult::HelpRequested);
}
