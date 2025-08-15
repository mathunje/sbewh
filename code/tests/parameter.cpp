#include "gtest/gtest.h"
#include "io/parameter/ParameterBool.h"
#include "io/parameter/ParameterNumeric.hpp"
#include "io/parameter/ParameterVector.hpp"
#include "io/parameter/ParameterString.h"

namespace {

TEST(Parameter, General){
    bool v = true;
    std::string n = "name";
    std::string d = "description";
    ParameterBool p(n, d, v);
    EXPECT_EQ(p.getName(), n);
    EXPECT_EQ(p.getDescription(), d);
    EXPECT_EQ(p.hasError(), false);
}

TEST(Parameter, String){
    std::string test = "";
    ParameterString p("", "", test);
    p.parse("testStr");
    EXPECT_EQ(p.hasError(), false);
    EXPECT_EQ(test.compare("testStr"), 0);
}

TEST(Parameter, Bool){
    bool v = true;
    ParameterBool p("", "", v, false);
    EXPECT_EQ(v, false);
    p.parse("true");
    EXPECT_EQ(p.hasError(), false);
    EXPECT_EQ(v, true);

    p.parse("false");
    EXPECT_EQ(p.hasError(), false);
    EXPECT_EQ(v, false);

    p.parse("1");
    EXPECT_EQ(p.hasError(), false);
    EXPECT_EQ(v, true);

    p.parse("0");
    EXPECT_EQ(p.hasError(), false);
    EXPECT_EQ(v, false);

    p.parse("invalid");
    EXPECT_EQ(p.hasError(), true);
}

TEST(Parameter, Numeric){
    int value = 42;
    ParameterNumeric<int> p("", "", value);
    p.parse("133");
    EXPECT_EQ(p.hasError(), false);
    EXPECT_EQ(value, 133);
    p.setMinValue(0);
    p.parse("-120");
    EXPECT_EQ(value, 133);
    EXPECT_EQ(p.hasError(), true);
    p.parse("0");
    EXPECT_EQ(value, 0);
    EXPECT_EQ(p.hasError(), false);

    p.setMaxValue(10);
    p.parse("12");
    EXPECT_EQ(p.hasError(), true);
    EXPECT_EQ(value, 0);

    ParameterNumeric<int> q("", "", value);
    q.setRange(-10, 10);
    q.parse("-23");
    EXPECT_EQ(q.hasError(), true);
    EXPECT_EQ(value, 0);

    q.parse("23");
    EXPECT_EQ(q.hasError(), true);
    EXPECT_EQ(value, 0);

    q.parse("-1");
    EXPECT_EQ(q.hasError(), false);
    EXPECT_EQ(value, -1);
}

TEST(Parameter, Int){
    int v = 0;
    ParameterNumeric<int> p("", "", v);
    p.parse("42");
    EXPECT_EQ(v, 42);
    p.parse("4.2");
    EXPECT_EQ(v, 4);
    p.parse("spam");
    EXPECT_EQ(p.hasError(), true);
}

TEST(Parameter, Uint){
    unsigned int v = 0;
    ParameterNumeric<unsigned int> p("", "", v);
    p.parse("42");
    EXPECT_EQ(v, 42);
    p.parse("-12");
    unsigned int res = -12;
    EXPECT_EQ(v, res);
    p.parse("spam");
    EXPECT_EQ(p.hasError(), true);
}

TEST(Parameter, LLint){
    long long int v = 0;
    ParameterNumeric<long long int> p("", "", v);
    p.parse("42");
    EXPECT_EQ(v, 42);
    p.parse("4.2");
    EXPECT_EQ(v, 4);
    p.parse("spam");
    EXPECT_EQ(p.hasError(), true);
}

TEST(Parameter, ULLint){
    unsigned long long int v = 0;
    ParameterNumeric<unsigned long long int> p("", "", v);
    p.parse("42");
    EXPECT_EQ(v, 42);
    p.parse("-12");
    unsigned long long int res = -12;
    EXPECT_EQ(v, res);
    p.parse("spam");
    EXPECT_EQ(p.hasError(), true);
}

TEST(Parameter, Float){
    float v = 0;
    ParameterNumeric<float> p("", "", v);
    p.parse("1.2");
    EXPECT_NEAR(v, 1.2, 1e-6);
    p.parse("spam");
    EXPECT_EQ(p.hasError(), true);
}

TEST(Parameter, Double){
    double v = 0;
    ParameterNumeric<double> p("", "", v);
    p.parse("1.2");
    EXPECT_NEAR(v, 1.2, 1e-12);
    p.parse("spam");
    EXPECT_EQ(p.hasError(), true);
}

TEST(Parameter, Vector){
    std::vector<int> vec;
    ParameterVector<int> p("", "", vec);
    p.parse("1, 5, -4");
    ASSERT_EQ(vec.size(), 3);
    EXPECT_EQ(vec[0], 1);
    EXPECT_EQ(vec[1], 5);
    EXPECT_EQ(vec[2], -4);
}


}
