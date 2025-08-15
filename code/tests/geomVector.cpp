#include "gtest/gtest.h"
#include "util/GeomVector.hpp"

namespace{


TEST(GeomVector, constructor){
    std::vector<double> v { 1, 2, 4};
    GeomVector3d g(v);
    for(unsigned i=0; i<v.size(); i++)
        EXPECT_EQ(v[i], g[i]);
}

TEST(GeomVector, norm){
    GeomVector3d a ({3, 4, 12});
    GeomVector3d b ({1, 2, 3});
    EXPECT_NEAR(a.norm(), 13, 1e-5);
    EXPECT_NEAR(b.norm2(), 14, 1e-5);
    a.normalize();
    EXPECT_NEAR(a.norm2(), 1, 1e-5);
}

TEST(GeomVector, vecSpecial){
    GeomVector3d a {1, 2, 3};
    GeomVector3d b {3, 2, 1};
    EXPECT_NEAR(a.dot(b), b.dot(a), 1e-5);
    EXPECT_NEAR(a.dot(b), 10, 1e-5);
    GeomVector3d c = a.cross(b);
    GeomVector3d cm = b.cross(a);
    for(unsigned i=0; i<3; i++)
        EXPECT_NEAR(c[i], -cm[i], 1e-5);
    EXPECT_NEAR(a.cross(b).dot(a), 0, 1e-5);
    EXPECT_NEAR(a.cross(b).dot(b), 0, 1e-5);
    EXPECT_NEAR(c[0], -4, 1e-5);
    EXPECT_NEAR(c[1], 8, 1e-5);
    EXPECT_NEAR(c[2], -4, 1e-5);

    GeomVector3d dir { 1, 1, 0};
    dir.normalize();
    a = {1, 1, 4};
    auto p = a.getParallel(dir);
    EXPECT_NEAR(p[0], 1, 1e-5);
    EXPECT_NEAR(p[1], 1, 1e-5);
    EXPECT_NEAR(p[2], 0, 1e-5);

}

TEST(GeomVector, scalarOp){
    GeomVector3d a {2, 4, 6};
    GeomVector3d a3m = 3.0 * a;
    GeomVector3d a3 = a * 3;
    GeomVector3d ap = a + 2;
    GeomVector3d apm = 2.0 + a;
    GeomVector3d am = a - 2.0;
    GeomVector3d amm = 2.0 - a;
    GeomVector3d aD2 = a / 2.0;
    for(unsigned i=0; i<3; i++){
        EXPECT_NEAR(a3m[i], a[i]*3, 1e-5);
        EXPECT_NEAR(a3[i], a[i]*3, 1e-5);
        EXPECT_NEAR(aD2[i], a[i]/2, 1e-5);
        EXPECT_NEAR(apm[i], a[i]+2, 1e-5);
        EXPECT_NEAR(ap[i], a[i]+2, 1e-5);
        EXPECT_NEAR(am[i], a[i]-2, 1e-5);
        EXPECT_NEAR(amm[i], 2-a[i], 1e-5);
    }
    a *= 3;
    for(unsigned i=0; i<3; i++)
        EXPECT_NEAR(a3[i], a[i], 1e-5);
    a /= 6;
    for(unsigned i=0; i<3; i++)
        EXPECT_NEAR(aD2[i], a[i], 1e-5);
    a *= 2;
    a += 2;
    for(unsigned i=0; i<3; i++)
        EXPECT_NEAR(ap[i], a[i], 1e-5);
}

TEST(GeomVector, vecArith){
    GeomVector3d a{2, 4, 6};
    GeomVector3d b{4, 7, 23};
    auto sum = a + b;
    auto diff = a - b;
    for(unsigned i=0; i<3; i++){
        EXPECT_NEAR(a[i]+b[i], sum[i], 1e-5);
        EXPECT_NEAR(a[i]-b[i], diff[i], 1e-5);
    }
}


}
