#include "gtest/gtest.h"

#include "util/unitConversion.h"

namespace{

namespace au = atomicUnits;

TEST(UnitConversion, BaseConversion){
    double eps = 1e-3;
    EXPECT_NEAR(1, 24.2 / au::to_as(1), eps);
    EXPECT_NEAR(1, 0.52917e-10 / au::to_m(1), eps);
    EXPECT_NEAR(1, 3.51e16 / au::to_W_cm2(1), eps);
    EXPECT_NEAR(1, 5.14e11/ au::to_V_m(1), eps);
    EXPECT_NEAR(1, 27.2 / au::to_eV(1), eps);
    EXPECT_NEAR(1, 137 / au::speedOfLight, eps);
    EXPECT_NEAR(1, au::from_K(315775.13), 1e-3);
}


TEST(UnitConversion, InversionTime){
    double v = 42;
    EXPECT_FLOAT_EQ(v, au::from_zs(au::to_zs(v)));
    EXPECT_FLOAT_EQ(v, au::from_as(au::to_as(v)));
    EXPECT_FLOAT_EQ(v, au::from_fs(au::to_fs(v)));
    EXPECT_FLOAT_EQ(v, au::from_us(au::to_us(v)));
    EXPECT_FLOAT_EQ(v, au::from_ms(au::to_ms(v)));
    EXPECT_FLOAT_EQ(v, au::from_s(au::to_s(v)));

    EXPECT_FLOAT_EQ(v*1e-21, au::to_s(au::from_zs(v)));
    EXPECT_FLOAT_EQ(v*1e-18, au::to_s(au::from_as(v)));
    EXPECT_FLOAT_EQ(v*1e-15, au::to_s(au::from_fs(v)));
    EXPECT_FLOAT_EQ(v*1e-12, au::to_s(au::from_ps(v)));
    EXPECT_FLOAT_EQ(v*1e-9, au::to_s(au::from_ns(v)));
    EXPECT_FLOAT_EQ(v*1e-6, au::to_s(au::from_us(v)));
    EXPECT_FLOAT_EQ(v*1e-3, au::to_s(au::from_ms(v)));
}

TEST(UnitConversion, InversionLength){
    double v = 42;
    EXPECT_FLOAT_EQ(v, au::from_A(au::to_A(v)));
    EXPECT_FLOAT_EQ(v, au::from_nm(au::to_nm(v)));
    EXPECT_FLOAT_EQ(v, au::from_um(au::to_um(v)));
    EXPECT_FLOAT_EQ(v, au::from_mm(au::to_mm(v)));
    EXPECT_FLOAT_EQ(v, au::from_m(au::to_m(v)));

    EXPECT_FLOAT_EQ(v*1e-1, au::to_nm(au::from_A(v)));
    EXPECT_FLOAT_EQ(v*1e-4, au::to_um(au::from_A(v)));
    EXPECT_FLOAT_EQ(v*1e-7, au::to_mm(au::from_A(v)));
    EXPECT_FLOAT_EQ(v*1e-10, au::to_m(au::from_A(v)));
}

TEST(UnitConversion, InversionIntensity){
    double v = 42;
    EXPECT_FLOAT_EQ(v, au::from_W_cm2(au::to_W_cm2(v)));
    EXPECT_FLOAT_EQ(v, au::from_W_m2(au::to_W_m2(v)));

    EXPECT_FLOAT_EQ(v*(100*100), au::to_W_m2(au::from_W_cm2(v)));
}

TEST(UnitConversion, InversionFieldStrength){
    double v = 42;
    EXPECT_FLOAT_EQ(v, au::from_V_nm(au::to_V_nm(v)));
    EXPECT_FLOAT_EQ(v, au::from_V_um(au::to_V_um(v)));
    EXPECT_FLOAT_EQ(v, au::from_V_mm(au::to_V_mm(v)));
    EXPECT_FLOAT_EQ(v, au::from_V_m(au::to_V_m(v)));

    EXPECT_FLOAT_EQ(v*1e3, au::to_V_um(au::from_V_nm(v)));
    EXPECT_FLOAT_EQ(v*1e6, au::to_V_mm(au::from_V_nm(v)));
    EXPECT_FLOAT_EQ(v*1e9, au::to_V_m(au::from_V_nm(v)));
}

TEST(UnitConversion, InversionEnergy){
    double v = 42;
    EXPECT_FLOAT_EQ(v, au::from_eV(au::to_eV(v)));
    EXPECT_FLOAT_EQ(v, au::from_meV(au::to_meV(v)));

    EXPECT_FLOAT_EQ(v*1e3, au::to_meV(au::from_eV(v)));
}

TEST(UnitConversion, InversionTemperature){
    double v = 42;
    EXPECT_FLOAT_EQ(v, au::to_K(au::from_K(v)));
}

}
