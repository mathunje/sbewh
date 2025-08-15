#include "gtest/gtest.h"
#include "math.h"

#include "pulse/PureGauss1D.h"
#include "pulse/GaussianPulse1D.h"
#include "pulse/Sin2Pulse1D.h"
#include "pulse/Sin2RampPulse1D.h"
#include "pulse/Pulse3D.h"

namespace{

TEST(Pulse, PureGauss){
    double E0 = atomicUnits::from_V_nm(1);
    double fwhm = atomicUnits::from_fs(20);
    PureGauss1D pulse(E0, fwhm, -5, 5);
    double timeDiff = pulse.getEndTime() - pulse.getStartTime();
    double err = sqrt(pulse.numericalIntegrationErrorE()) / (E0 * timeDiff);
    EXPECT_NEAR(err, 0, 1e-5);
    EXPECT_NEAR(pulse.getE(pulse.getStartTime()) / E0, 0, 1e-10);
    EXPECT_NEAR(pulse.getE(pulse.getEndTime()) / E0, 0, 1e-10);
    EXPECT_NEAR(pulse.getA0(), E0 * fwhm * sqrt(M_PI / 4 / log(2)) , 1e-8);
}

TEST(Pulse, GaussianPulse){
    double E0 = atomicUnits::from_V_nm(1);
    double fwhm = atomicUnits::from_fs(20);
    double omega = 2 * M_PI * atomicUnits::speedOfLight / atomicUnits::from_nm(5000);
    GaussianPulse1D l1(E0, omega, 0, fwhm, -5, 5);
    double timeDiff = l1.getEndTime() - l1.getStartTime();
    double err = sqrt(l1.numericalIntegrationErrorE()) / (E0 * timeDiff);
    EXPECT_NEAR(err, 0, 1e-5);
    EXPECT_NEAR(l1.getE(0) / E0, 1, 1e-8);
    EXPECT_NEAR(l1.getE(l1.getStartTime()) / E0, 0, 1e-10);
    EXPECT_NEAR(l1.getE(l1.getEndTime()) / E0, 0, 1e-10);
    EXPECT_NEAR(l1.numericalIntegrationE() / E0 / timeDiff, 0, 1e-6);
    EXPECT_NEAR(l1.getA0(), E0/ omega, 1e-8);
}

TEST(Pulse, Sin2Pulse){
    double E0 = atomicUnits::from_V_nm(1);
    double fwhm = atomicUnits::from_fs(20);
    double omega = 2 * M_PI * atomicUnits::speedOfLight / atomicUnits::from_nm(5000);
    Sin2Pulse1D l1(E0, omega, 0, fwhm, -2.5, 2.5);
    double timeDiff = l1.getEndTime() - l1.getStartTime();
    double err = sqrt(l1.numericalIntegrationErrorE()) / (E0 * timeDiff);
    EXPECT_NEAR(err, 0, 1e-5);
    EXPECT_NEAR(l1.getE(0) / E0, 1, 1e-8);
    EXPECT_NEAR(l1.getE(l1.getStartTime()) / E0, 0, 1e-10);
    EXPECT_NEAR(l1.getE(l1.getEndTime()) / E0, 0, 1e-10);
    EXPECT_NEAR(l1.numericalIntegrationE() / E0 / timeDiff, 0, 1e-6);
    EXPECT_NEAR(l1.getA0(), E0/omega, 1e-8);
}

TEST(Pulse, Sin2RampPulse){
    double E0 = atomicUnits::from_V_nm(1);
    double risingCycles = 2.3;
    double omega = 2 * M_PI * atomicUnits::speedOfLight / atomicUnits::from_nm(5000);
    Sin2RampPulse1D l1(E0, omega, risingCycles, 4);
    double timeDiff = l1.getEndTime() - l1.getStartTime();
    double err = sqrt(l1.numericalIntegrationErrorE()) / (E0 * timeDiff);
    EXPECT_NEAR(err, 0, 1e-4);
    EXPECT_NEAR(l1.getE(0) / E0, 1, 1e-8);
    EXPECT_NEAR(l1.getE(l1.getStartTime()) / E0, 0, 1e-10);
    EXPECT_NEAR(l1.getE(l1.getEndTime()) / E0, 0, 1e-10);
    EXPECT_NEAR(l1.numericalIntegrationE() / E0 / timeDiff, 0, 1e-6);
    EXPECT_NEAR(l1.getA0(), E0/omega, 1e-8);

    Sin2RampPulse1D l2(E0, omega, risingCycles, 5);
    err = sqrt(l2.numericalIntegrationErrorE()) / (E0 * timeDiff);
    EXPECT_NEAR(err, 0, 1e-4);
}

TEST(Pulse, Pulse3D){
    double E0 = atomicUnits::from_V_nm(1);
    double fwhm = atomicUnits::from_fs(20);
    double tShift = 100;
    double omega = 2 * M_PI * atomicUnits::speedOfLight / atomicUnits::from_nm(5000);
    auto l1 = std::shared_ptr<Pulse1D>(new Sin2Pulse1D (E0, omega, 0, fwhm, -2.5, 2.5));
    auto l2 = std::shared_ptr<Pulse1D>(new Sin2Pulse1D (E0, 2*omega, 0, 2*fwhm, -2.5, 2.5));
    auto l3 = std::shared_ptr<Pulse1D>(new Sin2Pulse1D (E0, 2*omega, 0, 3*fwhm, -2.5, 2.5));
    Pulse3D pulse;
    pulse.addPulse(GeomVector3d({0, 0, 2}), l1);
    pulse.addPulse(GeomVector3d({0, 1, 0}), l2);
    pulse.addPulse(GeomVector3d({1, 0, 0}), l3, -tShift);
    pulse.addPulse(GeomVector3d({1, 0, 0}), l1, tShift);

    pulse.print(stdout);

    double t = 122;
    GeomVector3d E = pulse.getE(t);
    EXPECT_EQ(E[0], l3->getE(t + tShift) + l1->getE(t - tShift) );
    EXPECT_EQ(E[1], l2->getE(t));
    EXPECT_EQ(E[2], l1->getE(t));

    double start = std::min({l1->getStartTime() + tShift, l1->getStartTime(), l2->getStartTime(), l3->getStartTime() - tShift});
    double end = std::max({l1->getEndTime() + tShift, l1->getEndTime(), l2->getEndTime(), l3->getEndTime() - tShift});
    EXPECT_EQ(start, pulse.getStartTime());
    EXPECT_EQ(end, pulse.getEndTime());
}


}
