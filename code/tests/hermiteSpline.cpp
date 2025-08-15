#include "gtest/gtest.h"
#include "math.h"

#include "util/HermiteSpline.hpp"

namespace{

TEST(Spline, Hermite1D_quadratic_boundary){
    unsigned N = 10;
    double a = 12, b = 125.5, c = -9;
    double x0 = 0.2, dx = 0.7;
    std::vector<double> data(N+2);
    for(int i=0; i<N+2; i++){
        double x = x0 + (i-1) * dx;
        data[i] = a * x *x + b * x + c;
    }
    RegularHermiteSpline1D<double> spline(&data[1], N, x0, dx, data[0], data[N+1]);
    for(double t: {-2.5, -1.0, 0.0, 0.5, 1.0, 2.0, 2.5, 2.8, 4.5, N-2.0, N-1.5,  N-1.0, N+10.0}){
        double x = x0 + t * dx;
        if (t <= -1){
            EXPECT_NEAR(spline.eval(x), data[0], 1e-10);
            continue;
        }
        if (t>=N){
            EXPECT_NEAR(spline.eval(x), data[N+1], 1e-10);
            continue;
        }
        EXPECT_NEAR(spline.eval(x), a * x *x + b*x + c, 1e-10);
    }
    FILE * f = fopen("spl.txt", "w");
    double t = x0 - 5 * dx;
    while ( t< x0+ (N+5)*dx){
        t += dx/100;
        fprintf(f, "%lg %lg\n", t, spline.eval(t));
    }
    fclose(f);
}

TEST(Spline, Hermite1D_cubic_boundary){
    unsigned N = 10;
    double a = 12, b = 125.5, c = -9, d = 0.13;
    double x0 = 0.2, dx = 0.7;
    std::vector<double> data(N+2);
    for(int i=0; i<N+2; i++){
        double x = x0 + (i-1) * dx;
        data[i] = a *x*x*x + b * x*x + c * x + d;
    }
    RegularHermiteSpline1D<double> spline(&data[1], N, x0, dx, data[0], data[N+1]);
    double eps = 1e-4;
    for(double t: { 0.0, 1-eps, 2-eps }){
        double x = x0 + t * dx;
        EXPECT_NEAR(spline.eval(x), a * x*x*x + b*x*x + c*x+d, eps * a);
    }
    FILE * f = fopen("spl2.txt", "w");
    double t = x0 - 5 * dx;
    while ( t< x0+ (N+5)*dx){
        t += dx/100;
        fprintf(f, "%lg %lg\n", t, spline.eval(t));
    }
    fclose(f);
}


TEST(Spline, Hermite1D_periodic){
    unsigned N = 10;
    std::vector<double> data(N);
    double dx = 2*M_PI / N;
    for(unsigned i=0; i<N; i++)
        data[i] = sin(dx*i) + 4.5 * cos(dx * (2*i));
    PeriodicHermiteSpline<double, 1> spline(&data[0], N, 0, dx);

    for(double t=-4; t<3; t+=0.01)
        EXPECT_NEAR(spline.eval(2*M_PI+t), spline.eval(t), 1e-10);

    FILE * f = fopen("splP.txt", "w");
    double t = -5 *dx;
    while ( t< (N+5)*dx){
        t += dx/10;
        fprintf(f, "%lg %lg\n", t, spline.eval(t));
    }
    fclose(f);
}

TEST(Spline, Hermite2D_periodic){

    unsigned N = 10, M = 19;

    std::vector<double> data(N*M);
    std::array<double, 2> x0 = {0.1, 0.4};
    std::array<double, 2> dx = { 2* M_PI/N, 2*M_PI/M};
    std::array<unsigned, 2> size = {N, M};

    for(unsigned x=0; x<N; x++)
        for(unsigned y=0; y<M; y++){
            double rx = x * dx[0] + x0[0];
            double ry = y * dx[1] + x0[1];
            data[x*M+y] = sin(rx) * cos(ry) + sin(ry);
        }
    PeriodicHermiteSpline<double, 2> spline(&data[0], size, x0, dx);

    FILE * f = fopen("splP2.txt", "w");
    for(double tx=-4; tx<4; tx+=0.01)
        for(double ty=-4; ty<4; ty+=0.01){
            fprintf(f, "%lg %lg %lg\n", tx, ty, spline.eval({tx, ty}));
            EXPECT_NEAR(sin(tx) * cos(ty) + sin(ty) , spline.eval({tx, ty}), 6e-3);
        }
    fclose(f);
}

TEST(Spline, Hermite3D_periodic){


    std::array<unsigned, 3> size = {31, 29, 23};
    std::array<double, 3> x0 = {0.1, 0.4, 0.6};
    std::array<double, 3> dx = { 2* M_PI/size[0], 2*M_PI/size[1], 2*M_PI/size[2]};
    std::vector<double> data(size[0]*size[1]*size[2]);

    for(unsigned x=0; x<size[0]; x++)
        for(unsigned y=0; y<size[1]; y++)
            for(unsigned z=0; z<size[2]; z++){
                double rx = x * dx[0] + x0[0];
                double ry = y * dx[1] + x0[1];
                double rz = z * dx[2] + x0[2];
                data[x*size[1]*size[2]+y*size[2]+z] = sin(rx) *  cos(ry) + sin(ry)*sin(2*rz) - sin(rz);
        }
    PeriodicHermiteSpline<double, 3> spline(&data[0], size, x0, dx);

    FILE * f = fopen("splP3.txt", "w");
    for(double tx=-4; tx<4; tx+=0.05)
        for(double ty=-4; ty<4; ty+=0.05)
          for(double tz=-4; tz<4; tz+=0.05){
                fprintf(f, "%lg %lg %lg %lg\n", tx, ty, tz, spline.eval({tx, ty, tz}));
                EXPECT_NEAR(sin(tx) * cos(ty) + sin(ty)*sin(2*tz) - sin(tz) , spline.eval({tx, ty, tz}), 6e-3);
        }
    fclose(f);
}


}
