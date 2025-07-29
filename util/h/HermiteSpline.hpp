#ifndef SBE_WH_HERMITE_SPLINE_H
#define SBE_WH_HERMITE_SPLINE_H

#include<array>


template<typename T>
class _Spline1Dbase {
public:
    static inline std::array<T, 4> getCoeffs(T y_m1, T y, T y_p1, T y_p2, double cdx) {
        return {
                y,
                (y_p1 - y_m1) / 2 / cdx,
                (2 * y_m1 - 5 * y + 4 * y_p1 - y_p2) / 2 / cdx / cdx,
                (- y_m1 + 3 * y - 3 * y_p1 + y_p2) / 2 / cdx / cdx / cdx
               };
    }

    static inline T evalDirect(T y_m1, T y, T y_p1, T y_p2, double cdx, double s){
        const std::array<T, 4> c = getCoeffs(y_m1, y, y_p1, y_p2, cdx);
        return  c[0] + s * ( c[1] + s * ( c[2] + s * c[3]) );
    }
};


template<typename T>
class RegularHermiteSpline1D
{
    double x0;
    double dx;
    unsigned size;
    std::vector< std::array<T, 4>  > regionCoeffs; // sum_i x^i c[i]

    inline std::array<T, 4> getCoeffs(T y_m1, T y, T y_p1, T y_p2) const {
        return _Spline1Dbase<T>::getCoeffs(y_m1, y, y_p1, y_p2, dx);
    }
public:
    RegularHermiteSpline1D(const T * data, unsigned size, double x0, double dx,
                           T yLeft3,  T yLeft2, T yLeft1,
                           T yRight1, T yRight2, T  yRight3) : size(size), x0(x0), dx(dx)
    {
        assert(size >= 3);
        regionCoeffs.push_back(getCoeffs(yLeft1, yLeft2, yLeft3, data[0]));
        regionCoeffs.push_back(getCoeffs(yLeft1, yLeft2, data[0], data[1]));
        regionCoeffs.push_back(getCoeffs(yLeft2, data[0], data[1], data[2]));
        for(unsigned i=1; i<size-2; i++)
            regionCoeffs.push_back(getCoeffs(data[i-1], data[i], data[i+1], data[i+2]));
        unsigned li = size - 1;
        regionCoeffs.push_back(getCoeffs(data[li-2], data[li-1], data[li], yRight1));
        regionCoeffs.push_back(getCoeffs(data[li-1], data[li], yRight1, yRight2));
        regionCoeffs.push_back(getCoeffs(data[li], yRight1, yRight2, yRight3));
    }

    // constant extrapolation!
    RegularHermiteSpline1D(const T * data, unsigned size, double x0, double dx, T yLeft, T yRight)
        : RegularHermiteSpline1D(data, size, x0, dx, yLeft, yLeft, yLeft, yRight, yRight, yRight) {}

    RegularHermiteSpline1D(const T *data, unsigned size, double x0, double dx) : size(size), x0(x0), dx(dx) {}

    T eval(double x){
        double r = x - x0;
        int ind = std::floor(r / dx) + 2;
        double s = r - (ind-2) * dx;
        if (ind < 0){
            ind = 0;
            s = 0;
        }
        if ( ind > size + 2 ){
            ind = size + 2;
            s = dx;
        }
        const std::array<T, 4> &c = regionCoeffs[ind];
        return  c[0] + s * ( c[1] + s * ( c[2] + s * c[3]) );
    }
};



template<typename T, unsigned dim>
class PeriodicHermiteSpline
{
    std::array<unsigned, dim> size;
    std::array<double, dim> x0;
    std::array<double, dim> dx;
    std::vector< PeriodicHermiteSpline<T,dim-1> > splines;

    inline T eval(double *x){
        double r = x[0] - x0[0];
        int ind = std::floor(r / dx[0]);
        double s = r - ind * dx[0];
        ind %=(int)size[0];
        if ( ind < 0 )
            ind += size[0];
        T y_m1 = splines[(ind+size[0]-1)%size[0]].eval(&x[1]);
        T y = splines[ind].eval(&x[1]);
        T y_p1 = splines[(ind+1)%size[0]].eval(&x[1]);
        T y_p2 = splines[(ind+2)%size[0]].eval(&x[1]);
        return _Spline1Dbase<T>::evalDirect(y_m1, y, y_p1, y_p2, dx[0], s);
    }

    template<typename Q>
    inline std::array<Q, dim-1> subA(const std::array<Q, dim> &inp)
    {
        std::array<Q, dim-1> res;
        for(unsigned i=1; i<dim; i++)
            res[i-1] = inp[i];
        return res;
    }

    friend class PeriodicHermiteSpline<T, dim+1>;

public:

    PeriodicHermiteSpline(const T * data, const std::array<unsigned, dim> size,
                          std::array<double, dim> x0, std::array<double, dim> dx)
        : size(size), x0(x0), dx(dx)
    {
        assert(size[0]>=2);
        unsigned subSize = size[1];
        for(unsigned i=2; i<dim; i++)
            subSize *= size[i];

        for(unsigned i=0; i<size[0]; i++){
            splines.push_back( PeriodicHermiteSpline<T, dim-1> (data + i * subSize, subA(size), subA(x0), subA(dx) ) );
        }
    }


    T eval(std::array<T, dim> x) { return eval(x.data()); }
};


template<typename T, unsigned dim> class PeriodicHermiteSpline;

template<typename T>
class PeriodicHermiteSpline<T, 1>
{
    const unsigned size;
    const double x0;
    const double dx;
    std::vector< std::array<T, 4>  > regionCoeffs; // sum_i x^i c[i]

    inline std::array<T, 4> getCoeffs(T y_m1, T y, T y_p1, T y_p2) const {
        return _Spline1Dbase<T>::getCoeffs(y_m1, y, y_p1, y_p2, dx);
    }

    inline T eval(const double *x) const {  return eval(*x); }

    friend class PeriodicHermiteSpline<T, 2>;
public:


    PeriodicHermiteSpline(const T * data, unsigned size, double x0, double dx)
        : size(size), x0(x0), dx(dx)
    {
        assert(size >= 3);
        unsigned li = size - 1;
        regionCoeffs.push_back(getCoeffs(data[li], data[0], data[1], data[2]));
        for(unsigned i=1; i<size-2; i++)
            regionCoeffs.push_back(getCoeffs(data[i-1], data[i], data[i+1], data[i+2]));
        regionCoeffs.push_back(getCoeffs(data[li-2], data[li-1], data[li], data[0]));
        regionCoeffs.push_back(getCoeffs(data[li-1], data[li], data[0], data[1]));
    }

    PeriodicHermiteSpline(const T * data, const unsigned *size, const double * x0, const double * dx)
        : PeriodicHermiteSpline(data, *size, *x0, *dx) {}

    PeriodicHermiteSpline(const T * data, const std::array<unsigned, 1> size, const std::array<double, 1> x0,
                          const std::array<double, 1> dx)
        : PeriodicHermiteSpline(data, size[0], x0[0], dx[0]) {}


    inline T eval(double x) const {
        double r = x - x0;
        int ind = std::floor(r / dx);
        double s = r - ind * dx;
        ind %= (int)(size);
        if ( ind < 0 ){
            ind += size;
        }
        const std::array<T, 4> &c = regionCoeffs[ind];
        return  c[0] + s * ( c[1] + s * ( c[2] + s * c[3]) );
    }

    inline T eval(const std::array<double, 1> x) const { return eval(x[0]); }

};

#endif
