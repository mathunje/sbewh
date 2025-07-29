#ifndef SBE_WH_GEOM_VECTOR_HPP
#define SBE_WH_GEOM_VECTOR_HPP

#include<initializer_list>
#include<array>
#include<vector>
#include<cassert>

#include<math.h>

template<typename T, unsigned N>
class GeomVector{
private:
    std::array<T, N> v;
public:

    GeomVector(){
        for(unsigned i=0; i<N; i++)
            v[i] = 0;
    }

    GeomVector(std::vector<T> data){
        assert(v.size() == N);
        for(unsigned i=0; i<N; i++)
            v[i] = data[i];
    }

    /*
    GeomVector(const std::array<T, N> data){
        v = data;
    } */

    GeomVector(std::initializer_list<T> data){
        assert(data.size() == N);
        unsigned i = 0;
        for(T d : data){
            v[i++] = d;
        }
    }

    GeomVector(const T data[N]){
        for(unsigned i=0; i<N; i++)
            v[i] = data[i];
    }

    GeomVector(const std::array<T, N> data){
        v = data;
    }

    /*GeomVector(const GeomVector<T, N> &gv){
        for(unsigned i=0; i<N; i++)
            v[i] = gv.at(i);
    }*/

    T & operator[](unsigned index) {
        return v[index];
    }
    T at(unsigned index) const{
        return v[index];
    }

    void operator+=(const GeomVector &b){
        for(unsigned i=0; i<N; i++)
            v[i] += b.at(i);
    }
    /*
    void operator+=(const GeomVector b){
        for(unsigned i=0; i<N; i++)
            v[i] += b.at(i);
    }*/
    void operator+=(const T c){
        for(unsigned i=0; i<N; i++)
            v[i] += c;
    }
    GeomVector operator+(const GeomVector &b) const{
        std::array<T, N> res;
        for(unsigned i=0; i<N; i++)
            res[i] = v[i] + b.at(i);
        return res;
    }
    /*
    GeomVector operator+(const GeomVector b){
        std::array<T, N> res;
        for(unsigned i=0; i<N; i++)
            res[i] += v[i] + b.at(i);
        return GeomVector<T, N>(res);
    } */
    GeomVector operator+(const T c) const{
        std::array<T, N> res;
        for(unsigned i=0; i<N; i++)
            res[i] = v[i] + c;
        return GeomVector<T, N>(res);
    }

    void operator-=(const GeomVector &b){
        for(unsigned i=0; i<N; i++)
            v[i] -= b.at(i);
    }
    /*
    void operator-=(const GeomVector b){
        for(unsigned i=0; i<N; i++)
            v[i] -= b.at(i);
    } */
    void operator-=(const T c){
        for(unsigned i=0; i<N; i++)
            v[i] -= c;
    }
    GeomVector operator-(const GeomVector &b) const{
        std::array<T, N> res;
        for(unsigned i=0; i<N; i++)
            res[i] = v[i] - b.at(i);
        return res;
    }
    /*
    GeomVector operator-(const GeomVector b){
        std::array<T, N> res;
        for(unsigned i=0; i<N; i++)
            res[i] += v[i] - b.at(i);
        return GeomVector(res);
    }*/

    GeomVector operator-() const{
        std::array<T, N> res;
        for(unsigned i=0; i<N; i++)
            res[i] = -v[i];
        return GeomVector(res);
    }

    GeomVector operator-(const T c) const{
        std::array<T, N> res;
        for(unsigned i=0; i<N; i++)
            res[i] = v[i] - c;
        return GeomVector(res);
    }

    GeomVector operator*(const T c) const{
        std::array<T, N> res;
        for(unsigned i=0; i<N; i++)
            res[i] = v[i] * c;
        return GeomVector(res);
    }
    void operator*=(const T c){
        for(unsigned i=0; i<N; i++)
            v[i] *= c;
    }

    GeomVector operator/(const T c) const{
        std::array<T, N> res;
        for(unsigned i=0; i<N; i++)
            res[i] = v[i] / c;
        return GeomVector(res);
    }
    void operator/=(const T c){
        for(unsigned i=0; i<N; i++)
            v[i] /= c;
    }

    T dot(const GeomVector &b) const{
        T sum = 0;
        for(unsigned i=0; i<N; i++)
            sum += v[i] * b.at(i);
        return sum;
    }

    T norm2() const {
        return dot(*this);
    }
    T norm() const {
        return sqrt(norm2());
    }

    GeomVector & normalize()
    {
        T n = norm();
        for(unsigned i=0; i<N; i++)
            v[i] /= n;
        return *this;
    }

    GeomVector cross(const GeomVector &b) const {
        assert(N == 3);
        return GeomVector { v[1] * b.at(2) - v[2] * b.at(1),
                            v[2] * b.at(0) - v[0] * b.at(2),
                            v[0] * b.at(1) - v[1] * b.at(0)};
    }

    GeomVector getParallel(const GeomVector &direction) const{
        return direction * dot(direction);
    }

    std::array<T, N> getArray() const{
        return v;
    }

    bool operator<(const GeomVector &b) const{
        for(unsigned i=0; i<N-1; i++)
            if (v[i] != b.at(i) )
                return v[i] < b.at(i);
        return v[N-1] < b.at(N-1);
    }

    bool operator==(const GeomVector &b) const{
        for(unsigned i=0; i<N; i++)
            if ( v[i] != b.at(i))
                return false;
        return true;
    }

    bool operator!=(const GeomVector &b) const{
        return !(*this == b);
    }

};


template<typename T, unsigned N>
GeomVector<T, N> operator-(const T c, const GeomVector<T, N> &v){
    return -(v - c);
}

template<typename T, unsigned N>
GeomVector<T, N> operator+(const T c, const GeomVector<T, N> &v){
    return v + c;
}

template<typename T, unsigned N>
GeomVector<T, N> operator*(const T c, const GeomVector<T, N> &v){
    return v * c;
}


typedef GeomVector<double, 3> GeomVector3d;


typedef GeomVector<int, 3> CellIndex;
template<>
struct std::hash<CellIndex>
{
    std::size_t operator()(const CellIndex & ci) const noexcept {
        std::size_t h0 = std::hash<int>{}(ci.at(0));
        std::size_t h1 = std::hash<int>{}(ci.at(1));
        std::size_t h2 = std::hash<int>{}(ci.at(2));
        return h0 ^ (h1 << 1) ^ (h2 << 2);
    }
};

#endif
