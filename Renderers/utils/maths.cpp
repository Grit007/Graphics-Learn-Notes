#include "maths.h"

template <> template <> vec<3,int>  ::vec(const vec<3,float> &v) : x(int(v.x+.5f)),y(int(v.y+.5f)),z(int(v.z+.5f)) {}
template <> template <> vec<3,float>::vec(const vec<3,int> &v)   : x(v.x),y(v.y),z(v.z) {}
template <> template <> vec<2,int>  ::vec(const vec<2,float> &v) : x(int(v.x+.5f)),y(int(v.y+.5f)) {}
template <> template <> vec<2,float>::vec(const vec<2,int> &v)   : x(v.x),y(v.y) {}

float Clamp(float x, float min, float max)
{
    x = (x>max?max:(x<min?min:x));
    return x;
//    if (x > max)
//        return max;
//    if (x < min)
//        return min;
//    return x;
}

Vec3f Mix(Vec3f a, Vec3f b, float weight){
    Vec3f res;
    res.x =  a.x * (1 - weight) + b.x * weight;
    res.y =  a.y * (1 - weight) + b.y * weight;
    res.z =  a.z * (1 - weight) + b.z * weight;
    return res;
}

float Length(Vec3f a){
    return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

Vec3f Vecmul(Vec3f a, Vec3f b){
    return Vec3f(a.x*b.x, a.y*b.y, a.z*b.z);
}