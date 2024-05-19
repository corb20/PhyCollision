#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>

using namespace Eigen;


double IsEqual(double a, double b);

class Segment
{
public:
    Vector3d start;
    Vector3d end;

    Segment() : start(Vector3d(0, 0, 0)), end(Vector3d(0, 0, 0)) {}
    Segment(Vector3d start, Vector3d end) : start(start), end(end) {}

    double GetClosestDist(Vector3d point);
};

class Ray
{
public:
    Vector3d org;
    Vector3d dir;

    Ray() : org(Vector3d(0, 0, 0)), dir(Vector3d(0, 0, 0)) {}
    Ray(Vector3d start, Vector3d direction) : org(start), dir(direction) {}
};

class Triangle
{
public:
    //默认方向为逆时针
    Vector3d a;
    Vector3d b;
    Vector3d c;

    Vector3d n;

    Triangle() : a(Vector3d(0, 0, 0)), b(Vector3d(0, 0, 0)), c(Vector3d(0, 0, 0)) {
        n=Vector3d(0,0,0);
    }
    Triangle(Vector3d a, Vector3d b, Vector3d c) : a(a), b(b), c(c) {
        normal();
    }

    Vector3d normal()
    {
        n = (b - a).cross(c - a).normalized();
        return n;
    }
    Vector3d normal2Org(){
        Vector3d n1 = normal();
        if(n1.dot(a) > 0){
            return -n1;
        }
        else{
            return n1;
        }
    }

    double GetClosestDist(Vector3d point);
    double GetClosestDistWithotRange(Vector3d point);
    Vector3d Barycentric(Vector3d p);
    bool IsInside(Vector3d p);

    bool IsValid();
};

//四面体
class Tetrahedron
{
public:
    Vector3d a;
    Vector3d b;
    Vector3d c;
    Vector3d d;

    Tetrahedron() : a(Vector3d(0, 0, 0)), b(Vector3d(0, 0, 0)), c(Vector3d(0, 0, 0)), d(Vector3d(0, 0, 0)) {}
    Tetrahedron(Vector3d a, Vector3d b, Vector3d c, Vector3d d) : a(a), b(b), c(c), d(d) {}
    //求重心坐标
    Vector4d Barycentric(Vector3d p);
    
    bool IsInside(Vector3d p);

    Triangle GetClosestTriangle(Vector3d p);

    //检测是否合法
    bool IsValid();

    bool IsAllZero(){
        return a.isZero() && b.isZero() && c.isZero() && d.isZero();
    }
};