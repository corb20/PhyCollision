#pragma once
#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Dense>

using namespace Eigen;

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

    Triangle() : a(Vector3d(0, 0, 0)), b(Vector3d(0, 0, 0)), c(Vector3d(0, 0, 0)) {}
    Triangle(Vector3d a, Vector3d b, Vector3d c) : a(a), b(b), c(c) {}

    Vector3d normal()
    {
        return (b - a).cross(c - a).normalized();
    }

    double GetClosestDist(Vector3d point);
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
};