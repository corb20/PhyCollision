#include "PhyBaseStruct.h"

double Segment::GetClosestDist(Vector3d point){
    Vector3d ab = end - start;
    Vector3d e_ab = ab.normalized();
    Vector3d ap = point - start;
    Vector3d ans = ap- e_ab.dot(ap)*e_ab;
    return ans.norm();
}

double Triangle::GetClosestDist(Vector3d point){
    Vector3d ab = b - a;
    Vector3d ac = c - a;
    Vector3d ap = point - a;
    Vector3d bp = point - b;
    Vector3d cp = point - c;
    Vector3d n = ab.cross(ac);
    Vector3d n1 = ab.cross(n);
    Vector3d n2 = n.cross(ac);
    if(n.dot(ap) > 0 && n1.dot(bp) > 0 && n2.dot(cp) > 0){
        return n.norm();
    }
    else{
        double d1 = Segment(a, b).GetClosestDist(point);
        double d2 = Segment(b, c).GetClosestDist(point);
        double d3 = Segment(c, a).GetClosestDist(point);
        return std::min(std::min(d1, d2), d3);
    }
}