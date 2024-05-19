#include "PhyBaseStruct.h"

double bias = 1e-6;

double IsEqual(double a, double b){
    return abs(a - b) < bias;
}

double Segment::GetClosestDist(Vector3d point){
    Vector3d ab = end - start;
    Vector3d e_ab = ab.normalized();
    Vector3d ap = point - start;
    Vector3d ans = ap- e_ab.dot(ap)*e_ab;
    return ans.norm();
}

bool Triangle::IsValid(){
    //检验三角形是否合法,先检测是否有重合点
    if((a - b).norm() < bias || (a - c).norm() < bias || (b - c).norm() < bias){
        return false;
    }
    //检测是否共线
    Vector3d ab = b - a;
    Vector3d ac = c - a;
    if(IsEqual(ab.cross(ac).norm(), 0)){
        return false;
    }
    return true;
}

Vector3d Triangle::Barycentric(Vector3d p){
    Vector3d ab = b - a;
    Vector3d ac = c - a;
    Vector3d ap = p - a;
    Vector3d n = this->normal();

    Vector2d map_ab = Vector2d(ab.x(),ab.y());
    Vector2d map_ac = Vector2d(ac.x(),ac.y());
    Vector2d map_ap = Vector2d(ap.x(),ap.y());
    if(IsEqual(n.z(), 0)){
        map_ab = Vector2d(ab.x(),ab.z());
        map_ac = Vector2d(ac.x(),ac.z());
        map_ap = Vector2d(ap.x(),ap.z());
    }
    double S = map_ab.x()*map_ac.y() - map_ab.y()*map_ac.x();
    //计算重心坐标
    double v = (map_ap.x()*map_ac.y() - map_ap.y()*map_ac.x())/S;
    double w = (map_ab.x()*map_ap.y() - map_ab.y()*map_ap.x())/S;
    double u = 1 - v - w;
    return Vector3d(u, v, w);
}

bool Triangle::IsInside(Vector3d p){
    Vector3d bary = Barycentric(p);
    return bary.x() >= 0 && bary.y() >= 0 && bary.z() >= 0;
}

double Triangle::GetClosestDist(Vector3d point){
    Vector3d ab = b - a;
    Vector3d ac = c - a;
    Vector3d ap = point - a;
    Vector3d n = ab.cross(ac).normalized();
    Vector3d pn = ap - n.dot(ap)*n;
    if(IsInside(pn)){
        return abs(n.dot(ap));
    }
    else{
        double d1 = Segment(a, b).GetClosestDist(point);
        double d2 = Segment(b, c).GetClosestDist(point);
        double d3 = Segment(c, a).GetClosestDist(point);
        return std::min(std::min(d1, d2), d3);
    }
}
double Triangle::GetClosestDistWithotRange(Vector3d point){
    Vector3d ab = b - a;
    Vector3d ac = c - a;
    Vector3d ap = point - a;
    Vector3d n = ab.cross(ac).normalized();
    return abs(n.dot(ap));
}

bool Tetrahedron::IsValid(){
    //检验四面体是否合法
    //检测abc是否共线
    Triangle abc = Triangle(a,b,c);
    if(abc.IsValid() == false){
        return false;
    }
    //检测d是否在abc所在的平面上
    Vector3d n = abc.normal();
    if(IsEqual(n.dot(d - a), 0)){
        return false;
    }
    return true;
}

Triangle Tetrahedron::GetClosestTriangle(Vector3d p){
    Triangle trList[] = {Triangle(a,b,c),Triangle(a,b,d),Triangle(a,c,d),Triangle(b,c,d)};
    double min_dist = trList[0].GetClosestDist(p);
    Triangle ans = trList[0];
    for(auto tr:trList){
        double dist = tr.GetClosestDist(p);
        if(dist < min_dist){
            min_dist = dist;
            ans = tr;
        }
    }
    return ans;
}

Vector4d Tetrahedron::Barycentric(Vector3d p){
    //计算重心坐标
    Vector3d ab = b - a;
    Vector3d ac = c - a;
    Vector3d ad = d - a;
    Vector3d ap = p - a;

    Matrix3d M;
    M << ab, ac, ad;

    Vector3d solve = M.inverse()*ap;
    Vector4d ans;
    ans << 1-solve.x()- solve.y()- solve.z(), solve.x(), solve.y(), solve.z();
    return ans;
}

bool Tetrahedron::IsInside(Vector3d p){
    Vector4d bary = Barycentric(p);
    return bary.x() >= 0 && bary.y() >= 0 && bary.z() >= 0 && bary.w() >= 0;
}