#include "PhyBaseStruct.h"
#include "iostream"

double bias = 1e-6;

double IsEqual(double a, double b, double _bias){
    return abs(a - b) < _bias;
}

//世界坐标转化为局部坐标
Vector3d Coord3d::ConversionFromWorld(Vector3d p){
    Vector3d ans;
    Vector3d k = i.cross(j);
    Matrix3d M;
    M << i,j,k;
    M= M.transpose();
    ans = M*p;
    return ans;
}

//局部坐标转化为世界坐标
Vector3d Coord3d::ConversionToWorld(Vector3d p){
    Vector3d ans;
    Vector3d k = i.cross(j);
    Matrix3d M;
    M << i,j,k;
    ans = M*p;
    return ans;

}

//其他坐标系的坐标转化为本坐标系
Vector3d Coord3d::ConversionFromOtherCoord(Vector3d p, Coord3d otherCoord){
    Vector3d p1 = otherCoord.ConversionToWorld(p);
    Vector3d ans = ConversionFromWorld(p1);
    return ans;
}

double Segment::GetClosestDist(Vector3d point,Vector3d& closestPoint){
    Vector3d ab = end - start;
    Vector3d e_ab = ab.normalized();
    Vector3d ap = point - start;
    double t = fmax(0, e_ab.dot(ap));
    t = fmin(t, ab.norm());
    closestPoint = start + t * e_ab;
    Vector3d ans = ap - t * e_ab;
    return ans.norm();
}

//第一个是自己身上的最近点，第二个是seg上的最近点
double Segment::GetClosestDistWithSegment(Segment seg, Vector3d& closestPoint1, Vector3d& closestPoint2){
    //可以用纯代数法解决 --- 通过求导数计算
    Vector3d ca = start - seg.start;
    Vector3d d1 = end - start;
    Vector3d d2 = seg.end - seg.start;
    double A = d1.dot(d1);
    double B = d2.dot(d2);
    double C = d1.dot(d2);

    double E = -d1.dot(ca);
    double F = -d2.dot(ca);

    //二元一次方程组如下
    // A*t1 - C*t2 = E
    // C*t1 - B*t2 = F
    //解出t1和t2
    double delt = -A*B + C*C;
    double t1 = (-B*E + C*F)/delt;
    double t2 = (-C*E + A*F)/delt;

    t1 = fmax(0, t1);
    t1 = fmin(t1, 1);
    t2 = fmax(0, t2);
    t2 = fmin(t2, 1);

    closestPoint1 = start + t1 * d1;
    closestPoint2 = seg.start + t2 * d2;

    Vector3d ans = closestPoint1 - closestPoint2;
    return ans.norm();
}

Vector3d Segment::GetRandomNormal(){
    Vector3d ab = end - start;
    Vector3d n = ab.cross(Vector3d(0, 0, 1));
    if(n.norm() < bias){
        n = ab.cross(Vector3d(0, 1, 0));
    }
    return n.normalized();
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
        Vector3d closestPoint;
        double d1 = Segment(a, b).GetClosestDist(point,closestPoint);
        double d2 = Segment(b, c).GetClosestDist(point,closestPoint);
        double d3 = Segment(c, a).GetClosestDist(point,closestPoint);
        return fmin(fmin(d1, d2), d3);
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
    double min_dist = trList[0].GetClosestDistWithotRange(p);
    Triangle ans = trList[0];
    for(auto tr:trList){
        double dist = tr.GetClosestDistWithotRange(p);
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