#include "PhyAlgrithom.h"
#include <vector>
#include <queue>

using namespace Eigen;
namespace PhyA
{
    const int MAX_ITERATE_GJK = 100;
    const int MAX_ITERATE_EPA = 100;

    struct TriangleMKFSJ
    {
        Vector3i ptIndex;
        //始终指向原点的法向量
        Vector3d normal;
        double minDist;

        TriangleMKFSJ(Vector3i _ptIndex, std::vector<Vector3d> &ptList){
            ptIndex = _ptIndex;
            Triangle t(ptList[ptIndex[0]], ptList[ptIndex[1]], ptList[ptIndex[2]]);
            normal = t.normal2Org();
            minDist = t.GetClosestDist(Vector3d(0, 0, 0));
        }

        bool operator<(const TriangleMKFSJ& t) const
        {
            return minDist < t.minDist;
        }
        bool operator>(const TriangleMKFSJ& t) const
        {
            return minDist > t.minDist;
        }
    };

    IntesectData MKFSJCollision(PhyObject* obj1, PhyObject* obj2, double bias)
    {
        //GJK算法第一次迭代
        Vector3d s0 = obj1->SupportPosition(Vector3d(1, 0, 0)) - obj2->SupportPosition(Vector3d(-1, 0, 0));
        Vector3d d0 = -s0;
        Vector3d s1 = obj1->SupportPosition(d0) - obj2->SupportPosition(-d0);
        Vector3d d1 = GetIterateDir(s0, s1);
        Vector3d s2 = obj1->SupportPosition(d1) - obj2->SupportPosition(-d1);
        Vector3d d2 = GetIterateDir(s0, s1, s2);
        Vector3d s3 = obj1->SupportPosition(d2) - obj2->SupportPosition(-d2);
        //GJK算法开始跌代
        Tetrahedron simplex(s0, s1, s2, s3);
        int iterate = 0;
        while (simplex.IsInside(Vector3d(0,0,0)) == false)
        {
            //取得距离原点最近的平面的3个点以及朝向原点方向的法向量
            Triangle t = simplex.GetClosestTriangle(Vector3d(0, 0, 0));
            Vector3d di = t.normal2Org();
            Vector3d si = obj1->SupportPosition(di) - obj2->SupportPosition(-di); 
            if (si.dot(t.normal()) < 0)
            {
                return {false, Vector3d(0, 0, 0), 0};
            }
            if(iterate > MAX_ITERATE_GJK)
            {
                //迭代次数过多,认为没有碰撞
                return {false, Vector3d(0, 0, 0), 0};
            }
            simplex = Tetrahedron(si, t.a, t.b, t.c);
            iterate++;
        }
        //迭代结束

        // 随后使用EPA算法求具体的碰撞信息-- 但跌代结束的条件是什么呢？
        std::vector<Vector3d> points;
        points.push_back(simplex.a);
        points.push_back(simplex.b);
        points.push_back(simplex.c);
        points.push_back(simplex.d);
        //使用小顶堆存储三角形
        std::priority_queue<TriangleMKFSJ, std::vector<TriangleMKFSJ>, std::greater<TriangleMKFSJ> > pq;
        pq.push(TriangleMKFSJ(Vector3i(0, 1, 2), points));
        pq.push(TriangleMKFSJ(Vector3i(0, 1, 3), points));
        pq.push(TriangleMKFSJ(Vector3i(0, 2, 3), points));
        pq.push(TriangleMKFSJ(Vector3i(1, 2, 3), points));

        iterate = 0;
        bool isOver = false;
        IntesectData ans = IntesectData();
        while (!isOver && iterate <= MAX_ITERATE_EPA)
        {
            TriangleMKFSJ t = pq.top();
            pq.pop();
            Vector3d n = t.normal;
            Vector3d s = obj1->SupportPosition(n) - obj2->SupportPosition(-n);
            double dist = s.norm();
            
            for (int i = 0; i < points.size(); i++)
            {
                if (points[i].isApprox(s, bias))
                {
                    isOver = true;
                    break;
                }
            }
            if(ans.isIntersecting == false || dist < ans.penetration)
            {
                ans.isIntersecting = true;
                ans.normal = s.normalized();
                ans.penetration = dist;
            }

            pq.push(TriangleMKFSJ(Vector3i(t.ptIndex[0], t.ptIndex[1], points.size()), points));
            pq.push(TriangleMKFSJ(Vector3i(t.ptIndex[0], t.ptIndex[2], points.size()), points));
            pq.push(TriangleMKFSJ(Vector3i(t.ptIndex[1], t.ptIndex[2], points.size()), points));
            points.push_back(s);
        }
        return ans;
    }

    Tetrahedron GetRandomTetrahedron(Vector3d p1,Vector3d p2,PhyObject* obj1, PhyObject* obj2)
    {
        Triangle t;
        for(int i = 0;i<3;i++){
            Vector3d dir = Vector3d(0,0,0);
            dir[i] = 1;
            Vector3d s1 = obj1->SupportPosition(dir) - obj2->SupportPosition(-dir);
            t = Triangle(p1,p2,s1);
            if(t.IsValid()){
                break;
            }
        }
        //t在此时一定是合法的
        Vector3d n = t.normal2Org();
        Vector3d s = obj1->SupportPosition(n) - obj2->SupportPosition(-n);
        //但最后的四面体不一定合法
        return Tetrahedron(s,p1,p2,t.a);
    }

    Tetrahedron GJKCollision(PhyObject* obj1, PhyObject* obj2)
    {
        //GJK算法第一次迭代
        Vector3d s0 = obj1->SupportPosition(Vector3d(1, 0, 0)) - obj2->SupportPosition(Vector3d(-1, 0, 0));
        Vector3d d0 = -s0;
        Vector3d s1 = obj1->SupportPosition(d0) - obj2->SupportPosition(-d0);
        Vector3d d1 = GetIterateDir(s0, s1);
        if(d1.isApprox(Vector3d(0,0,0)))
        {
            //迭代方向为0,说明点在该直线上，有碰撞
            return GetRandomTetrahedron(s0,s1,obj1,obj2);
        }
        Vector3d s2 = obj1->SupportPosition(d1) - obj2->SupportPosition(-d1);
        Vector3d d2 = GetIterateDir(s0, s1, s2);
        Vector3d s3 = obj1->SupportPosition(d2) - obj2->SupportPosition(-d2);
        //GJK算法开始跌代
        Tetrahedron simplex(s0, s1, s2, s3);
        if(simplex.IsValid() == false)
        {
            return simplex;
        }
        int iterate = 0;
        while (simplex.IsInside(Vector3d(0,0,0)) == false)
        {
            //取得距离原点最近的平面的3个点以及朝向原点方向的法向量
            Triangle t = simplex.GetClosestTriangle(Vector3d(0, 0, 0));
            Vector3d di = t.normal2Org();
            Vector3d si = obj1->SupportPosition(di) - obj2->SupportPosition(-di); 
            if (si.dot(t.normal()) < 0)
            {
                return Tetrahedron();
            }
            if(iterate > MAX_ITERATE_GJK)
            {
                //迭代次数过多,认为没有碰撞
                return Tetrahedron();
            }
            simplex = Tetrahedron(si, t.a, t.b, t.c);
            iterate++;
        }
        return simplex;
    }

    IntesectData EPACollision(Tetrahedron& simplex, PhyObject* obj1, PhyObject* obj2, double bias)
    {
        if(simplex.IsAllZero())
        {
            return {false, Vector3d(0, 0, 0), 0};
        }

        if(simplex.IsValid() == false)
        {
            //两者恰好接触
            return {true, Vector3d(0, 0, 0), 0};
        }
        // 随后使用EPA算法求具体的碰撞信息-- 但跌代结束的条件是什么呢？
        std::vector<Vector3d> points;
        points.push_back(simplex.a);
        points.push_back(simplex.b);
        points.push_back(simplex.c);
        points.push_back(simplex.d);
        //使用小顶堆存储三角形
        std::priority_queue<TriangleMKFSJ, std::vector<TriangleMKFSJ>, std::greater<TriangleMKFSJ> > pq;
        pq.push(TriangleMKFSJ(Vector3i(0, 1, 2), points));
        pq.push(TriangleMKFSJ(Vector3i(0, 1, 3), points));
        pq.push(TriangleMKFSJ(Vector3i(0, 2, 3), points));
        pq.push(TriangleMKFSJ(Vector3i(1, 2, 3), points));

        int iterate = 0;
        bool isOver = false;
        IntesectData ans = IntesectData();
        while (!isOver && iterate <= MAX_ITERATE_EPA)
        {
            TriangleMKFSJ t = pq.top();
            pq.pop();
            Vector3d n = t.normal;
            Vector3d s = obj1->SupportPosition(n) - obj2->SupportPosition(-n);
            double dist = s.norm();
            
            for (int i = 0; i < points.size(); i++)
            {
                if (points[i].isApprox(s, bias))
                {
                    isOver = true;
                    break;
                }
            }
            if(ans.isIntersecting == false || dist < ans.penetration)
            {
                ans.isIntersecting = true;
                ans.normal = s.normalized();
                ans.penetration = dist;
            }

            pq.push(TriangleMKFSJ(Vector3i(t.ptIndex[0], t.ptIndex[1], points.size()), points));
            pq.push(TriangleMKFSJ(Vector3i(t.ptIndex[0], t.ptIndex[2], points.size()), points));
            pq.push(TriangleMKFSJ(Vector3i(t.ptIndex[1], t.ptIndex[2], points.size()), points));
            points.push_back(s);
        }
        return ans;
    }

    Vector3d GetIterateDir(Vector3d s0, Vector3d s1, Vector3d s2)
    {
        Vector3d v01 = s1 - s0;
        Vector3d v02 = s2 - s0;

        Vector3d n = v01.cross(v02).normalized();
        if (n.dot(s0) > 0)
        {
            return -n;
        }
        else
        {
            return n;
        }
    }

    Vector3d GetIterateDir(Triangle& t)
    {
        Vector3d n = t.normal();
        if (n.dot(t.a) > 0)
        {
            return -n;
        }
        else
        {
            return n;
        }
    }

    Vector3d GetIterateDir(Vector3d s0, Vector3d s1)
    {
        Vector3d v01 = s1 - s0;
        Vector3d e_v01 = v01.normalized();

        return -s0 - (-s0).dot(e_v01) * e_v01;
    }
}