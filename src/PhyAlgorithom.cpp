#include "PhyAlgrithom.h"
#include <vector>
#include <queue>
#include <unordered_set>

using namespace Eigen;

#define IsSupportZero(s) if((s).isApprox(Vector3d(0,0,0),bias)){return true;} 

namespace PhyA
{
    const int MAX_ITERATE_GJK = 100;
    const int MAX_ITERATE_EPA = 100;

    struct TriangleEPA
    {
        Vector3i ptIndex;
        //始终指向原点的法向量
        Vector3d normal;
        double minDist;

        TriangleEPA(Vector3i _ptIndex, std::vector<Vector3d> &ptList){
            ptIndex = _ptIndex;
            Triangle t(ptList[ptIndex[0]], ptList[ptIndex[1]], ptList[ptIndex[2]]);
            normal = t.normal();
            normalOutOrg(ptList);
            minDist = t.GetClosestDistWithotRange(Vector3d(0, 0, 0));
            //如果出现原点在三角形上面的情况，需要利用其他点来确定法向量
            if(IsEqual(minDist, 0)){
                normalConvex(ptList);
            }
        }

        //从圆心指向外的法向量
        void normalOutOrg(std::vector<Vector3d> &ptList){
            Vector3d a = ptList[ptIndex[0]];
            if(a.dot(normal) < 0){
                normal = -normal;
            }
        }

        //输入另一个点，获取凸多边形的法向量
        void normalConvex(std::vector<Vector3d> &ptList){
            Vector3d otherPoint;
            for(int i = 0;i<ptList.size();i++){
                if(i != ptIndex[0] && i != ptIndex[1] && i != ptIndex[2]){
                    otherPoint = ptList[i];
                    break;
                }
            }
            Vector3d ap = otherPoint - ptList[ptIndex[0]];
            if(normal.dot(ap) > 0){
                normal = -normal;
            }
        }

        bool operator<(const TriangleEPA& t) const
        {
            return minDist < t.minDist;
        }
        bool operator>(const TriangleEPA& t) const
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
        std::priority_queue<TriangleEPA, std::vector<TriangleEPA>, std::greater<TriangleEPA> > pq;
        pq.push(TriangleEPA(Vector3i(0, 1, 2), points));
        pq.push(TriangleEPA(Vector3i(0, 1, 3), points));
        pq.push(TriangleEPA(Vector3i(0, 2, 3), points));
        pq.push(TriangleEPA(Vector3i(1, 2, 3), points));

        iterate = 0;
        bool isOver = false;
        IntesectData ans = IntesectData();
        while (!isOver && iterate <= MAX_ITERATE_EPA)
        {
            TriangleEPA t = pq.top();
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

            pq.push(TriangleEPA(Vector3i(t.ptIndex[0], t.ptIndex[1], points.size()), points));
            pq.push(TriangleEPA(Vector3i(t.ptIndex[0], t.ptIndex[2], points.size()), points));
            pq.push(TriangleEPA(Vector3i(t.ptIndex[1], t.ptIndex[2], points.size()), points));
            points.push_back(s);
        }
        return ans;
    }

    Tetrahedron GetRandomTetrahedron(Vector3d p1,Vector3d p2,PhyObject* obj1, PhyObject* obj2)
    {
        //前提是p1和p2不重合
        if(p1.isApprox(p2)){
            return Tetrahedron();
        }
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

    bool GJKCollision(PhyObject *obj1, PhyObject *obj2, Tetrahedron& simplex, double bias)
    {
        simplex = Tetrahedron();
        //GJK算法第一次迭代
        Vector3d s0 = obj1->SupportPosition(Vector3d(1, 0, 0)) - obj2->SupportPosition(Vector3d(-1, 0, 0));
        IsSupportZero(s0);
        Vector3d d0 = -s0;
        Vector3d s1 = obj1->SupportPosition(d0) - obj2->SupportPosition(-d0);
        IsSupportZero(s1);
        Vector3d d1 = GetIterateDir(s0, s1);
        if(d1.isApprox(Vector3d(0,0,0)))
        {
            //迭代方向为0,说明点在该直线上，有碰撞
            simplex = GetRandomTetrahedron(s0,s1,obj1,obj2);
            return true;
        }
        Vector3d s2 = obj1->SupportPosition(d1) - obj2->SupportPosition(-d1);
        IsSupportZero(s2);
        Vector3d d2 = GetIterateDir(s0, s1, s2);
        Vector3d s3 = obj1->SupportPosition(d2) - obj2->SupportPosition(-d2);
        IsSupportZero(s3);
        //GJK算法开始跌代
        simplex = Tetrahedron(s0, s1, s2, s3);
        //说明存在共面的情况
        if(simplex.IsValid() == false)
        {
            return true;
        }
        int iterate = 0;
        while (simplex.IsInside(Vector3d(0,0,0)) == false)
        {
            //取得距离原点最近的平面的3个点以及朝向原点方向的法向量
            Triangle t = simplex.GetClosestTriangle(Vector3d(0, 0, 0));
            Vector3d di = t.normal2Org();
            Vector3d si = obj1->SupportPosition(di) - obj2->SupportPosition(-di); 
            if(si.isApprox(Vector3d(0,0,0)))
            {
                //如果支撑点位0，那么说明两者恰好接触
                simplex = Tetrahedron();
                return true;
            }
            if (si.dot(t.normal()) < 0)
            {
                return false;
            }
            if(iterate > MAX_ITERATE_GJK)
            {
                //迭代次数过多,认为没有碰撞
                return false;
            }
            simplex = Tetrahedron(si, t.a, t.b, t.c);
            iterate++;
        }
        return true;
    }

    //EPA算法,只有在GJK算法返回true的情况下才会调用
    IntesectData EPACollision(Tetrahedron& simplex, PhyObject* obj1, PhyObject* obj2, double bias)
    {
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
        std::vector<TriangleEPA> TriList;
        TriList.push_back(TriangleEPA(Vector3i(0, 1, 2), points));
        TriList.push_back(TriangleEPA(Vector3i(0, 1, 3), points));
        TriList.push_back(TriangleEPA(Vector3i(0, 2, 3), points));
        TriList.push_back(TriangleEPA(Vector3i(1, 2, 3), points));

        int iterate = 0;
        IntesectData ans = IntesectData();
        int minTriIndex = GetMinTriIndex(TriList);

        while (iterate <= MAX_ITERATE_EPA)
        {
            //根据最近三角形的法向求新的支撑点
            TriangleEPA t = TriList[minTriIndex];
            Vector3d n = t.normal;
            Vector3d s = obj1->SupportPosition(n) - obj2->SupportPosition(-n);
            double newMinDist = s.dot(n);
            //缩小容差
            if(IsEqual(newMinDist, t.minDist,1e-3))
            {
                //说明已经找到了最近的三角形
                break;
            }

            points.push_back(s);
            //更新三角形列表，并求一个新的最小三角形，这里要删除掉所有的沿着n方向的三角形
            int pct = points.size();
            std::unordered_set<int> UniqueEdge;
            std::vector<TriangleEPA> delTriList;
            for(int i = 0;i<TriList.size();i++){
                TriangleEPA tri = TriList[i];
                //常用的一边检测一边删除的O(1)的算法，总共需要O(n)的时间复杂度
                //通过判断点是否在三角面的外部来判断是否要删除该三角面
                Vector3d a = points[tri.ptIndex[0]];
                Vector3d as = s - a;
                if(tri.normal.dot(as) > 0){
                    //删除掉所有沿着n方向的三角形
                    delTriList.push_back(tri);
                    //并且删掉原来的三角形
                    TriList[i] = TriList.back();
                    TriList.pop_back();
                    i--;
                }
            }
            //更新边缘列表
            GetUniqueEdges(delTriList, UniqueEdge, pct);
            //重构三角形列表
            for(auto edge:UniqueEdge){
                int a = edge/pct;
                int b = edge%pct;
                TriList.push_back(TriangleEPA(Vector3i(a,b,pct-1),points));
            }
            minTriIndex = GetMinTriIndex(TriList);
        }
        return ans;
    }

    int GetMinTriIndex(std::vector<TriangleEPA> &TriList){
        int minIndex = -1;
        double minTriDist = MAXFLOAT;
        for(int i = 0;i<TriList.size();i++){
            if(TriList[i].minDist < minTriDist){
                minTriDist = TriList[i].minDist;
                minIndex = i;
            }
        }
        return minIndex;
    }

    void GetUniqueEdges(std::vector<TriangleEPA> &TriList, std::unordered_set<int> &UniqueEdge, int pct){
        for(auto tri:TriList){
            for(int i = 0;i<3;i++){
                int a = tri.ptIndex[i];
                int b = tri.ptIndex[(i+1)%3];
                //一定是小的在前面
                if(a>b)
                    std::swap(a,b);
                int edge = a*pct+b;
                if(UniqueEdge.find(edge) == UniqueEdge.end()){
                    UniqueEdge.insert(edge);
                }
                else{
                    UniqueEdge.erase(edge);
                }
            }
        }
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