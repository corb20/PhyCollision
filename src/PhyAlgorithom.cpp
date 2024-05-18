// #include "PhyAlgrithom.h"
// #include "PhyBaseStruct.h"
// #include <vector>

// using namespace Eigen;
// namespace PhyA
// {
//     IntesectData MKFSJCollision(PhyObject* obj1, PhyObject* obj2, double bias)
//     {
//         //GJK算法第一次迭代
//         Vector3d s0 = obj1->SupportPosition(Vector3d(1, 0, 0)) - obj2->SupportPosition(Vector3d(-1, 0, 0));
//         Vector3d d0 = -s0;
//         Vector3d s1 = obj1->SupportPosition(d0) - obj2->SupportPosition(-d0);
//         Vector3d d1 = GetIterateDir(s0, s1);
//         Vector3d s2 = obj1->SupportPosition(d1) - obj2->SupportPosition(-d1);
//         Vector3d d2 = GetIterateDir(s0, s1, s2);
//         Vector3d s3 = obj1->SupportPosition(d2) - obj2->SupportPosition(-d2);
//         //将采访过的点都放入一个列表中
//         std::vector<Vector3d> visited = {s0, s1, s2, s3};

//         Tetrahedron simplex(s0, s1, s2, s3);
//         while (simplex.IsInside(Vector3d(0,0,0)) == false)
//         {
//             //取得距离原点最近的平面的3个点以及朝向原点方向的法向量
            
//         }
        
//         // 得到第一次的迭代结果
//     }

//     Vector3d GetIterateDir(Vector3d s0, Vector3d s1, Vector3d s2)
//     {
//         Vector3d v01 = s1 - s0;
//         Vector3d v02 = s2 - s0;

//         Vector3d n = v01.cross(v02);
//         if (n.dot(s0) > 0)
//         {
//             return n;
//         }
//         else
//         {
//             return -n;
//         }
//     }

//     Vector3d GetIterateDir(Vector3d s0, Vector3d s1)
//     {
//         Vector3d v01 = s1 - s0;
//         Vector3d e_v01 = v01.normalized();
//         return -s0 - s0.dot(e_v01) * e_v01;
//     }
// }