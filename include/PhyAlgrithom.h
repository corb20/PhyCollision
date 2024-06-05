#pragma once
#include "PhyObject.h"
#include "PhyBaseStruct.h"

namespace PhyA
{
    IntesectData MKFSJCollision(PhyObject* obj1, PhyObject* obj2, double bias = 1e-6);

    /// @brief gjk算法
    /// @param obj1 in:需要计算求交的凸物体1
    /// @param obj2 in:需要计算求交的凸物体2
    /// @param simplex [out] 算法结束后的最简单四面体
    /// @param bias in:误差
    /// @return 是否相交
    bool GJKCollision(PhyObject *obj1, PhyObject *obj2, Tetrahedron& simplex, double bias= 1e-6);

    IntesectData EPACollision(Tetrahedron &simplex, PhyObject *obj1, PhyObject *obj2, double bias);

    Vector3d GetIterateDir(Vector3d s0, Vector3d s1, Vector3d s2);

    Vector3d GetIterateDir(Triangle& t);

    Vector3d GetIterateDir(Vector3d s0, Vector3d s1);
}