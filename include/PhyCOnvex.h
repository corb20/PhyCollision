#pragma once
#include "PhyObject.h"

class PhyConvex : public PhyObject
{
public:
    //只需要存顶点的位置——凸包
    std::vector<Eigen::Vector3d> vertices;

    PhyConvex(std::vector<Eigen::Vector3d> &_vertices) : vertices(_vertices) {}

    Eigen::Vector3d SupportPosition(const Eigen::Vector3d& direction) override;
    IntesectData CollisionWith(PhyObject* other) override;

    //可以用EPA算法去求一个点集的唯一凸包
};