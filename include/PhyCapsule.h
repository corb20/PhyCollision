#pragma once
#include "PhyObject.h"
#include "PhyBaseStruct.h"

class PhyCapsule : public PhyObject
{
public:
    Segment segment;
    double radius;

    PhyCapsule(){
        //默认生成一个竖直向上的胶囊
        Segment seg(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 1));
        radius = 1;
    }
    PhyCapsule(Segment& _segment, double _radius) : segment(_segment), radius(_radius) {}

    Eigen::Vector3d SupportPosition(const Eigen::Vector3d& direction) override;
    IntesectData CollisionWith(PhyObject* other) override;
};