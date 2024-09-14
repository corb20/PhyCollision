#pragma once
#include "PhyObject.h"
#include <vector>

class PhyBox : public PhyObject
{
public:
    //坐标原点为长方体的左后下角
    Coord3d coord;
    Eigen::Vector3d minP;
    Eigen::Vector3d maxP;

    PhyBox() : coord(Coord3d()), minP(Eigen::Vector3d(-1, -1, -1)), maxP(Eigen::Vector3d(1, 1, 1)) {}
    PhyBox(Vector3d _minP, Vector3d _maxP) : coord(Coord3d()), minP(_minP), maxP(_maxP) {}
    PhyBox(Coord3d _coord, Eigen::Vector3d _minP, Eigen::Vector3d _maxP) : coord(_coord), minP(_minP), maxP(_maxP) {}
    
    std::vector<Eigen::Vector3d> GetLocalVertices();
    std::vector<Eigen::Vector3d> GetWorldVertices();

    double GetClosestDist(const Eigen::Vector3d& point, Eigen::Vector3d& closestPoint);
    Eigen::Vector3d GetRandomNormal(Eigen::Vector3d p);

    Eigen::Vector3d SupportPosition(const Eigen::Vector3d& direction) override;
    IntesectData CollisionWith(PhyObject* other) override;
};