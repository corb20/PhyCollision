#include "PhyConvex.h"
#include "PhyAlgrithom.h"

Eigen::Vector3d PhyConvex::SupportPosition(const Eigen::Vector3d& direction){
    double max = -1e9;
    Eigen::Vector3d result;
    for(auto &v: vertices){
        double dot = v.dot(direction);
        if(dot > max){
            max = dot;
            result = v;
        }
    }
    return result;
}

IntesectData PhyConvex::CollisionWith(PhyObject* other){
    IntesectData data = PhyA::MKFSJCollision(this,other);
    return data;
}