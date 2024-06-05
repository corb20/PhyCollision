
#include "PhyCapsule.h"
#include "PhySphere.h"
#include "PhyBox.h"
#include "PhyConvex.h"
#include "PhyAlgrithom.h"
    
Eigen::Vector3d PhyCapsule::SupportPosition(const Eigen::Vector3d& direction){
    double p1DotDir = segment.start.dot(direction);
    double p2DotDir = segment.end.dot(direction);

    Vector3d selectP = p1DotDir > p2DotDir ? segment.start : segment.end;

    return selectP + radius * direction.normalized();
}
IntesectData PhyCapsule::CollisionWith(PhyObject* other){
    IntesectData data;
    data.isIntersecting = false;
    data.normal = Eigen::Vector3d(0, 0, 0);
    data.penetration = 0;

    PhySphere* sphere = dynamic_cast<PhySphere*>(other);
    if (sphere)
    {
        data = sphere->CollisionWith(this);
        data.normal = -data.normal;
        return data;
    }

    // 检测是否是Capsule
    PhyCapsule* capsule = dynamic_cast<PhyCapsule*>(other);
    if (capsule)
    {
        //求两根线段的最近距离
    }

}