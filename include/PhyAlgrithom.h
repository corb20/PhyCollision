#include "PhyObject.h"
#include "PhyBaseStruct.h"

namespace PhyA
{
    IntesectData MKFSJCollision(PhyObject* obj1, PhyObject* obj2, double bias = 1e-6);
    Vector3d GetIterateDir(Vector3d s0, Vector3d s1, Vector3d s2);
    Vector3d GetIterateDir(Triangle& t);
    Vector3d GetIterateDir(Vector3d s0, Vector3d s1);
}