#include "include/test.h"
#include "PhyAlgrithom.h"
#include "PhySphere.h"
#include <iostream>

int main(){
    testPrint();
    PhySphere sphere1(Eigen::Vector3d(0, 0, 0), 5);
    PhySphere sphere2(Eigen::Vector3d(5, 0, 0), 5);
    IntesectData data = PhyA::MKFSJCollision(&sphere1, &sphere2);
    std::cout << "isIntersecting: " << data.isIntersecting<<" penetration: "<< data.penetration << std::endl;
    return 0;
}