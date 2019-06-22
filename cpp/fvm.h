//
// Created by klaus on 17.04.19.
//

#ifndef ULB_FVM_H
#define ULB_FVM_H

#include "mesh.h"

class FVM
{
public:
   FVM(IMesh* mesh);

   Mesh<2, 2>* getDualMesh() const;

private:
   std::unique_ptr<Mesh<2, 2>> voronoi;

};


#endif //ULB_FVM_H
