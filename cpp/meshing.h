//
// Created by klaus on 08.03.19.
//

#ifndef ULB_MESHING_H
#define ULB_MESHING_H

#include "mesh.h"

template<uint Dim, uint TopDim>
class Meshing
{
public:
   Meshing();

   void generate(Mesh<Dim, TopDim>& mesh);

};

#endif //ULB_MESHING_H
