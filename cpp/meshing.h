//
// Created by klaus on 08.03.19.
//

#ifndef ULB_MESHING_H
#define ULB_MESHING_H

#include "mesh.h"

class Meshing
{
public:
   Meshing();

   template<uint Dim, uint TopDim>
   void generate(Mesh<Dim, TopDim>& mesh);

};

#endif //ULB_MESHING_H
