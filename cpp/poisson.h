//
// Created by klaus on 02.03.19.
//

#ifndef ULB_POISSON_H
#define ULB_POISSON_H

#include "eigen.h"
#include "mesh.h"

template<uint Dim, uint TopDim>
class Poisson
{
   static_assert(Dim >= TopDim && TopDim >= 1, "Dimension combination not supported");

public:
   Poisson(IMesh* mesh);

private:
   Eigen::SparseMatrix<double> K;

};


#endif //ULB_POISSON_H
