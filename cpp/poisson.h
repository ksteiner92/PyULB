//
// Created by klaus on 02.03.19.
//

#ifndef ULB_POISSON_H
#define ULB_POISSON_H

#include <vector>

#include "eigen.h"
#include "mesh.h"

template<uint Dim, uint TopDim>
class Poisson
{
   static_assert(Dim >= TopDim && TopDim >= 1, "Dimension combination not supported");

public:
   Poisson(IMesh* mesh);

   void solve(const std::string &phi,
              const std::string &rho,
              const std::vector<std::size_t> &dirichlet);

   void grad(const std::string& u, const std::string& E) const;

private:
   IMesh* mesh;
   Eigen::SparseMatrix<double, Eigen::RowMajor> K;
   std::vector<long long int> interior;

};


#endif //ULB_POISSON_H
