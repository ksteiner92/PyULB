//
// Created by klaus on 02.03.19.
//

#include <vector>
#include <array>

#include "poisson.h"
#include "eigen.h"

using namespace std;
using namespace Eigen;

template<uint Dim, uint TopDim>
Poisson<Dim, TopDim>::Poisson(IMesh* mesh)
{
   const size_t nvertices = mesh->getNumVertices();
   vector<long long int> interior(nvertices);
   IMesh* hull = mesh->getHull();
   size_t refcount = 0;
   for (size_t iv = 0; iv < nvertices; iv++) {
      if (hull->getVertexByID(iv) == nullptr) {
         interior[iv] = refcount;
         refcount++;
      } else
         interior[iv] = -1;
   }
   K.resize(refcount, refcount);
   Mesh<Dim, 0>* mesh0D = dynamic_cast<Mesh<Dim, 0>*>(mesh);
   for (size_t ic = 0; ic < mesh->getNumBodies(); ic++) {
      MeshElement* body = mesh->getBody(ic);
      array<Matrix<double, Dim, 1>, 3> p;
      for (uint8_t i = 0; i < 3; i++)
         p[i] = mesh0D->getPoint((*body)[i]);
      array<double, 3> b;
      for (uint8_t i = 0; i < 3; i++)
         b[i] = p[(i + 1) % 3](1) - p[(i + 2) % 3](1);
      array<double, 3> c;
      for (uint8_t i = 0; i < 3; i++)
         c[i] = p[(i + 2) % 3](0) - p[(i + 1) % 3](0);
      const double Omega = 2 * (p[1](0) * p[2](1) + p[0](0) * p[1](1) + p[0](1) * p[2](0) -
              p[1](0) * p[0](1) - p[2](0) * p[1](1) - p[2](1) * p[0](0));
      for (uint8_t i = 0; i < 3; i++) {
         const size_t gi = interior[(*body)[i]];
         if (gi < 0)
            continue;
         for (uint8_t j = i; j < 3; j++) {
            const size_t gj = interior[(*body)[i]];
            if (gj >= 0)
               K.coeffRef(gi, gj) += (b[j] * b[i] + c[j] * c[i]) / Omega;
         }
      }
   }
}

template class Poisson<1, 1>;
template class Poisson<2, 1>;
template class Poisson<2, 2>;
template class Poisson<3, 1>;
template class Poisson<3, 2>;
template class Poisson<3, 3>;