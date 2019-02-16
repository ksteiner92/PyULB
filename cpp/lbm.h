//
// Created by klaus on 20.01.19.
//

#ifndef LBM_LBM_H
#define LBM_LBM_H

#include <array>
#include <unordered_set>

#include "mesh.h"
#include "utils.h"
#include "eigen.h"

class LatticeBolzmann2D
{
public:
   LatticeBolzmann2D(Mesh<2, 2>* mesh);

   void calculate();

private:
   Mesh<2, 2>* mesh;
   LineHash linhash;
   Attribute<std::vector<Eigen::Vector2d>>* E_ij;
   Attribute<std::vector<Eigen::Matrix<double, 2, 1>>>* C_ij;
   Attribute<std::vector<Eigen::Vector2d>>* AttrN_K;
   Attribute<std::vector<int>>* AttrP_K;
   Attribute<Eigen::Vector2d>* C_i;
   Attribute<Eigen::Vector2d>* E_i;
   Attribute<std::vector<std::array<double, 9>>>* f_ti;
   std::vector<Eigen::Matrix<double, 9, Eigen::Dynamic>> S;
   std::vector<std::vector<double>> C;
   std::vector<std::vector<int>> P_K;
   std::array<Eigen::Vector2d, 9> c_i;
   double dt;
   double tau;

   static double calcPivot(const Eigen::Vector2d &P,
                           const Eigen::Vector2d &P1,
                           const Eigen::Vector2d &P2,
                           const Eigen::Vector2d &C,
                           Eigen::Vector2d &N1,
                           Eigen::Vector2d &N2);

   inline static double triangleArea(const Eigen::Vector2d &A,
                                     const Eigen::Vector2d &B,
                                     const Eigen::Vector2d &C);

   template<class T, uint SimplexDim>
   static T interpolateVertexAttributeOnSimplex(Simplex<2, SimplexDim>* simplex, Attribute<T>* attr);
};


#endif //LBM_LBM_H
