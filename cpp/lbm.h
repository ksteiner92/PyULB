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
#include "logger.h"

class LatticeBoltzmann2D
{
public:
   LatticeBoltzmann2D(Mesh<2, 2>* mesh, const std::vector<Mesh<2, 1>*>& boundaries, double tau, double c_s);

   void calculateNextTime(double dt);

   void update();

private:
   static constexpr double epsilon = 1.0e-10;
   Mesh<2, 2>* mesh;
   ListAttribute<std::vector<Eigen::Matrix<double, 9, 1>>>* f;
   ListAttribute<std::vector<Eigen::Matrix<double, 9, 1>>>* fe;
   ListAttribute<std::vector<Eigen::Vector2d>>* u;
   ListAttribute<std::vector<double>>* rho;

   std::vector<Eigen::Matrix<double, 9, Eigen::Dynamic>> S;
   std::vector<std::vector<double>> C;
   std::vector<std::vector<int>> P_K;
   Eigen::Matrix<double, 2, 9> c_i;
   Eigen::Matrix<double, 9, 2> betac_i;
   Eigen::Matrix<double, 9, 1> w_i;
   std::unordered_map<std::size_t, uint8_t> infacing_boundary_c;
   double c_s;
   double beta;
   double tau;

   void init();

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
