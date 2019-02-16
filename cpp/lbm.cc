//
// Created by klaus on 20.01.19.
//

#include "eigen.h"
#include "lbm.h"

using namespace std;
using namespace Eigen;

double LatticeBolzmann2D::triangleArea(const Vector2d& A, const Vector2d& B, const Vector2d& C)
{
   return abs(A[0] * (B[1] - C[1]) + B[0] * (C[1] - A[1]) + C[0] * (A[1] - B[1])) * 0.5;
}

template<class T, uint SimplexDim>
T LatticeBolzmann2D::interpolateVertexAttributeOnSimplex(Simplex<2, SimplexDim>* simplex,
        Attribute<T>* attr)
{
   T val = 0;
   for (int i = 0; i < (SimplexDim + 1); i++)
      val += attr->getValue((*simplex)[i]);
   return val / ((T) (SimplexDim + 1));
}

void LatticeBolzmann2D::calculate()
{
   const double dttau = dt / tau;
   for (size_t iP = 0; iP < P_K.size(); iP++) {
      array<double, 9>& f_P = f_ti->getValue(iP).back();
      const auto& S_ik = S[iP];
      const vector<double>& C_k = C[iP];
      for (int iK = 0; iK < C_k.size(); iK++) {
         array<double, 9>& f_PK = f_ti->getValue(P_K[iP][iK]).back();
         for (int i = 0; i < 9; i++) {
            double f_i = f_P[i];
            f_i += S_ik(i,iK) * f_PK[i];
            /*if (i == iK) {
               f_i -= dttau * C[iK] * ()
            }*/
         }
      }
   }
}

double LatticeBolzmann2D::calcPivot(const Vector2d &P,
                                    const Vector2d &E1,
                                    const Vector2d &E2,
                                    const Vector2d &C,
                                    Vector2d &N1,
                                    Vector2d &N2)
{
   const Vector2d CE1 = E1 - C;
   const Vector2d CE2 = E2 - C;
   N1 << CE1(1), CE1(0);
   N2 << CE2(1), CE2(0);
   if (((E1(0) - P(0)) * (C(1) - P(1)) - (E1(1) - P(1)) * (C(0) - P(0))) > 0) {
      N1(0) = -N1(0);
      N2(1) = -N2(1);
   } else {
      N1(1) = -N1(1);
      N2(0) = -N2(0);
   }
   return triangleArea(P, E1, C) + triangleArea(P, E2, C);
}

LatticeBolzmann2D::LatticeBolzmann2D(Mesh<2, 2>* mesh)
   : mesh(mesh)
{
   E_ij = mesh->getOrCreateAttributeOnVertex<vector<Vector2d>>("E_ij");
   C_ij = mesh->getOrCreateAttributeOnVertex<vector<Matrix<double, 2, 1>>>("C_ij");
   C_i = mesh->getOrCreateAttributeOnFace<Vector2d>("C_i");
   E_i = mesh->getOrCreateAttributeOnEdge<Vector2d>("E_i");
   AttrN_K = mesh->getOrCreateAttributeOnVertex<vector<Vector2d>>("N_K");
   AttrP_K = mesh->getOrCreateAttributeOnVertex<vector<int>>("P_K");
   f_ti = mesh->getOrCreateAttributeOnVertex<vector<array<double, 9>>>("f_ti");

   P_K.resize(mesh->getNumVertices());
   vector<vector<Vector2d>> N_K(mesh->getNumVertices());
   vector<vector<double>> V_K(mesh->getNumVertices());
   vector<double> V_P(mesh->getNumVertices(), 0.0);
   array<Vector2d, 3> E;
   array<Vector2d, 3> P;
   array<int, 3> Pidx;
   for (size_t i = 0; i < mesh->getNumFaces(); i++) {
      // Calculate face center point
      const Face<2>* face = mesh->getFace(i);
      Pidx[0] = (*face)[0];
      P[0] = mesh->getPoint(Pidx[0]);
      Pidx[1] = (*face)[1];
      P[1] = mesh->getPoint(Pidx[1]);
      Pidx[2] = (*face)[2];
      P[2] = mesh->getPoint(Pidx[2]);
      E[0] = P[0] + (P[1] - P[0]) * 0.5;
      const Vector2d C = P[2] + (E[0] - P[2]) * 2.0 / 3.0;
      C_i->setValue(i, C);
      C_ij->getValue(Pidx[0]).push_back(C);
      C_ij->getValue(Pidx[1]).push_back(C);
      C_ij->getValue(Pidx[2]).push_back(C);

      // Calculate edge center points
      E_i->setValue(mesh->getEdge(Pidx[0], Pidx[1])->getID(), E[0]);
      E_ij->getValue(Pidx[0]).push_back(E[0]);
      E_ij->getValue(Pidx[1]).push_back(E[0]);
      E[1] = P[0] + (P[2] - P[0]) * 0.5;
      E_i->setValue(mesh->getEdge(Pidx[0], Pidx[2])->getID(), E[1]);
      E_ij->getValue(Pidx[0]).push_back(E[1]);
      E_ij->getValue(Pidx[2]).push_back(E[1]);
      E[2] = P[1] + (P[2] - P[1]) * 0.5;
      E_i->setValue(mesh->getEdge(Pidx[1], Pidx[2])->getID(), E[2]);
      E_ij->getValue(Pidx[1]).push_back(E[2]);
      E_ij->getValue(Pidx[2]).push_back(E[2]);

      const array<int, 3> iE1 = {0, 0, 2};
      const array<int, 3> iE2 = {1, 2, 1};
      const array<int, 6> iP_K = {1, 2, 0, 2, 0, 1};
      Vector2d N1, N2;
      for (int iP = 0; iP < 3; iP++) {
         const double V = calcPivot(P[iP], E[iE1[iP]], E[iE2[iP]], C, N1, N2);
         V_P[Pidx[iP]] += V;
         Vector2d N = 5.0 / 12.0 * N1 + 2.0 / 12.0 * N2;
         for (int j = 0; j < 2; j++) {
            const auto it = find(P_K[Pidx[iP]].begin(), P_K[Pidx[iP]].end(), Pidx[iP_K[2 * iP + j]]);
            if (it == P_K[Pidx[iP]].end()) {
               P_K[Pidx[iP]].push_back(Pidx[iP_K[2 * iP  + j]]);
               AttrP_K->getValue(Pidx[iP]).push_back(Pidx[iP_K[2 * iP  + j]]);
               N_K[Pidx[iP]].push_back(N);
               AttrN_K->getValue(Pidx[iP]).push_back(N);
               V_K[Pidx[iP]].push_back(V);
            } else {
               const size_t idx = distance(P_K[Pidx[iP]].begin(), it);
               N_K[Pidx[iP]][idx] += N;
               AttrN_K->getValue(Pidx[iP])[idx] = N_K[Pidx[iP]][idx];
               V_K[Pidx[iP]][idx] += V;
            }
            N = 5.0 / 12.0 * N2 + 2.0 / 12.0 * N1;
         }
      }
   }

   S.resize(mesh->getNumVertices());
   C.resize(mesh->getNumVertices());
   const double c = 1.0 / sqrt(2.0);
   c_i[0] << 0.0, 0.0;
   c_i[1] << 1.0, 0.0;
   c_i[2] << 0.0, 1.0;
   c_i[3] << -1.0, 0.0;
   c_i[4] << 0.0, -1.0;
   c_i[5] << c, c;
   c_i[6] << -c, c;
   c_i[7] << -c, -c;
   c_i[8] << c, -c;
   for (size_t ivertex = 0; ivertex < mesh->getNumVertices(); ivertex++) {
      const int n = P_K[ivertex].size() + 1;
      Matrix<double, 9, Dynamic> S_ik(9, n);
      C[ivertex].resize(n);
      C[ivertex][0] = 1.0 / 3.0;
      const double V = V_P[ivertex];
      for (int i = 0; i < 9; i++) {
         S_ik(i, 0) = 0.0;
         for (int k = 1; k < n; k++) {
            S_ik(i, k) =  c_i[i].dot(N_K[ivertex][k - 1]) / V;
            if (i == 0)
               C[ivertex][k] = V_K[ivertex][k - 1] / (3.0 * V);
         }
      }
      S[ivertex] = S_ik;
   }
}