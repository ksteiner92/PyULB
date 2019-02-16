//
// Created by klaus on 06.01.19.
//

#ifndef LBM_DELAUNAY_H
#define LBM_DELAUNAY_H

#include <array>
#include <vector>
#include <unordered_map>

#include "mesh.h"
#include "utils.h"
#include "eigen.h"

class Circle
{
public:
   Circle(const Eigen::Vector2d& center, double rad);

   Eigen::Vector2d getCenter() const;

   double getRadius() const;

   static Circle circumCircle(const Eigen::Vector2d& a, const Eigen::Vector2d& b, const Eigen::Vector2d& c);

private:
   Eigen::Vector2d center;
   double rad;
};

class Triangle
{
public:
   Triangle() : pts(nullptr), pidx({-1}) {}

   Triangle(std::vector<Eigen::Vector2d>* pts, const std::array<int, 3>& pidx)
      : pts(pts), pidx(pidx) {}

   virtual ~Triangle() {}

   int operator[](int idx) const
   {
      if (idx < 0 || idx >= 3) {
         std::cout << "Triangle: Index out of bounds" << std::endl;
         std::cout.flush();
         throw std::out_of_range("Index out of bounds");
      }
      return pidx[idx];
   }

private:
   std::array<int, 3> pidx;
   std::vector<Eigen::Vector2d>* pts;
};

class Delaunay2D : public IMesher<2, 2>
{
public:
   Delaunay2D();

   void generate(Mesh<2, 2>& mesh) override;

private:
   struct TriangleTreeElement
   {
      TriangleTreeElement() : stat(1), d({-1}), p({0}) {}
      int stat;
      std::array<int, 3> d;
      std::array<int, 3> p;
   };

   static constexpr double fuzz = 1.0e-6;
   static constexpr double epsilon = 1.0e-10;
   static constexpr double bigscale = 1000.0;
   typedef unsigned long long int Ullong;

   std::vector<Eigen::Vector2d> pts;
   std::unordered_map<Ullong, int, NullHash> linmap;
   std::unordered_map<Ullong, int, NullHash> trimap;
   std::vector<int> hull;
   std::vector<int> pidx;
   LineHash linhash;
   TriangleHash trihash;
   std::vector<TriangleTreeElement> triangles;
   std::vector<int> perm;
   size_t t;
   int ntree;
   int ntri;
   int mostleftbottom;
   double delx;
   double dely;
   static unsigned int jran;
   int option;
   int npts;
   int ntreemax;

   void calcHull();

   static double orientation(const Eigen::Vector2d &p,
                             const Eigen::Vector2d &q,
                             const Eigen::Vector2d &i);

   inline int pointInTriangle(const TriangleTreeElement& trielm, const Eigen::Vector2d& p) const;

   void insertPoint(int r);

   int whichContainsPoint(const Eigen::Vector2d& p, int strict = 0);

   int storeTriangle(int a, int b, int c);

   void eraseTriangle(int a, int b, int c, int d0, int d1, int d2);

   static inline double inCircle(const Eigen::Vector2d &d,
                                 const Eigen::Vector2d &a,
                                 const Eigen::Vector2d &b,
                                 const Eigen::Vector2d &c);

};


#endif //LBM_DELAUNAY_H
