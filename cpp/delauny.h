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
   Eigen::Vector2d centeroid;
   std::array<int, 3> pidx;
   std::vector<Eigen::Vector2d>* pts;
};

class Delaunay2D : public MeshGenerator<2>
{
public:
   Delaunay2D();

   //double interpolate(const Eigen::Vector2d* p, const std::vector<double>& fnvals, double defval = 0.0) const;

   void getFaces(std::vector<Triangle> &f) const;

   void getHull(std::vector<int>& hull) const;

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

   inline Edge<2>* getOrCreateEdge(Mesh<2> &mesh, int a, int b) const
   {
      auto key = linhash(std::make_tuple(a, b));
      auto it = mesh.edge2edge.find(key);
      if (it == mesh.edge2edge.end()) {
         key = linhash(std::make_tuple(b, a));
         it = mesh.edge2edge.find(key);
      }
      if (it == mesh.edge2edge.end()) {
         mesh.edges.push_back(std::make_unique<Edge<2>>(Edge<2>(&mesh,
                 {mesh.vertices[a].get(), mesh.vertices[b].get()}, mesh.edges.size())));
         Edge<2>* edge = mesh.edges.back().get();
         mesh.edge2edge[key] = edge;
         mesh.point2edge.insert({a, edge});
         mesh.point2edge.insert({b, edge});
         return edge;
      }
      return it->second;
   }

   inline Face<2>* getOrCreateFace(Mesh<2> &mesh, Edge<2>* e1, Edge<2>* e2, Edge<2>* e3) const
   {
      const auto key = trihash(std::make_tuple(e1->getID(), e2->getID(), e3->getID()));
      const auto it = mesh.face2face.find(key);
      if (it == mesh.face2face.end()) {
         mesh.faces.push_back(std::make_shared<Face<2>>(Face<2>(&mesh, {e1, e2, e3}, mesh.faces.size())));
         Face<2>* face = mesh.faces.back().get();
         mesh.face2face[key] = face;
         for (uint i = 0; i < 3; i++)
            mesh.point2face.insert(std::make_pair((*face)[i], face));
         //mesh.point2edge.insert({key, e1});
         //mesh.point2edge.insert({b, e1});
         return face;
      }
      return it->second;
   }

   static inline double inCircle(const Eigen::Vector2d& d, const Eigen::Vector2d& a, const Eigen::Vector2d& b, const Eigen::Vector2d& c);

};


#endif //LBM_DELAUNAY_H
