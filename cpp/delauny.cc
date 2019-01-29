//
// Created by klaus on 06.01.19.
//

#include <iostream>
#include <memory>

#include "delauny.h"
#include "utils.h"

using namespace std;
using namespace Eigen;

Circle::Circle(const Vector2d& center, double rad)
   : center(center), rad(rad)
{}

Vector2d Circle::getCenter() const
{
   return center;
}

double Circle::getRadius() const
{
   return rad;
}

Circle Circle::circumCircle(const Vector2d& a, const Vector2d& b, const Vector2d& c)
{
   const double a0 = a[0] - b[0];
   const double a1 = a[1] - b[1];
   const double c0 = c[0] - b[0];
   const double c1 = c[1] - b[1];
   double det = a0 * c1 - c0 * a1;
   if (det == 0.0) {
      std::cout << "Circle: Index out of bounds" << std::endl;
      std::cout.flush();
      throw ("No circle through colinear points");
   }
   det = 0.5 / det;
   const double asq = a0 * a0 + a1 * a1;
   const double csq = c0 * c0 + c1 * c1;
   const double ctr0 = det * (asq * c1 - csq * a1);
   const double ctr1 = det * (csq * a0 - asq * c0);
   const double rad2 = ctr0 * ctr0 + ctr1 * ctr1;
   return Circle({ctr0 + b[0], ctr1 + b[1]}, sqrt(rad2));
}
Delaunay2D::Delaunay2D()
   : npts(0), ntri(0), ntree(0), ntreemax(0)
{}

void Delaunay2D::generate(Mesh<2> &mesh)
{
   npts = mesh.points->size();
   ntri = 0;
   ntree = 0;
   ntreemax = 10 * npts + 1000;
   triangles.resize(ntreemax);
   pts.resize(npts + 3);
   perm.resize(npts);
   mesh.vertices.reserve(npts);
   linmap.reserve(6 * npts + 12);
   trimap.reserve(2 * npts + 6);
   double x1 = (*mesh.points)[0][0];
   double xh = x1;
   double y1 = (*mesh.points)[0][1];
   double yh = y1;

   for (int i = 0; i < npts; i++) {
      pts[i] = (*mesh.points)[i];
      mesh.vertices.push_back(make_unique<Vertex<2>>(Vertex<2>(&mesh, i, i)));
      perm[i] = i;
      x1 = min(x1, (*mesh.points)[i][0]);
      xh = max(xh, (*mesh.points)[i][0]);
      y1 = min(y1, (*mesh.points)[i][1]);
      yh = max(yh, (*mesh.points)[i][1]);
   }

   delx = xh - x1;
   dely = yh - y1;
   pts[npts + 0] = {0.5 * (x1 + xh), xh + bigscale * dely};
   pts[npts + 1] = {x1 - 0.5 * bigscale * delx, y1 - 0.5 * bigscale * dely};
   pts[npts + 2] = {xh + 0.5 * bigscale * delx, y1 - 0.5 * bigscale * dely};
   storeTriangle(npts, npts + 1, npts + 2);
   for (int i = npts; i > 0; i--)
      swap(perm[i - 1], perm[RandomHash::int64(jran++) % i]);
   for (int i = 0; i < npts; i++)
      insertPoint(perm[i]);

   mesh.edges.reserve(linmap.size());
   mesh.faces.reserve(trimap.size());
   for (int i = 0; i < ntree; i++) {
      if (triangles[i].stat > 0) {
         if (triangles[i].p[0] >= npts || triangles[i].p[1] >= npts || triangles[i].p[2] >= npts) {
            triangles[i].stat = -1;
            ntri--;
         } else {
            Edge<2>* e1 = getOrCreateEdge(mesh, triangles[i].p[1], triangles[i].p[2]);
            Edge<2>* e2 = getOrCreateEdge(mesh, triangles[i].p[2], triangles[i].p[0]);
            Edge<2>* e3 = getOrCreateEdge(mesh, triangles[i].p[0], triangles[i].p[1]);
            getOrCreateFace(mesh, e1, e2, e3);
         }
      }
   }
}

/*double Delaunay2D::interpolate(const Vector2d* p, const std::vector<double>& fnvals, double defval) const
{

}*/

void Delaunay2D::getFaces(vector<Triangle> &f) const
{
   /*f.resize(faces.size());
   copy(faces.begin(), faces.end(), f.begin());*/
}

void Delaunay2D::insertPoint(int r)
{
   int i;
   int tidx = -1;
   for (i = 0; i < 3; i++) {
      tidx = whichContainsPoint(pts[r], 1);
      if (tidx >= 0)
         break;
      pts[r][0] += fuzz * delx * (RandomHash::doub(jran++) - 0.5);
      pts[r][1] += fuzz * dely * (RandomHash::doub(jran++) - 0.5);
   }
   if (i == 3) {
      cout << "Points degenerated even after fuzzing" << endl;
      cout.flush();
      throw ("Points degenerated even after fuzzing");
   }
   int ntask = 0;
   i = triangles[tidx].p[0];
   int j = triangles[tidx].p[1];
   int k = triangles[tidx].p[2];
   if (option & 2 && i < npts && j < npts && k < npts)
      return;
   array<int, 50> tasks;
   array<int, 50> taski;
   array<int, 50> taskj;
   int d0 = storeTriangle(r, i, j);
   tasks[++ntask] = r;
   taski[ntask] = i;
   taskj[ntask] = j;
   int d1 = storeTriangle(r, j, k);
   tasks[++ntask] = r;
   taski[ntask] = j;
   taskj[ntask] = k;
   int d2 = storeTriangle(r, k, i);
   tasks[++ntask] = r;
   taski[ntask] = k;
   taskj[ntask] = i;
   eraseTriangle(i, j, k, d0, d1, d2);
   while (ntask) {
      const int s = tasks[ntask];
      i = taski[ntask];
      j = taskj[ntask--];
      unsigned long long int key = linhash(make_tuple(j, i));
      const auto lit = linmap.find(key);
      if (lit == linmap.end())
         continue;
      const int l = lit->second;
      if (inCircle(pts[l], pts[j], pts[s], pts[i]) > 0.0) {
         d0 = storeTriangle(s, l, j);
         d1 = storeTriangle(s, i, l);
         eraseTriangle(s, i, j, d0, d1, -1);
         eraseTriangle(l, j, i, d0, d1, -1);
         key = linhash(make_tuple(i, j));
         linmap.erase(key);
         key = 0 - key;
         linmap.erase(key);
         tasks[++ntask] = s;
         taski[ntask] = l;
         taskj[ntask] = j;
         tasks[++ntask] = s;
         taski[ntask] = i;
         taskj[ntask] = l;
      }
   }
}

int Delaunay2D::pointInTriangle(const TriangleTreeElement& trielm, const Vector2d& pt) const
{
   int res = 0;
   for (int i = 0; i < 3; i++) {
      const int j = (i + 1) % 3;
      const double d = (pts[trielm.p[j]][0] - pts[trielm.p[i]][0]) * (pt[1] - pts[trielm.p[i]][1]) -
                       (pts[trielm.p[j]][1] - pts[trielm.p[i]][1]) * (pt[0] - pts[trielm.p[i]][0]);
      if (d < 0.0)
         return -1;
      if (d == 0.0)
         res = 1;
   }
   return (res ? 0 : 1);
}

double Delaunay2D::inCircle(const Vector2d& d, const Vector2d& a, const Vector2d& b, const Vector2d& c)
{
   const Circle cc = Circle::circumCircle(a, b, c);
   const double radd = sqr(d[0] - cc.getCenter()[0]) + sqr(d[1] - cc.getCenter()[1]);
   return (sqr(cc.getRadius()) - radd);
}


int Delaunay2D::whichContainsPoint(const Vector2d& p, int strict)
{
   int k = 0;
   int j;
   int i;
   while (triangles[k].stat <= 0) {
      for (i = 0; i < 3; i++) {
         if ((j = triangles[k].d[i]) < 0)
            continue;
         if (strict) {
            if (pointInTriangle(triangles[j], p) > 0)
               break;
         } else {
            if (pointInTriangle(triangles[j], p) >= 0)
               break;
         }
      }
      if (i == 3)
         return -1;
      k = j;
   }
   return k;
}

int Delaunay2D::storeTriangle(int a, int b, int c)
{
   triangles[ntree].p = {a, b, c};
   trimap[trihash(make_tuple(a, b, c))] = ntree;
   linmap[linhash(make_tuple(b, c))] = a;
   linmap[linhash(make_tuple(c, a))] = b;
   linmap[linhash(make_tuple(a, b))] = c;
   if (++ntree == ntreemax) {
      cout << "Faces list size is too small" << endl;
      cout.flush();
      throw ("Faces list size is too small");
   }
   ntri++;
   return (ntree - 1);
}

void Delaunay2D::eraseTriangle(int a, int b, int c, int d0, int d1, int d2)
{
   const auto it = trimap.find(trihash(make_tuple(a, b, c)));
   if (it == trimap.end()) {
      cout << "Non existent triangle" << endl;
      cout.flush();
      throw ("Non existent triangle");
   }
   const int i = it->second;
   triangles[i].d[0] = d0;
   triangles[i].d[1] = d1;
   triangles[i].d[2] = d2;
   triangles[i].stat = 0;
   ntri--;
   trimap.erase(it);
}

unsigned int Delaunay2D::jran = 14921620;
