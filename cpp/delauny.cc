//
// Created by klaus on 06.01.19.
//

#include <iostream>
#include <memory>

#include "delauny.h"
#include "utils.h"
#include "logger.h"

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
   if (det == 0.0)
      throw logic_error("No circle through colinear points");
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

double Delaunay2D::orientation(const Vector2d &p,
                               const Vector2d &q,
                               const Vector2d &i)
{
   const Vector2d pq = q - p;
   const Vector2d pi = i - p;
   const double o = pq(0) * pi(1) - pq(1) * pi(0);
   if (abs(o) > Delaunay2D::epsilon)
      return o;
   return 0.0;
}

void Delaunay2D::calcHull()
{
   cout << "Initialize hull calculation ..." << endl;
   pidx.resize(npts);
   iota(pidx.begin(), pidx.end(), 0);
   swap(pidx[0], pidx[mostleftbottom]);
   cout << "Sorting points ..." << endl;
   Vector2d p = pts[mostleftbottom];
   size_t pstart = 0;
   double maxcolin = 0.0;
   sort(pidx.begin() + 1, pidx.end(), [&p, &maxcolin, &pstart, this](int a, int b) {
      const double o = Delaunay2D::orientation(p, pts[a], pts[b]);
      if (o == 0.0) {
         const double pa = (pts[a] - p).squaredNorm();
         const double pb = (pts[b] - p).squaredNorm();
         const double m = max(pa, pb);
         const bool res = pa < pb;
         if (m > maxcolin) {
            maxcolin = m;
            pstart = res ? b : a;
         }
         return res;
      }
      return o > 0.0;
   });
   cout << "Calculating hull ..." << endl;
   size_t istart = 0;
   hull.push_back(mostleftbottom);
   if (maxcolin > 0) {
      const auto itstart = find(pidx.begin(), pidx.end(), pstart);
      cout << "it found" << endl;
      cout.flush();
      hull.insert(hull.begin() + 1, pidx.begin() + 1, itstart + 1);
      cout << "inserted" << endl;
      cout.flush();
      istart = distance(pidx.begin(), itstart);
   }
   /*p = pts[pidx[istart]];
   Vector2d q = pts[pidx[istart + 1]];
   Vector2d i = pts[pidx[istart + 2]];*/
   cout << "istart: " << istart << endl;
   cout.flush();
   hull.push_back(pidx[istart + 1]);
   hull.push_back(pidx[istart + 2]);
   double o;
   int top;
   for (size_t i = 3; i < pidx.size(); i++) {
      top = hull.back();
      hull.pop_back();
      while ((o = orientation(pts[hull.back()], pts[top], pts[pidx[i]])) < 0)   {
         top = hull.back();
         hull.pop_back();
      }
      if (o == 0.0) {
         const size_t chkpoint = hull.size();
         do {
            hull.push_back(top);
            hull.push_back(pidx[i]);
            top = hull.back();
         } while ((o = orientation(pts[hull.back()], pts[top], pts[pidx[i]])) == 0);
         if (o > 0.0)
            hull.resize(chkpoint);
      }
      hull.push_back(top);
      hull.push_back(pidx[i]);
   }
   /*for (size_t ip = istart + 1; ip < (pidx.size() - 2); ip++) {
      const double o = orientation(p, q, i);
      if (o > 0) {
         hull.push_back(pidx[ip]);
         p = q;
         q = i;
         if (ip == (pidx.size() - 3))
            break;
         i = pts[ip + 2];
         cout << "#" << ip << " is boundary" << endl;
      } else if (o == 0.0) {
         cout << "#" << ip << " is colinear" << endl;
      } else {
         cout << "#" << ip << " is not boundary" << endl;
         q = p;
         p = pts[hull[hull.size() - 2]];
         //ip--;
      }
   }*/
   hull.push_back(pidx[istart]);
   cout << "Hull contains " << (hull.size() - 1) << " points" << endl;
}

void Delaunay2D::generate(Mesh<2, 2> &mesh)
{
   npts = mesh.getNumPoints();
   if (npts < 3)
      throw logic_error("There is no grid for only two points");
   ntri = 0;
   ntree = 0;
   ntreemax = 10 * npts + 1000;
   triangles.resize(ntreemax);
   pts.resize(npts + 3);
   perm.resize(npts);
   linmap.reserve(6 * npts + 12);
   trimap.reserve(2 * npts + 6);
   const Vector2d& p = mesh.getPoint(0);
   double x1 = p[0];
   double xh = x1;
   double y1 = p[1];
   double yh = y1;

   mostleftbottom = 0;
   for (int i = 0; i < npts; i++) {
      const Vector2d& p = mesh.getPoint(i);
      mesh.getOrCreateVertexByID(i);
      pts[i] = p;
      perm[i] = i;
      const double xdiff = p[0] - x1;
      const bool ylower = p[1] < y1;
      if (ylower)
         y1 = p[1];
      if (abs(xdiff) <= epsilon) {
         if (ylower)
            mostleftbottom = i;
      } else if (xdiff < 0.0) {
         x1 = p[0];
         mostleftbottom = i;
      }
      xh = max(xh, p[0]);
      yh = max(yh, p[1]);
   }

   //calcHull();

   delx = xh - x1;
   dely = yh - y1;
   pts[npts + 0] = {0.5 * (x1 + xh), xh + bigscale * dely};
   pts[npts + 1] = {x1 - 0.5 * bigscale * delx, y1 - 0.5 * bigscale * dely};
   pts[npts + 2] = {xh + 0.5 * bigscale * delx, y1 - 0.5 * bigscale * dely};
   storeTriangle(npts, npts + 1, npts + 2);
   for (int i = npts; i > 0; i--)
      swap(perm[i - 1], perm[RandomHash::int64(jran++) % i]);
   LOG_T(INFO) << "Meshing ..." << LogFlags::ENDL;
   ProgressBar progress;
   progress.start(npts);
   for (int i = 0; i < npts; i++) {
      insertPoint(perm[i]);
      progress.update(i);
   }
   progress.stop();

   LOG_T(INFO) << "Storing mesh information ..." << LogFlags::ENDL;
   progress.start(ntree);
   vector<uint8_t> nneibours;
   for (int i = 0; i < ntree; i++) {
      if (triangles[i].stat > 0) {
         if (triangles[i].p[0] >= npts || triangles[i].p[1] >= npts || triangles[i].p[2] >= npts) {
            triangles[i].stat = -1;
            ntri--;
         } else {
            Edge<2>* e1 = mesh.getOrCreateEdge(triangles[i].p[1], triangles[i].p[2]);
            Edge<2>* e2 = mesh.getOrCreateEdge(triangles[i].p[2], triangles[i].p[0]);
            Edge<2>* e3 = mesh.getOrCreateEdge(triangles[i].p[0], triangles[i].p[1]);
            mesh.getOrCreateFace(e1, e2, e3);
            const size_t maxid = max(e1->getID(), max(e2->getID(), e3->getID()));
            if (maxid >= nneibours.size())
               nneibours.resize(maxid + 1, 0);
            nneibours[e1->getID()]++;
            nneibours[e2->getID()]++;
            nneibours[e3->getID()]++;
         }
      }
      progress.update(i);
   }
   progress.stop();
   LOG_T(INFO) << "Storing hull information ..." << LogFlags::ENDL;
   Mesh<2, 1>* hull = mesh.getHull();
   for (size_t i = 0; i < nneibours.size(); i++)
      if (nneibours[i] == 1)
         hull->getOrCreateEdgeByID(i);
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
   if (i == 3)
      throw invalid_argument("Points degenerated even after fuzzing");
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
   if (++ntree == ntreemax)
      throw overflow_error("Faces list size is too small");
   ntri++;
   return (ntree - 1);
}

void Delaunay2D::eraseTriangle(int a, int b, int c, int d0, int d1, int d2)
{
   const auto it = trimap.find(trihash(make_tuple(a, b, c)));
   if (it == trimap.end())
      throw range_error("Non existent triangle");
   const int i = it->second;
   triangles[i].d[0] = d0;
   triangles[i].d[1] = d1;
   triangles[i].d[2] = d2;
   triangles[i].stat = 0;
   ntri--;
   trimap.erase(it);
}

unsigned int Delaunay2D::jran = 14921620;
