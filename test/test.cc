//
// Created by klaus on 16.02.19.
//

#include <iostream>
#include <vector>
#include <set>
#include <tuple>

#include "eigen.h"
#include "mesh.h"
#include "delauny.h"
#include "lbm.h"
#include "logger.h"

#include "tetgen/tetgen.h"
#define VOID void
#include "triangle/triangle.h"

using namespace std;
using namespace Eigen;

int main(int argc, char** argv)
{
   const size_t nx = 50;
   const size_t ny = 50;
   const double lx = 1.0;
   const double ly = 1.0;
   const double dx = lx / (double) (nx - 1);
   const double dy = ly / (double) (ny - 1);

   vector<Vector2d> pts(nx * ny);

   vector<double> points(nx * ny * 2, 0.0);


   for (size_t x = 0; x < nx; x++) {
      for (size_t y = 0; y < ny; y++) {
         pts[y * nx + x] << x * dx, y * dy;
         points[y * nx * 2 + x * 2 + 0] = pts[y * nx + x](0);
         points[y * nx * 2 + x * 2 + 1] = pts[y * nx + x](1);
      }
   }

   triangulateio triio;
   triio.pointlist = &points[0];
   triio.numberofpoints = nx * ny;
   triio.numberofsegments = 0;
   vector<int> edgemarks(nx * ny, 0);
   triio.edgemarkerlist = &edgemarks[0];
   vector<int> pointmarks(nx * ny, 0);
   triio.pointmarkerlist = &pointmarks[0];
   triio.numberofpointattributes = 0;
   triio.pointmarkerlist = nullptr;
   triio.numberofholes = 0;
   triio.numberofregions = 0;

   triangulateio triout;

   triout.pointlist = nullptr;            /* Not needed if -N switch used. */
   /* Not needed if -N switch used or number of point attributes is zero: */
   triout.pointattributelist = nullptr;
   triout.pointmarkerlist = nullptr; /* Not needed if -N or -B switch used. */
   triout.trianglelist = nullptr;          /* Not needed if -E switch used. */
   /* Not needed if -E switch used or number of triangle attributes is zero: */
   triout.triangleattributelist = nullptr;
   triout.neighborlist = nullptr;         /* Needed only if -n switch used. */
   /* Needed only if segments are output (-p or -c) and -P not used: */
   triout.segmentlist = nullptr;
   /* Needed only if segments are output (-p or -c) and -P and -B not used: */
   triout.segmentmarkerlist = nullptr;
   triout.edgelist = nullptr;             /* Needed only if -e switch used. */
   triout.edgemarkerlist = nullptr;   /* Needed if -e used and -B not used. */

   triangulate((char *) "peczq", &triio, &triout, nullptr);

   Mesh<2, 2> mesh(&pts);
   Mesh<2, 1>* hull = dynamic_cast<Mesh<2, 1>*>(mesh.getHull());
   size_t hullcount = 0;
   for (size_t iv = 0; iv < triout.numberofpoints; iv++) {
      mesh.getOrCreateVertexByID(iv);
      if (triout.pointmarkerlist[iv] > 0)
         hull->getOrCreateVertexByID(iv);
   }
   for (size_t ie = 0; ie < triout.numberofedges; ie++) {
      Edge<2>* edge = mesh.getOrCreateEdge(triout.edgelist[ie * 2 + 0], triout.edgelist[ie * 2 + 1]);
      if (triout.edgemarkerlist[ie] > 0) {
         hull->getOrCreateEdgeByID(edge->getID());
      }
   }
   for (size_t it = 0; it < triout.numberoftriangles; it++) {
      array<Edge<2>*, 3> edges  = {
         mesh.getEdge(triout.trianglelist[it * 3 + 0], triout.trianglelist[it * 3 + 1]),
         mesh.getEdge(triout.trianglelist[it * 3 + 0], triout.trianglelist[it * 3 + 2]),
         mesh.getEdge(triout.trianglelist[it * 3 + 1], triout.trianglelist[it * 3 + 2])
      };
      mesh.getOrCreateFace(edges[0], edges[1], edges[2]);
   }

   /*Mesh<2, 2> mesh(&points);
   Delaunay2D d;
   d.generate(mesh);

   Edge<2>* edge = mesh.getEdge(5);
   cout << "Edge 5: " << (*edge)[0] << ", " << (*edge)[1] << endl;
   set<tuple<size_t, size_t, size_t>> facesid;
   for (size_t i = 0; i < mesh.getNumFaces(); i++) {
      Face<2>* face = mesh.getFace(i);
      facesid.insert(make_tuple((*face)[0], (*face)[1], (*face)[2]));
   }
   cout << "Num Faces: " << facesid.size() << endl;*/
   //LatticeBoltzmann2D lbm(&mesh, {dynamic_cast<Mesh<2, 1>*>(mesh.getHull())}, 1.0, 1.0 / sqrt(3.0));
   //lbm.calculate(0.1, 10);



   return 0;
}