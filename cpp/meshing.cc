//
// Created by klaus on 08.03.19.
//

#include "meshing.h"
#include "logger.h"
#include "eigen.h"
#define VOID void
#define REAL double
#include "triangle/triangle.h"

using namespace std;
using namespace Eigen;

template<uint Dim, uint TopDim>
Meshing<Dim, TopDim>::Meshing()
{
}

template<uint Dim, uint TopDim>
void Meshing<Dim, TopDim>::generate(Mesh<Dim, TopDim>& mesh)
{
}

template<>
Meshing<2, 2>::Meshing()
{

}

template<>
void Meshing<2, 2>::generate(Mesh<2, 2>& mesh)
{
   const size_t npts = mesh.getNumPoints();
   if (npts < 3)
      throw logic_error("There is no grid for only two points");

   vector<double> coords(npts * 2, 0.0);
   for (int i = 0; i < npts; i++) {
      const Vector2d &p = mesh.getPoint(i);
      coords[i * 2 + 0] = p(0);
      coords[i * 2 + 1] = p(1);
   }

   triangulateio triio;
   triio.pointlist = &coords[0];
   triio.numberofpoints = npts;
   triio.numberofsegments = 0;
   vector<int> edgemarks(npts, 0);
   triio.edgemarkerlist = &edgemarks[0];
   vector<int> pointmarks(npts, 0);
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

   triangulate((char *) "peczqQ", &triio, &triout, nullptr);

   Mesh<2, 1>* hull = dynamic_cast<Mesh<2, 1>*>(mesh.getHull());
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
      array<Edge<2>*, 3> edges = {
              mesh.getEdge(triout.trianglelist[it * 3 + 0], triout.trianglelist[it * 3 + 1]),
              mesh.getEdge(triout.trianglelist[it * 3 + 0], triout.trianglelist[it * 3 + 2]),
              mesh.getEdge(triout.trianglelist[it * 3 + 1], triout.trianglelist[it * 3 + 2])
      };
      mesh.getOrCreateFace(edges[0], edges[1], edges[2]);
   }
}

template<>
Meshing<3, 3>::Meshing()
{

}

template<>
void Meshing<3, 3>::generate(Mesh<3, 3>& mesh)
{

}
