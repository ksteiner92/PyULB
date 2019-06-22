//
// Created by klaus on 17.04.19.
//

#include <array>
#include <vector>

#include "fvm.h"
#include "logger.h"
#define VOID void
#define REAL double
#include "triangle/triangle.h"

using namespace std;
using namespace Eigen;

FVM::FVM(IMesh *mesh)
{
   const size_t npts = mesh->getNumVertices();
   cout << "npts: " << npts << endl;
   triangulateio triio;
   triio.pointlist = &mesh->getCoords()[0];
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

   triangulateio vout;
   vout.pointlist = nullptr;            /* Not needed if -N switch used. */
   /* Not needed if -N switch used or number of point attributes is zero: */
   vout.pointattributelist = nullptr;
   vout.pointmarkerlist = nullptr; /* Not needed if -N or -B switch used. */
   vout.trianglelist = nullptr;          /* Not needed if -E switch used. */
   /* Not needed if -E switch used or number of triangle attributes is zero: */
   vout.triangleattributelist = nullptr;
   vout.neighborlist = nullptr;         /* Needed only if -n switch used. */
   /* Needed only if segments are output (-p or -c) and -P not used: */
   vout.segmentlist = nullptr;
   /* Needed only if segments are output (-p or -c) and -P and -B not used: */
   vout.segmentmarkerlist = nullptr;
   vout.edgelist = nullptr;             /* Needed only if -e switch used. */
   vout.edgemarkerlist = nullptr;   /* Needed if -e used and -B not used. */

   triangulate((char *) "veczq", &triio, &triout, &vout);
   cout << "Triangulation done" << endl;
   cout.flush();
   vector<double> voronoi_coords(vout.numberofpoints * 2);
   copy(vout.pointlist, vout.pointlist + vout.numberofpoints * 2, voronoi_coords.begin());
   voronoi = make_unique<Mesh<2, 2>>(voronoi_coords);
   for (size_t iv = 0; iv < vout.numberofpoints; iv++) {
      voronoi->getOrCreateVertexByID(iv);
   }
   for (size_t ie = 0; ie < vout.numberofedges; ie++) {
      if (vout.edgelist[ie * 2 + 0] > 0 && vout.edgelist[ie * 2 + 1] > 0)
         Edge<2>* edge = voronoi->getOrCreateEdge(vout.edgelist[ie * 2 + 0], vout.edgelist[ie * 2 + 1]);
   }

}

Mesh<2, 2>* FVM::getDualMesh() const
{
   return voronoi.get();
}