//
// Created by klaus on 16.02.19.
//

#include <iostream>
#include <vector>

#include "eigen.h"
#include "mesh.h"
#include "delauny.h"

using namespace std;
using namespace Eigen;

int main(int argc, char** argv)
{
   const size_t nx = 5;
   const size_t ny = 5;
   const double lx = 40.0;
   const double ly = 40.0;
   const double dx = lx / (double) (nx - 1);
   const double dy = ly / (double) (ny - 1);

   vector<Vector2d> points(nx * ny);
   for (size_t x = 0; x < nx; x++)
      for (size_t y = 0; y < ny; y++)
         points[x * ny + y] << x * dx, y * dy;
   Mesh<2, 2> mesh(&points);
   Delaunay2D d;
   d.generate(mesh);
   Mesh<2, 1>* hull = mesh.getHull();
   return 0;
}