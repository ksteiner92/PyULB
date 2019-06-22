%module pyulb
%{
#define SWIG_FILE_WITH_INIT
#include <Python.h>
#include "lbm.h"
//#include "delauny.h"
%}

%include <typemaps.i>
%include <std_vector.i>
%include <std_string.i>
%include <std_array.i>
%include "eigen.i"

%eigen_typemaps(Eigen::Vector2d)
%eigen_typemaps(Eigen::Vector3d)
%eigen_typemaps(Eigen::Matrix<double, 1, 1>)
%eigen_typemaps(Eigen::Matrix<double, 2, 1>)
%eigen_typemaps(Eigen::Matrix<double, 3, 1>)

%template(Point1DList) std::vector<Eigen::Matrix<double, 1, 1>>;
%template(Point2DList) std::vector<Eigen::Matrix<double, 2, 1>>;
%template(Point3DList) std::vector<Eigen::Matrix<double, 3, 1>>;

//%include "delauny.h"
%include "mesh.i"
%include "poisson.i"

%template(VectorDouble) std::vector<double>;
%template(VectorInt) std::vector<int>;

%template(ArrayInt1) std::array<int, 1>;
%template(ArrayInt2) std::array<int, 2>;
%template(ArrayInt3) std::array<int, 3>;
%template(ArrayInt4) std::array<int, 4>;

//%ignore Circle;
//%ignore Triangle;
//%nodefaultctor Delaunay2D;

%include "lbm.h"