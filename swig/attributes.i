%{
#include "attributes.h"
%}

%include <std_vector.i>
%include <cpointer.i>
%include <stdint.i>
%include "attributes.h"

%template(AttributeDouble) Attribute<double>;
%template(AttributeInt) Attribute<int>;
%template(AttributeUInt8) Attribute<uint8_t>;
%template(AttributeIntList) Attribute<std::vector<int>>;
%template(AttributePoint1D) Attribute<Eigen::Matrix<double, 1, 1>>;
%template(AttributePoint2D) Attribute<Eigen::Matrix<double, 2, 1>>;
%template(AttributePoint3D) Attribute<Eigen::Matrix<double, 3, 1>>;
%template(AttributePoint1DList) Attribute<std::vector<Eigen::Matrix<double, 1, 1>>>;
%template(AttributePoint2DList) Attribute<std::vector<Eigen::Matrix<double, 2, 1>>>;
%template(AttributePoint3DList) Attribute<std::vector<Eigen::Matrix<double, 3, 1>>>;
%pointer_class(double, doublep);
%pointer_class(uint8_t, uint8_tp);