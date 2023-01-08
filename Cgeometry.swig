%module Cgeometry;

%{
#include "geometry2D.h"

%}

%include "cinder-lite/Vector.h"

%template(Vec2) cinderlite::Vec2<double>;

%include "geometry2D.h"


/*
%extend Vec2 {
  std::string __repr__() const {
    return "toto";
  }
}
*/

%pythoncode %{

def Vec2__repr__(self):
   return "Vec2(%f, %f)" % (self.x, self.y)

Vec2.__repr__ = Vec2__repr__




%}

