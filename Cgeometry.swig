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

%inline %{

size_t swig_ptr_as_int(const void *x) {
    return (size_t)x;
}

Circle* downcast_Circle(Shape2D *shape) {
   return dynamic_cast<Circle*>(shape);
}


%}






%pythoncode %{

def Vec2__repr__(self):
   return "Vec2(%f, %f)" % (self.x, self.y)

Vec2.__repr__ = Vec2__repr__

def BBox__repr__(self):
   return "Bbox(%s, %s)" % (self.Cmin, self.Cmax)

BBox.__repr__ = BBox__repr__


# to convert a C++ style iterator to one that can be used
# directly in python

def py_style_iterator(self):
    while self.has_next():
        yield self.next()       

IntersectingLeavesIterator.__iter__ = py_style_iterator
ShapeVectorIterator.__iter__ = py_style_iterator
LeafIterator.__iter__ = py_style_iterator
TwoLeavesIterator.__iter__ = py_style_iterator

%}

