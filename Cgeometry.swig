%module Cgeometry;

%{
#include "geometry2D.h"

%}

%include "cinder-lite/Vector.h"

%template(Vec2) cinderlite::Vec2<double>;

%include "geometry2D.h"



