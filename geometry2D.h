

#pragma once


#include "cinder-lite/Vector.h"


using Vec2 = cinderlite::Vec2<double>;


int contact_3circle(Vec2 c1, double r1,
		    Vec2 c2, double r2,
		    double r3,
		    Vec2 &c31, Vec2 &c32,
		    bool inside_c2=false);
