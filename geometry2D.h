

#pragma once


#include "cinder-lite/Vector.h"


using Vec2 = cinderlite::Vec2<double>;


int contact_3circle(Vec2 c1, double r1,
		    Vec2 c2, double r2,
		    double r3,
		    Vec2 &c31, Vec2 &c32,
		    bool inside_c2=false);

struct BBox {
    Vec2 Cmin, Cmax; 
};

struct Object2D {

    virtual bool intersects(const BBox & bbox) const = 0;

    virtual ~Object2D() {}

};
    

struct Circle: Object2D {
    Vec2 c;
    double r;
    int id;

    Circle(Vec2 c, double r, int id=-1):
    c(c), r(r), id(id) {}

    Circle(double x, double y, double r, int id=-1):
    c(x, y), r(r), id(id) {}

    bool intersects(const BBox & bbox) const final;

    bool intersects(const Circle & other) const; 
    
    virtual ~Circle() {}

};
