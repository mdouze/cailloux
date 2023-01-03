
#include "geometry2D.h"


#include <cmath>

int solve_2nd_degree(double a, double b, double c,
		     double &x0, double &x1) {
    double delta = b * b - 4 * a * c;
    if (delta < 0 || a == 0) {
       return 0;
    } else {	
	double sdelta = sqrt(delta); 
	x0 = (-b + sdelta) / (2 * a);
	x1 = (-b - sdelta) / (2 * a);
	return 2; 
    }
}


// find one point on the line defined by a.x = b
Vec2 find_origin(Vec2 a, double b) {
    if (fabs(a.x) > fabs(a.y)) {
	return Vec2(b / a.x, 0); 
    } else {
	return Vec2(0, b / a.y); 
    }
}


double sqr(double x) {
    return x * x;
}

// Find c3 st. (c1, r1), (c2, r2) and (c3, r3) are tangent    
int contact_3circle(Vec2 c1, double r1,
		    Vec2 c2, double r2,
		    double r3,
		    Vec2 &c31, Vec2 &c32,
		    bool inside_c2) {


    double s_r3 = inside_c2 ? -r3 : r3;
    
    double b = (sqr(r1 + r3) - sqr(r2 + s_r3) - c1.lengthSquared() + c2.lengthSquared()) / 2;
    Vec2 a = c2 - c1;
    // a.c3 = b
    Vec2 aorth = Vec2(a.y, -a.x); // vector orthogonal to a
    Vec2 o = find_origin(a, b);
        
    // c3 = o + t * aorth, find t
    Vec2 oprime = o - c1;
    double t1, t2; 
        
    int nres = solve_2nd_degree(
	    aorth.lengthSquared(),
	    2 * oprime.dot(aorth),
	    oprime.lengthSquared() - sqr(r1 + r3),
	    t1, t2
    );
				
    if (nres >= 1)  {
	c31 = o + t1 * aorth;
    }

    if (nres >= 2) {
	c32 = o + t2 * aorth; 
    }
	
    return nres;     
}

bool intersect_range(double amin, double amax, double bmin, double bmax) {
    return !(bmax < amin || amax < bmin);
}

bool Circle::intersects(const BBox & bbox) const {
    double x = c.x, y = c.y;
    double xmin = bbox.Cmin.x, ymin = bbox.Cmin.y;
    double xmax = bbox.Cmax.x, ymax = bbox.Cmax.y;

    if (x + r < xmin) {
	return false;
    } else if (x < xmin) {
        double yo = sqrt(r * r - sqr(xmin - x));
        return intersect_range(y - yo, y + yo, ymin, ymax);
    } else if (x < xmax) {
        return intersect_range(y - r, y + r, ymin, ymax);
    } else if (x - r < xmax) {
        double yo = sqrt(r * r - sqr(x - xmax));
        return intersect_range(y - yo, y + yo, ymin, ymax);
    } else {
	return false; 
    }
}


bool Circle::intersects(const Circle  & other) const {
    return (c - other.c).lengthSquared() < sqr(other.r + r); 
}


