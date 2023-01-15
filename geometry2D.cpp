
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


/********************************************
 * KDTree
 ********************************************/

namespace {

void rec_delete_nodes(Node *root) {
    if (!root->is_leaf()) {
	rec_delete_nodes(root->child1);
	rec_delete_nodes(root->child2);	
    }
    delete root;
}
    /*
bool rec_check_consistent(Node * root, const std::vector<Shape2D*> flat) {
    if (node->is_leaf()) {
	std::unordered_set<Shape2D*> sub_shapes;
	for (auto x : flat) {
	    if (x.intersects(bbox)) {
		sub_shapes.add(x); 
	    }	   
	}

	for (auto y : shapes) {
	    if (sub_sahpes
	    
	


    } else {
	rec_check_consistent(root->child1, flat); 
	rec_check_consistent(root->child2, flat); 
    }
}
    */
    
};

KDTree::~KDTree() {    
    rec_delete_nodes(root);
    for (auto i: flat_shapes) {
	delete i;
    }
}

namespace {

double range_distance(double amin, double amax, double bmin, double bmax) {
    if (amax < bmin) 
        return bmin - amax;
    if (bmax < amin) 
        return amin - bmax;
    // the ranges intersect
    return 0;
}

}

Vec2 BBox::nearest_corner(Vec2 p) const {
    return Vec2(p.x < (Cmin.x  + Cmax.x) / 2 ? Cmin.x : Cmax.x, 
		p.y < (Cmin.y  + Cmax.y) / 2 ? Cmin.y : Cmax.y);
}

// shortest distance between the two
double BBox::distance(const BBox & o) const {
    double dy = range_distance(Cmin.y, Cmax.y, o.Cmin.y, o.Cmax.y); 
    double dx = range_distance(Cmin.x, Cmax.x, o.Cmin.x, o.Cmax.x); 
    if (dx != 0 && dy != 0) {
	// it is a corner-wise distance
	// which corner?
	Vec2 c1 = nearest_corner(o.Cmin);
	Vec2 c2 = o.nearest_corner(Cmin);
	return c1.distance(c2);	
    }
    return std::max(dx, dy); 
}


void Node::split() {
    Vec2 v = bbox.Cmax - bbox.Cmin;
    BBox bbox1, bbox2; 
    if (v.x > v.y) { // split vertically
	double xmid = (bbox.Cmax.x + bbox.Cmin.x) / 2;
	bbox1 = BBox(bbox.Cmin, Vec2(xmid, bbox.Cmax.y));
	bbox2 = BBox(Vec2(xmid, bbox.Cmin.y), bbox.Cmax);		            
    } else { // split horizontally
	double ymid = (bbox.Cmax.y + bbox.Cmin.y) / 2;
	bbox1 = BBox(bbox.Cmin, Vec2(bbox.Cmax.x, ymid));
	bbox2 = BBox(Vec2(bbox.Cmin.x, ymid), bbox.Cmax);		            	
    }
    child1 = new Node(bbox1);
    child2 = new Node(bbox2);

    for(auto shape: shapes) {
	if (shape->intersects(bbox1)) {
	    child1->shapes.push_back(shape); 
	}
	if (shape->intersects(bbox2)) {
	    child2->shapes.push_back(shape); 
	}	
    }
    shapes.resize(0);       
}

void KDTree::add_shape(const Shape2D &shape_in) {
    Shape2D *shape = shape_in.clone();    
    flat_shapes.push_back(shape);
    root->add_shape(shape, max_per_leaf);
}
    
void Node::add_shape(Shape2D *shape, int max_per_leaf) {
    if (!shape->intersects(bbox)) {
	return; 
    }
    if (is_leaf()) {
	shapes.push_back(shape);
	if (shapes.size() > max_per_leaf) {
	    split(); 
	}
    } else {
	child1->add_shape(shape, max_per_leaf);
	child2->add_shape(shape, max_per_leaf);	
    }   
	
}
