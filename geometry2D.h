

#pragma once

#include <vector>

#include "cinder-lite/Vector.h"


using Vec2 = cinderlite::Vec2<double>;


int contact_3circle(Vec2 c1, double r1,
		    Vec2 c2, double r2,
		    double r3,
		    Vec2 &c31, Vec2 &c32,
		    bool inside_c2=false);

struct BBox {
    Vec2 Cmin, Cmax;
    double surface() const {
        Vec2 v = Cmax - Cmin; 
	return v.x * v.y;
    }
};

struct Shape2D {
    int id;
    Shape2D(int id = -1): id(id) {}
    
    virtual bool intersects(const BBox & bbox) const = 0;
    virtual Shape2D * clone() const = 0; 
    virtual ~Shape2D() {}
    
};
    

struct Circle: Shape2D {
    Vec2 c;
    double r;

    Circle(Vec2 c, double r, int id=-1):
	Shape2D(id), c(c), r(r) {}

    Circle(double x, double y, double r, int id=-1):
	Shape2D(id), c(x, y), r(r) {}

    bool intersects(const BBox & bbox) const final;

    bool intersects(const Circle & other) const; 

    Shape2D * clone() const {return new Circle(*this); }
        
    virtual ~Circle() {}
  
};


/********************************************
 * KDTree
 ********************************************/

struct Node {
    BBox bbox;
    
    std::vector<Shape2D *> shapes;
    Node * child1 = nullptr;
    Node * child2 = nullptr; 
    

    bool is_leaf() const {
	return child1 == nullptr;
    }

    Node(const BBox & bbox): bbox(bbox) {}

    void split(); 
    
    void add_shape(Shape2D *shape, int max_per_leaf); 

};


struct KDTree {

    Node * root;
    int max_per_leaf = 10; 
    std::vector<Shape2D*> flat_shapes; 
    
    KDTree(const BBox & bbox) {
	root = new Node(bbox);	
    }

    ~KDTree();

    void add_shape(const Shape2D &shape); 
    

};
 
