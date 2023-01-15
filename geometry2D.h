

#pragma once

#include <vector>
#include <cassert>

#include "cinder-lite/Vector.h"


using Vec2 = cinderlite::Vec2<double>;


int contact_3circle(Vec2 c1, double r1,
		    Vec2 c2, double r2,
		    double r3,
		    Vec2 &c31, Vec2 &c32,
		    bool inside_c2=false);

struct BBox {
    Vec2 Cmin, Cmax;

    BBox(double xmin = 0, double ymin = 0,
	 double xmax = 0, double ymax = 0):
	Cmin(xmin, ymin), Cmax(xmax, ymax) {}

    BBox(Vec2 Cmin, Vec2 Cmax): Cmin(Cmin), Cmax(Cmax) {}
   
    double surface() const {
        Vec2 v = Cmax - Cmin; 
	return v.x * v.y;
    }

    // shortest distance between the two
    double distance(const BBox & other) const;

    // nearest bbox corner w.r.t. p
    Vec2 nearest_corner(Vec2 p) const; 
	
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
    
    ///    bool check_consistent() const; 
};
 
/********************************************
 * Iterators
 ********************************************/

struct ShapeVectorIterator {
    const std::vector<Shape2D*> &shapes;     
    size_t i = 0; 
    
    ShapeVectorIterator(const std::vector<Shape2D*> &shapes):
	shapes(shapes) {}

    bool has_next() {
	return i < shapes.size();
    }

    const Shape2D * next() {
	if(!has_next()) return nullptr; 
	return shapes[i++]; 
    }   

};


struct LeafIterator {
    std::vector<const Node*> stack;

    LeafIterator(const KDTree & tree) {
	// reach leftmost leaf
	stack.push_back(tree.root);
    }

    bool has_next() {
	return !stack.empty(); 
    }

    const Node * next() {
	if (!has_next()) return nullptr;
	// find next leaf
	while(!stack.back()->is_leaf()) {
	    const Node * n = stack.back();
	    stack.pop_back(); 
	    stack.push_back(n->child2);
	    stack.push_back(n->child1);	    
	}
	// we have our leaf
	const Node *ret = stack.back();
	stack.pop_back();
	return ret;    
    }       
    
};


struct IntersectingLeavesIterator {
    std::vector<const Node*> stack;
    const Shape2D &shape;

    IntersectingLeavesIterator(const KDTree & tree, const Shape2D &shape):shape(shape) {
	// reach leftmost leaf
	if (shape.intersects(tree.root->bbox))
	    stack.push_back(tree.root);
    }

    bool has_next() {
	return !stack.empty(); 
    }

    const Node * next() {
	if (!has_next()) return nullptr;
	// find next leaf
	while(!stack.back()->is_leaf()) {
	    const Node * n = stack.back();
	    stack.pop_back();
	    if (shape.intersects(n->child2->bbox))
		stack.push_back(n->child2);
	    if (shape.intersects(n->child1->bbox))
		stack.push_back(n->child1);
	    if (!has_next()) return nullptr;
	}
	// we have our leaf
	const Node *ret = stack.back();
	stack.pop_back();
	return ret;    
    }       
    
};


struct TwoNodes {
    const Node* node1;
    const Node *node2;

    TwoNodes(const Node* node1 = nullptr, 
	     const Node *node2 = nullptr):
	node1(node1), node2(node2) {}

};

struct TwoLeavesIterator {
    double distance;
    std::vector<TwoNodes> stack;

    TwoLeavesIterator(const KDTree & tree, double distance):
	distance(distance) {
	push_single(tree.root); 
	while(!is_returnable()) {
	    step_stack();
	}
    }


    bool is_returnable() const {
	if (stack.empty())
	    return true;
	TwoNodes top = stack.back();
	if(!top.node2)
	    return false;
	return top.node1->is_leaf() && top.node2->is_leaf(); 	    
    }

    void push_single(const Node *n) {
	if (!n->is_leaf()) {
	    stack.push_back(TwoNodes(n));
	}
    }
		     
    void push_pair(const Node *node1, const Node *node2) {
	if (node1->bbox.distance(node2->bbox) < distance) {
	    stack.push_back(TwoNodes(node1, node2)); 
	}
    }
    
    void step_stack() {
	assert(!is_returnable()); 
	const TwoNodes & top = stack.back();
	stack.pop_back();
	const Node *root1 = top.node1;	    
	const Node *root2 = top.node2;
	if(!top.node2) {
	    if(root1->is_leaf()) {
		// nothing
	    } else {
		push_pair(root1->child1, root1->child2);
		push_single(root1->child1);
		push_single(root1->child2);		
	    }	    
	} else {
	    if (root1->is_leaf()) {
		assert(!root2->is_leaf());
		push_pair(root1, root2->child1);
		push_pair(root1, root2->child2);
	    } else if(root2->is_leaf()) {
		push_pair(root1->child1, root2);
		push_pair(root1->child2, root2);		
	    } else {
		push_pair(root1->child1, root2->child1);
		push_pair(root1->child2, root2->child1);		
		push_pair(root1->child1, root2->child2);
		push_pair(root1->child2, root2->child2);		
	    }  	
	}	
    }
    
    bool has_next() {
	return !stack.empty(); 
    }

    TwoNodes next() {
	assert(is_returnable()); 
	if (!has_next()) {
	    return TwoNodes();
	}
	TwoNodes top = stack.back();
	stack.pop_back();
	while(!is_returnable()) {
	    step_stack();
	}
	return top;
    }    
    

};
