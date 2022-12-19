import numpy as np


##########################################################
# 3 tangent circles 
##########################################################

dot = np.dot

def solve_2nd_degree(a, b, c): 
    delta = b ** 2 - 4 * a * c
    if delta < 0: 
        return []
    return (
        (-b + np.sqrt(delta)) / (2 * a),
        (-b - np.sqrt(delta)) / (2 * a)
    )

def find_origin(a, b): 
    """ find one point for the line defined by a.x = b"""
    ax, ay = a
    if np.abs(ax) > np.abs(ay): 
        o = np.array([b/ax, 0])
    else: 
        o = np.array([0, b/ay])
    # print(np.dot(o, a) - b)
    return o
    
def norm2(x): 
    return np.dot(x, x)

def norm(x): 
    return np.sqrt(norm2(x))

def contact_3cricle(
        c1, r1, 
        c2, r2,
        r3): 
    """ 
    Find c3 st. (c1, r1), (c2, r2) and (c3, r3) are tangent    
    """
    b = ((r1 + r3)**2 - (r2 + r3)**2 - norm2(c1) + norm2(c2)) / 2
    a = c2 - c1 
    # a.c3 = b
    aorth = np.array([a[1], -a[0]])
    
    o = find_origin(a, b)
        
    # c3 = o + t * aorth, find t
    oprime = o - c1
        
    ts = solve_2nd_degree(
        norm2(aorth), 
        2 * dot(oprime, aorth), 
        norm2(oprime) - (r1 + r3)**2
    )
    
    return [o + t * aorth for t in ts]

def contact_3circle_inside(
        c1, r1, 
        c2, r2, 
        r3): 
    """ 
    Find c3 st. (c1, r1), (c2, r2) and (c3, r3) are tangent
    but c1 and c3 are inside c2, r2.
    """
    b = ((r1 + r3)**2 - (r2 - r3)**2 - norm2(c1) + norm2(c2)) / 2
    a = c2 - c1 
    # a.c3 = b
    aorth = np.array([a[1], -a[0]])
    
    o = find_origin(a, b)
        
    # c3 = o + t * aorth, find t
    oprime = o - c1
        
    ts = solve_2nd_degree(
        norm2(aorth), 
        2 * dot(oprime, aorth), 
        norm2(oprime) - (r1 + r3)**2
    )
    
    return [o + t * aorth for t in ts]



##########################################################
# KDTree for circles
##########################################################



def intersect_range(amin, amax, bmin, bmax): 
    return not(bmax < amin or amax < bmin)


def circle_intersects_bbox(circle, bbox): 
    c, r = circle
    x, y = c
    xmin, ymin, xmax, ymax = bbox
    # horizontal intersect 
    if x + r < xmin: 
        return False
    elif x < xmin: 
        yo = np.sqrt(r ** 2 - (xmin - x) ** 2)
        return intersect_range(y - yo, y + yo, ymin, ymax)
    elif x < xmax: 
        return intersect_range(y - r, y + r, ymin, ymax)
    elif x - r < xmax: 
        yo = np.sqrt(r ** 2 - (x - xmax) ** 2)
        return intersect_range(y - yo, y + yo, ymin, ymax)
    else: 
        return False        


def filter_with_bbox(circles, bbox):                     
    return list(filter(
        lambda circle: circle_intersects_bbox((circle.c, circle.r), bbox), 
        circles
    ))

class Circle:

    def __init__(self, c, r):
        self.c = c
        self.r = r

    def intersects_bbox(self, bbox):
        return circle_intersects_bbox((self.c, self.r), bbox)


max_per_leaf = 10

class Node: 
    
    def __init__(self, bbox, path, circles=None): 
        self.bbox = bbox
        self.path = path
        self.is_leaf = True
        self.circles = circles if circles else []

    def surface(self): 
        xmin, ymin, xmax, ymax = self.bbox
        return (xmax - xmin) * (ymax - ymin)
        
    def split(self): 
        xmin, ymin, xmax, ymax = self.bbox
        if xmax - xmin > ymax - ymin: 
            # split vertically
            xmid = (xmin + xmax) / 2 
            bbox1 = (xmin, ymin, xmid, ymax)
            bbox2 = (xmid, ymin, xmax, ymax)
        else: 
            # split horizontally
            ymid = (ymin + ymax) / 2 
            bbox1 = (xmin, ymin, xmax, ymid)
            bbox2 = (xmin, ymid, xmax, ymax)
        circles = self.circles 
        del self.circles 
        
        self.child1 = Node(bbox1, self.path + "1", filter_with_bbox(circles, bbox1))
        self.child2 = Node(bbox2, self.path + "2", filter_with_bbox(circles, bbox2))
        
        self.is_leaf = False

    def add_circle(self, circle): 
        if not circle.intersects_bbox(self.bbox): 
            return 
        if self.is_leaf: 
            self.circles.append(circle)
            if len(self.circles) > max_per_leaf: 
                self.split()
        else: 
            self.child1.add_circle(circle)
            self.child2.add_circle(circle)

            
def enumerate_intersecting_leaves(root, circle):
    if not circle.intersects_bbox(root.bbox):
        return 
    if root.is_leaf:
        yield root
    else:
        enumerate_intersecting_leaves(root.child1, circle)
        enumerate_intersecting_leaves(root.child2, circle)    

            
##########################################################
# Explore KDTree 
##########################################################
            
def enumerate_leaves(root):
    """ enumerate all leaves in the KTree """
    if root.is_leaf: 
        yield root
    else: 
        yield from enumerate_leaves(root.child1)
        yield from enumerate_leaves(root.child2)
        


def range_distance(amin, amax, bmin, bmax): 
    if amax < bmin: 
        return bmin - amax
    if bmax < amin:
        return amin - bmax
    # the ranges intersect
    return 0

def nearest_corner(bbox, xy): 
    """ find nearest corner of the bbox w.r.t. xy"""
    x, y = xy
    xmin, ymin, xmax, ymax = bbox
    return (
        xmin if x < (xmin + xmax) / 2 else xmax, 
        ymin if y < (ymin + ymax) / 2 else ymax
    )

def bbox_distance(bbox1, bbox2): 
    """ minimum distance between contents of 2 bounding boxes """
    xmin1, ymin1, xmax1, ymax1 = bbox1
    xmin2, ymin2, xmax2, ymax2 = bbox2

    dy = range_distance(ymin1, ymax1, ymin2, ymax2)
    dx = range_distance(xmin1, xmax1, xmin2, xmax2)

    if dx != 0 and dy != 0: 
        # then it is a corner-wise distance 
        # find which corner... 
        x1, y1 = nearest_corner(bbox1, (xmin2, ymin2))
        x2, y2 = nearest_corner(bbox2, (xmin1, ymin1))
        return np.hypot(x1 - x2, y1 - y2)        
    else: 
        return max(dx, dy)



def enumerate_pairs_2(root1, root2, distance): 
    
    # check if the nodes are close enough
    d = bbox_distance(root1.bbox, root2.bbox)
    
    if d > distance: 
        return 
    #print(pref, end='\r', flush=True)
    #time.sleep(0.02)
    
    if root1.is_leaf and root2.is_leaf: 
        yield root1, root2
    elif root1.is_leaf: 
        yield from enumerate_pairs_2(root1, root2.child1, distance)
        yield from enumerate_pairs_2(root1, root2.child2, distance)
    elif root2.is_leaf:
        yield from enumerate_pairs_2(root1.child1, root2, distance)
        yield from enumerate_pairs_2(root1.child2, root2, distance)
    else: 
        yield from enumerate_pairs_2(root1.child1, root2.child1, distance)
        yield from enumerate_pairs_2(root1.child1, root2.child2, distance)
        yield from enumerate_pairs_2(root1.child2, root2.child1, distance)
        yield from enumerate_pairs_2(root1.child2, root2.child2, distance)
        
            
def enumerate_pairs(root, dis, pref=""):
    " yield pairs of leaves such that the min distance between leave bbox is < distance"

    if root.is_leaf: 
        return 
    else:         
        yield from enumerate_pairs_2(root.child1, root.child2, dis)
        yield from enumerate_pairs(root.child1, dis)
        yield from enumerate_pairs(root.child2, dis)
            
            

        
