

from geometry import Node, Circle, enumerate_intersecting_leaves, enumerate_pairs, \
    contact_3circle_inside, contact_3circle, norm2, norm, enumerate_leaves

import Cgeometry
import numpy as np



######################################################################
# Reference implementation
######################################################################



def generate_circles_gravity(c0, r0, c1, r1, radiuses):
    """
    Make a field with a set of circles inside (c0, r0)
    The first circle is (c1, r1)
    The other ones have radiuses in the list of radiuses, with priority to the ones
    with lowest y. 
    """
    
    circles = [(c1, r1)]

    tot1 = tot2 = 0
    for i, r3 in enumerate(radiuses): 
        nc = len(circles)
        c3s = []
        
        # check contact with great circle 
        for k in range(nc): 
            c2, r2 = circles[k]
            for c3 in contact_3circle_inside(c2, r2, c0, r0, r3):                
                for l in range(nc): 
                    tot1 += 1
                    if l == k: 
                        continue
                    cl, rl = circles[l]
                    if norm2(c3 - cl) < (rl + r3) ** 2: 
                        break 
                else: 
                    c3s.append(c3)
                    tot2 += 1

        for j in range(nc):
            c1, r1 = circles[j]
            for k in range(j + 1, nc): 
                c2, r2 = circles[k]
                for c3 in contact_3circle(c1, r1, c2, r2, r3): 
                    # check if intersection with other circles 
                    for l in range(nc): 
                        tot1 += 1
                        if l == k or l == j: 
                            continue
                        cl, rl = circles[l]
                        if norm2(c3 - cl) < (rl + r3) ** 2:
                            break 
                    else: 
                        if norm(c0 - c3) + r3 < r0: 
                            c3s.append(c3)
                            tot2 += 1
                        
                
        # impossible to put this circle in the field
        # print(f"{i=:} REF c3s=", c3s)
        if len(c3s) == 0: 
            continue
        if True: 
            # pick the c3 that has lowest y
            c3s.sort(key=lambda x: x[1])
        else: 
            # pick the c3 that is closest to 0
            c3s.sort(key=lambda x: norm2(x))
        c3 = c3s[0]
        # print("REF c3=", c3)
        circles.append((c3, r3))
        print(f"{i=:} nb circles: {len(circles)} "
              f"nb c3: {len(c3s)} {tot1=:} {tot2=:}", end="\r", flush=True)

    return circles



######################################################################
# Construction with kdtree 
######################################################################



def circles_intersect(cir1, cir2): 
    return norm(cir1.c - cir2.c) < cir1.r + cir2.r

def any_intersection(root, circle, exclude=None): 
    if exclude is None: 
        seen = set()
    else: 
        seen = set(exclude)
    for leaf in enumerate_intersecting_leaves(root, circle):
        for circle2 in leaf.circles: 
            # print(f"KDTREE {leaf.path} {circle} check {circle2}")

            if circle2 in seen: 
                continue
            seen.add(circle2)
            if circles_intersect(circle, circle2): 
                # print("    INTER")
                return True
    return False

def any_intersection_flat(circles, circle, exclude=None): 
    if exclude is None: 
        seen = set()
    else: 
        seen = set(exclude)
    for circle2 in circles: 
        if circle2 in seen: 
            continue
        # print(f"FLAT {circle} check {circle2}")
        if circles_intersect(circle, circle2): 
            # print("    INTER")
            return True
    return False


def generate_circles_gravity_kdtree(c0, r0, c1, r1, radiuses): 
    circles = [Circle(c1, r1, name=-1)]

    root = Node(
        (c0[0] - r0, c0[1] - r0, c0[0] + r0, c0[1] + r0),
        ""
    )
    root.add_circle(circles[0])
    
    tot1 = tot2 = 0
    for rno, r3 in enumerate(radiuses): 
        nc = len(circles)
        c3s = []
        
        # check contact with great circle 
        for cir2 in circles: 
            for c3 in contact_3circle_inside(cir2.c, cir2.r, c0, r0, r3): 
                if not any_intersection(root, Circle(c3, r3), exclude=[cir2]):
                    c3s.append(c3)
                    tot2 += 1
        
        seen = set()
        
        def handle_circle_pair(cir1, cir2): 
            if (cir1, cir2) in seen or (cir2, cir1) in seen: 
                return 
            seen.add((cir1, cir2))
            for c3 in contact_3circle(cir1.c, cir1.r, cir2.c, cir2.r, r3): 
                # check if intersection with other circles 
                if not any_intersection(
                    root, Circle(c3, r3), 
                    exclude=[cir1, cir2]):
                    if norm(c0 - c3) + r3 < r0: 
                        c3s.append(c3)

        if True: 

            # handle pairs of circles within a leaf 
            for leaf in enumerate_leaves(root): 
                nc = len(leaf.circles)
                for i in range(nc): 
                    for j in range(i + 1, nc):
                        handle_circle_pair(leaf.circles[i], leaf.circles[j]) 
                        
            # handle pairs of circles within 2 leaves
            for leaf1, leaf2 in enumerate_pairs(root, 2 * r3): 
                for ci in leaf1.circles: 
                    for cj in leaf2.circles: 
                        if ci != cj: 
                            handle_circle_pair(ci, cj)
        else:
            # brute force
            nc = len(circles)
            for i in range(nc): 
                for j in range(i + 1, nc):
                    handle_circle_pair(circles[i], circles[j]) 


        if len(c3s) == 0: 
            continue
        # pick the c3 that has lowest y
        c3s.sort(key=lambda x: x[1])
        c3 = c3s[0]
        cir3 = Circle(c3, r3, name=rno)
        circles.append(cir3)
        root.add_circle(cir3)
        print(f"iter {rno} nb circles: {len(circles)} "
              f"nb c3: {len(c3s)} {tot1=:} {tot2=:}", end="\r", flush=True)
    print()
    # root.display()
    return circles


######################################################################
# Construction with C version
######################################################################

def contact_3circle_C(cir2, cir0, r3, inside_c2=False): 
    c31 = Cgeometry.Vec2()
    c32 = Cgeometry.Vec2()
    nres = Cgeometry.contact_3circle(cir2.c, cir2.r, cir0.c, cir0.r, r3, c31, c32, inside_c2)
    if nres == 0: 
        return []
    assert nres == 2
    return [c31, c32]


def intersects_any_circle_C(cir, circles, exclude): 
    for cir1 in circles: 
        if cir1.id in exclude: 
            continue
        if cir1.intersects(cir): 
            return True
    return False


def make_circle_C(c, r, id=-1): 
    if type(c) == np.ndarray: 
        c = Cgeometry.Vec2(float(c[0]), float(c[1]))
    return Cgeometry.Circle(c.x, c.y, r, id)


def generate_circles_gravity_C(c0, r0, c1, r1, radiuses): 
    cir0 = make_circle_C(c0, r0)
    circles = [make_circle_C(c1, r1, -1)]
    
    tot1 = tot2 = 0
    for i, r3 in enumerate(radiuses): 
        nc = len(circles)
        c3s = []
        
        # check contact with great circle 
        for k in range(nc): 
            cir2 = circles[k]
            for c3 in contact_3circle_C(cir2, cir0, r3, inside_c2=True): 
                cir3 = make_circle_C(c3, r3)
                if not intersects_any_circle_C(cir3, circles, exclude=[cir2.id]): 
                    c3s.append(c3)
                    tot2 += 1

        for j in range(nc):
            cir1 = circles[j]
            for k in range(j + 1, nc): 
                cir2 = circles[k]
                for c3 in contact_3circle_C(cir1, cir2, r3): 
                    cir3 = make_circle_C(c3, r3)
                    if not intersects_any_circle_C(cir3, circles, exclude=[cir1.id, cir2.id]):
                        if cir0.c.distance(c3) + r3 < r0: 
                            c3s.append(c3)
                            tot2 += 1
                        
        # print("c3s=", c3s)     
        if len(c3s) == 0: 
            continue
        # pick the c3 that has lowest y
        c3s.sort(key=lambda c3: c3.y)
        c3 = c3s[0]
        # print(f"{c3=:}")
        circles.append(make_circle_C(c3, r3, i))
        print(f"{i=:} nb circles: {len(circles)} "
              f"nb c3: {len(c3s)} {tot1=:} {tot2=:}", end="\r", flush=True)

    return circles



def intersects_any_circle_C_kdtree(kdtree, cir, exclude):
    seen = set(exclude)
    for node in Cgeometry.IntersectingLeavesIterator(kdtree, cir):
        for shape in Cgeometry.ShapeVectorIterator(node.shapes):
            if shape.id in seen:
                continue
            seen.add(shape.id)
            cir2 = Cgeometry.downcast_Circle(shape)
            if cir2.intersects(cir):
                return True
    return False
    

    

def generate_circles_gravity_C_kdtree(c0, r0, c1, r1, radiuses): 
    cir0 = make_circle_C(c0, r0)
    circles = [make_circle_C(c1, r1, -1)]
    kdtree = Cgeometry.KDTree(Cgeometry.BBox(
        cir0.c.x - cir0.r, cir0.c.y - cir0.r,
        cir0.c.x + cir0.r, cir0.c.y + cir0.r
    ))
    kdtree.add_shape(circles[0])        
    
    tot1 = tot2 = 0
    for it, r3 in enumerate(radiuses): 
        nc = len(circles)
        c3s = []
        
        # check contact with great circle 
        for k in range(nc): 
            cir2 = circles[k]
            for c3 in contact_3circle_C(cir2, cir0, r3, inside_c2=True): 
                cir3 = make_circle_C(c3, r3)
                if not intersects_any_circle_C_kdtree(kdtree, cir3, exclude=[cir2.id]): 
                    c3s.append(c3)
                    tot2 += 1

        seen = set()

        def handle_circle_pair(cir1, cir2):            
            k = tuple(sorted([cir1.id, cir2.id]))
            # print("TRY", k)
            if k in seen:
                return
            seen.add(k)                
            for c3 in contact_3circle_C(cir1, cir2, r3):
                #print("TRY", c3)
                cir3 = make_circle_C(c3, r3)
                if not intersects_any_circle_C_kdtree(kdtree, cir3, exclude=[cir1.id, cir2.id]):
                    if cir0.c.distance(c3) + r3 < r0:
                        c3s.append(c3)

        # handle pairs of circles within a leaf 
        for leaf in Cgeometry.LeafIterator(kdtree):
            leaf_circles = [
                Cgeometry.downcast_Circle(shape) for shape in Cgeometry.ShapeVectorIterator(leaf.shapes)
            ]
            # print(leaf_circles)
            nc = len(leaf_circles)
            for i in range(nc): 
                for j in range(i + 1, nc):
                    handle_circle_pair(leaf_circles[i], leaf_circles[j]) 

        # handle pairs of circles between 2 leaves        
        for twoleaves in Cgeometry.TwoLeavesIterator(kdtree, r3):
            for cir1 in Cgeometry.ShapeVectorIterator(twoleaves.node1.shapes):
                cir1 = Cgeometry.downcast_Circle(cir1)
                for cir2 in Cgeometry.ShapeVectorIterator(twoleaves.node2.shapes):
                    cir2 = Cgeometry.downcast_Circle(cir2)
                    handle_circle_pair(cir1, cir2)
                        
        # print(f"{it=:} NEW c3s=", c3s)     
        if len(c3s) == 0: 
            continue
        # pick the c3 that has lowest y
        c3s.sort(key=lambda c3: c3.y)
        c3 = c3s[0]
        # print(f"NEW {c3=:}")
        circles.append(make_circle_C(c3, r3, it))
        kdtree.add_shape(make_circle_C(c3, r3, it))
        print(f"{i=:} nb circles: {len(circles)} "
              f"nb c3: {len(c3s)} {tot1=:} {tot2=:}", end="\r", flush=True)

    return circles


