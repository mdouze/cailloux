

from geometry import Node, Circle, enumerate_intersecting_leaves, enumerate_pairs, \
    contact_3circle_inside, contact_3circle, norm2, norm, enumerate_leaves

import Cgeometry
import numpy as np




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
        # r3 = 0.25 * rs.rand() ** 3 + 0.01
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
        if len(c3s) == 0: 
            continue
        if True: 
            # pick the c3 that has lowest y
            c3s.sort(key=lambda x: x[1])
        else: 
            # pick the c3 that is closest to 0
            c3s.sort(key=lambda x: norm2(x))
        c3 = c3s[0]
        #print(f"{c3=:}")
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





