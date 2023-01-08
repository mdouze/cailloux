import unittest
import geometry
import fields
import Cgeometry
import numpy as np


class TestCircle(unittest.TestCase):

    def test2(self):
        c1 = np.array([0.4, 0.3]); r1 = 0.10
        c2 = np.array([0.7, 0.5]); r2 = 0.25
        
        r3 = 0.15
        c3s = geometry.contact_3circle(c1, r1, c2, r2, r3)

        for c3 in c3s:
            np.testing.assert_almost_equal(geometry.norm(c1 - c3), r1 + r3)
            np.testing.assert_almost_equal(geometry.norm(c2 - c3), r2 + r3)


    def test3C(self):
        def to_vec2(x):
            return Cgeometry.Vec2(x[0], x[1])
        
        c1 = np.array([0.4, 0.3]); r1 = 0.10
        c2 = np.array([0.7, 0.5]); r2 = 0.25
        
        r3 = 0.15
        c31 = Cgeometry.Vec2()
        c32 = Cgeometry.Vec2()
        
        nv = Cgeometry.contact_3circle(to_vec2(c1), r1, to_vec2(c2), r2, r3, c31, c32)
        self.assertEqual(nv, 2)
        c3s = [
            np.array([c31.x, c31.y]),
            np.array([c32.x, c32.y])            
        ]
        for c3 in c3s:
            np.testing.assert_almost_equal(geometry.norm(c1 - c3), r1 + r3)
            np.testing.assert_almost_equal(geometry.norm(c2 - c3), r2 + r3)

              

def random_bbox(rs): 
    xmin, xmax, ymin, ymax = rs.rand(4) 
    if xmax < xmin: 
        xmin, xmax = xmax, xmin
    if ymax < ymin: 
        ymin, ymax = ymax, ymin
    return (xmin, ymin, xmax + 0.1, ymax + 0.1)        


def plot_bbox(bbox): 
    xmin, ymin, xmax, ymax = bbox
    pyplot.gca().add_patch(
        pyplot.Rectangle((xmin, ymin), xmax - xmin, ymax - ymin, color=np.random.rand(3))    
    )

class TestBBoxDistance(unittest.TestCase):

    def test_symmetry(self):
        rs = np.random.RandomState(456)        
        for a in range(50): 
            bbox1 = random_bbox(rs)
            bbox2 = random_bbox(rs)
            d12 = geometry.bbox_distance(bbox1, bbox2)
            d21 = geometry.bbox_distance(bbox2, bbox1)
            # print(bbox1, bbox2, d12, d21)
            np.testing.assert_almost_equal(d12, d21)

def check_consistent(root, circles):
    
    if root.is_leaf:
        set1 = set()        
        for circle in root.circles: 
            if circle.intersects_bbox(root.bbox): 
                set1.add(str(circle))
        set2 = set()
        for circle in circles: 
            if circle.intersects_bbox(root.bbox): 
                set2.add(str(circle))
        assert set1 == set2
    else: 
        check_consistent(root.child1, circles)
        check_consistent(root.child2, circles)



class TestKDTree(unittest.TestCase):


    def test_build(self):
        rs = np.random.RandomState(456)        
        circles = [
            geometry.Circle(rs.rand(2), rs.rand(1) ** 5)
            for _ in range(50)
        ]

        root = geometry.Node((-2, -2, 2, 2), path="")
        for circle in circles:
            root.add_circle(circle)

        check_consistent(root, circles)

        
    def test_enumerate_pairs(self):     
        rs = np.random.RandomState(123)        
        circles = [
            geometry.Circle(rs.rand(2), rs.rand(1) ** 5)
            for _ in range(50)
        ]

        root = geometry.Node((-2, -2, 2, 2), path="")
        dis = 0.012

        ref = {}

        print("REF")
        for leaf1 in geometry.enumerate_leaves(root):
            for leaf2 in geometry.enumerate_leaves(root):       
                if leaf1 != leaf2 and bbox_distance(leaf1.bbox, leaf2.bbox) < dis: 
                    x = [str(leaf1), str(leaf2)]
                    x.sort()
                    ref[tuple(x)] = (leaf1, leaf2)

        new = {}
        print("NEW")

        for leaf1, leaf2 in geometry.enumerate_pairs(root, dis): 
            x = [str(leaf1), str(leaf2)]
            x.sort()
            new[tuple(x)] = (leaf1, leaf2)

        self.assertEqual(ref.keys(), new.keys())
    
            
class TestField(unittest.TestCase):

    def test_kdtree(self):
        nc = 25
        rs = np.random.RandomState(345)

        radiuses = [0.5 * rs.rand() ** 3 for _ in range(nc)]
        # radiuses.sort(reverse=True)
        r0 = radiuses[0]
        
        circles_ref = fields.generate_circles_gravity(
            np.array([0, 0]), 1,
            np.array([0, -1 + r0]), r0,
            radiuses[1:]
        )

        circles_new = fields.generate_circles_gravity_kdtree(
            np.array([0, 0]), 1,
            np.array([0, -1 + r0]), r0,
            radiuses[1:]
        )

        def rd(x): 
            m = 1e6
            return np.floor(x * m) / m

        circles_ref_s = set((rd(c[0]), rd(c[1]), rd(r)) for (c, r) in circles_ref)
        circles_new_s = set((rd(cir.c[0]), rd(cir.c[1]), rd(cir.r)) for cir in circles_new)
        self.assertEqual(circles_ref_s, circles_new_s)

    def test_C(self):
        rs = np.random.RandomState(345)

        # times = OrderedDict()
        nc = 5

        radiuses = [0.5 * rs.rand() ** 3 for _ in range(nc)]
        # radiuses.sort(reverse=True)
        r0 = radiuses[0]

        circles_ref = fields.generate_circles_gravity(
            np.array([0, 0]), 1,
            np.array([0, -1 + r0]), r0,
            radiuses[1:]
        )
        print()

        circles_new = fields.generate_circles_gravity_C(
            np.array([0, 0]), 1,
            np.array([0, -1 + r0]), r0,
            radiuses[1:]
        )

        def rd(x): 
            m = 1e6
            return np.floor(x * m) / m

        circles_ref_s = set((rd(c[0]), rd(c[1]), rd(r)) for (c, r) in circles_ref)
        circles_new_s = set((rd(cir.c.x), rd(cir.c.y), rd(cir.r)) for cir in circles_new)
        self.assertEqual(circles_ref_s, circles_new_s)


