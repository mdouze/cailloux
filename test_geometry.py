
import unittest
import geometry
import numpy as np



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
            if geometry.circle_intersects_bbox(circle, root.bbox): 
                set1.add(str(circle))
        set2 = set()
        for circle in circles: 
            if geometry.circle_intersects_bbox(circle, root.bbox): 
                set2.add(str(circle))
        assert set1 == set2
    else: 
        check_consistent(root.child1, circles)
        check_consistent(root.child2, circles)



class TestKDTree(unittest.TestCase):


    def test_build(self):
        rs = np.random.RandomState(456)        
        circles = [
            (rs.rand(2), rs.rand(1) ** 5)
            for _ in range(50)
        ]

        root = geometry.Node((-2, -2, 2, 2), path="")
        for circle in circles:
            root.add_circle(circle)

        check_consistent(root, circles)

        
    def test_enumerate_pairs(self):     
        rs = np.random.RandomState(123)        
        circles = [
            (rs.rand(2), rs.rand(1) ** 5)
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
    
            

