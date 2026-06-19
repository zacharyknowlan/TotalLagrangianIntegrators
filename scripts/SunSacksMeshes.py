import gmsh
from numpy import sqrt

L = 0.025
cx, cy = L/2., L/2.
r = L/5.
NE = 1024
tol = 1e-6 # For locating things in bounding box (do not make smaller)

def main():

    gmsh.initialize()
    gmsh.model.add("Sun-Sacks-Square")

    square = gmsh.model.occ.addRectangle(0., 0., 0., L, L)
    p_left   = gmsh.model.occ.addPoint(0., cy, 0.)
    p_right  = gmsh.model.occ.addPoint(L, cy, 0.)
    p_bottom = gmsh.model.occ.addPoint(cx, 0., 0.)
    p_top    = gmsh.model.occ.addPoint(cx, L, 0.)
    hline = gmsh.model.occ.addLine(p_left, p_right)
    vline = gmsh.model.occ.addLine(p_bottom, p_top)
    gmsh.model.occ.fragment([(2, square)], [(1, hline), (1, vline)])
    gmsh.model.occ.synchronize()

    edges = gmsh.model.getEntities(1)
    surfaces = gmsh.model.getEntities(2)

    left, right, bottom, top, hsplit, vsplit = [], [], [], [], [], []

    for dim, tag in edges:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)

        if abs(xmin) < tol and abs(xmax) < tol:
            left.append(tag)
        elif abs(xmin - L) < tol and abs(xmax - L) < tol:
            right.append(tag)
        elif abs(ymin) < tol and abs(ymax) < tol:
            bottom.append(tag)
        elif abs(ymin - L) < tol and abs(ymax - L) < tol:
            top.append(tag)
        elif abs(ymin - cy) < tol and abs(ymax - cy) < tol:
            hsplit.append(tag)
        elif abs(xmin - cx) < tol and abs(xmax - cx) < tol:
            vsplit.append(tag)

    gmsh.model.addPhysicalGroup(1, bottom, tag=1, name="Bottom")
    gmsh.model.addPhysicalGroup(1, right, tag=2, name="Right")
    gmsh.model.addPhysicalGroup(1, top, tag=3, name="Top")
    gmsh.model.addPhysicalGroup(1, left, tag=4, name="Left")
    gmsh.model.addPhysicalGroup(1, hsplit, tag=5, name="hSplit")
    gmsh.model.addPhysicalGroup(1, vsplit, tag=6, name="vSplit")
    gmsh.model.addPhysicalGroup(2, [s[1] for s in surfaces], tag=8, name="Domain")

    # Set mesh size at all mesh points
    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), L/sqrt(NE))

    # Make the mesh quadrilateral
    for dim, tag in surfaces:
        gmsh.model.mesh.setTransfiniteSurface(tag)
    
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.ElementOrder", 1)
    gmsh.model.mesh.generate(dim=2)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write("../meshes/SunSacksSquare.msh")
    gmsh.clear()
    gmsh.finalize()

    gmsh.initialize()
    gmsh.model.add("Sun-Sacks-Square-With-Inclusion")

    square = gmsh.model.occ.addRectangle(0., 0., 0., L, L)
    circle = gmsh.model.occ.addDisk(cx, cy, 0., r, r)
    cut, _ = gmsh.model.occ.cut( [(2, square)], [(2, circle)], 
                                    removeObject=True, removeTool=True)
    p_left = gmsh.model.occ.addPoint(0., cx, 0.)
    p_left_inc = gmsh.model.occ.addPoint(cx-r, cx, 0.)
    p_right_inc = gmsh.model.occ.addPoint(cx+r, cx, 0.)
    p_right  = gmsh.model.occ.addPoint(L, cx, 0.)
    line_left = gmsh.model.occ.addLine(p_left, p_left_inc)
    line_right = gmsh.model.occ.addLine(p_right_inc, p_right)
    gmsh.model.occ.fragment(cut, [(1, line_left), (1, line_right)])
    gmsh.model.occ.synchronize()

    edges = gmsh.model.getEntities(1)
    surfaces = gmsh.model.getEntities(2)

    left, right, bottom, top, centerline, inc = [], [], [], [], [], []

    for dim, tag in edges:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)

        com = gmsh.model.occ.getCenterOfMass(dim, tag)
        dist = sqrt((com[0] - cx)**2 + (com[1] - cy)**2)

        if abs(xmin) < tol and abs(xmax) < tol:
            left.append(tag)
        elif abs(xmin - L) < tol and abs(xmax - L) < tol:
            right.append(tag)
        elif abs(ymin) < tol and abs(ymax) < tol:
            bottom.append(tag)
        elif abs(ymin - L) < tol and abs(ymax - L) < tol:
            top.append(tag)
        elif abs(ymin - cy) < tol and abs(ymax - cy) < tol and (xmax > (cx + r) or xmin < (cx - r)):
            centerline.append(tag)
        elif abs(dist - r) < tol:
            inc.append(tag)

    gmsh.model.addPhysicalGroup(1, bottom, tag=1, name="Bottom")
    gmsh.model.addPhysicalGroup(1, right, tag=2, name="Right")
    gmsh.model.addPhysicalGroup(1, top, tag=3, name="Top")
    gmsh.model.addPhysicalGroup(1, left, tag=4, name="Left")
    gmsh.model.addPhysicalGroup(1, centerline, tag=5, name="Centerline")
    gmsh.model.addPhysicalGroup(1, inc, tag=6, name="Inclusion")
    gmsh.model.addPhysicalGroup(2, [s[1] for s in surfaces], tag=7, name="Domain")

    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", L/(2.*sqrt(NE)))
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", L/sqrt(NE))
    gmsh.option.setNumber("Mesh.ElementOrder", 1)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.model.mesh.generate(dim=2)
    gmsh.write("../meshes/SunSacksSquareWithInclusion.msh")
    gmsh.clear()
    gmsh.finalize()

if __name__ == "__main__":
    main()
