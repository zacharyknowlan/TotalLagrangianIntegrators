import gmsh

l = 33 # Ligament length (mm)
w = 5.5 # Half of the ligament width (mm)

nx = 30 # Elements along length
ny = 5 # Elements along half-width

def main():

    gmsh.initialize()
    gmsh.model.add("HalfLigament")

    # Define geometry points
    p1 = gmsh.model.occ.addPoint(0, 0, 0)
    p2 = gmsh.model.occ.addPoint(l, 0, 0)
    p3 = gmsh.model.occ.addPoint(l, w, 0)
    p4 = gmsh.model.occ.addPoint(0, w, 0)

    # Define boundary from points
    bottom = gmsh.model.occ.addLine(p1, p2)
    right = gmsh.model.occ.addLine(p2, p3)
    top = gmsh.model.occ.addLine(p3, p4)
    left = gmsh.model.occ.addLine(p4, p1)
    loop = gmsh.model.occ.addCurveLoop([bottom, right, top, left])

    # Create surface
    surface = gmsh.model.occ.addPlaneSurface([loop])
    gmsh.model.occ.synchronize()

    # Make all elements linear quads
    gmsh.model.mesh.setTransfiniteCurve(bottom, nx+1)
    gmsh.model.mesh.setTransfiniteCurve(right, ny+1)
    gmsh.model.mesh.setTransfiniteCurve(top, nx+1)
    gmsh.model.mesh.setTransfiniteCurve(left, ny+1)
    gmsh.model.mesh.setTransfiniteSurface(surface, cornerTags=[p1, p2, p3, p4])
    gmsh.model.mesh.setRecombine(2, surface)

    # Add physical groups to interface with MFEM
    gmsh.model.addPhysicalGroup(1, [bottom], name="Bottom")
    gmsh.model.addPhysicalGroup(1, [right], name="Right")
    gmsh.model.addPhysicalGroup(1, [top], name="Top")
    gmsh.model.addPhysicalGroup(1, [left], name="Left")
    gmsh.model.addPhysicalGroup(2, [surface], name="Domain")

    # Generate and write mesh
    gmsh.model.mesh.generate(2)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write("../meshes/HalfLigament.msh")

    gmsh.finalize()


if __name__ == "__main__":
    main()
