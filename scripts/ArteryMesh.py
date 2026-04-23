import gmsh

r = 2. # Inner radius (mm)
t = 0.2 # Thickness (mm)

nr = 10 # Elements through artery thickness
nt = 30 # Elements along arc

def main():

    gmsh.initialize()
    gmsh.model.add("QuarterArtery")

    # Define geometry points
    c = gmsh.model.occ.addPoint(0, 0, 0)
    p1 = gmsh.model.occ.addPoint(r, 0, 0)
    p2 = gmsh.model.occ.addPoint(0, r, 0)
    p3 = gmsh.model.occ.addPoint(r+t, 0, 0)
    p4 = gmsh.model.occ.addPoint(0, r+t, 0)

    # Define boundary from points
    inner_arc = gmsh.model.occ.addCircleArc(p1, c, p2)
    outer_arc = gmsh.model.occ.addCircleArc(p3, c, p4)
    horizontal = gmsh.model.occ.addLine(p1, p3)
    vertical = gmsh.model.occ.addLine(p2, p4)
    loop = gmsh.model.occ.addCurveLoop([horizontal, outer_arc, vertical, inner_arc])

    # Create surface
    surface = gmsh.model.occ.addPlaneSurface([loop])
    gmsh.model.occ.synchronize()

    # Make all elements linear quads
    gmsh.model.mesh.setTransfiniteCurve(inner_arc, nt+1)
    gmsh.model.mesh.setTransfiniteCurve(outer_arc, nt+1)
    gmsh.model.mesh.setTransfiniteCurve(horizontal, nr+1)
    gmsh.model.mesh.setTransfiniteCurve(vertical, nr+1)
    gmsh.model.mesh.setTransfiniteSurface(surface)
    gmsh.model.mesh.setRecombine(2, surface)

    # Add physical groups to interface with MFEM
    gmsh.model.addPhysicalGroup(1, [horizontal], tag=1, name="HorizontalLine")
    gmsh.model.addPhysicalGroup(1, [outer_arc], tag=2, name="OuterArc")
    gmsh.model.addPhysicalGroup(1, [vertical], tag=3, name="VerticalLine")
    gmsh.model.addPhysicalGroup(1, [inner_arc], tag=4, name="InnerArc")
    gmsh.model.addPhysicalGroup(2, [surface], tag=5, name="Surface")

    # Generate and write mesh
    gmsh.model.mesh.generate(2)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write("../meshes/QuarterArtery.msh")

    gmsh.finalize()

if __name__ == "__main__":
    main()
