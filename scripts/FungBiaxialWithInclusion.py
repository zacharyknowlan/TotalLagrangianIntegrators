import subprocess
from FungBiaxialExtension import a, A1, A2, A3, A4, A5, A6, u_x, u_y, increments

def main():

    CommandLineInput = ["../build/problems/FBE"]
    CommandLineInput.extend(["--MeshFile", "../meshes/SunSacksSquareWithInclusion.msh"])
    CommandLineInput.extend(["--ResultFile", "../results/BiaxialExtension/SquareWithInclusion.vtk"])
    CommandLineInput.extend(["-a", str(a)])
    CommandLineInput.extend(["-A1", str(A1)])
    CommandLineInput.extend(["-A2", str(A2)])
    CommandLineInput.extend(["-A3", str(A3)])
    CommandLineInput.extend(["-A4", str(A4)])
    CommandLineInput.extend(["-A5", str(A5)])
    CommandLineInput.extend(["-A6", str(A6)])
    CommandLineInput.extend(["-u_x", str(u_x[0])])
    CommandLineInput.extend(["-u_y", str(u_y[0])])
    CommandLineInput.extend(["-i", str(increments[0])])

    _ = subprocess.run(CommandLineInput)

if __name__ == "__main__":
    main()
