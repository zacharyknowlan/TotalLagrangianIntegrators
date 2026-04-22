import subprocess

Lambda = 2.
mu = 1.

pressures = [(0.012 + 0.0005*i) for i in range(0,9)] # in MPa

CommandLineInput = ["../build/example_problems/HyperElasticArtery"]
CommandLineInput.extend(["--MeshFile", "../meshes/QuarterArtery.msh"])
CommandLineInput.extend(["--ResultFile", " "])
CommandLineInput.extend(["--lambda", str(Lambda)])
CommandLineInput.extend(["--mu", str(mu)])
CommandLineInput.extend(["--p", " "])

for i in range(0,len(pressures)):
    CommandLineInput[4] = str("../results/HyperElasticArtery_" + str(i) + ".vtk")
    CommandLineInput[10] = str(pressures[i])
    _ = subprocess.run(CommandLineInput)
