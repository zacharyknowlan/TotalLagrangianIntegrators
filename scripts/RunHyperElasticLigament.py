import subprocess
from LigamentMesh import l

Lambda = 2.
mu = 1.

strains = [0.01*i for i in range(0,7)]
displacements = [l*strain for strain in strains]

CommandLineInput = ["../build/example_problems/HyperElasticLigament"]
CommandLineInput.extend(["--MeshFile", "../meshes/HalfLigament.msh"])
CommandLineInput.extend(["--ResultFile", " "])
CommandLineInput.extend(["--lambda", str(Lambda)])
CommandLineInput.extend(["--mu", str(mu)])
CommandLineInput.extend(["--u", " "])

for i in range(0,len(displacements)):
    CommandLineInput[4] = str("../results/HyperElasticLigament_" + str(i) + ".vtk")
    CommandLineInput[10] = str(displacements[i])
    _ = subprocess.run(CommandLineInput)
