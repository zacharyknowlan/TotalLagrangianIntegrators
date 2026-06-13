import subprocess
from pathlib import Path

a = 5.12
A1 = 60.12
A2 = 86.34
A3 = 2.00
A4 = 203.16
A5 = 43.05
A6 = 42.14
strains = [0.01*i for i in range(1,21)]

def main():

    CommandLineInput = ["../build/problems/FUE"]
    CommandLineInput.extend(["--MeshFile", "../meshes/Square.msh"])
    CommandLineInput.extend(["--ResultFile", " "])
    CommandLineInput.extend(["-a", str(a)])
    CommandLineInput.extend(["-A1", str(A1)])
    CommandLineInput.extend(["-A2", str(A2)])
    CommandLineInput.extend(["-A3", str(A3)])
    CommandLineInput.extend(["-A4", str(A4)])
    CommandLineInput.extend(["-A5", str(A5)])
    CommandLineInput.extend(["-A6", str(A6)])
    CommandLineInput.extend(["-u_x", str(strains[0])])

    DOF_filepath = Path("../results/DOF.csv")
    DOF_filepath.unlink(missing_ok=True)

    solve_time_filepath = Path("../results/solve_time.csv")
    solve_time_filepath.unlink(missing_ok=True)

    free_energy_filepath = Path("../results/free_energy.csv")
    free_energy_filepath.unlink(missing_ok=True)

    processes = []
    for i in range(0, len(strains)):
        CommandLineInput[4] = str("../results/Fung_Uniaxial_Extension_Result_" + str(i) + ".vtk")
        CommandLineInput[20] = str(strains[i])
        process = subprocess.Popen(CommandLineInput)
        processes.append(process)

    for process in processes:
        returncode = process.wait()
        assert int(returncode) == 0, "\033[31mProcess Failed\033[0m"

if __name__ == "__main__":
    main()
