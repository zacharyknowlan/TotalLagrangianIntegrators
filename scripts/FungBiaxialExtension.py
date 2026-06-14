import subprocess
from FungUniaxialExtension import a, A1, A2, A3, A4, A5, A6, E11
from pathlib import Path

ResultDirectory = "../results/BiaxialExtension/"

E22 = [((strain + 0.04)/1.79) for strain in E11]

def main():

    CommandLineInput = ["../build/problems/FBE"]
    CommandLineInput.extend(["--MeshFile", "../meshes/Square.msh"])
    CommandLineInput.extend(["--ResultFile", " "])
    CommandLineInput.extend(["-a", str(a)])
    CommandLineInput.extend(["-A1", str(A1)])
    CommandLineInput.extend(["-A2", str(A2)])
    CommandLineInput.extend(["-A3", str(A3)])
    CommandLineInput.extend(["-A4", str(A4)])
    CommandLineInput.extend(["-A5", str(A5)])
    CommandLineInput.extend(["-A6", str(A6)])
    CommandLineInput.extend(["-u_x", " "])
    CommandLineInput.extend(["-u_y", " "])
    CommandLineInput.extend(["-o", "1"])

    DOF_filepath = Path(ResultDirectory + "Fung_DOF.csv")
    DOF_filepath.unlink(missing_ok=True)

    solve_time_filepath = Path(ResultDirectory + "Fung_solve_time.csv")
    solve_time_filepath.unlink(missing_ok=True)

    free_energy_filepath = Path(ResultDirectory + "Fung_free_energy.csv")
    free_energy_filepath.unlink(missing_ok=True)

    processes = []
    for i in range(0, len(E11)):
        CommandLineInput[4] = str(ResultDirectory + "Fung_Biaxial_Result_" + str(i) + ".vtk")
        CommandLineInput[20] = str(E11[i])
        CommandLineInput[22] = str(E22[i])
        process = subprocess.Popen(CommandLineInput)
        processes.append(process)

    for process in processes:
        returncode = process.wait()
        assert int(returncode) == 0, "\033[31mProcess Failed\033[0m"

if __name__ == "__main__":
    main()
