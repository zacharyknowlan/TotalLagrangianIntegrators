import subprocess
from pathlib import Path

ResultDirectory = "../results/BiaxialExtension/"

a = 4.80e3
A1 = 64.99
A2 = 43.55
A3 = 1.67
A4 = 2.47
A5 = -5.25
A6 = 2.46

T = [50e3*i for i in range(1,21)]
increments = [int(10*i) for i in range(1,len(T)+1)]

def main():

    CommandLineInput = ["../build/problems/FBE"]
    CommandLineInput.extend(["--MeshFile", "../meshes/SunSacksSquare.msh"])
    CommandLineInput.extend(["--ResultFile", " "])
    CommandLineInput.extend(["-a", str(a)])
    CommandLineInput.extend(["-A1", str(A1)])
    CommandLineInput.extend(["-A2", str(A2)])
    CommandLineInput.extend(["-A3", str(A3)])
    CommandLineInput.extend(["-A4", str(A4)])
    CommandLineInput.extend(["-A5", str(A5)])
    CommandLineInput.extend(["-A6", str(A6)])
    CommandLineInput.extend(["-T", " "])
    CommandLineInput.extend(["-i", " "])
    CommandLineInput.extend(["-o", "1"])

    DOF_filepath = Path(ResultDirectory + "Fung_DOF.csv")
    DOF_filepath.unlink(missing_ok=True)

    solve_time_filepath = Path(ResultDirectory + "Fung_solve_time.csv")
    solve_time_filepath.unlink(missing_ok=True)

    free_energy_filepath = Path(ResultDirectory + "Fung_free_energy.csv")
    free_energy_filepath.unlink(missing_ok=True)

    processes = []
    for i in range(0, len(T)):
        CommandLineInput[4] = str(ResultDirectory + "Fung_Biaxial_Result_" + str(i) + ".vtk")
        CommandLineInput[20] = str(T[i])
        CommandLineInput[22] = str(increments[i])
        process = subprocess.Popen(CommandLineInput)
        processes.append(process)

    for process in processes:
        returncode = process.wait()
        assert int(returncode) == 0, "\033[31mProcess Failed\033[0m"

if __name__ == "__main__":
    main()
