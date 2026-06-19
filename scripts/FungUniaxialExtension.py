import subprocess
from pathlib import Path
from SunSacksMeshes import L

ResultDirectory = "../results/UniaxialExtension/"

a = 5.12e3 # Pa
A1 = 60.12
A2 = 86.34
A3 = 2.00
A4 = 203.16
A5 = 43.05
A6 = 42.14
E11 = [0.01*i for i in range(1,21)]
u_x = [L*strain for strain in E11]
increments = [int(10.*u/(L/100.)) for u in u_x]

def main():

    CommandLineInput = ["../build/problems/FUE"]
    CommandLineInput.extend(["--MeshFile", "../meshes/SunSacksSquare.msh"])
    CommandLineInput.extend(["--ResultFile", " "])
    CommandLineInput.extend(["-a", str(a)])
    CommandLineInput.extend(["-A1", str(A1)])
    CommandLineInput.extend(["-A2", str(A2)])
    CommandLineInput.extend(["-A3", str(A3)])
    CommandLineInput.extend(["-A4", str(A4)])
    CommandLineInput.extend(["-A5", str(A5)])
    CommandLineInput.extend(["-A6", str(A6)])
    CommandLineInput.extend(["-u", str(E11[0])])
    CommandLineInput.extend(["-o", "1"])

    DOF_filepath = Path(ResultDirectory + "Fung_DOF.csv")
    DOF_filepath.unlink(missing_ok=True)

    solve_time_filepath = Path(ResultDirectory + "Fung_solve_time.csv")
    solve_time_filepath.unlink(missing_ok=True)

    free_energy_filepath = Path(ResultDirectory + "Fung_free_energy.csv")
    free_energy_filepath.unlink(missing_ok=True)

    processes = []
    for i in range(0, len(E11)):
        CommandLineInput[4] = str(ResultDirectory + "Fung_Uniaxial_Result_" + str(i) + ".vtk")
        CommandLineInput[20] = str(u_x[i])
        CommandLineInput[22] = str(increments[i])
        process = subprocess.Popen(CommandLineInput)
        processes.append(process)

    for process in processes:
        returncode = process.wait()
        assert int(returncode) == 0, "\033[31mProcess Failed\033[0m"

if __name__ == "__main__":
    main()
