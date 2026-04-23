import matplotlib.pyplot as plt
from LigamentAnalysis import GetArray

def GetRadialVMStress(nodes, sigma_VM):
    
    r, sigma_VM_r = [], []

    for i in range(0,nodes.shape[0]):
        if abs(nodes[i,1]) < 1e-8:
            r.append(nodes[i,0])
            sigma_VM_r.append(sigma_VM[i]) 
    
    # Due to formatting of vtk cell data there are duplicates that must be removed
    unique_r = []
    for val in r:
        if val not in unique_r:
            unique_r.append(val)
    
    unique_sigma_VM_r = []
    for val in sigma_VM_r:
        if val not in unique_sigma_VM_r:
            unique_sigma_VM_r.append(val)

    # Use cell center positions for stress locations
    unique_cell_center_r = [(unique_r[i] + unique_r[i+1]) for i in range(len(unique_r) - 1)]

    return unique_cell_center_r, unique_sigma_VM_r


def GetRadialDisplacements(nodes, u):

    r, u_r = [], []

    for i in range(0,nodes.shape[0]):
        if abs(nodes[i,1]) < 1e-8:
            r.append(nodes[i,0])
            u_r.append(u[i]) 

    return r, u_r

def main():

    nodes = GetArray("../results/HyperElasticArtery_0.vtk", "nodes")

    lp_he_sigma_vm = GetArray("../results/HyperElasticArtery_0.vtk", "sigma_VM")
    r, lp_he_sigma_vm_r = GetRadialVMStress(nodes, lp_he_sigma_vm) # Low pressure hyper elastic stress

    hp_he_sigma_vm = GetArray("../results/HyperElasticArtery_8.vtk", "sigma_VM")
    r, hp_he_sigma_vm_r = GetRadialVMStress(nodes, hp_he_sigma_vm) # High pressure hyper elastic stress

    lp_fung_sigma_vm = GetArray("../results/FungArtery_0.vtk", "sigma_VM")
    r, lp_fung_sigma_vm_r = GetRadialVMStress(nodes, lp_fung_sigma_vm) # Low pressure fung stress

    hp_fung_sigma_vm = GetArray("../results/FungArtery_8.vtk", "sigma_VM")
    r, hp_fung_sigma_vm_r = GetRadialVMStress(nodes, hp_fung_sigma_vm) # High pressure fung stress

    fig1 = plt.figure(figsize=[8,4])
    if fig1:
        plt.plot(r, lp_he_sigma_vm_r, color="darkorange", marker="s",
                    linewidth=3, markersize=8, label=r"Hyper Elastic, $p$ = 16 kPa")
        plt.plot(r, hp_he_sigma_vm_r, color="royalblue", marker="^", 
                    linewidth=3, markersize=8, label=r"Hyper Elastic, $p$ = 24 kPa")
        plt.plot(r, lp_fung_sigma_vm_r, color="mediumpurple", marker="8",
                    linewidth=3, markersize=8, label=r"Fung, $p$ = 16 kPa")
        plt.plot(r, hp_fung_sigma_vm_r, color="firebrick", marker="o", 
                    linewidth=3, markersize=8, label=r"Fung, $p$ = 24 kPa")
        plt.legend(fontsize=13)
        plt.xlabel("Radial Position (mm)", fontsize=14)
        plt.xticks(fontsize=14)
        plt.ylabel("Von Mises Stress (MPa)", fontsize=14)
        plt.yticks(fontsize=14)
        plt.tight_layout(pad=0.8)
        plt.savefig("../results/RadialVMStressComparison.png", dpi=300)
        plt.close(fig1)

    lp_he_u = GetArray("../results/HyperElasticArtery_0.vtk", "u")[:,0]
    r, lp_he_u_r = GetRadialDisplacements(nodes, lp_he_u) # Low pressure hyper elastic displacements

    hp_he_u = GetArray("../results/HyperElasticArtery_8.vtk", "u")[:,0]
    r, hp_he_u_r = GetRadialDisplacements(nodes, hp_he_u) # High pressure hyper elastic displacements

    lp_fung_u = GetArray("../results/FungArtery_0.vtk", "u")[:,0]
    r, lp_fung_u_r = GetRadialDisplacements(nodes, lp_fung_u) # Low pressure Fung displacements

    hp_fung_u = GetArray("../results/FungArtery_8.vtk", "u")[:,0]
    r, hp_fung_u_r = GetRadialDisplacements(nodes, hp_fung_u) # High pressure Fung displacements

    fig2 = plt.figure(figsize=[8,4])
    if fig2:
        plt.plot(r, lp_he_u_r, color="darkorange", marker="s",
                    linewidth=3, markersize=8, label=r"Hyper Elastic, $p$ = 16 kPa")
        plt.plot(r, hp_he_u_r, color="royalblue", marker="^", 
                    linewidth=3, markersize=8, label=r"Hyper Elastic, $p$ = 24 kPa")
        plt.plot(r, lp_fung_u_r, color="mediumpurple", marker="8",
                    linewidth=3, markersize=8, label=r"Fung, $p$ = 16 kPa")
        plt.plot(r, hp_fung_u_r, color="firebrick", marker="o", 
                    linewidth=3, markersize=8, label=r"Fung, $p$ = 24 kPa")
        plt.legend(fontsize=13)
        plt.xlabel("Radial Position (mm)", fontsize=14)
        plt.xticks(fontsize=14)
        plt.ylabel("Radial Displacement (mm)", fontsize=14)
        plt.yticks(fontsize=14)
        plt.ylim([0.2, 0.5])
        plt.tight_layout(pad=0.8)
        plt.savefig("../results/RadialDisplacementComparison.png", dpi=300)
        plt.close(fig2)

if __name__ == "__main__":
    main()
