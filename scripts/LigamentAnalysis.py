import matplotlib.pyplot as plt
import vtk
import vtk.util.numpy_support as vns #type: ignore
from RunHyperElasticLigament import strains
from LigamentMesh import l

def GetArray(filename, Array):
    
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()
    Output = reader.GetOutput()
    
    if Array == "nodes":
        ArrayData = vns.vtk_to_numpy(Output.GetPoints().GetData())
    else:
        ArrayData = vns.vtk_to_numpy(Output.GetPointData().GetArray(Array))
    
    return ArrayData

def GetRightEdgeAvgStress(nodes, sigma_xx, edge_x=l):
    
    node_count = 0
    field_value_sum = 0.

    for i in range(0,nodes.shape[0]):
        if abs(nodes[i,0] - edge_x) < 1e-8:
            field_value_sum += abs(sigma_xx[i]) 
            node_count += 1

    return (field_value_sum/node_count)

def main():

    nodes = GetArray("../results/HyperElasticLigament_0.vtk", "nodes")
    hyper_elastic_stresses, fung_stresses = [], []
    
    for i in range(0, len(strains)):

        hyper_elastic_filename = str("../results/HyperElasticLigament_" + str(i) + ".vtk")
        hyper_elastic_sigma_xx = GetArray(hyper_elastic_filename, "sigma0")
        hyper_elastic_stresses.append(GetRightEdgeAvgStress(nodes, hyper_elastic_sigma_xx))

        fung_filename = str("../results/FungLigament_" + str(i) + ".vtk")
        fung_sigma_xx = GetArray(fung_filename, "sigma0")
        fung_stresses.append(GetRightEdgeAvgStress(nodes, fung_sigma_xx))

    fig1 = plt.figure(figsize=[8,4])
    if fig1:
        plt.plot([0.] + strains, [0.] + hyper_elastic_stresses, color="darkorange", marker="s",
                 linewidth=3, markersize=8, label="Hyper Elastic")
        plt.plot([0.] + strains, [0.] + fung_stresses, color="royalblue", marker="^", 
                 linewidth=3, markersize=8, label="Fung Exponential")
        plt.legend(fontsize=14)
        plt.xlabel("Strain", fontsize=14)
        plt.xticks(fontsize=14)
        plt.ylabel("Stress (MPa)", fontsize=14)
        plt.yticks(fontsize=14)
        plt.tight_layout(pad=0.8)
        plt.savefig("../results/StressStrainComparison.png", dpi=300)
        plt.close(fig1)

if __name__ == "__main__":
    main()
