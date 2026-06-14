#include "mfem.hpp"
#include "affine.hpp"
#include <chrono>
#include <filesystem>

int main(int argc, char** argv)
{
    // Material parameters
    double a, A1, A2, A3, A4, A5, A6;

    // Top and right displacement
    double u_x = 0., u_y = 0.;
    
    // I/O parameters
    std::string MeshFile, ResultFile;
    int output_csv;

    mfem::OptionsParser args(argc, argv);
    args.AddOption(&MeshFile, "-mf", "--MeshFile", "Mesh File");
    args.AddOption(&ResultFile, "-rf", "--ResultFile", "Result File");
    args.AddOption(&a, "-a", "--a", "a");
    args.AddOption(&A1, "-A1", "--A1", "A1");
    args.AddOption(&A2, "-A2", "--A2", "A2");
    args.AddOption(&A3, "-A3", "--A3", "A3");
    args.AddOption(&A4, "-A4", "--A4", "A4");
    args.AddOption(&A5, "-A5", "--A5", "A5");
    args.AddOption(&A6, "-A6", "--A6", "A6");
    args.AddOption(&u_x, "-u_x", "--u_x", "Right Displacement");
    args.AddOption(&u_y, "-u_y", "--u_y", "Top Displacement");
    args.AddOption(&output_csv, "-o", "--output_csv", "Ouput the CSV Data Files");
    args.Parse();
    if (!args.Good())
    {
       args.PrintUsage(std::cout);
       return 1;
    }

    auto start = std::chrono::system_clock::now();

    mfem::Mesh mesh = mfem::Mesh(MeshFile.c_str(), 1, 1);
    int dim = mesh.Dimension();

    auto u_ec = mfem::H1_FECollection(2, dim, mfem::BasisType::GaussLobatto); 
    auto u_space = mfem::FiniteElementSpace(&mesh, &u_ec, dim);

    mfem::GridFunction u(&u_space);
    u = 0.;
    mfem::GridFunction f(&u_space);
    f = 0.;

    auto a_coeff = mfem::ConstantCoefficient(a);
    auto A1_coeff = mfem::ConstantCoefficient(A1);
    auto A2_coeff = mfem::ConstantCoefficient(A2);
    auto A3_coeff = mfem::ConstantCoefficient(A3);
    auto A4_coeff = mfem::ConstantCoefficient(A4);
    auto A5_coeff = mfem::ConstantCoefficient(A5);
    auto A6_coeff = mfem::ConstantCoefficient(A6);

    mfem::Array<int> ess_tdofs, tmp_tdofs, right_tdofs, top_tdofs;

    mfem::Array<int> bottom({1, 0, 0, 0, 0, 0});
    mfem::Array<int> right({0, 1, 0, 0, 0, 0});
    mfem::Array<int> top({0, 0, 1, 0, 0, 0});
    mfem::Array<int> left({0, 0, 0, 1, 0, 0});

    u_space.GetEssentialTrueDofs(bottom, tmp_tdofs, 1);
    ess_tdofs.Append(tmp_tdofs);
    
    u_space.GetEssentialTrueDofs(right, right_tdofs, 0);
    ess_tdofs.Append(right_tdofs);

    u_space.GetEssentialTrueDofs(top, top_tdofs, 1);
    ess_tdofs.Append(top_tdofs);
    
    u_space.GetEssentialTrueDofs(left, tmp_tdofs, 0);
    ess_tdofs.Append(tmp_tdofs);

    auto B = mfem::NonlinearForm(&u_space);
    B.AddDomainIntegrator(new FungExponentialIntegrator(a_coeff, A1_coeff, A2_coeff, A3_coeff, 
                                                        A4_coeff, A5_coeff, A6_coeff));
    B.SetEssentialTrueDofs(ess_tdofs);

    auto prec = mfem::UMFPackSolver();
    auto ns = mfem::NewtonSolver();
    ns.SetOperator(B);
    ns.SetPreconditioner(prec);
    ns.SetRelTol(1e-12);
    ns.SetAbsTol(1e-8);
    ns.SetMaxIter(40);
    ns.SetPrintLevel(0);

    // Number of increments depends on strain applied
    double max_u = std::max(u_x, u_y);
    int N_increments = static_cast<int>(std::abs(max_u/0.0005)); 
    for (int i=0; i<N_increments; i++)
    {
        //mfem::out << "Solving increment " << (i+1) << " out of " << N_increments << " \n";
        u.SetSubVector(right_tdofs, (static_cast<double>(i+1)/N_increments)*u_x);
        u.SetSubVector(top_tdofs, (static_cast<double>(i+1)/N_increments)*u_y);
        ns.Mult(f, u);
    }

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    auto dg_ec = mfem::DG_FECollection(0, dim, mfem::BasisType::GaussLegendre);
    auto dg_tensor_space = mfem::FiniteElementSpace(&mesh, &dg_ec, dim*dim);

    auto E = mfem::GridFunction(&dg_tensor_space);
    CalcGreenLagrangeStrain(u, E);

    auto sigma = mfem::GridFunction(&dg_tensor_space);
    CalcFungCauchyStress(u, E, a_coeff, A1_coeff, A2_coeff, A3_coeff, 
                            A4_coeff, A5_coeff, A6_coeff, sigma);
    
    auto P = mfem::GridFunction(&dg_tensor_space);
    CalcFungPK1Stress(u, E, a_coeff, A1_coeff, A2_coeff, A3_coeff, 
                        A4_coeff, A5_coeff, A6_coeff, P);
    
    double free_energy = GetFungHelmholtzEnergy(E, a_coeff, A1_coeff, A2_coeff, A3_coeff, 
                                                A4_coeff, A5_coeff, A6_coeff);

    std::ofstream file(ResultFile);
    file.precision(16);
    mesh.PrintVTK(file, 0);
    u.SaveVTK(file, "u", 0);
    E.SaveVTK(file, "E", 0);
    sigma.SaveVTK(file, "sigma", 0);
    P.SaveVTK(file, "P", 0);
    file.close();

    if (output_csv)
    {
        std::string DOF_file_name = "../results/BiaxialExtension/Fung_DOF.csv";
        if (std::filesystem::exists(DOF_file_name)) 
        {
            std::ofstream DOF_file(DOF_file_name, std::ios::app);
            DOF_file << u_x << ", " << u_y << ", " << u.Size() << "\n";
            DOF_file.close();
        } 
        else 
        {
            std::ofstream DOF_file(DOF_file_name);
            DOF_file << "# DOF counts\n";
            DOF_file << u_x << ", " << u_y << ", " << u.Size() << "\n";
            DOF_file.close();
        }

        std::string solve_time_file_name = "../results/BiaxialExtension/Fung_solve_time.csv";
        if (std::filesystem::exists(solve_time_file_name)) 
        {
            std::ofstream solve_time_file(solve_time_file_name, std::ios::app);
            solve_time_file << u_x << ", " << u_y << ", " << elapsed.count() << "\n";
            solve_time_file.close();
        } 
        else 
        {
            std::ofstream solve_time_file(solve_time_file_name);
            solve_time_file << "# Solve times for each strain\n";
            solve_time_file << u_x << ", " << u_y << ", " << elapsed.count() << "\n";
            solve_time_file.close();
        }
        
        std::string free_energy_file_name = "../results/BiaxialExtension/Fung_free_energy.csv";
        if (std::filesystem::exists(free_energy_file_name))
        {
            std::ofstream free_energy_file(free_energy_file_name, std::ios::app);
            free_energy_file << u_x << ", " << u_y << ", " << free_energy << "\n";
            free_energy_file.close();
        } 
        else 
        {
            std::ofstream free_energy_file(free_energy_file_name);
            free_energy_file << "# Total Helmholtz free energy for each strain\n";
            free_energy_file << u_x << ", " << u_y << ", " << free_energy << "\n";
            free_energy_file.close();
        }
    }

    return 0;
}
