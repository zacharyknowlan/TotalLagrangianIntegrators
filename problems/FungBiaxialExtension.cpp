#include "mfem.hpp"
#include "affine.hpp"
#include <chrono>
#include <filesystem>

double T_mag = 1e6;

void T_bottom_func(const mfem::Vector& x, double pseudo_time, mfem::Vector& T) 
{
    T.SetSize(2);
    T[0] = 0.;
    T[1] = -T_mag * pseudo_time;
}

void T_right_func(const mfem::Vector& x, double pseudo_time, mfem::Vector& T) 
{
    T.SetSize(2);
    T[0] = T_mag * pseudo_time;
    T[1] = 0.;
}

void T_top_func(const mfem::Vector& x, double pseudo_time, mfem::Vector& T) 
{
    T.SetSize(2);
    T[0] = 0.;
    T[1] = T_mag * pseudo_time;
}

void T_left_func(const mfem::Vector& x, double pseudo_time, mfem::Vector& T) 
{
    T.SetSize(2);
    T[0] = -T_mag * pseudo_time;
    T[1] = 0.;
}

int main(int argc, char** argv)
{
    // Material parameters
    double a, A1, A2, A3, A4, A5, A6;
    
    // I/O parameters
    std::string MeshFile, ResultFile;
    int output_csv;
    
    int num_increments = 1; // Number of load increments

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
    args.AddOption(&T_mag, "-T", "--T_mag", "TractionMagnitude");
    args.AddOption(&num_increments, "-i", "--increments", "Number of Load Increments");
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

    mfem::Array<int> ess_tdofs, vdofs;
    mfem::DenseMatrix coords;
    double tol = 1e-6, L = 0.025;
    for(int bdr_el=0; bdr_el<u_space.GetNBE(); bdr_el++)
    {
        const mfem::FiniteElement *BE = u_space.GetBE(bdr_el);
        int dof = BE->GetDof();
        
        coords.SetSize(dim, dof);
        coords = 0.;
        
        mfem::ElementTransformation *Tr = u_space.GetBdrElementTransformation(bdr_el);
        const mfem::IntegrationRule nodes = BE->GetNodes(); 
        Tr->Transform(nodes, coords);

        vdofs.SetSize(dim*dof);
        u_space.GetBdrElementVDofs(bdr_el, vdofs);

        for(int i=0; i<dof; i++)
        { 
            if (abs(coords(0, i)) < tol && abs(coords(1, i)) < tol)
            {
                // Origin is fixed in both x and y
                ess_tdofs.Append(vdofs[i]);
                ess_tdofs.Append(vdofs[i+dof]);
            }
            else if (abs(coords(0, i) - L) < tol && abs(coords(1, i)) < tol)
            {
                // Fix bottom right corner y dof
                ess_tdofs.Append(vdofs[i+dof]);
            }
            else if (abs(coords(0, i)) < tol && abs(coords(1, i) - L) < tol)
            {
                // Fix top left corner x dof
                ess_tdofs.Append(vdofs[i]); 
            }
        }
    }

    mfem::Array<int> bottom({1, 0, 0, 0, 0, 0});
    mfem::Array<int> right({0, 1, 0, 0, 0, 0});
    mfem::Array<int> top({0, 0, 1, 0, 0, 0});
    mfem::Array<int> left({0, 0, 0, 1, 0, 0});

    auto T_bottom = mfem::VectorFunctionCoefficient(2, T_bottom_func);
    auto T_right = mfem::VectorFunctionCoefficient(2, T_right_func);
    auto T_top = mfem::VectorFunctionCoefficient(2, T_top_func);
    auto T_left = mfem::VectorFunctionCoefficient(2, T_left_func);

    auto B = mfem::NonlinearForm(&u_space);
    B.AddDomainIntegrator(new FungExponentialIntegrator(a_coeff, A1_coeff, A2_coeff, A3_coeff, 
                                                        A4_coeff, A5_coeff, A6_coeff));
    B.AddBoundaryIntegrator(new PK1TractionIntegrator(T_bottom), bottom);
    B.AddBoundaryIntegrator(new PK1TractionIntegrator(T_right), right);
    B.AddBoundaryIntegrator(new PK1TractionIntegrator(T_top), top);
    B.AddBoundaryIntegrator(new PK1TractionIntegrator(T_left), left);
    B.SetEssentialTrueDofs(ess_tdofs);

    auto prec = mfem::UMFPackSolver();
    auto ns = mfem::NewtonSolver();
    ns.SetOperator(B);
    ns.SetPreconditioner(prec);
    ns.SetRelTol(1e-12);
    ns.SetAbsTol(1e-8);
    ns.SetMaxIter(40);
    ns.SetPrintLevel(0);

    for (int i=0; i<num_increments; i++)
    {
        //mfem::out << "Solving increment " << (i+1) << " out of " << N_increments << " \n";
        double load_fraction = static_cast<double>(i+1)/num_increments;
        T_bottom.SetTime(load_fraction);
        T_right.SetTime(load_fraction);
        T_top.SetTime(load_fraction);
        T_left.SetTime(load_fraction);
        ns.Mult(f, u);
    }

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    auto dg_ec = mfem::DG_FECollection(0, dim, mfem::BasisType::GaussLegendre);
    auto dg_tensor_space = mfem::FiniteElementSpace(&mesh, &dg_ec, dim*dim);

    auto E = mfem::GridFunction(&dg_tensor_space);
    CalcGreenLagrangeStrain(u, E);

    auto S = mfem::GridFunction(&dg_tensor_space);
    CalcFungPK2Stress(u, a_coeff, A1_coeff, A2_coeff, A3_coeff, 
                        A4_coeff, A5_coeff, A6_coeff, S);

    auto P = mfem::GridFunction(&dg_tensor_space);
    CalcFungPK1Stress(u, a_coeff, A1_coeff, A2_coeff, A3_coeff, 
                        A4_coeff, A5_coeff, A6_coeff, P);
                        
    auto sigma = mfem::GridFunction(&dg_tensor_space);
    CalcFungCauchyStress(u, a_coeff, A1_coeff, A2_coeff, A3_coeff, 
                            A4_coeff, A5_coeff, A6_coeff, sigma);
    
    double free_energy = GetFungHelmholtzEnergy(u, a_coeff, A1_coeff, A2_coeff, A3_coeff, 
                                                A4_coeff, A5_coeff, A6_coeff);

    std::ofstream file(ResultFile);
    file.precision(16);
    mesh.PrintVTK(file, 0);
    u.SaveVTK(file, "u", 0);
    E.SaveVTK(file, "E", 0);
    S.SaveVTK(file, "S", 0);
    P.SaveVTK(file, "P", 0);
    sigma.SaveVTK(file, "sigma", 0);
    file.close();

    if (output_csv)
    {
        std::string DOF_file_name = "../results/BiaxialExtension/Fung_DOF.csv";
        if (std::filesystem::exists(DOF_file_name)) 
        {
            std::ofstream DOF_file(DOF_file_name, std::ios::app);
            DOF_file << T_mag << ", " << u.Size() << "\n";
            DOF_file.close();
        } 
        else 
        {
            std::ofstream DOF_file(DOF_file_name);
            DOF_file << "# DOF counts\n";
            DOF_file << T_mag << ", " << u.Size() << "\n";
            DOF_file.close();
        }

        std::string solve_time_file_name = "../results/BiaxialExtension/Fung_solve_time.csv";
        if (std::filesystem::exists(solve_time_file_name)) 
        {
            std::ofstream solve_time_file(solve_time_file_name, std::ios::app);
            solve_time_file << T_mag << ", " << elapsed.count() << "\n";
            solve_time_file.close();
        } 
        else 
        {
            std::ofstream solve_time_file(solve_time_file_name);
            solve_time_file << "# Solve times for each strain\n";
            solve_time_file << T_mag << ", " << elapsed.count() << "\n";
            solve_time_file.close();
        }
        
        std::string free_energy_file_name = "../results/BiaxialExtension/Fung_free_energy.csv";
        if (std::filesystem::exists(free_energy_file_name))
        {
            std::ofstream free_energy_file(free_energy_file_name, std::ios::app);
            free_energy_file << T_mag << ", " << free_energy << "\n";
            free_energy_file.close();
        } 
        else 
        {
            std::ofstream free_energy_file(free_energy_file_name);
            free_energy_file << "# Total Helmholtz free energy for each strain\n";
            free_energy_file << T_mag << ", " << free_energy << "\n";
            free_energy_file.close();
        }
    }

    return 0;
}
