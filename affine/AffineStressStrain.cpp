#include "AffineStressStrain.hpp"

void CalcGreenLagrangeStrain(const mfem::GridFunction& u, mfem::GridFunction& E)
{
    int NE = u.FESpace()->GetNE();
    int dim = u.FESpace()->GetMesh()->SpaceDimension();

    mfem::DenseMatrix dNdeta, dNdx, F(dim), E_(dim);
    mfem::Array<int> vdofs;

    E = 0.;
    for(int el=0; el<NE; el++)
    {
        const mfem::FiniteElement *FE = u.FESpace()->GetFE(el);
        int dof = FE->GetDof();

        dNdeta.SetSize(dof, dim);
        dNdx.SetSize(dof, dim);
        vdofs.SetSize(dim*dof);

        u.FESpace()->GetElementVDofs(el,vdofs);

        const mfem::IntegrationRule *int_rule = &mfem::IntRules.Get(FE->GetGeomType(), 2*FE->GetOrder());
        mfem::ElementTransformation *Tr = u.FESpace()->GetElementTransformation(el);

        for(int ip=0; ip<int_rule->GetNPoints(); ip++)
        {
            const mfem::IntegrationPoint &int_point = int_rule->IntPoint(ip);
            Tr->SetIntPoint(&int_point);
            FE->CalcDShape(int_point, dNdeta);
            mfem::Mult(dNdeta, Tr->InverseJacobian(), dNdx);
            
            F = 0.;

            for(int i=0; i<dim; i++)
            {
                for(int j=0; j<dim; j++)
                {
                    for(int A=0; A<dof; A++) {F(i,j) += dNdx(A,j) * u[vdofs[A+i*dof]];}
                    if (i == j) {F(i,j) += 1.;}
                }
            }

            mfem::MultAtB(F, F, E_);

            for (int i=0; i<dim; i++)
            {
                for (int j=0; j<dim; j++)
                {
                    if (i == j) {E_(i,j) -= 1.;}
                    E_(i,j) /= 2.;
                }
            }

            for (int i=0; i<dim; i++)
            {
                for (int j=0; j<dim; j++)
                {
                    E[el+(dim*i+j)*NE] += E_(i,j) * int_point.weight;
                }
            }
        }
    }
}

void CalcHyperElasticCauchyStress(const mfem::GridFunction& u, mfem::Coefficient& mu, 
                                    mfem::Coefficient& lambda, mfem::GridFunction& sigma)
{
    int NE = u.FESpace()->GetNE();
    int dim = u.FESpace()->GetMesh()->SpaceDimension();
    
    double lambda_val, mu_val, J;
    mfem::DenseMatrix F(dim), S(dim), E_(dim), dNdeta, dNdx, P(dim), sigma_(dim);
    mfem::Array<int> vdofs;
    
    sigma = 0.;
    for(int el=0; el<NE; el++)
    {
        const mfem::FiniteElement *FE = u.FESpace()->GetFE(el);
        int dof = FE->GetDof();

        dNdeta.SetSize(dof, dim);
        dNdx.SetSize(dof, dim);
        vdofs.SetSize(dim*dof);

        u.FESpace()->GetElementVDofs(el, vdofs);

        const mfem::IntegrationRule *int_rule = &mfem::IntRules.Get(FE->GetGeomType(), 2*FE->GetOrder());
        mfem::ElementTransformation *Tr = u.FESpace()->GetElementTransformation(el);

        for(int ip=0; ip<int_rule->GetNPoints(); ip++)
        {
            const mfem::IntegrationPoint &int_point = int_rule->IntPoint(ip);
            Tr->SetIntPoint(&int_point);
            FE->CalcDShape(int_point, dNdeta);
            mfem::Mult(dNdeta, Tr->InverseJacobian(), dNdx);

            lambda_val = lambda.Eval(*Tr, int_point);
            mu_val = mu.Eval(*Tr, int_point);

            F = 0.; S = 0.;

            for(int i=0; i<dim; i++)
            {
                for(int j=0; j<dim; j++)
                {
                    for (int A=0; A<dof; A++) {F(i,j) += dNdx(A,j) * u[vdofs[A+i*dof]];}
                    if (i == j) {F(i,j) += 1.;}
                }
            }

            mfem::MultAtB(F, F, E_);

            for (int i=0; i<dim; i++)
            {
                for (int j=0; j<dim; j++)
                {
                    if (i == j) {E_(i,j) -= 1.;}
                    E_(i,j) /= 2.;
                }
            }

            for (int i=0; i<dim; i++)
            {
                for (int j=0; j<dim; j++)
                {
                    S(i,j) += 2.*mu_val*E_(i,j);
                    if (i == j)
                    {
                        for (int k=0; k<dim; k++)
                        {
                            S(i,j) += lambda_val*E_(k,k);
                        }
                    }
                }
            }

            J = F.Det();
            mfem::MultABt(S, F, P);
            mfem::Mult(F, P, sigma_);

            for (int i=0; i<dim; i++)
            {
                for (int j=0; j<dim; j++)
                {
                    sigma[NE*(dim*i+j) + el] += sigma_(i,j)*int_point.weight/J;
                }
            }
        }
    }
}

void CalcFungPK2Stress(const mfem::GridFunction& u, mfem::Coefficient& a, mfem::Coefficient& A1, 
                        mfem::Coefficient& A2, mfem::Coefficient& A3, mfem::Coefficient& A4, 
                        mfem::Coefficient& A5, mfem::Coefficient& A6, mfem::GridFunction& S)
{
    int NE = u.FESpace()->GetNE();
    int dim = u.FESpace()->GetMesh()->SpaceDimension();

    MFEM_VERIFY(dim == 2, "PK2 Fung stress recovery only implemented for 2-D.");

    double a_val, A1_val, A2_val, A3_val, A4_val, A5_val, A6_val, Lambda, J;
    mfem::DenseMatrix F(dim), E_(dim), S_(dim), dNdeta, dNdx;
    mfem::Array<int> vdofs;
    
    S = 0.;
    for(int el=0; el<NE; el++)
    {
        const mfem::FiniteElement *FE = u.FESpace()->GetFE(el);
        int dof = FE->GetDof();

        dNdeta.SetSize(dof, dim);
        dNdx.SetSize(dof, dim);
        vdofs.SetSize(dim*dof);

        u.FESpace()->GetElementVDofs(el, vdofs);

        const mfem::IntegrationRule *int_rule = &mfem::IntRules.Get(FE->GetGeomType(), 2*FE->GetOrder());
        mfem::ElementTransformation *Tr = u.FESpace()->GetElementTransformation(el);

        for(int ip=0; ip<int_rule->GetNPoints(); ip++)
        {
            const mfem::IntegrationPoint &int_point = int_rule->IntPoint(ip);
            Tr->SetIntPoint(&int_point);
            FE->CalcDShape(int_point, dNdeta);
            mfem::Mult(dNdeta, Tr->InverseJacobian(), dNdx);

            a_val = a.Eval(*Tr, int_point);
            A1_val = A1.Eval(*Tr, int_point);
            A2_val = A2.Eval(*Tr, int_point);
            A3_val = A3.Eval(*Tr, int_point);
            A4_val = A4.Eval(*Tr, int_point);
            A5_val = A5.Eval(*Tr, int_point);
            A6_val = A6.Eval(*Tr, int_point);

            F = 0.;

            for(int i=0; i<dim; i++)
            {
                for(int j=0; j<dim; j++)
                {
                    for (int A=0; A<dof; A++) {F(i,j) += dNdx(A,j) * u[vdofs[A+i*dof]];}
                    if (i == j) {F(i,j) += 1.;}
                }
            }

            mfem::MultAtB(F, F, E_);

            for (int i=0; i<dim; i++)
            {
                for (int j=0; j<dim; j++)
                {
                    if (i == j) {E_(i,j) -= 1.;}
                    E_(i,j) /= 2.;
                }
            }

            Lambda = A1_val*std::pow(E_(0,0), 2);
            Lambda += A2_val*std::pow(E_(1,1), 2);
            Lambda += 2.*A3_val*E_(0,0)*E_(1,1);
            Lambda += A4_val*std::pow(E_(0,1), 2);
            Lambda += 2.*A5_val*E_(0,1)*E_(0,0);
            Lambda += 2.*A6_val*E_(0,1)*E_(1,1);

            S_(0,0) = a_val*(A1_val*E_(0,0) + A3_val*E_(1,1) + A5_val*E_(0,1))*std::exp(Lambda);
            S_(0,1) = a_val*(A4_val*E_(0,1) + A5_val*E_(0,0) + A6_val*E_(1,1))*std::exp(Lambda);
            S_(1,0) = a_val*(A4_val*E_(0,1) + A5_val*E_(0,0) + A6_val*E_(1,1))*std::exp(Lambda);
            S_(1,1) = a_val*(A2_val*E_(1,1) + A3_val*E_(0,0) + A6_val*E_(0,1))*std::exp(Lambda);

            for (int i=0; i<dim; i++)
            {
                for (int j=0; j<dim; j++)
                {
                    S[NE*(dim*i+j) + el] += S_(i,j)*int_point.weight;
                }
            }
        }
    }    
}

void CalcFungPK1Stress(const mfem::GridFunction& u, mfem::Coefficient& a, mfem::Coefficient& A1, 
                        mfem::Coefficient& A2, mfem::Coefficient& A3, mfem::Coefficient& A4, 
                        mfem::Coefficient& A5, mfem::Coefficient& A6, mfem::GridFunction& P)
{
    int NE = u.FESpace()->GetNE();
    int dim = u.FESpace()->GetMesh()->SpaceDimension();
    
    double a_val, A1_val, A2_val, A3_val, A4_val, A5_val, A6_val, Lambda, J;
    mfem::DenseMatrix F(dim), E_(dim), S(dim), dNdeta, dNdx, P_(dim);
    mfem::Array<int> vdofs;
    
    P = 0.;
    for(int el=0; el<NE; el++)
    {
        const mfem::FiniteElement *FE = u.FESpace()->GetFE(el);
        int dof = FE->GetDof();

        dNdeta.SetSize(dof, dim);
        dNdx.SetSize(dof, dim);
        vdofs.SetSize(dim*dof);

        u.FESpace()->GetElementVDofs(el, vdofs);

        const mfem::IntegrationRule *int_rule = &mfem::IntRules.Get(FE->GetGeomType(), 2*FE->GetOrder());
        mfem::ElementTransformation *Tr = u.FESpace()->GetElementTransformation(el);
        
        for(int ip=0; ip<int_rule->GetNPoints(); ip++)
        {
            const mfem::IntegrationPoint &int_point = int_rule->IntPoint(ip);
            Tr->SetIntPoint(&int_point);
            FE->CalcDShape(int_point, dNdeta);
            mfem::Mult(dNdeta, Tr->InverseJacobian(), dNdx);

            a_val = a.Eval(*Tr, int_point);
            A1_val = A1.Eval(*Tr, int_point);
            A2_val = A2.Eval(*Tr, int_point);
            A3_val = A3.Eval(*Tr, int_point);
            A4_val = A4.Eval(*Tr, int_point);
            A5_val = A5.Eval(*Tr, int_point);
            A6_val = A6.Eval(*Tr, int_point);

            F = 0.;

            for(int i=0; i<dim; i++)
            {
                for(int j=0; j<dim; j++)
                {
                    for (int A=0; A<dof; A++) {F(i,j) += dNdx(A,j) * u[vdofs[A+i*dof]];}
                    if (i == j) {F(i,j) += 1.;}
                }
            }

            mfem::MultAtB(F, F, E_);

            for (int i=0; i<dim; i++)
            {
                for (int j=0; j<dim; j++)
                {
                    if (i == j) {E_(i,j) -= 1.;}
                    E_(i,j) /= 2.;
                }
            }

            Lambda = A1_val*std::pow(E_(0,0), 2);
            Lambda += A2_val*std::pow(E_(1,1), 2);
            Lambda += 2.*A3_val*E_(0,0)*E_(1,1);
            Lambda += A4_val*std::pow(E_(0,1), 2);
            Lambda += 2.*A5_val*E_(0,1)*E_(0,0);
            Lambda += 2.*A6_val*E_(0,1)*E_(1,1);

            S(0,0) = a_val*(A1_val*E_(0,0) + A3_val*E_(1,1) + A5_val*E_(0,1))*std::exp(Lambda);
            S(0,1) = a_val*(A4_val*E_(0,1) + A5_val*E_(0,0) + A6_val*E_(1,1))*std::exp(Lambda);
            S(1,0) = a_val*(A4_val*E_(0,1) + A5_val*E_(0,0) + A6_val*E_(1,1))*std::exp(Lambda);
            S(1,1) = a_val*(A2_val*E_(1,1) + A3_val*E_(0,0) + A6_val*E_(0,1))*std::exp(Lambda);

            mfem::MultABt(S, F, P_);

            for (int i=0; i<dim; i++)
            {
                for (int j=0; j<dim; j++)
                {
                    P[NE*(dim*i+j) + el] += P_(i,j)*int_point.weight;
                }
            }
        }
    }    
}

void CalcFungCauchyStress(const mfem::GridFunction& u, mfem::Coefficient& a, mfem::Coefficient& A1, 
                            mfem::Coefficient& A2, mfem::Coefficient& A3, mfem::Coefficient& A4, 
                            mfem::Coefficient& A5, mfem::Coefficient& A6, mfem::GridFunction& sigma)
{
    int NE = u.FESpace()->GetNE();
    int dim = u.FESpace()->GetMesh()->SpaceDimension();
    
    double a_val, A1_val, A2_val, A3_val, A4_val, A5_val, A6_val, Lambda, J;
    mfem::DenseMatrix F(dim), E_(dim), S(dim), dNdeta, dNdx, P(dim), sigma_(dim);
    mfem::Array<int> vdofs;
    
    sigma = 0.;
    for(int el=0; el<NE; el++)
    {
        const mfem::FiniteElement *FE = u.FESpace()->GetFE(el);
        int dof = FE->GetDof();

        dNdeta.SetSize(dof, dim);
        dNdx.SetSize(dof, dim);
        vdofs.SetSize(dim*dof);

        u.FESpace()->GetElementVDofs(el, vdofs);

        const mfem::IntegrationRule *int_rule = &mfem::IntRules.Get(FE->GetGeomType(), 2*FE->GetOrder());
        mfem::ElementTransformation *Tr = u.FESpace()->GetElementTransformation(el);

        for(int ip=0; ip<int_rule->GetNPoints(); ip++)
        {
            const mfem::IntegrationPoint &int_point = int_rule->IntPoint(ip);
            Tr->SetIntPoint(&int_point);
            FE->CalcDShape(int_point, dNdeta);
            mfem::Mult(dNdeta, Tr->InverseJacobian(), dNdx);

            a_val = a.Eval(*Tr, int_point);
            A1_val = A1.Eval(*Tr, int_point);
            A2_val = A2.Eval(*Tr, int_point);
            A3_val = A3.Eval(*Tr, int_point);
            A4_val = A4.Eval(*Tr, int_point);
            A5_val = A5.Eval(*Tr, int_point);
            A6_val = A6.Eval(*Tr, int_point);

            F = 0.;

            for(int i=0; i<dim; i++)
            {
                for(int j=0; j<dim; j++)
                {
                    for (int A=0; A<dof; A++) {F(i,j) += dNdx(A,j) * u[vdofs[A+i*dof]];}
                    if (i == j) {F(i,j) += 1.;}
                }
            }

            mfem::MultAtB(F, F, E_);

            for (int i=0; i<dim; i++)
            {
                for (int j=0; j<dim; j++)
                {
                    if (i == j) {E_(i,j) -= 1.;}
                    E_(i,j) /= 2.;
                }
            }

            Lambda = A1_val*std::pow(E_(0,0), 2);
            Lambda += A2_val*std::pow(E_(1,1), 2);
            Lambda += 2.*A3_val*E_(0,0)*E_(1,1);
            Lambda += A4_val*std::pow(E_(0,1), 2);
            Lambda += 2.*A5_val*E_(0,1)*E_(0,0);
            Lambda += 2.*A6_val*E_(0,1)*E_(1,1);

            S(0,0) = a_val*(A1_val*E_(0,0) + A3_val*E_(1,1) + A5_val*E_(0,1))*std::exp(Lambda);
            S(0,1) = a_val*(A4_val*E_(0,1) + A5_val*E_(0,0) + A6_val*E_(1,1))*std::exp(Lambda);
            S(1,0) = a_val*(A4_val*E_(0,1) + A5_val*E_(0,0) + A6_val*E_(1,1))*std::exp(Lambda);
            S(1,1) = a_val*(A2_val*E_(1,1) + A3_val*E_(0,0) + A6_val*E_(0,1))*std::exp(Lambda);

            J = F.Det();
            mfem::MultABt(S, F, P);
            mfem::Mult(F, P, sigma_);

            for (int i=0; i<dim; i++)
            {
                for (int j=0; j<dim; j++)
                {
                    sigma[NE*(dim*i+j) + el] += sigma_(i,j) * int_point.weight / J;
                }
            }
        }
    }    
}

void CalcVonMisesStress(mfem::GridFunction& sigma, mfem::GridFunction& VMStress)
{       
    int NE = sigma.FESpace()->GetNE();
    int dim = sigma.FESpace()->GetMesh()->SpaceDimension();

    double InnerProd, AvgTr;
    
    VMStress = 0.;
    for (int el=0; el<NE; el++)
    {
        InnerProd = 0., 
        AvgTr = 0.;
        for (int k=0; k<dim; k++)
        {
            AvgTr += sigma[NE*(dim*k+k)+el]/dim;
        }
        for (int k=0; k<dim; k++)
        {
            sigma[NE*(dim*k+k)+el] -= AvgTr;
        }
        for (int i=0; i<dim; i++)
        {
            for (int j=0; j<dim; j++)
            {
                InnerProd += std::pow(sigma[NE*(dim*i+j)+el], 2);
            }
        }
        VMStress[el] = std::sqrt(3./2.*InnerProd);
    }
}

double GetFungHelmholtzEnergy(const mfem::GridFunction& u, mfem::Coefficient& a, 
                                mfem::Coefficient& A1, mfem::Coefficient& A2, 
                                mfem::Coefficient& A3, mfem::Coefficient& A4, 
                                mfem::Coefficient& A5, mfem::Coefficient& A6)
{
    int NE = u.FESpace()->GetNE();
    int dim = u.FESpace()->GetMesh()->SpaceDimension();
    
    double a_val, A1_val, A2_val, A3_val, A4_val, A5_val, A6_val, Lambda, free_energy = 0.;
    mfem::DenseMatrix F(dim), E(dim), dNdeta, dNdx;
    mfem::Array<int> vdofs;
    
    for(int el=0; el<NE; el++)
    {
        const mfem::FiniteElement *FE = u.FESpace()->GetFE(el);
        int dof = FE->GetDof();

        dNdeta.SetSize(dof, dim);
        dNdx.SetSize(dof, dim);
        vdofs.SetSize(dim*dof);

        u.FESpace()->GetElementVDofs(el, vdofs);

        const mfem::IntegrationRule *int_rule = &mfem::IntRules.Get(FE->GetGeomType(), 2*FE->GetOrder());
        mfem::ElementTransformation *Tr = u.FESpace()->GetElementTransformation(el);

        for(int ip=0; ip<int_rule->GetNPoints(); ip++)
        {
            const mfem::IntegrationPoint &int_point = int_rule->IntPoint(ip);
            Tr->SetIntPoint(&int_point);
            FE->CalcDShape(int_point, dNdeta);
            mfem::Mult(dNdeta, Tr->InverseJacobian(), dNdx);

            a_val = a.Eval(*Tr, int_point);
            A1_val = A1.Eval(*Tr, int_point);
            A2_val = A2.Eval(*Tr, int_point);
            A3_val = A3.Eval(*Tr, int_point);
            A4_val = A4.Eval(*Tr, int_point);
            A5_val = A5.Eval(*Tr, int_point);
            A6_val = A6.Eval(*Tr, int_point);

            F = 0.;
            for(int i=0; i<dim; i++)
            {
                for(int j=0; j<dim; j++)
                {
                    for (int A=0; A<dof; A++) {F(i,j) += dNdx(A,j) * u[vdofs[A+i*dof]];}
                    if (i == j) {F(i,j) += 1.;}
                }
            }

            mfem::MultAtB(F, F, E);

            for (int i=0; i<dim; i++)
            {
                for (int j=0; j<dim; j++)
                {
                    if (i == j) {E(i,j) -= 1.;}
                    E(i,j) /= 2.;
                }
            } 

            Lambda = A1_val*std::pow(E(0,0), 2);
            Lambda += A2_val*std::pow(E(1,1), 2);
            Lambda += 2.*A3_val*E(0,0)*E(1,1);
            Lambda += A4_val*std::pow(E(0,1), 2);
            Lambda += 2.*A5_val*E(0,1)*E(0,0);
            Lambda += 2.*A6_val*E(0,1)*E(1,1);

            free_energy += (a_val/2.) * (std::exp(Lambda) - 1.) * int_point.weight * Tr->Weight();
        }
    }

    return free_energy;
}
