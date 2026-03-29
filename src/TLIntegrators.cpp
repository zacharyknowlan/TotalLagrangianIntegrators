#include "TLIntegrators.hpp"

using FEI = FungExponentialIntegrator;

FEI::FungExponentialIntegrator(mfem::Coefficient& a, mfem::Coefficient& A1, mfem::Coefficient& A2,
                                mfem::Coefficient& A3, mfem::Coefficient& A4, mfem::Coefficient& A5,
                                mfem::Coefficient& A6)
{
    mat_props_.SetSize(7);
    mat_props_[0] = &a;
    mat_props_[1] = &A1;
    mat_props_[2] = &A2;
    mat_props_[3] = &A3;
    mat_props_[4] = &A4;
    mat_props_[5] = &A5;
    mat_props_[6] = &A6;
}

int FEI::MapVoigt2D(const int i, const int j)
{
        if ((i == 0 && j == 1) || (i == 1 && j == 0)) {return 2;} // 12
        return i; // For 11 and 22
}

void FEI::ComputeE()
{   
    tmp1_.SetSize(Gradu_.Height());
    mfem::MultAtB(Gradu_, Gradu_, tmp1_);
    for (int i=0; i<Gradu_.Height(); i++)
    {
        for (int j=0; j<Gradu_.Width(); j++)
        {
            E_(i,j) = (Gradu_(i,j) + Gradu_(j,i) + tmp1_(i,j))/2.;
        }
    }
}

void FEI::ComputeLambda()
{   
    Lambda_ = A1_*std::pow(E_(0,0), 2);
    Lambda_ += A2_*std::pow(E_(1,1), 2);
    Lambda_ += 2.*A3_*E_(0,0)*E_(1,1);
    Lambda_ += A4_*std::pow(E_(0,1), 2);
    Lambda_ += 2.*A5_*E_(0,1)*E_(0,0);
    Lambda_ += 2.*A6_*E_(0,1)*E_(1,1);
}

void FEI::ComputedLambdadE()
{
    dLambdadE_(0,0) = 2.*(A1_*E_(0,0) + A3_*E_(1,1) + A5_*E_(0,1));
    dLambdadE_(0,1) = 2.*(A4_*E_(0,1) + A5_*E_(0,0) + A6_*E_(1,1));
    dLambdadE_(1,0) = dLambdadE_(0,1);
    dLambdadE_(1,1) = 2.*(A2_*E_(1,1) + A3_*E_(0,0) + A6_*E_(0,1));
}

void FEI::ComputeA()
{
    A_(0,0) = 2.*A1_;
    A_(1,1) = 2.*A2_;
    A_(0,1) = A_(1,0) = 2.*A3_;
    A_(2,2) = 2.*A4_;
    A_(2,0) = A_(0,2) = 2.*A5_;
    A_(2,1) = A_(1,2) = 2.*A6_;
}

void FEI::AssembleElementVector(const mfem::FiniteElement& el, mfem::ElementTransformation& Tr,
                                const mfem::Vector& elfun, mfem::Vector& elvect)
{
    int I; // Aux. index;
    int dof = el.GetDof();
    const int dim = Tr.GetSpaceDim();
    MFEM_VERIFY(dim == 2, "Fung integrator only supported for dim == 2.");

    dNdeta_.SetSize(dof, dim);
    dNdx_.SetSize(dof, dim);  
    Gradu_.SetSize(dim);
    F_.SetSize(dim);
    S_.SetSize(dim);
    tmp1_.SetSize(dim);
    dLambdadE_.SetSize(dim);
    E_.SetSize(dim);
    elvect.SetSize(dof*dim);

    const mfem::IntegrationRule *ir= &mfem::IntRules.Get(el.GetGeomType(), 3*el.GetOrder());

    elvect = 0.; // Initialize element residual vector to zero
    for (int ip=0; ip<ir->GetNPoints(); ip++)
    {
        const mfem::IntegrationPoint &int_point = ir->IntPoint(ip);
        Tr.SetIntPoint(&int_point);
        el.CalcDShape(int_point, dNdeta_);
        mfem::Mult(dNdeta_, Tr.InverseJacobian(), dNdx_);

        a_ = mat_props_[0]->Eval(Tr, int_point);
        A1_ = mat_props_[1]->Eval(Tr, int_point);
        A2_ = mat_props_[2]->Eval(Tr, int_point);
        A3_ = mat_props_[3]->Eval(Tr, int_point);
        A4_ = mat_props_[4]->Eval(Tr, int_point);
        A5_ = mat_props_[5]->Eval(Tr, int_point);
        A6_ = mat_props_[6]->Eval(Tr, int_point);

        // Compute displacement gradient and deformation gradient
        Gradu_ = 0.;
        for (int i=0; i<dim; i++)
        {
            for (int j=0; j<dim; j++)
            {
                for (int C=0; C<dof; C++)
                {
                    Gradu_(i,j) += dNdx_(C,j)*elfun[i*dof+C];
                }
                F_(i,j) = Gradu_(i,j) + static_cast<double>(i==j);
            }
        }

        ComputeE(); // Green-Lagrange strain stored in E_
        ComputeLambda(); // Exponential energy
        ComputedLambdadE(); // Exponential energy derivative stored in dLambdadE_

        // Compute stress measure S
        S_ = 0.;
        for (int k=0; k<dim; k++)
        {
            for (int m=0; m<dim; m++)
            {
                S_(k,m) += (a_/2.)*dLambdadE_(k,m)*std::exp(Lambda_);
            }
        }

        // Transform S
        mfem::MultABt(S_, F_, tmp1_);

        // Compute the residual vector
        for (int i=0; i<dim; i++)
        {
            for (int B=0; B<dof; B++)
            {
                I = i*dof+B;
                for (int j=0; j<dim; j++)
                {
                    elvect[I] += dNdx_(B,j) * tmp1_(j,i) * int_point.weight * Tr.Weight(); // Transposed
                }
            }
        }
    }
}



void FEI::AssembleElementGrad(const mfem::FiniteElement& el, mfem::ElementTransformation& Tr,
                                const mfem::Vector& elfun, mfem::DenseMatrix& elmat)
    {
    int I, J, K, L; // Aux. indices
    int dof = el.GetDof();
    const int dim = Tr.GetSpaceDim();
    MFEM_VERIFY(dim == 2, "Fung integrator only supported for dim == 2.");

    dNdeta_.SetSize(dof, dim);
    dNdx_.SetSize(dof, dim);  
    Gradu_.SetSize(dim);
    F_.SetSize(dim);
    S_.SetSize(dim);
    dSdu_.SetSize(dim*dim, dof*dim);
    tmp1_.SetSize(dim*dim, dof*dim);
    E_.SetSize(dim);
    dLambdadE_.SetSize(dim);
    A_.SetSize(dim*(dim+1)/2);
    elmat.SetSize(dof*dim, dof*dim);

    const mfem::IntegrationRule *ir= &mfem::IntRules.Get(el.GetGeomType(), 3*el.GetOrder());

    elmat = 0.; // Set element gradient matrix to zero
    for (int ip=0; ip<ir->GetNPoints(); ip++)
    {
        const mfem::IntegrationPoint &int_point = ir->IntPoint(ip);
        Tr.SetIntPoint(&int_point);
        el.CalcDShape(int_point, dNdeta_);
        mfem::Mult(dNdeta_, Tr.InverseJacobian(), dNdx_);

        a_ = mat_props_[0]->Eval(Tr, int_point);
        A1_ = mat_props_[1]->Eval(Tr, int_point);
        A2_ = mat_props_[2]->Eval(Tr, int_point);
        A3_ = mat_props_[3]->Eval(Tr, int_point);
        A4_ = mat_props_[4]->Eval(Tr, int_point);
        A5_ = mat_props_[5]->Eval(Tr, int_point);
        A6_ = mat_props_[6]->Eval(Tr, int_point);
        
        // Compute displacement gradient and deformation gradient
        Gradu_ = 0.;
        for (int i=0; i<dim; i++)
        {
            for (int j=0; j<dim; j++)
            {
                for (int C=0; C<dof; C++)
                {
                    Gradu_(i,j) += dNdx_(C,j)*elfun[i*dof+C];
                }
                F_(i,j) = Gradu_(i,j) + static_cast<double>(i==j);
            }
        }

        ComputeE(); // Green-Lagrange strain stored in E_
        ComputeLambda(); // Exponential energy
        ComputedLambdadE(); // Exponential energy derivative stored in dLambdadE_
        ComputeA(); // Material tensor A_{klmn} stored in A_

        // Compute stress measure S and its derivatives dSdu
        S_ = 0.;
        dSdu_ = 0.;
        for (int k=0; k<dim; k++)
        {
            for (int m=0; m<dim; m++)
            {
                S_(k,m) += (a_/2.)*dLambdadE_(k,m)*std::exp(Lambda_);

                I = k*dim+m;
                K = MapVoigt2D(k,m);

                for (int v=0; v<dim; v++)
                {
                    for (int H=0; H<dof; H++)
                    {
                        J = v*dof+H;
                        for (int n=0; n<dim; n++)
                        {
                            for (int p=0; p<=n; p++) // Prevents double-counting
                            {
                                L = MapVoigt2D(n,p);
                                
                                dSdu_(I,J) += (a_/4.)*A_(K,L)*dNdx_(H,p)*(n==v)*std::exp(Lambda_);
                                dSdu_(I,J) += (a_/4.)*A_(K,L)*dNdx_(H,n)*(p==v)*std::exp(Lambda_);
                                dSdu_(I,J) += (a_/4.)*A_(K,L)*dNdx_(H,n)*Gradu_(v,p)*std::exp(Lambda_);
                                dSdu_(I,J) += (a_/4.)*A_(K,L)*dNdx_(H,p)*Gradu_(v,n)*std::exp(Lambda_);
                                dSdu_(I,J) += (a_/4.)*dLambdadE_(k,m)*dLambdadE_(n,p)*dNdx_(H,p)*(n==v)*std::exp(Lambda_);
                                dSdu_(I,J) += (a_/4.)*dLambdadE_(k,m)*dLambdadE_(n,p)*dNdx_(H,n)*(p==v)*std::exp(Lambda_);
                                dSdu_(I,J) += (a_/4.)*dLambdadE_(k,m)*dLambdadE_(n,p)*dNdx_(H,n)*Gradu_(v,p)*std::exp(Lambda_);
                                dSdu_(I,J) += (a_/4.)*dLambdadE_(k,m)*dLambdadE_(n,p)*dNdx_(H,p)*Gradu_(v,n)*std::exp(Lambda_);     
                            }
                        }
                    }
                }
            }
        }
        
        tmp1_ = 0.;
        for (int i=0; i<dim; i++)
        {
            for (int j=0; j<dim; j++)
            {
                for (int vH=0; vH<(dim*dof); vH++)
                {
                    for (int m=0; m<dim; m++)
                    {
                        tmp1_((j*dim+i),vH) += dSdu_((j*dim+m),vH)*F_(i,m);
                    }
                }
            }
        }

        for (int i=0; i<dim; i++)
        {
            for (int B=0; B<dof; B++)
            {
                I = i*dof+B;
                for (int v=0; v<dim; v++)
                {
                    for (int H=0; H<dof; H++)
                    {
                        J = v*dof+H;
                        for (int j=0; j<dim; j++)
                        {
                            // Note that both tmp1 and S are transposed
                            elmat(I,J) += dNdx_(B,j) * tmp1_((j*dim+i),(v*dof+H)) * int_point.weight * Tr.Weight();
                            for (int m=0; m<dim; m++)
                            {
                                elmat(I,J) += dNdx_(B,j) * S_(j,m) * dNdx_(H,m) * (i==v) * int_point.weight * Tr.Weight();
                            }
                        }
                    }
                }
            }
        }
    }
}
