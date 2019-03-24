#include "conjgrad.hh"
#include <vector>

void conjGrad(std::vector<double> &x,
                    std::vector<std::vector<double> >&A,
                            std::vector<double>&b,
                            size_t Nx , size_t NxNy , size_t N_dof){

    //std::vector<double> residual;
    std::vector<double> r1(N_dof, 0.);
    std::vector<double> r2(N_dof, 0.);
    std::vector<double> p(N_dof, 0.);

    r1 = b;
    p = r1;
    double eps = 1.;
    size_t n_iter = 1;

    std::vector<double> A_p(N_dof , 0.);

    double alpha;
    double beta;

     while(0.01 < eps){

        multiplyDiag(A_p, A, p, Nx , NxNy , N_dof);

        alpha = MultiplyVV(r1,r1)/MultiplyVV(A_p,p);

        for(size_t i = 0 ; i< x.size() ;++i){
            x[i]+= alpha*p[i];
        }

        for(size_t i = 0 ; i < r2.size(); i++){
            r2[i] = r1[i] - alpha * A_p[i];
        }

        beta = MultiplyVV(r2,r2)/MultiplyVV(r1,r1);

        for(size_t i = 0 ; i < p.size(); i++){
            p[i] = r2[i] + beta*p[i];
        }

        r1 = r2;

        eps = NormV(r2)/NormV(x);

        //residual.push_back(eps);
        ++n_iter;

     }
}

void multiplyDiag(std::vector<double> &y, std::vector<std::vector<double> >&A,std::vector<double >&x,
                        size_t Nx , size_t NxNy , size_t N_dof){

    y[0] = A[0][3]*x[0] + A[0][4]*x[1];
    y[N_dof-1] = A[N_dof-1][3]*x[N_dof-1] + A[N_dof-1][2]*x[N_dof-2];
    for(size_t i = 1; i < N_dof - 2;++i){
        y[i] = A[i][2] * x[i-1] + A[i][3] * x[i] + A[i][4] * x[i+1];
    }

    for(size_t i = NxNy ; i < N_dof ; ++i){
            y[i] += A[i][0] * x[i-NxNy];
    }

    for(size_t i = Nx; i < N_dof ; ++i){
        y[i] += A[i][1] * x[i-Nx];
    }

    for(size_t i = 0; i < N_dof-Nx;++i){
        y[i] += A[i][5]*x[i+Nx];
    }


   for(size_t i = 0 ; i<N_dof - NxNy;++i ){
       y[i] += A[i][6]*x[i+NxNy];
   }
}

double MultiplyVV(std::vector<double>&a, std::vector<double>&b){
    if( a.size()!=b.size() ){

    }
    double c = 0.;
    for(size_t i = 0; i < a.size(); i++){
        c=c+ a[i]*b[i];
    }
    return c;
}

double NormV(std::vector<double>&x){
    double a = 0.;
    for(size_t i = 0 ; i< x.size(); ++i){
        a = a + x[i]*x[i];
    }
    return sqrt(a);
}
