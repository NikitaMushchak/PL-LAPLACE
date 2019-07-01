#include "ailibrary/ai.hh"
#include "planar3D.hh"

void multiplyFiveDiagonalsFast(
    std::vector<double> &vector,
    std::vector<double> &result,
    const double coeff,
    const std::size_t Nx,
    const std::size_t Ny
){
    const std::size_t NxNy = Nx * Ny;

    result.resize(NxNy);

    result[0] = coeff * vector[0];
    result[NxNy - 1] = coeff * vector[NxNy - 1];

    for(std::size_t i = 1; i + 1 < NxNy; ++i){
        result[i] = coeff * vector[i] 
            - vector[i - 1] - vector[i + 1];
    }
    for(std::size_t i = 0; i < Ny; ++i){
        result[i] -= 2 * vector[i + Ny];
    }
    for(std::size_t i = Ny; i < NxNy - Ny; ++i){
        result[i] -= vector[i + Ny] + vector[i - Ny];
    }
    for(std::size_t i = NxNy - Ny; i < NxNy; ++i){
        result[i] -= vector[i - Ny];
    }
}
double MultiplyVV(std::vector<double>&a, std::vector<double>&b){

    double c = 0.;
    for(size_t i = 0; i < a.size(); i++){
        c+= a[i]*b[i];
    }
    return c;
}
double NormV(std::vector<double>&x){
    double a = 0.;
    for(size_t i = 0 ; i< x.size(); ++i){
        a = a+x[i]*x[i];
    }
    return a;
}
void conjugateGradient(
    std::vector<double> &b,
    std::vector<double> &x,
    const double coeff,
    const std::size_t Nx,
    const std::size_t Ny
){
    const std::size_t N = b.size();

    x.resize(N);
    std::fill(x.begin(), x.end(), 0.);

    std::vector<double> r1 = b;
    std::vector<double> r2 = r1;
    std::vector<double> p = r1;

    std::vector<double> A_p;

    double error = 1.;
    double alpha;
    double beta;

    std::size_t counter = 0;

    while(0.00001 < std::abs(error)){
        ++counter;
        multiplyFiveDiagonalsFast(p, A_p, coeff, Nx, Ny);

        alpha = MultiplyVV(r1,r1)/MultiplyVV(A_p,p);
        for(size_t i = 0 ; i < x.size() ;++i){
            x[i]+= alpha*p[i];
            r2[i] = r1[i] - alpha * A_p[i];
        }

        beta = MultiplyVV(r2,r2)/MultiplyVV(r1,r1);
        for(size_t i = 0 ; i < p.size(); i++){
            p[i] = r2[i] + beta*p[i];
            r1[i] = r2[i];
        }
        error = NormV(r2)/NormV(x);
     }
}