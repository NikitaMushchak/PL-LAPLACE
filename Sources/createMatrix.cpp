#include <vector>


void createMatrixDiag(std::vector<std::vector<double> >& A,
                            size_t& N_dof, size_t& Nx, size_t& Ny, size_t& Nz){
    size_t NxNy = Nx * Ny;
    for(size_t i = 0; i < N_dof ;++i){
        A[i][0] = 1.;
        A[i][1] = 1.;
        A[i][2] = 1.;
        A[i][3] = -6.;
        A[i][4] = 1.;
        A[i][5] = 1.;
        A[i][6] = 1.;
    }
   double n;
   double p ;
   for(size_t j = 0 ; j < Ny; ++j){
       for(size_t k = 0 ; k < Nz; ++k){
           n = j * Nx + k * NxNy;
           A[n][2] = 0.;
           A[n][4] = 2.;

           p = n + Nx - 1;
           A[p][4] = 0.;
           A[p][3] = -5.;
       }
   }
   for(size_t i = 0; i < Nx; ++i ){
      for(size_t k = 0; k < Nz; ++k){
          n = i  + k * NxNy;
          A[n][1] = 0.;
          A[n][3] = -5.;

          p = n + (Ny-1)*Nx;
          A[p][5] = 0.;
          A[p][3] = -5.;
      }
  }
  
  for(size_t i = 0; i<Nx ;++i){
      for(size_t j = 0; j < Ny; ++j){
          n = i  + j * Nx;
          A[i + j*Nx][0] = 0.;

          A[n + (Nz-1)*NxNy][5] = 0.;
          A[n + (Nz-1)*NxNy][3] = -5.;
      }
  }
}
