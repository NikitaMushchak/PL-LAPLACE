#ifdef OPENMP
    #include <omp.h>
#endif

#include "ailibrary/ai.hh"
#include "planar3D.hh"
#include "conjgrad.hh"

/*!
 \details Функция умножает матрицу на вектор и возвращает результат по ссылке

 \param[in] matrix - исходная матрица (N*N)
 \param[in] vector - исходный вектор (N)
 \param[out] result - результирующий вектор (N)
*/
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

     //while(0.01 < eps){
        //ai::saveVector("p", p);
        multiplyDiag(A_p, A, p, Nx , NxNy , N_dof);
        //ai::saveVector("A_p ", A_p);
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

    // }
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



inline void multiply(
    std::vector< std::vector<double> > &matrix,
    std::vector<double> &vector,
    std::vector<double> &result
){
    int size = (int) vector.size();

    result.resize(size);

    #pragma omp parallel for
    for(int i = 0; i < size; ++i){
          int j = 0;

        for(; j <= size - 16; j += 16){
            result[i] += matrix[i][j] * vector[j]
                + matrix[i][j + 1] * vector[j + 1]
                + matrix[i][j + 2] * vector[j + 2]
                + matrix[i][j + 3] * vector[j + 3]
                + matrix[i][j + 4] * vector[j + 4]
                + matrix[i][j + 5] * vector[j + 5]
                + matrix[i][j + 6] * vector[j + 6]
                + matrix[i][j + 7] * vector[j + 7]
                + matrix[i][j + 8] * vector[j + 8]
                + matrix[i][j + 9] * vector[j + 9]
                + matrix[i][j + 10] * vector[j + 10]
                + matrix[i][j + 11] * vector[j + 11]
                + matrix[i][j + 12] * vector[j + 12]
                + matrix[i][j + 13] * vector[j + 13]
                + matrix[i][j + 14] * vector[j + 14]
                + matrix[i][j + 15] * vector[j + 15];
        }

        for(; j < size; ++j){
            result[i] += matrix[i][j] * vector[j];
        }
    }
}



/*!
 \details Функция рассчитывает давление в активных элементах с помощью матрицы
 коэффициентов влияния, раскрытия и заданного контраста напряжений

 \param[out] pressure – вектор давлений
 \param[in] index - матрица индексов
 \param[in] activeElements - вектор активных элементов.
 \param[in] influenceMatrix - матрица коэффициентов влияния
 \param[in] opening - вектор раскрытий
 \param[in] stress - вектор сжимающих напряжений в разных слоях
*/
void calculatePressure(
    std::vector<double> &pressure,
    std::vector< std::vector<std::size_t> > &index,
    std::vector< std::vector<std::size_t> > &activeElements,
    std::vector< std::vector<double> > &influenceMatrix,
    std::vector<double> &opening,
    std::vector<double> &stress,
    size_t N_dof,
    size_t Nx,
    size_t Ny
){
    std::fill(pressure.begin(), pressure.end(), 0.);
    //ai::printMarker();
    //std::vector<double> partialPressure(N_dof, 0.);;
    std::vector<double> T(N_dof, 0.);  // unknowns in  finite difference discretization of Laplace equation (AT = b)

    std::vector<double> b(N_dof , 0.);

    for(size_t i = 0 ; i < Nx; ++i){
        for(size_t j = 0 ; j < Ny; ++j){
            b[i+Nx*j] = opening[j+Ny*i];
        }
    }


    std::cout<<"N_dof = "<<N_dof<<" Nx = "<<Nx<<"  Ny = "<<Ny<<std::endl;
    ai::saveMatrix("inf",influenceMatrix);
    ai::saveVector("b", b);
    conjGrad(T, influenceMatrix, b, Nx , Nx*Ny , N_dof);

    ai::saveVector("T", T);
    //std::cout<<"act Elements size = "<<activeElements.size()<<std::endl;
    //ai::saveMatrix("actEl", activeElements);
    std::vector<double> press(N_dof, 0.);
    for(size_t i =0 ; i < Nx;++i){
        for(size_t j = 0 ; j < Ny;++j){
            press[i+Nx*j] = (0.5*1./(1.- 0.25*0.25)) * (-b[i+Nx*j] - T[i + j*Nx] )/dx;
        }
    }

    ai::saveVector("press", press);
    for(std::size_t k = 0; k < activeElements.size(); ++k){
        const size_t i = activeElements[k][0];
        const size_t j = activeElements[k][1];
        // pressure[index[i][j]] = partialPressure[k] / dx + stress[j];
        pressure[index[i][j]] = -press[i+Nx*j]/dx ;//+stress[j];
    }
    ai::saveVector("pr", pressure);
}

/*!
 \details Функция рассчитывает скорость фронта в заданном режиме развития

 \param[out] velocities - матрица скоростей фронта
 \param[in] leakOff - вектор коэффициентов Картера в разных слоях
 \param[in] index - матрица индексов
 \param[in] ribbons - вектор граничных элементов
 \param[in] distances - матрица расстояний до фронта трещины
 \param[in] opening - вектор раскрытий
 \param[in] n - индекс реологии жидкости
*/
void calculateFrontVelocity(
    std::vector< std::vector<double> > &velocities,
    std::vector<double> &leakOff,
    std::vector< std::vector<size_t> > &index,
    std::vector<Ribbon> &ribbons,
    std::vector< std::vector<double> > &distances,
    std::vector<double> &opening,
    double n
){
    velocities.resize(distances.size());
    std::vector<double> zeroVector(distances[0].size(), 0);
    std::fill(velocities.begin(), velocities.end(), zeroVector);

    if(1 == regime){
        //viscosity dominated regime
        const double D3 = 1. / (1. - alpha);

        for(size_t i = 0; i < ribbons.size(); ++i){
            const size_t iRibbon = ribbons[i].i;
            const size_t jRibbon = ribbons[i].j;

            const double distance = distances[iRibbon][jRibbon];

            if(epsilon <= distance){
                velocities[iRibbon][jRibbon] =
                    std::pow(opening[index[iRibbon][jRibbon]]
                        / (Amu * std::pow(distance, alpha)), D3);
            }else{
                velocities[iRibbon][jRibbon] = 0.;
            }
        }
    }else{
        if(-1 == regime){
            //leak-off dominated regime
            const double alphaL = (n + 4.) / (4. * n + 4.);
            const double D3 = 3. / (1. - alphaL);
            alpha = 2. / (n + 2.);
            double BAlpha = 0.25 * alphaL * tan(0.5 * M_PI - M_PI * (1 - alphaL));
            double Ainf = std::pow(BAlpha * (1 - alphaL), -0.5 * alphaL);
            Ainf = std::pow(Ainf, (1. + 2. * alphaL) / 3.);

            for(size_t i = 0; i < ribbons.size(); ++i){
                const size_t iRibbon = ribbons[i].i;
                const size_t jRibbon = ribbons[i].j;

                const double distance = distances[iRibbon][jRibbon];

                if(epsilon <= distance){
                    velocities[iRibbon][jRibbon] =
                        std::pow(opening[index[iRibbon][jRibbon]]
                            / (Ainf * std::pow(distance, alphaL)), D3)
                        / std::pow(4. * leakOff[jRibbon], 2);
                }else{
                    velocities[iRibbon][jRibbon] = 0.;
                }
            }
        }
    }
}

/*!
 \details Функция рассчитывает производные раскрытия и концентрации пропанта по
 времени с учетом давления и раскрытия в соседних элементах, а также утечки
 жидкости в пласт

 \param[in] dWdt - вектор скоростей изменения раскрытия
 \param[in] dCdt - вектор скоростей изменения раскрытия
 \param[in] opening - вектор раскрытий
 \param[in] concentration - вектор концентраций пропанта
 \param[in] pressure - вектор давлений
 \param[in] leakOff - вектор коэффициентов Картера в разных слоях
 \param[in] index - матрица индексов
 \param[in] activeElements - вектор активных элементов.
 \param[in] elementIsActive - матрица активных элементов
 \param[in] activationTime - время активации элемента
 \param[in] currentTime - текущее расчетное время
 \param[in] fluidDensity - плотность жидкости
 \param[in] proppantDensity - плотность пропанта
 \param[in] n - индекс реологии жидкости
*/
void calculateOpeningAndConcentrationSpeeds(
    std::vector<double> &dWdt,
    std::vector<double> &dCdt,
    std::vector<double> &opening,
    std::vector<double> &concentration,
    std::vector<double> &pressure,
    std::vector<double> &leakOff,
    std::vector< std::vector<size_t> > &index,
    std::vector< std::vector<size_t> > &activeElements,
    std::vector< std::vector<bool> > &elementIsActive,
    std::vector<double> &activationTime,
    double currentTime,
    double fluidDensity,
    double proppantDensity,
    double n
){
    dWdt.resize(opening.size());
    std::fill(dWdt.begin(), dWdt.end(), 0.);

    dCdt.resize(concentration.size());
    std::fill(dCdt.begin(), dCdt.end(), 0.);

    const double WPow = (2. * n + 1.) / n;
    const double PPow = 1. / n;

    /// Параметры модели распространения пропанта

    const double maximumConcentration = 0.585;
    const double CPow = -2.5;

    /// Диаметр пропанта

    const double proppantDiameter = 2. * std::pow(10., -3);

    /// Скорость оседания

    const double settlingVelocity = 476.7 * 9.81
        * (proppantDensity - fluidDensity)
        * std::pow(proppantDiameter, n + 1.) / (18. * std::pow(3., n - 1.));

    double fluidFlowRight;
    double fluidFlowBottom;
    double proppantFlowRight;
    double proppantFlowBottom;
    double settlingFlow;
    double pressureIJ;
    double pressureDrop;
    double openingOnTheBorder;
    double concentrationOnTheBorder;

    for(std::size_t k = 0; k < activeElements.size(); ++k){
        const size_t i = activeElements[k][0];
        const size_t j = activeElements[k][1];

        fluidFlowRight = 0.;
        fluidFlowBottom = 0.;
        proppantFlowRight = 0.;
        proppantFlowBottom = 0.;
        settlingFlow = 0.;
        pressureIJ = pressure[index[i][j]];

        /// Расчёт потоков жидкости и пропанта

        if(elementIsActive[i + 1][j]){
            pressureDrop = pressure[index[i + 1][j]] - pressureIJ;
            openingOnTheBorder = 0.5
                * (opening[index[i][j]] + opening[index[i + 1][j]]);
            concentrationOnTheBorder = 0.5
                * (concentration[index[i][j]] + concentration[index[i + 1][j]]);
            fluidFlowRight = ai::sign(pressureDrop)
                * std::pow(openingOnTheBorder, WPow)
                * std::pow(std::abs(pressureDrop) / dx, PPow)
                / std::pow(
                    1. - concentrationOnTheBorder / 0.6, // / maximumConcentration,
                    CPow
                );

            if(
                openingOnTheBorder > epsilon
                && !(
                    (concentration[index[i][j]] >= maximumConcentration && pressureDrop > 0)
                    || (concentration[index[i + 1][j]] >= maximumConcentration && pressureDrop < 0)
                )
            ){
                proppantFlowRight = fluidFlowRight / openingOnTheBorder;

                if(0. < proppantFlowRight){
                    proppantFlowRight *= concentration[index[i + 1][j]]
                        * opening[index[i + 1][j]];
                }else{
                    proppantFlowRight *= concentration[index[i][j]]
                        * opening[index[i][j]];
                }
            }
        }
        if(elementIsActive[i][j + 1]){
            pressureDrop = pressure[index[i][j + 1]] - pressureIJ;
            openingOnTheBorder = 0.5
                * (opening[index[i][j]] + opening[index[i][j + 1]]);
            concentrationOnTheBorder = 0.5
                * (concentration[index[i][j]] + concentration[index[i][j + 1]]);
            fluidFlowBottom = ai::sign(pressureDrop)
                * std::pow(openingOnTheBorder, WPow)
                * std::pow(std::abs(pressureDrop) / dy, PPow)
                / std::pow(
                    1. - concentrationOnTheBorder / 0.6,// / maximumConcentration,
                    CPow
                );

            if(
                true
            ){
                if(
                    openingOnTheBorder > epsilon
                    && !(
                        (concentration[index[i][j]] >= maximumConcentration && pressureDrop > 0)
                        || (concentration[index[i][j + 1]] >= maximumConcentration && pressureDrop < 0)
                    )
                ){
                    proppantFlowBottom = fluidFlowBottom / openingOnTheBorder;
                }


                settlingFlow = settlingVelocity * (
                    1. - concentrationOnTheBorder / 0.6
                );

                if(0. < proppantFlowBottom){
                    proppantFlowBottom *= concentration[index[i][j + 1]]
                        * opening[index[i][j + 1]];
                    settlingFlow *= concentration[index[i][j + 1]]
                        * opening[index[i][j + 1]];
                }else{
                    proppantFlowBottom *= concentration[index[i][j]]
                        * opening[index[i][j]];
                    settlingFlow *= concentration[index[i][j]]
                        * opening[index[i][j]];
                }
            }
        }

        dWdt[index[i][j]] += (fluidFlowRight + fluidFlowBottom) / dx
            - leakOff[j] / std::sqrt(currentTime - activationTime[k]);
        dWdt[index[i + 1][j]] -= fluidFlowRight / dx;
        dWdt[index[i][j + 1]] -= fluidFlowBottom / dy;

        dCdt[index[i][j]] += (proppantFlowRight + proppantFlowBottom + settlingFlow) / dx;
        dCdt[index[i + 1][j]] -= proppantFlowRight / dx;
        dCdt[index[i][j + 1]] -= (proppantFlowBottom + settlingFlow) / dy;

        if(i00 == i){
            dWdt[index[i][j]] += fluidFlowRight / dx;
            dCdt[index[i][j]] += proppantFlowRight / dx;
        }
    }

    dWdt[index[i00][j00]] += fluidInjection / (dx * dy);
    dCdt[index[i00][j00]] += proppantInjection * fluidInjection
        / (proppantDensity * dx * dy);
        // * maximumConcentration / (proppantDensity * dx * dy);
}
