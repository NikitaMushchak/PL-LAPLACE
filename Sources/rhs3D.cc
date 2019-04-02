#ifdef OPENMP
    #include <omp.h>
#endif

#include "ailibrary/ai.hh"
#include "planar3D.hh"

/*!
 \details Функция умножает матрицу на вектор и возвращает результат по ссылке

 \param[in] matrix - исходная матрица (N*N)
 \param[in] vector - исходный вектор (N)
 \param[out] result - результирующий вектор (N)
*/
void multiplyDiag(std::vector<double> &y,
                        //std::vector<double>& A0,
                        std::vector<double>& A1,
                        std::vector<double>& A2,
                        std::vector<double>& A3,
                        std::vector<double>& A4,
                        std::vector<double>& A5,
                        //std::vector<double>& A6,
                        std::vector<double >&x,
                        size_t Nx,
                        size_t NxNy,
                        size_t N_dof){
    //std::vector<double> y(N_dof, 0.);
    //% multiply by 3rd, 4th, 5th diagonals
    y[0] = A3[0]*x[0] + A4[0]*x[1];
    y[N_dof-1] = A3[N_dof-1]*x[N_dof-1] + A2[N_dof-1]*x[N_dof-2];

    for(size_t i = 1; i < N_dof - 2;++i){
        y[i] = A2[i] * x[i-1] + A3[i] * x[i] + A4[i] * x[i+1];
    }
    //% multiply by the 1st diagonal
        for(size_t i = NxNy ; i < N_dof ; ++i){
            y[i] += x[i-NxNy];
    }
    //% multiply by the 2nd diagonal
    // for(size_t i = Nx; i < N_dof ; ++i){
    //     y[i] += A1[i] * x[i-Nx];
    // }
    for(size_t i = 0; i < N_dof - Nx ; ++i){
        y[i+Nx] += A1[i] * x[i];
    }
    //     % multiply by the 6th diagonal
    for(size_t i = 0; i < N_dof-Nx;++i){
        y[i] += A5[i]*x[i+Nx];
    }

    //% multiply by the 7th diagonal
   for(size_t i = 0 ; i < N_dof - NxNy;++i ){
       y[i] += x[i + NxNy];
   }
}

double MultiplyVV(std::vector<double>&a, std::vector<double>&b){

    double c = 0.;
    for(size_t i = 0; i < a.size(); i++){
        c += a[i]*b[i];
    }
    return c;
}

double NormV(std::vector<double>&x){
    double a = 0.;
    for(size_t i = 0 ; i< x.size(); ++i){
        a += x[i]*x[i];
    }
    return a;
}



void conjGrad(
                            std::vector<double> &x,                   // выход функции
                            std::vector<double>& A1,
                            std::vector<double>& A2,
                            std::vector<double>& A3,
                            std::vector<double>& A4,
                            std::vector<double>& A5,
                            //std::vector<double>& A6,// семидиагональная матрица
                            std::vector<double>& b,       // вектор раскрытий
                            std::vector<double>& r1,
                            std::vector<double>& r2,
                            std::vector<double>& p,
                            std::vector<double>& A_p,
                            size_t Nx , size_t NxNy , size_t N_dof){

for(size_t i = 0 ; i < N_dof ; ++i){
        r2[i] = 0.;
        p[i]  = 0.;
        A_p[i]= 0.;
        r1[i] = b[i];
        p[i] = r1[i];
    }
    double eps = 1.;
    double alpha;
    double beta;

    // while(0.0001 < eps){
        multiplyDiag(A_p, A1, A2, A3, A4, A5, p, Nx , NxNy , N_dof);
        // ai::printMarker();
		// ai::saveVector("p",p);
		// ai::saveVector("A_p", A_p);
        alpha = MultiplyVV(r1,r1)/MultiplyVV(A_p,p);
        for(size_t i = 0 ; i < x.size() ;++i){
            x[i]+= alpha*p[i];
            r2[i] = r1[i] - alpha * A_p[i];
        }
        // ai::saveVector("x",x);
        // ai::saveVector("r2",r2);
        beta = MultiplyVV(r2,r2)/MultiplyVV(r1,r1);
         //std::cout<<"beta = "<<std::fixed<<std::setprecision(15)<<beta<<std::endl;
        for(size_t i = 0 ; i < p.size(); i++){
            p[i] = r2[i] + beta*p[i];
            r1[i] = r2[i];
        }
         //ai::saveVector("p",p);
        // r1 = r2;
        eps = NormV(r2)/NormV(x);
        // std::cout<<"eps = "<<eps<<std::endl;
     // }
}
//

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
    std::vector<std::vector<std::size_t> > &index,
    std::vector<std::vector<std::size_t> > &activeElements,
    std::vector<double>& A1,
    std::vector<double>& A2,
    std::vector<double>& A3,
    std::vector<double>& A4,
    std::vector<double>& A5,
    std::vector<double> &opening,
    std::vector<double> &stress,
    std::vector<double> &T,
    std::vector<double> &b,
    std::vector<double> &r1,
    std::vector<double> &r2,
    std::vector<double> &p,
    std::vector<double> &A_p,
    size_t N_dof,
    size_t Nx,
    size_t Ny
){

    std::fill(pressure.begin(), pressure.end(), 0.);
    //ai::printMarker();
    //std::vector<double> partialPressure(N_dof, 0.);;
    //std::vector<double> T(N_dof, 0.);  // unknowns in  finite difference discretization of Laplace equation (AT = b)
    //std::fill(T.begin(), T.end(), 0.);
    //std::vector<double> b(N_dof , 0.);
    //std::fill(b.begin(), b.end(), 0.);
    // ai::saveVector("cop", opening);
    for(size_t i = 0 ;i < N_dof; ++i){
        T[i] = 0.;
        b[i] = 0.;
    }

    for(std::size_t k = 0; k < activeElements.size(); ++k){
        const size_t i = activeElements[k][0];
        const size_t j = activeElements[k][1];

        b[i+(Nx)*j] = - opening[k]; ////0.5!!!!
    }

    // std::cout<<"N_dof = "<<N_dof<<" Nx = "<<Nx<<"  Ny = "<<Ny<<std::endl;
    // ai::saveMatrix("inf",influenceMatrix);
    ai::saveVector("b", b);


    conjGrad(
             T,
             A1,
             A2,
             A3,
             A4,
             A5,
             b,
             r1,
             r2,
             p,
             A_p,
             Nx,
             Nx*Ny,
             N_dof
             );
    ai::saveVector("T", T);
    //std::cout<<"act Elements size = "<<activeElements.size()<<std::endl;
    //ai::saveMatrix("actEl", activeElements);
    // std::vector<double> press(N_dof, 0.);
    // for(size_t i =0 ; i < Nx; ++i){
    //     for(size_t j = 0 ; j < Ny;++j){
    //         press[i+Nx*j] = (0.5*1./(1.- 0.25*0.25)) * (-b[i+Nx*j] - T[i + Nx*j] )/dx;
    //     }
    // }
    std::vector<std::vector<double> > pre(Nx);
    for(size_t i =0 ; i< Nx;++i){
        pre[i].resize(Ny);
    }
    // ai::saveVector("press", press);
    for(std::size_t k = 0; k < activeElements.size(); ++k){
        const size_t i = activeElements[k][0];
        const size_t j = activeElements[k][1];
                                    // E        nu * nu
        pressure[ index[i][j] ] = (0.5*1./(1.- 0.25*0.25)) * (- b[i+Nx*j] - T[i + Nx*j] )/dx + stress[j];
        pre[i][j] = pressure[ index[i][j] ];
    }
    ai::saveVector("pr", pressure);
    ai::saveMatrix("davl", pre);
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
