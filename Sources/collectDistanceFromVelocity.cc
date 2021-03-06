#include <algorithm>

#include "planar3D.hh"

/*!
 \details Функция расчитывает скорость граничного элемента с учётом
 соседних граничных элементов и значения раскрытия и определяет расстояние 
 до фронта в соответствии с выбранным режимом распространения

 \param[in] i - индекс элемента
 \param[in] j - индекс элемента
 \param[in] n - индекс реологии жидкости
 \param[in] opening - значение раскрытия в текущем элементе
 \param[in] leakOff - вектор коэффициентов Картера в разных слоях
 \param[out] ribbons - вектор граничных элементов трещины
 \param[in] velocities - матрица скоростей фронта
 \return Расстояние до фронта
*/
double collectDistanceFromVelocity(
    const size_t i,
    const size_t j,
    double n,
    double opening,
    std::vector<double> &leakOff,
    std::vector<Ribbon> &ribbons,
    std::vector< std::vector<double> > &velocities
){
    double averageVelocity = 0;
    
    size_t counter = 1;
    
    auto collectVelocity = [&averageVelocity, &counter, &ribbons, &velocities]
        (const size_t i, const size_t j) -> void
    {
        const Ribbon element(i, j);
        
        if(std::find(ribbons.begin(), ribbons.end(), element) != ribbons.end()){
            averageVelocity += velocities[i][j];
            ++counter;
        }
    };
    
    collectVelocity(i - 1, j);
    collectVelocity(i + 1, j);
    collectVelocity(i, j + 1);
    collectVelocity(i, j - 1);
    
    collectVelocity(i - 1, j - 1);
    collectVelocity(i + 1, j + 1);
    collectVelocity(i - 1, j + 1);
    collectVelocity(i + 1, j - 1);
    
    averageVelocity /= (double) counter;
    
    if(-1 == regime){
        const double alphaL = (n + 4.) / (4. * n + 4.);
        const double Ainf = std::pow(Amu, (1. + 2. * alphaL) / 3.);
        
        return std::pow(
            std::pow(
                averageVelocity * std::pow(4. * leakOff[j], 2), 
                (alphaL - 1.) / 3.
            ) * opening / Ainf, 
            1. / alphaL
        );
    }
    
    return std::pow(
        std::pow(averageVelocity, alpha - 1.) * opening / Amu, 
        1. / alpha
    );
}

