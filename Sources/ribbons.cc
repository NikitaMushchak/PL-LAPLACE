#include <algorithm>

#include "planar3D.hh"
#include "collectDistanceFromVelocity.hh"

/*!
 \details Функция определяет начальные граничные элементы из списка
 активных элементов с учётом заданной сетки

 \param[in] mesh - сеточная матрица
 \param[in] activeElements - вектор активных элементов.
 \param[in] distances - матрица расстояний до фронта трещины
 \param[in] testDistance - проверочное растояние до фронта трещины
 \param[in] dMin - параметр определения принадлежности ячейки к граничным
 элементам
 \param[in] dMax  - параметр определения принадлежности ячейки к граничным
 элементам
 \return Вектор граничных элементов
*/
std::vector<Ribbon> findRibbons(
    std::vector< std::vector<Cell> > &mesh,
    const std::vector< std::vector<size_t> > activeElements,
    std::vector< std::vector<double> > &distances,
    double testDistance,
    const double dMin,
    const double dMax
){
    std::vector<double> zeroVector(mesh[0].size(), 0.);

    distances.resize(mesh.size());
    std::fill(distances.begin(), distances.end(), zeroVector);

    testDistance += epsilon;

    const double testDistanceMin = std::pow(testDistance - dMax, 2);
    const double testDistanceMax = std::pow(testDistance - dMin, 2);

    std::vector<Ribbon> ribbons;

    for(size_t k = 0; k < activeElements.size(); ++k){
        const size_t i = activeElements[k][0];
        const size_t j = activeElements[k][1];

        const double d = std::pow(mesh[i][j].x, 2)
            + std::pow(mesh[i][j].y, 2);

        if(testDistanceMin < d && testDistanceMax > d){
            ribbons.push_back(Ribbon(i, j));
            distances[i][j] = testDistance - std::sqrt(d);
            mesh[i][j].type = RIBBON;
        }
    }

    return ribbons;
}

/*!
 \details Функция определяет новые граничные элементы из списка активных
 элементов, рассчитывает скорость фронта в них и записывает время активации
 (для учёта утечек в пласт)

 \param[in] i - индекс элемента
 \param[in] j - индекс элемента
 \param[in] d - растояние до фронта для текущего элемента
 \param[in] dMin - параметр определения принадлежности ячейки к граничным
 элементам
 \param[in] dCenter  - параметр определения принадлежности ячейки к граничным
 элементам
 \param[in] n - индекс реологии жидкости
 \param[in] opening - вектор раскрытий
 \param[in] leakOff - вектор коэффициентов Картера в разных слоях
 \param[out] mesh - сеточная матрица
 \param[in] index - матрица индексов
 \param[out] activeElements - вектор активных элементов.
 \param[out] elementIsActive - матрица активных элементов
 \param[out] activationTime - время активации элемента
 \param[out] ribbons - вектор новых граничных элементов трещины
 \param[in] ribbonsOld - вектор исходных граничных элементов трещины
 \param[in] distances - матрица расстояний до фронта трещины
 \param[in] velocities - матрица скоростей фронта
 \param[in] currentTime - текущее расчетное время
 \param[in] openingAtPoint - значение раскрытия в текущем элементе
*/
void findNewRibbons(
    const size_t i,
    const size_t j,
    const double d,
    const double dMin,
    const double dCenter,
    double n,
    std::vector<double> &opening,
    std::vector<double> &leakOff,
    std::vector< std::vector<Cell> > &mesh,
    std::vector< std::vector<size_t> > &index,
    std::vector< std::vector<size_t> > &activeElements,
    std::vector< std::vector<bool> > &elementIsActive,
    std::vector<double> &activationTime,
    std::vector<Ribbon> &ribbons,
    std::vector<Ribbon> &ribbonsOld,
    std::vector< std::vector<double> > &distances,
    std::vector< std::vector<double> > &velocities,
    const double currentTime,
    const double openingAtPoint
){
    const std::vector<size_t> element{i, j};

    const bool elementIsAMember = (
        std::find(activeElements.begin(), activeElements.end(), element)
        != activeElements.end()
    );

    if(elementIsAMember && d > dMin && 1 > mesh[i][j].type){
        ribbons.push_back(Ribbon(i, j));

        mesh[i][j].type = RIBBON;

        if(0 == regime){
            distances[i][j] = 0.25 * sqrt(0.5 * M_PI) * E * openingAtPoint
                / Kic;
            distances[i][j] *= distances[i][j];
        }else{
            distances[i][j] = collectDistanceFromVelocity(
                i,
                j,
                n,
                opening[index[i][j]],
                leakOff,
                ribbonsOld,
                velocities
            );
        }

        if(epsilon > distances[i][j]){
            distances[i][j] = epsilon;
        }
    }

    if(!elementIsAMember && d >= dCenter){
        activeElements.push_back(std::vector<size_t>{i, j});
        elementIsActive[i][j] = true;
        activationTime.push_back(currentTime);
    }
}
