#include "planar3D.hh"

/*!
 \details Функция определяет активные элементы в пределах
 указанного расстояния, опираясь на заданную сетку

 \param[out] activeElements - вектор активных элементов.
 \param[out] elementIsActive - матрица активных элементов
 \param[in] x - вектор координат центров ячеек по х
 \param[in] y - вектор координат центров ячеек по y
 \param[in] testDistance - расстояние до фронта трещины
*/
void findActiveElements(
    std::vector< std::vector<std::size_t> > &activeElements,
    std::vector< std::vector<bool> > &elementIsActive,
    std::vector<double> &x,
    std::vector<double> &y,
    double testDistance
){
    const std::size_t xSize = x.size();
    const std::size_t ySize = y.size();

    testDistance = std::pow(testDistance, 2);

    elementIsActive.resize(xSize);

    activeElements.clear();

    for(std::size_t i = 0; i < xSize; ++i){
        elementIsActive[i].resize(ySize);

        for(std::size_t j = 0; j < ySize; ++j){
            if(
                std::pow(std::abs(x[i]) + 0.25 * dx, 2)
                + std::pow(std::abs(y[j]) + 0.25 * dy, 2)
                < testDistance
            ){
                elementIsActive[i][j] = true;

                activeElements.push_back(std::vector<std::size_t>{i, j});
            }else{
                elementIsActive[i][j] = false;
            }
        }
    }
}
