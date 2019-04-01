#pragma once

/*!
 \brief Нахождение активных элементов
*/
void findActiveElements(
    std::vector< std::vector<std::size_t> > &activeElements,
    std::vector< std::vector<bool> > &elementIsActive,
    std::vector<double> &x,
    std::vector<double> &y,
    double testDistance
);