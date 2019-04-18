#pragma once

/*!
 \brief Заполнение матрицы коэффициентов влияния
*/
void buildInfluenceMatrix(
    std::vector< std::vector<double> > &influenceMatrix,
    const size_t xSize,
    const size_t ySize
);

/*!
 \brief Создание укороченной мартицы влияния (между активными элементами)
*/
void buildPartialInfluenceMatrix(
    std::vector< std::vector<double> > &influenceMatrix,
    std::vector< std::vector<size_t> > &activeElements,
    std::vector<double> &opening,
    std::vector<double> &partialOpening,
    std::vector< std::vector<double> > &partialInfluenceMatrix,
    std::vector< std::vector<size_t> > &index
);
