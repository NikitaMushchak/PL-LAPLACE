#pragma once

#if defined(USE_EIGEN)
    #include "eigen/Dense"
    #include "eigen/Sparse"
    #include "eigen/Core"
#endif

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
    #if defined(USE_EIGEN)
    Eigen::VectorXd &partialOpening,
    Eigen::MatrixXd &partialInfluenceMatrix,
    #else
    std::vector<double> &partialOpening,
    std::vector< std::vector<double> > &partialInfluenceMatrix,
    #endif
    std::vector< std::vector<size_t> > &index
);
