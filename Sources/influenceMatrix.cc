#include "ailibrary/ai.hh"
#include "planar3D.hh"

#if defined(USE_EIGEN)
    #include "eigen/Dense"
    #include "eigen/Sparse"
    #include "eigen/Core"
#endif

/*!
 \details Функция заполняет матрицу коэффициентами влияния с учётом 
симметрии вдоль одного из орт

 \param[out] influenceMatrix - матрица влияния
 \param[in] xSize - размер сетки по i
 \param[in] ySize - размер сетки по j
*/
void buildInfluenceMatrix(
    std::vector< std::vector<double> > &influenceMatrix,
    const size_t xSize,
    const size_t ySize
){
    const double coeff = -1. / (8. * M_PI);
    /// Вектор границ интегрирования по х
    std::vector<double> xIntegrationLimits;
    /// Вектор границ интегрирования по y
    std::vector<double> yIntegrationLimits;
    std::vector<double> xCollocationPoints;
    std::vector<double> yCollocationPoints;
    
    for(size_t i = 0; i < xSize; ++i){
        xCollocationPoints.push_back((double) i - i00);
    }
    // В предыдущих версиях – 0.25
    xCollocationPoints[0] = 0.01;
    
    for(size_t i = 0; i < ySize; ++i){
        yCollocationPoints.push_back((double) i + 1. - j00);
    }
    
    for(size_t i = 0; i <= xSize; ++i){
        xIntegrationLimits.push_back((double) i - i00 - 0.5);
    }
    xIntegrationLimits[0] = 0;
    
    for(size_t i = 0; i <= ySize; ++i){
        yIntegrationLimits.push_back((double) i - j00 + 0.5);
    }
    
    for(size_t k = 0; k < xSize; ++k){
        for(size_t m = 0; m < ySize; ++m){
            const double x1 = xCollocationPoints[k];
            const double x2 = yCollocationPoints[m];
            
            for(size_t i = 0; i < xSize; ++i){
                for(size_t j = 0; j < ySize; ++j){
                    const double a = xIntegrationLimits[i];
                    const double b = xIntegrationLimits[i + 1];
                    const double c = yIntegrationLimits[j];
                    const double d = yIntegrationLimits[j + 1];
                    
                    const double term1 = std::sqrt(std::pow(x1 - b, 2)
                        + std::pow(x2 - c, 2)) / (x1 - b)
                        - std::sqrt(std::pow(x1 - a, 2)
                        + std::pow(x2 - c, 2)) / (x1 - a);
                    const double term2 = std::sqrt(std::pow(x1 - a, 2)
                        + std::pow(x2 - d, 2)) / (x1 - a)
                        - std::sqrt(std::pow(x1 - b, 2)
                        + std::pow(x2 - d, 2)) / (x1 - b);
                    const double term3 = std::sqrt(std::pow(x1 + b, 2)
                        + std::pow(c - x2, 2)) / (x1 + b)
                        - std::sqrt(std::pow(x1 + a, 2)
                        + std::pow(c - x2, 2)) / (x1 + a);
                    const double term4 = std::sqrt(std::pow(x1 + a, 2)
                        + std::pow(d - x2, 2)) / (x1 + a)
                        - std::sqrt(std::pow(x1 + b, 2)
                        + std::pow(d - x2, 2)) / (x1 + b);
                    
                    influenceMatrix[k * ySize + m][i * ySize + j] = coeff 
                        * (
                            term1 / (x2 - c) + term2 / (x2 - d)
                            + term3 / (c - x2) + term4 / (d - x2)
                        );
                }
            }
        }
    }
}

/*!
\details Функция заполняет матрицу значениями из общей матрицы 
 коэффициентов влияния, также заполняя значениями укороченные столбцы 
 раскрытий и напряжений

 \param[in] influenceMatrix - полная матрица влияния
 \param[in] activeElements - матрица активных элементов
 \param[in] opening - вектор раскрытий
 \param[out] partialOpening - укороченный вектор раскрытий
 \param[out] partialInfluenceMatrix - укороченная матрица влияния 
 (только между активными элементами)
 \param[in] index - 
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
){
    #if defined(USE_EIGEN)
    partialOpening.resize(activeElements.size());
    partialOpening.setZero();
    partialInfluenceMatrix.resize(activeElements.size(), activeElements.size());
    partialInfluenceMatrix.setZero();
    #else
    std::vector<double> zeroVector(activeElements.size(), 0);
    partialOpening = zeroVector;
    partialInfluenceMatrix.resize(activeElements.size());
    std::fill(
        partialInfluenceMatrix.begin(),
        partialInfluenceMatrix.end(),
        zeroVector
    );
    #endif

    for(size_t j = 0; j < activeElements.size(); ++j){
        const size_t ai = activeElements[j][0];
        const size_t aj = activeElements[j][1];
        const size_t pc1 = index[ai][aj];
        
        #if defined(USE_EIGEN)
        partialOpening(j) = opening[pc1];
        #else
        partialOpening[j] = opening[pc1];
        #endif
        
        for(size_t i = 0; i < activeElements.size(); ++i){
            const size_t ai2 = activeElements[i][0];
            const size_t aj2 = activeElements[i][1];
            const size_t pc2 = index[ai2][aj2];
            
            #if defined(USE_EIGEN)
            partialInfluenceMatrix(i, j) = influenceMatrix[pc2][pc1];
            #else
            partialInfluenceMatrix[i][j] = influenceMatrix[pc2][pc1];
            #endif
        }
    }
}