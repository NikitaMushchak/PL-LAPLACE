#include <exception>
#include <vector>
#include <cmath>

/*!
 \brief Обращение матрицы
 \details Функция заменяет переданную матрицу на обратную
*/
void lo_ober(
    std::vector< std::vector<double> > &a
){
    const int n = a.size();
    double max;
    int k, index;
    const double eps = 0.0000001;  // точность
    std::vector<double> x(n);
    std::vector< std::vector<double> > y(n);

    k = 0;

    for (int i = 0; i < n; i++)
    {
        y[i].resize(n);
        for (int j = 0; j < n; j++)
        {
            if (i == j)
                y[i][j] = 1;
            else
                y[i][j] = 0;
        }
    }

    while (k < n)
    {
        // Поиск строки с максимальным a[i][k]
        max = std::abs(a[k][k]);
        index = k;
        for (int i = k + 1; i < n; i++)
        {
            if (std::abs(a[i][k]) > max)
            {
                max = std::abs(a[i][k]);
                index = i;
            }
        }
        // Перестановка строк
        if (max < eps)
        {
            // нет ненулевых диагональных элементов
            throw std::runtime_error("Cannot calculate matrix inverse due to zero diagonal.");
        }
        for (int j = 0; j < n; j++)
        {
            double temp = a[k][j];
            a[k][j] = a[index][j];
            a[index][j] = temp;

            double temp1 = y[k][j];
            y[k][j] = y[index][j];
            y[index][j] = temp1;
        }
        // Нормализация уравнений
        for (int i = k; i < n; i++)
        {
            double temp = a[i][k];
            if (std::abs(temp) < eps) continue; // для нулевого коэффициента пропустить
            for (int j = 0; j < n; j++)
            {
                a[i][j] = a[i][j] / temp;
                y[i][j] = y[i][j] / temp;
            }

            if (i == k)  continue; // уравнение не вычитать само из себя
            for (int j = 0; j < n; j++)
            {
                a[i][j] = a[i][j] - a[k][j];
                y[i][j] = y[i][j] - y[k][j];
            }

        }
        k++;
    }
    // обратная подстановка
    for (k = n - 1; k >= 0; k--)
    {
        for (int i = k; i > 0; i--)
        {
            double temp = a[i-1][k];
            for (int j = 0; j < n; j++)
            {
                y[i - 1][j] = y[i - 1][j] - y[k][j] * temp;
                a[i - 1][j] = a[i - 1][j] - a[k][j] * temp;
            }
        }
    }
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            a[i][j] = y[i][j];
        }
    }
    return;
}

/*!
 \brief Умножение матрицы на вектор
 \details Функция умножает матрицу на вектор, заменяя переданный вектор на
 результат умножения
*/
void lo_multiply(
    std::vector< std::vector<double> > &a,
    std::vector<double> &q
)
{
    const int n = q.size();
    std::vector<double> rez(n);

    for (int i = 0; i < n; i++)
    {
        rez[i] = 0;
        for (int j = 0; j < n; j++)
        {
            rez[i] = rez[i] + a[i][j] * q[j];
        }
    }
    for (int i = 0; i < n; i++)
    {
        q[i] = rez[i];
    }
    return;
}

/*!
 \brief Вычисление коэффициентов линейного оператора
 \details Функция вычисляет коэффициентов линейного оператора
*/
void createLinearOperator(
    std::vector<double> &Z,
    std::vector<double> &Q
){
    Z.insert(std::begin(Z), 0.);
    const int n = Q.size();
    std::vector< std::vector<double> > M(n);
    for (int i = 0; i < n; i++)
    {
        M[i].resize(n);
        M[i][i] = 0;
    }
    
    for (int i = 0; i < n; i++)
    {
        for (int j = i+1; j < n; j++)
        {
            M[i][j]= Z[j+1]- Z[j]+ M[i][j-1];
            M[j][i] = M[i][j];
        }
    }
    
    lo_ober(M);
    lo_multiply(M, Q);
}

/*!
 \brief Вычисление значений линейного оператора
 \details Функция вычисляет значения линейного оператора для заданных
 коэффициентов
*/
double calculateValueWithLinearOperator(
    const double z,
    const std::vector<double> &Z,
    const std::vector<double> &Q
){
    double q = 0.;
    
    for(std::size_t i = 0; i < Q.size(); ++i){
        q += Q[i] * std::abs(z - Z[i + 1]);
    }
    
    return q;
}