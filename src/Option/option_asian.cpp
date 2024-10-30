#include "option.hpp"
#include <algorithm>

OptionAsian::OptionAsian(double T, double K, OptionType type, double D, PnlVect *coeff)
    : Option(T, K, type, D, coeff)

{
}

double OptionAsian::payOff(PnlMat *matrix)
{
    // ============== méthode avec allocation : ============
    // PnlVect *ligne = pnl_vect_create(option_size);
    // double sum_d = 0.0;

    // for (int i = 0; i < matrix->m; i++)
    // {
    //     pnl_mat_get_row(ligne, matrix, i);
    //     sum_d += pnl_vect_scalar_prod(this->payoff_coeffcients, ligne);
    // }

    // sum_d = sum_d / (double)matrix->m - this->strike;
    // pnl_vect_free(&ligne);

    // return std::max(sum_d, 0.0);

    // ================== méthode sans allocation ==================
    double res = 0.0;
    for (int d = 0; d < matrix->n; d++)
    {
        double res2 = 0;
        for (int i = 0; i < matrix->m; i++)
        {
            res2 += pnl_mat_get(matrix, i, d);
        }
        res += GET(payoff_coeffcients, d) * res2;
    }
    return std::max(res / matrix->m - strike, 0.0);
}
OptionAsian::~OptionAsian() {}
