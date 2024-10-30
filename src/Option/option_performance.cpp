#include "option.hpp"
#include <algorithm>

OptionPerformance::OptionPerformance(double T, double K, OptionType type, double D, PnlVect *coeff)
    : Option(T, K, type, D, coeff)

{
}

double OptionPerformance::payOff(PnlMat *matrix)
{
    // ================== méthode avce allocation : =====================
    // PnlVect *S_ti = pnl_vect_create(option_size);
    // double sum_d_1 = 0.0;
    // double sum_d_2 = 0.0;
    // double res = 0.0;

    // for (int i = 1; i < matrix->m; i++)
    // {
    //     pnl_mat_get_row(S_ti, matrix, i - 1);
    //     sum_d_1 = pnl_vect_scalar_prod(this->payoff_coeffcients, S_ti);

    //     pnl_mat_get_row(S_ti, matrix, i);
    //     sum_d_2 = pnl_vect_scalar_prod(this->payoff_coeffcients, S_ti);

    //     res += std::max(sum_d_2 / sum_d_1 - 1.0, 0.0);
    // }

    // pnl_vect_free(&S_ti);

    // ================== méthode sans allocation : =====================
    double res = 0;
    for (int i = 1; i < matrix->m; i++)
    {
        double sum1 = 0;
        double sum2 = 0;
        double lamda_d;
        for (int d = 0; d < matrix->n; d++)
        {
            lamda_d = GET(payoff_coeffcients, d);
            sum1 += lamda_d * pnl_mat_get(matrix, i, d);
            sum2 += lamda_d * pnl_mat_get(matrix, i - 1, d);
        }

        res += std::max(sum1 / sum2 - 1, 0.0);
    }

    return 1 + res;
}

OptionPerformance::~OptionPerformance() {}
