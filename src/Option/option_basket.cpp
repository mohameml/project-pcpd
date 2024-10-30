#include "option.hpp"
#include <algorithm>
OptionBasket::OptionBasket(double T, double K, OptionType type, double D, PnlVect *coeff)
    : Option(T, K, type, D, coeff)

{
}

double OptionBasket::payOff(PnlMat *matrix)
{
    // ====================== méthode avec allocation : ==============
    // PnlVect *S_T = pnl_vect_create(option_size);
    // pnl_mat_get_row(S_T, matrix, matrix->m - 1);
    // double res = pnl_vect_scalar_prod(S_T, this->payoff_coeffcients) - this->strike;
    // pnl_vect_free(&S_T);
    // return std::max(res, 0.0);

    // ====================== méthode sans allocation : ==============
    double res = 0;
    for (int d = 0; d < this->option_size; d++)
    {
        res += MGET(matrix, matrix->m - 1, d) * GET(this->payoff_coeffcients, d);
    }

    return std::max(res - strike, 0.0);
}

OptionBasket::~OptionBasket()
{
}
