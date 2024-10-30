#include "option.hpp"
#include <algorithm>
OptionBasket::OptionBasket(double T, double K, OptionType type, double D, PnlVect *coeff)
    : Option(T, K, type, D, coeff)

{
}

double OptionBasket::payOff(PnlMat *matrix)
{
    int D = this->option_size;

    PnlVect *S_T = pnl_vect_create(D);
    pnl_mat_get_row(S_T, matrix, matrix->m - 1);

    // double res = 0;
    // for (int d = 0; d < this->option_size; d++)
    // {
    //     res += MGET(matrix, matrix->m - 1, d) * GET(this->payoff_coeffcients, d);
    // }

    double res = pnl_vect_scalar_prod(S_T, this->payoff_coeffcients) - this->strike;
    // free
    pnl_vect_free(&S_T);

    return std::max(0.0, res);
}

OptionBasket::~OptionBasket()
{
}
