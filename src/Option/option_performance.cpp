#include "option.hpp"
#include <algorithm>

OptionPerformance::OptionPerformance(double T, double K, OptionType type, double D, PnlVect *coeff)
    : Option(T, K, type, D, coeff)

{
}

double OptionPerformance::payOff(PnlMat *matrix)
{
    int rows = matrix->m; // N+1
    int cols = matrix->n; // D
    PnlVect *S_ti = pnl_vect_create(cols);
    // PnlVect *S2_ti = pnl_vect_create(cols);
    double sum_d_1 = 0.0;
    double sum_d_2 = 0.0;
    double res = 0.0;


    for(int i=1; i<rows; i++)
    {   
        pnl_mat_get_row(S_ti,matrix,i-1);    //
        sum_d_1 =  pnl_vect_scalar_prod(this->payoff_coeffcients,S_ti);

        pnl_mat_get_row(S_ti, matrix, i); //
        sum_d_2 = pnl_vect_scalar_prod(this->payoff_coeffcients, S_ti);

        res += std::max(sum_d_2/sum_d_1 - 1.0, 0.0);
    }

    // free
    pnl_vect_free(&S_ti);
    // pnl_vect_free(&S2_ti);

    return 1 + res;
}

OptionPerformance::~OptionPerformance() {}

// double res = 0;
// int rows = matrix->m; // D
// int cols = matrix->n; // N+1

// for (int j = 1; j < cols; j++)
// {
//     double sum1 = 0;
//     double sum2 = 0;
//     double lamda_d;
//     for (int d = 0; d < rows; d++)
//     {
//         lamda_d = GET(payoff_coeffcients, d);
//         sum1 += lamda_d * pnl_mat_get(matrix, d, j);
//         sum2 += lamda_d * pnl_mat_get(matrix, d, j - 1);
//     }

//     res += std::max(sum1 / sum2 - 1, 0.0);
// }