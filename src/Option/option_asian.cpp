#include "option.hpp"
#include <algorithm>

OptionAsian::OptionAsian(double T, double K, OptionType type, double D, PnlVect *coeff)
    : Option(T, K, type, D, coeff)

{
}

double OptionAsian::payOff(PnlMat *matrix)
{
    int rows = matrix->m; // N+1
    int cols = matrix->n; // D
    PnlVect *ligne = pnl_vect_create(cols);
    double sum_d = 0.0;

    for(int i=0; i<rows; i++)
    {   
        pnl_mat_get_row(ligne,matrix,i);
        sum_d +=  pnl_vect_scalar_prod(this->payoff_coeffcients,ligne);
    }

    sum_d = sum_d/(double)rows - this->strike;
    double zero = 0.0;
    // free
    pnl_vect_free(&ligne);
    return std::max(zero, sum_d);
}
OptionAsian::~OptionAsian() {}



    // for (int i = 0; d < cols; d++)
    // {
    //     double lamda_d = GET(payoff_coeffcients, d);

        
        
    //     res += lamda_d * pnl_vect_sum();
        
        
        
    //     double res2 = 0;
    //     for (int j = 0; j < cols; j++)
    //     {
    //         res2 += pnl_mat_get(matrix, d, j);
    //     }

    // }
    // res = pnl_mat_sum(matrix); 