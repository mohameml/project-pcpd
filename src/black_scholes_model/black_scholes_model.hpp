#ifndef BLACK_SCHOLES_MODEL_HPP
#define BLACK_SCHOLES_MODEL_HPP
#include "pnl/pnl_matvect.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_random.h"
class BlackScholesModel
{
public:
    double interest_rate;
    PnlVect *volatility; // pnl_vect_free
    PnlVect *spots;      // pnl_vect_free
    double correlation;
    int model_size;
    double time_step;
    PnlMat *L;

public:
    BlackScholesModel();
    BlackScholesModel(double rate, PnlVect *vol, PnlVect *spots, double corr, double time_step);
    ~BlackScholesModel();
    // void asset(double t, PnlVect *Dates, PnlMat *mat_asset);
    void asset(const PnlMat *past, double t, PnlMat *path, PnlRng *rng);
    void asset(PnlMat *path, PnlRng *rng);
    void get_matrix_Cholesky_corr(PnlMat *matrix);
    void shift_asset(int d, double h, PnlMat *original_paths);
    void shift_asset(int d, double t, double h, PnlMat *original_paths);
};
#endif