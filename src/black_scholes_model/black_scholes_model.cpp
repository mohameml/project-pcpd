#include "black_scholes_model.hpp"
#include <cmath>
#include <random>
#include "../Utils/utils.hpp"
#include <iostream>
#include <cassert>


BlackScholesModel::BlackScholesModel()
{
    this->volatility = pnl_vect_new();
    this->spots = pnl_vect_new();
}

BlackScholesModel::BlackScholesModel(double rate, PnlVect *vol, PnlVect *spots, double corr, double time_step)
    : interest_rate(rate),
      volatility(vol),
      spots(spots),
      correlation(corr),
      time_step(time_step)
{
    this->model_size = spots->size;
    this->L = pnl_mat_new();
    get_matrix_Cholesky_corr(L);
}

BlackScholesModel::~BlackScholesModel()
{
    pnl_vect_free(&volatility);
    pnl_vect_free(&spots);
    pnl_mat_free(&L);
}

void BlackScholesModel::asset(const PnlMat *past, double t, PnlMat *path, PnlRng *rng)
{
    int n = path->m; // n = N+1
    int D = this->model_size;
    // matrice L :

    PnlVect *L_d = pnl_vect_create(D);
    PnlVect *G = pnl_vect_create(D);
    // t>0
    // col contient le vecteur s_t
    PnlVect *s_t = pnl_vect_create(D);
    pnl_mat_get_row(s_t, past, past->m - 1);

    int last_index_t = past->m - 2;

    if(past->m >path->m)
    {
        PnlMat *sub_block = pnl_mat_new();
        pnl_mat_extract_subblock(sub_block, past, 0, last_index_t + 1, 0, D);
        pnl_mat_set_subblock(path, sub_block, 0, 0);
        pnl_mat_free(&sub_block);

    }
    else{
        pnl_mat_set_subblock(path, past, 0, 0);
    }

    // pnl_mat_set_subblock(path, past, 0, 0);
    // calcul de St_i+1
    double r = this->interest_rate;

    PnlVect *calcul = pnl_vect_create(D);


    pnl_vect_rng_normal(G, D, rng);

    // timstep0 = 
    double timstep0 = (time_step)*(last_index_t+1) -t;
    for (int d = 0; d < D; d++)
    {
        double s_t_d = pnl_vect_get(s_t, d);
        double sigma_d = pnl_vect_get(this->volatility, d);
        pnl_mat_get_row(L_d, L, d);

        double x =  s_t_d * exp((r - pow(sigma_d, 2) / 2.0) * (timstep0) + sigma_d * sqrt(timstep0) * pnl_vect_scalar_prod(L_d, G)); // t{i+1} - t{i} = T/N
        pnl_vect_set(calcul, d, x);

    }

    if (last_index_t + 1 != n ) {

        pnl_mat_set_row(path, calcul, last_index_t+1);
    }




    for (int ligne = last_index_t + 2; ligne < n; ligne++)
    {
        pnl_vect_rng_normal(G, D, rng);


        for (int d = 0; d < D; d++)
        {
            double s_t_d = pnl_vect_get(calcul, d);
            double sigma_d = pnl_vect_get(this->volatility, d);
            pnl_mat_get_row(L_d, L, d);

            double x =  s_t_d * exp((r - pow(sigma_d, 2) / 2.0) * (time_step) + sigma_d * sqrt(time_step) * pnl_vect_scalar_prod(L_d, G)); // t{i+1} - t{i} = T/N
            pnl_vect_set(calcul, d, x);

        }

        pnl_mat_set_row(path, calcul, ligne);
    }

    // free
    pnl_vect_free(&calcul);
    pnl_vect_free(&s_t);
    pnl_vect_free(&L_d);
    pnl_vect_free(&G);
}

// pour le cas t = 0
void BlackScholesModel::asset(PnlMat *path, PnlRng *rng)
{
    int n = path->m; // n = N+1
    int D = this->model_size;
    // matrice L :
    PnlMat *L = pnl_mat_new();
    get_matrix_Cholesky_corr(L);

    PnlVect *L_d = pnl_vect_create(D);
    PnlVect *G = pnl_vect_create(D);
    // t>0
    // col contient le vecteur s_0
    PnlVect *col = pnl_vect_create(D);
    pnl_vect_clone(col, spots);

    pnl_mat_set_row(path, spots, 0);

    // calcul de St_i+1

    double r = this->interest_rate;

    for (int ligne = 1; ligne < n; ligne++)
    {

        pnl_vect_rng_normal(G, D, rng);

        for (int d = 0; d < D; d++)
        {
            double s_t_i = pnl_vect_get(col, d);
            double sigma_d = pnl_vect_get(this->volatility, d);
            pnl_mat_get_row(L_d, L, d);

            double x = s_t_i * exp((r - pow(sigma_d, 2) / 2.0) * (time_step) + sigma_d * sqrt(time_step) * pnl_vect_scalar_prod(L_d, G)); // t{i+1} - t{i} = T/N
            pnl_vect_set(col, d, x);
        }

        pnl_mat_set_row(path, col, ligne);
    }

    // free
    pnl_vect_free(&col);
    pnl_mat_free(&L);
    pnl_vect_free(&L_d);
    pnl_vect_free(&G);
}

void BlackScholesModel::get_matrix_Cholesky_corr(PnlMat *matrix_chol)
{
    /*
    To use this function there are some steps to do :
    PnlMat* matrix = pnl_mat_new();
    get_matrix_Cholesky_corr(matrix);
    pnl_mat_free(&matrix);

    */
    pnl_mat_resize(matrix_chol, this->model_size, this->model_size);
    pnl_mat_set_all(matrix_chol, this->correlation);

    for (int i = 0; i < this->model_size; i++)
        pnl_mat_set_diag(matrix_chol, 1.0, i);

    pnl_mat_chol(matrix_chol);
}

void BlackScholesModel::shift_asset(int d, double h, PnlMat *original_paths)
{
    int nb_lines = original_paths->m;

    for (int i = 1; i < nb_lines; i++)
    {
        pnl_mat_set(original_paths, i, d, pnl_mat_get(original_paths, i, d) * (1 + h));
    };
}

void BlackScholesModel::shift_asset(int d, double t, double h, PnlMat *original_paths)
{
    int nb_lines = original_paths->m;
    double T = this->time_step * (nb_lines-1);
    int index = compute_last_index(t, T, nb_lines -1);
    for (int i = index + 1; i < nb_lines; i++)
    {
        pnl_mat_set(original_paths, i, d, pnl_mat_get(original_paths, i, d) * (1 + h));
    };
}
