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
    this->L = pnl_mat_create(model_size, model_size);
    pnl_mat_set_all(L, this->correlation);

    for (int i = 0; i < this->model_size; i++)
        pnl_mat_set_diag(L, 1.0, i);

    pnl_mat_chol(L);
}

BlackScholesModel::~BlackScholesModel()
{
    pnl_vect_free(&volatility);
    pnl_vect_free(&spots);
    pnl_mat_free(&L);
}

void BlackScholesModel::asset(PnlMat *path, PnlRng *rng)
{
    int D = this->model_size;
    double r = this->interest_rate;

    PnlVect *L_d = pnl_vect_create(D);
    PnlVect *G = pnl_vect_create(D);

    pnl_mat_set_row(path, spots, 0);

    for (int i = 1; i < path->m; i++)
    {

        pnl_vect_rng_normal(G, D, rng);

        for (int d = 0; d < D; d++)
        {
            double sigma_d = pnl_vect_get(this->volatility, d);
            pnl_mat_get_row(L_d, this->L, d);
            double s_t_i = MGET(path, i - 1, d);
            MLET(path, i, d) = s_t_i * exp((r - pow(sigma_d, 2) / 2.0) * (time_step) + sigma_d * sqrt(time_step) * pnl_vect_scalar_prod(L_d, G));
        }
    }

    pnl_vect_free(&L_d);
    pnl_vect_free(&G);
}

void BlackScholesModel::asset(const PnlMat *past, double t, double T, PnlMat *path, PnlRng *rng)
{
    int last_index = compute_last_index(t, T, path->m - 1);
    int D = this->model_size;
    double r = this->interest_rate;

    bool is_t_equals_t_i = past->m == last_index + 1;

    pnl_mat_set_subblock(path, past, 0, 0);

    PnlVect *L_d = pnl_vect_create(D);
    PnlVect *G = pnl_vect_create(D);

    for (int i = last_index + 1; i < path->m; i++)
    {
        pnl_vect_rng_normal(G, D, rng);

        for (int d = 0; d < D; d++)
        {
            double s_t_d = i == last_index + 1 && !is_t_equals_t_i ? MGET(past, last_index + 1, d) : MGET(path, i - 1, d);
            double dt = i == last_index + 1 && !is_t_equals_t_i ? (last_index + 1) * time_step - t : time_step;
            double sigma_d = pnl_vect_get(this->volatility, d);
            pnl_mat_get_row(L_d, L, d);
            MLET(path, i, d) = s_t_d * exp((r - pow(sigma_d, 2) / 2.0) * dt + sigma_d * sqrt(dt) * pnl_vect_scalar_prod(L_d, G));
        }
    }

    pnl_vect_free(&L_d);
    pnl_vect_free(&G);
}

void BlackScholesModel::shift_asset(int d, double h, PnlMat *original_paths)
{
    for (int i = 1; i < original_paths->m; i++)
    {
        pnl_mat_set(original_paths, i, d, pnl_mat_get(original_paths, i, d) * h);
    };
}

void BlackScholesModel::shift_asset(int d, double t, double h, PnlMat *original_paths)
{
    int nb_lines = original_paths->m;
    double T = this->time_step * (nb_lines - 1);
    int index = compute_last_index(t, T, nb_lines - 1);
    for (int i = index + 1; i < nb_lines; i++)
    {
        pnl_mat_set(original_paths, i, d, pnl_mat_get(original_paths, i, d) * h);
    };
}
