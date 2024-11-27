#include <cmath>
#include "MonteCarlo.hpp"
#include "pnl/pnl_cdf.h"
#include "compute_last_index.hpp"

MonteCarlo::MonteCarlo()
{
}

MonteCarlo::MonteCarlo(const nlohmann::json json)
{
    json.at("fixing dates number").get_to(fixing_dates_number);
    json.at("sample number").get_to(sample_number);
    json.at("fd step").get_to(fd_step);
    json.at("hedging dates number").get_to(hedging_dates_number);

    model = new BlackScholesModel(json);
    option = instance_option(json);
    rng = pnl_rng_create(PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, time(NULL));
}

MonteCarlo::~MonteCarlo()
{
    delete option;
    delete model;
    pnl_rng_free(&rng);
}

void MonteCarlo::price(double &price, double &price_std)
{
    int D = this->option->size;
    int M = this->sample_number;
    int N = this->fixing_dates_number;
    double r = this->model->interest_rate;
    double T = this->option->maturity;

    double v_t = 0.0;
    double price_std_dev = 0.0;

    PnlMat *matrix = pnl_mat_create(N + 1, D);

    for (int i = 0; i < M; i++)
    {
        this->model->asset(matrix, this->rng);
        double phi_j = this->option->payOff(matrix);
        v_t += phi_j;
        price_std_dev += pow(phi_j, 2);
    }

    double inv_M = 1.0 / (double)M;
    price = std::exp(-r * T) * inv_M * v_t;
    price_std = sqrt(exp(-2 * r * T) * (inv_M * price_std_dev - pow(inv_M * v_t, 2))) / sqrt(M);

    pnl_mat_free(&matrix);
}

void MonteCarlo::price_delta(double &price, double &price_std, PnlVect *deltas_vect, PnlVect *stddev_deltas_vect)
{

    int D = this->option->size;
    int M = this->sample_number;
    int N = this->fixing_dates_number;
    double h = this->fd_step;
    double r = this->model->interest_rate;
    double T = this->option->maturity;

    double v_t = 0.0;
    double price_std_dev = 0.0;

    PnlMat *matrix = pnl_mat_create(N + 1, D);

    for (int i = 0; i < M; i++)
    {
        this->model->asset(matrix, this->rng);
        double phi_j = this->option->payOff(matrix);
        v_t += phi_j;
        price_std_dev += pow(phi_j, 2);

        for (int d = 0; d < D; d++)
        {
            this->model->shift_asset(d, 1 + h, matrix);
            double payoff_plus = this->option->payOff(matrix);
            this->model->shift_asset(d, (1.0 - h) / (1.0 + h), matrix);
            double payoff_minus = this->option->payOff(matrix);
            this->model->shift_asset(d, 1.0 / (1.0 - h), matrix);
            double delta_j = payoff_plus - payoff_minus;
            LET(deltas_vect, d) = GET(deltas_vect, d) + delta_j;
            LET(stddev_deltas_vect, d) = GET(stddev_deltas_vect, d) + delta_j * delta_j;
        }
    }

    double inv_M = 1.0 / (double)M;
    price = std::exp(-r * T) * inv_M * v_t;
    price_std = sqrt(exp(-2 * r * T) * (inv_M * price_std_dev - pow(inv_M * v_t, 2))) / sqrt(M);

    for (int d = 0; d < D; d++)
    {
        double delta_sum = GET(deltas_vect, d) / (double)M;
        LET(deltas_vect, d) = delta_sum * exp(-r * T) / (2.0 * h * GET(model->spots, d));
        double var_delta = (GET(stddev_deltas_vect, d) / M - delta_sum * delta_sum) * exp(-2 * r * T) / pow(2.0 * h * GET(model->spots, d), 2);
        LET(stddev_deltas_vect, d) = sqrt(var_delta) / sqrt(M);
    }

    pnl_mat_free(&matrix);
}

void MonteCarlo::price_delta(double t, double &price, double &price_std, const PnlMat *Past, PnlVect *deltas_vect, PnlVect *stddev_deltas_vect)
{

    int D = this->option->size;
    int M = this->sample_number;
    int N = this->fixing_dates_number;
    double h = this->fd_step;
    double r = this->model->interest_rate;
    double T = this->option->maturity;

    double v_t = 0.0;
    double price_std_dev = 0.0;

    PnlMat *matrix = pnl_mat_create(N + 1, D);

    for (int i = 0; i < M; i++)
    {
        this->model->asset(Past, t, T, matrix, this->rng);
        double phi_j = this->option->payOff(matrix);
        v_t += phi_j;
        price_std_dev += pow(phi_j, 2);

        for (int d = 0; d < D; d++)
        {
            this->model->shift_asset(d, t, 1 + h, matrix);
            double payoff_plus = this->option->payOff(matrix);
            this->model->shift_asset(d, t, (1.0 - h) / (1.0 + h), matrix);
            double payoff_minus = this->option->payOff(matrix);
            this->model->shift_asset(d, t, 1.0 / (1.0 - h), matrix);
            double delta_j = payoff_plus - payoff_minus;
            LET(deltas_vect, d) = GET(deltas_vect, d) + delta_j;
            LET(stddev_deltas_vect, d) = GET(stddev_deltas_vect, d) + delta_j * delta_j;
        }
    }

    double inv_M = 1.0 / (double)M;
    price = std::exp(-r * (T - t)) * inv_M * v_t;
    price_std = sqrt(exp(-2 * r * (T - t)) * (inv_M * price_std_dev - pow(inv_M * v_t, 2))) / sqrt(M);

    for (int d = 0; d < D; d++)
    {
        double delta_sum = GET(deltas_vect, d) / (double)M;
        double s_t_d = MGET(Past, Past->m - 1, d);
        LET(deltas_vect, d) = delta_sum * exp(-r * (T - t)) / (2.0 * h * s_t_d);
        double var_delta = (GET(stddev_deltas_vect, d) / M - delta_sum * delta_sum) * exp(-2 * r * (T - t)) / pow(2.0 * h * s_t_d, 2);
        LET(stddev_deltas_vect, d) = sqrt(var_delta) / sqrt(M);
    }

    pnl_mat_free(&matrix);
}

void MonteCarlo::calculPAndL(PnlMat *market_data, double &p_and_l)
{
    double price = 0.0;
    double price_stdev = 0.0;
    PnlVect *delta = pnl_vect_create_from_zero(model->model_size);       // deltai
    PnlVect *delta_prec = pnl_vect_create_from_zero(model->model_size);  // delta{i-1}
    PnlVect *delta_stdev = pnl_vect_create_from_zero(model->model_size); // std_dev
    double step = option->maturity / (double)hedging_dates_number;       // T/H
    int H_N = hedging_dates_number / fixing_dates_number;                // H/N
    double r = model->interest_rate;
    double riskFreePortfolio = 0.0;

    PnlMat *past = pnl_mat_create(this->fixing_dates_number, option->size);
    pnl_mat_resize(past, 2, option->size);
    pnl_mat_set_row(past, model->spots, 0);

    // traitement de V0
    price_delta(price, price_stdev, delta, delta_stdev);
    pnl_vect_clone(delta_prec, delta);
    riskFreePortfolio = price - pnl_vect_scalar_prod(delta, model->spots);

    for (int i = 1; i < this->hedging_dates_number + 1; i++)
    {
        double tau_i = i * step;
        PnlVect Stau_i = pnl_vect_wrap_mat_row(market_data, i);
        if (i != 0 && (i - 1) % H_N == 0)
        {
            pnl_mat_resize(past, past->m + 1, past->n);
        }

        pnl_mat_set_row(past, &Stau_i, past->m - 1);
        price_delta(tau_i, price, price_stdev, past, delta, delta_stdev);
        riskFreePortfolio = riskFreePortfolio * exp(r * step) - pnl_vect_scalar_prod(delta, &Stau_i) + pnl_vect_scalar_prod(delta_prec, &Stau_i);
        pnl_vect_clone(delta_prec, delta);
    }

    PnlVect ST = pnl_vect_wrap_mat_row(market_data, market_data->m - 1);
    p_and_l = riskFreePortfolio - +pnl_vect_scalar_prod(delta, &ST) - option->payOff(past);

    pnl_mat_free(&past);
    pnl_vect_free(&delta);
    pnl_vect_free(&delta_prec);
    pnl_vect_free(&delta_stdev);
}
