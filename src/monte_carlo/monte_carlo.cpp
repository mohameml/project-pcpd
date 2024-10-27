#include "monte_carlo.hpp"
#include <cmath>
#include "../Utils/convert.hpp"
// #include "../Utils/construct_append_mat.hpp"
#include "../Utils/utils.hpp"
#include "pnl/pnl_cdf.h"

MonteCarlo::MonteCarlo(Option *option, BlackScholesModel *model, int N, int M, double H, double h)
    : option(option),
      model(model),
      fixing_dates_number(N),
      sample_number(M),
      hedging_dates_number(H),
      fd_step(h)
{
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

    int D = this->option->option_size;
    double r = this->model->interest_rate;
    int M = this->sample_number;
    int N = this->fixing_dates_number;
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

void MonteCarlo::price(double t, double &price, double &price_std, const PnlMat *Past)
{

    int D = this->option->option_size;
    double r = this->model->interest_rate;
    int M = this->sample_number;
    int N = this->fixing_dates_number;
    double T = this->option->maturity;

    double v_t = 0.0;
    double price_std_dev = 0.0;

    PnlMat *matrix = pnl_mat_create(N + 1, D);

    for (int i = 0; i < M; i++)
    {
        this->model->asset(Past, t, matrix, this->rng);
        double phi_j = this->option->payOff(matrix);
        v_t += phi_j;
        price_std_dev += pow(phi_j, 2);
    }

    double inv_M = 1.0 / (double)M;

    price = std::exp(-r * (T - t)) * inv_M * v_t;

    price_std = sqrt(exp(-2 * r * T) * (inv_M * price_std_dev - pow(inv_M * v_t, 2))) / sqrt(M);
    pnl_mat_free(&matrix);
}

void MonteCarlo::calculPAndL(PnlMat *market_data, double &p_and_l)
{
    double step = option->maturity / (double)this->hedging_dates_number; // step = T/H
    double r = model->interest_rate;
    double p;
    double price_stdev;
    double v;

    PnlMat *past = pnl_mat_new();
    PnlVect *St_i = pnl_vect_new();
    PnlVect *delta = pnl_vect_create(model->model_size);       // deltai
    PnlVect *delta_stdev = pnl_vect_create(model->model_size); // std_dev
    PnlVect *delta_1 = pnl_vect_create(model->model_size);     // delta{i+1}

    //
    // traitement de ti = 0
    price(p, price_stdev);
    // deltaCall(PnlVect *St_i, double t, PnlVect *deltan)
    this->delta(delta, delta_stdev);
    // pnl_vect_print(delta);

    v = p - pnl_vect_scalar_prod(delta, model->spots);

    for (int i = 1; i < this->hedging_dates_number + 1; i++)
    {
        double t = i * step;
        pnl_mat_get_row(St_i, market_data, i); // sti

        // delta_i
        get_cotations(t, past, market_data);
        // std::cout<<"------------------------------------------------------------------------"<< i <<std::endl;
        // deltaCall(St_i,t, delta_1);
        // std::cout<<"----------------------------delta call-----------------"<<std::endl;
        // pnl_vect_print(delta_1);

        this->delta(past, delta_1, delta_stdev, t);
        // std::cout<<"-----------------delta calcul-----------------"<<std::endl;
        // pnl_vect_print(delta_1);
        pnl_vect_minus_vect(delta, delta_1); // delat <- delta{i-1} - delta{i}
        v = v * exp(r * step) + pnl_vect_scalar_prod(delta, St_i);
        pnl_vect_clone(delta, delta_1);
    }

    // l'erreur de couverture
    int N = this->fixing_dates_number;
    int D = this->option->option_size;
    PnlMat *matrix = pnl_mat_create(N + 1, D);
    pnl_mat_extract_subblock(matrix, past, 0, N + 1, 0, D);
    
    p_and_l = v + pnl_vect_scalar_prod(delta_1, St_i) - option->payOff(matrix);

    // free
    pnl_mat_free(&past);
    pnl_vect_free(&St_i);
    pnl_vect_free(&delta);
    pnl_vect_free(&delta_1);
    pnl_vect_free(&delta_stdev);
}

void MonteCarlo::delta(PnlVect *deltas_vect, PnlVect *stddev_deltas_vect)
{
    // delta of the option at time t = 0
    int D = this->option->option_size;
    int M = this->sample_number;
    int N = this->fixing_dates_number;
    double h = this->fd_step;
    double r = this->model->interest_rate;
    double T = this->option->maturity;

    PnlMat *mat_asset_plus = pnl_mat_create(N + 1, D);
    PnlMat *mat_asset_minus = pnl_mat_create(N + 1, D);

    // Loop over assets
    for (int d = 0; d < D; d++)
    {
        double delta_sum = 0.0;
        double delta_square_sum = 0.0;
        // Loop over Monte Carlo simulations
        for (int j = 0; j < M; j++)
        {
            this->model->asset(mat_asset_plus, this->rng);
            pnl_mat_clone(mat_asset_minus, mat_asset_plus);
            // pnl_mat_get_row(st_0, mat_asset_plus, 0);
            this->model->shift_asset(d, h, mat_asset_plus);
            this->model->shift_asset(d, -h, mat_asset_minus);

            double payoff_plus = this->option->payOff(mat_asset_plus);
            double payoff_minus = this->option->payOff(mat_asset_minus);

            double delta_j = (payoff_plus - payoff_minus);

            // Accumulate delta and squared delta for variance calculation
            delta_sum += delta_j;
            delta_square_sum += delta_j * delta_j;
        };
        // pnl_mat_get_row(st_0, mat_asset_plus, 0);
        double delta_mean = delta_sum * exp(-r * T) / ((2.0 * h) * M * pnl_vect_get(model->spots, d));
        pnl_vect_set(deltas_vect, d, delta_mean);

        // Calculate the standard deviation of delta
        double delta_var = (delta_square_sum / M - pow(delta_sum / M, 2)) * exp(-2 * r * T) / pow(2.0 * h * pnl_vect_get(model->spots, d), 2);
        double delta_std_dev = sqrt(delta_var) / sqrt(M);
        pnl_vect_set(stddev_deltas_vect, d, delta_std_dev);
    };

    // Free allocated memory
    pnl_mat_free(&mat_asset_plus);
    pnl_mat_free(&mat_asset_minus);
}


void MonteCarlo::deltaCall(PnlVect *St_i, double t, PnlVect *delta)
{
    
    double inf = (GET(model->volatility,0) * sqrt(option->maturity - t));
    double d1 = (log(GET(St_i,0)/option->strike) + (model->interest_rate + pow(GET(model->volatility,0),2)/2)*(option->maturity - t))/inf;
    pnl_vect_set(delta,0,pnl_cdfnor(d1));
    
}

void MonteCarlo::delta(PnlMat *past, PnlVect *deltas_vect, PnlVect *stddev_deltas_vect, double t)
{
    /*
    Usage of the function
    BEFORE CALLING THE FUNCTION
    PnlVect* deltas_vect = pnl_vect_new();
    PnlVect* stddev_deltas_vect = pnl_vect_new();

    delta(deltas_vect);

    AFTER CALLING THE FUNCTION
    pnl_vect_free(&deltas_vect);
    pnl_vect_free(&stddev_deltas_vect);
    */
    int D = this->option->option_size; // Nb assets
    int M = this->sample_number;       // Nb Monte Carlo simulations
    int N = this->fixing_dates_number; // Nb fixing dates
    double h = this->fd_step;          // Finite difference step size
    double r = this->model->interest_rate;
    double T = this->option->maturity;

    // last index based on time t
    // int index = compute_last_index(t, T, N);

    // Create matrices for shifted assets
    PnlMat *mat_asset_plus = pnl_mat_create(N + 1, D);
    PnlMat *mat_asset_minus = pnl_mat_create(N + 1, D);
    PnlVect *st_i = pnl_vect_create(D); // Holds the asset values at time t

    // Temporary vectors to store sums and squared sums for each delta

    pnl_mat_get_row(st_i, past, past->m -1);

    // Loop over assets
    for (int d = 0; d < D; d++)
    {

        double delta_sum = 0.0;
        double delta_square_sum = 0.0;
        // Loop over Monte Carlo simulations
        for (int j = 0; j < M; j++)
        {
            this->model->asset(past, t, mat_asset_plus, this->rng);
            pnl_mat_clone(mat_asset_minus, mat_asset_plus);
            // pnl_mat_get_row(st_i, mat_asset_plus, index);
            // Shift the asset for finite difference calculation
            this->model->shift_asset(d, t, h, mat_asset_plus);
            this->model->shift_asset(d, t, -h, mat_asset_minus);

            double payoff_plus = this->option->payOff(mat_asset_plus);
            double payoff_minus = this->option->payOff(mat_asset_minus);

            // Calculate the delta for the current simulation
            double delta_j = (payoff_plus - payoff_minus);

            // Accumulate the sum of deltas and squared deltas
            delta_sum += delta_j;
            delta_square_sum += delta_j * delta_j;
        }

        // mean delta
        double delta_mean = delta_sum * exp(-r * (T - t)) / (2.0 * h * M * pnl_vect_get(st_i, d));
        pnl_vect_set(deltas_vect, d, delta_mean);

        // variance and standard deviation
        double delta_var = (delta_square_sum / M - pow(delta_sum / M, 2)) * exp(-2 * r * (T - t)) / pow(2.0 * h * pnl_vect_get(st_i, d), 2);
        double delta_std_dev = sqrt(delta_var) / sqrt((double)M);
        pnl_vect_set(stddev_deltas_vect, d, delta_std_dev);
    }


    // Free memory
    pnl_mat_free(&mat_asset_plus);
    pnl_mat_free(&mat_asset_minus);
    pnl_vect_free(&st_i);

}

void MonteCarlo::get_cotations(double t, PnlMat *cots, PnlMat *market_data)
{
    /*
    To use this function there are some steps to do :

    Step 1 :
    PnlVect* s_t = pnl_vect_new();
    PnlMat* cots = pnl_mat_new();

    get_cotations(t , cots , s_t);

    pnl_vect_free(&s_t);
    pnl_mat_free(&cots);
  */

    if (market_data == NULL)
    {
        throw std::invalid_argument("argument data.txt non fourni");
        exit(1);
    }

    int H = this->hedging_dates_number;
    int N = fixing_dates_number;
    double T = option->maturity;
    int D = option->option_size;

    int i = compute_last_index(t, T, N);

    pnl_mat_resize(cots, i + 2, D);

    PnlVect *col = pnl_vect_create(D);

    for (int j = 0; j < i + 1; j++)
    {
        pnl_mat_get_row(col, market_data, j * H / N);
        pnl_mat_set_row(cots, col, j);
    }

    // s_t
    pnl_mat_get_row(col, market_data, t * H / T);
    pnl_mat_set_row(cots, col, i + 1);

    // free
    pnl_vect_free(&col);
}
