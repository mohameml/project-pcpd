#ifndef OPTION_HPP
#define OPTION_HPP

#include <iostream>
#include "option_type.hpp"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"

class Option
{
public:
    double maturity;
    double strike;
    OptionType type;
    double option_size;
    PnlVect *payoff_coeffcients;

public:
    virtual double payOff(PnlMat *matrix) = 0;
    Option(double T, double K, OptionType type, double D, PnlVect *coeff);
    Option(const Option &autre);
    virtual ~Option();
    static Option *GetOption(double T, double K, OptionType type, double D, PnlVect *coeff);
};

class OptionAsian : public Option
{

public:
    OptionAsian(double T, double K, OptionType type, double D, PnlVect *coeff);
    double payOff(PnlMat *matrix) override;
    ~OptionAsian() override;
};

class OptionBasket : public Option
{

public:
    OptionBasket(double T, double K, OptionType type, double D, PnlVect *coeff);
    double payOff(PnlMat *matrix) override;
    ~OptionBasket() override;
};

class OptionPerformance : public Option
{

public:
    OptionPerformance(double T, double K, OptionType type, double D, PnlVect *coeff);
    double payOff(PnlMat *matrix) override;
    ~OptionPerformance() override;
};

#endif