#include "option.hpp"
#include <cmath>

Option::Option(double T, double K, OptionType type, double D, PnlVect *coeff)
    : maturity(T),
      strike(K),
      type(type),
      option_size(D),
      payoff_coeffcients(coeff)
{
}

Option::Option(const Option &autre)
{
    this->maturity = autre.maturity;
    this->strike = autre.strike;
    this->type = autre.type;
    this->option_size = autre.option_size;
    this->payoff_coeffcients = pnl_vect_copy(autre.payoff_coeffcients);
}

Option::~Option()
{
    pnl_vect_free(&this->payoff_coeffcients);
}

Option *Option::GetOption(double T, double K, OptionType type, double D, PnlVect *coeff)
{
    Option *option = nullptr;
    switch (type)
    {
    case OptionType::Basket:
        option = new OptionBasket(T, K, type, D, coeff);
        break;
    case OptionType::Asian:
        option = new OptionAsian(T, K, type, D, coeff);
        break;
    case OptionType::Performance:
        option = new OptionPerformance(T, K, type, D, coeff);
        break;
    default:
        throw std::invalid_argument("Type d'option non valide : " + type);
        break;
    }

    return option;
}
