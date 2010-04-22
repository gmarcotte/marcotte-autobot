
#ifndef ENVIRONMENT_HPP
#define ENVIRONMENT_HPP

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>

#include <stdexcept>

using namespace std;


// A constant source target, always providing the same reward
class ConstantTarget
{
private:
	double _value;

public:
	ConstantTarget(double val) {_value = val;}
	double Payout() {return _value;};
};
/*************************************/

// A discrete target, providing several possible rewards with different probabilities
class DiscreteTarget
{
private:
	gsl_rng* _randGen;
	gsl_ran_discrete_t* _discretePdf;
	gsl_vector* _values;

public:
	DiscreteTarget(int N, double* probs, double* rewards);
	~DiscreteTarget();
	double Payout();
};
/********************************************************************************/

#endif // ENVIRONMENT_HPP