#include <cmath>
#include <ctime>

#include <environment.hpp>
#include <utils.hpp>

using namespace std;


// Implementation of the discrete target
DiscreteTarget::DiscreteTarget(int N, double* probs, double* rewards)
{
	for (int i=0; i<N; i++)
		if (probs[i] <= 0.0)
			throw Exception("Error: all probabilities must be positive");

	// Set up the random number generator
//	const gsl_rng_type* T = gsl_rng_default;
	_randGen = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(_randGen, (unsigned long)time(NULL));

	// Set up the discrete distribution
	_discretePdf = gsl_ran_discrete_preproc(N, probs);

	// Set up the rewards	
	_values = gsl_vector_alloc(N);
	for (int i=0; i<N; i++)
		gsl_vector_set(_values, i, rewards[i]);
}

DiscreteTarget::~DiscreteTarget()
{
	gsl_rng_free(_randGen);
	gsl_ran_discrete_free(_discretePdf);
	gsl_vector_free(_values);
}

double DiscreteTarget::Payout()
{
	size_t index = gsl_ran_discrete(_randGen, _discretePdf);
	return gsl_vector_get(_values, index);
}
/*************************************************************/