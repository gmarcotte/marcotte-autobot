#include <gsl/gsl_rng.h>

#include <controller.hpp>
#include <utils.hpp>

#include <ctime>
#include <cmath>

using namespace std;



/* Implementation of a controller that changes directions randomly */
RandomController::RandomController(double prob)
{
	if (prob < 0.0 || prob > 1.0)
		throw Exception("Error: Probability must be in [0,1]");

	_prob = prob;
}

RandomController::~RandomController()
{
}

bool RandomController::ChangeDirection()
{
	double rnd = gsl_rng_uniform(get_rng());
	return (rnd < _prob);
}
/**********************************************************************/



/* Implementation of a basic NN controller as described in Niv et al. */

// Default constructor initializes randomly
NivController::NivController()
{
	// Synapse presence
	_useBlue = rand_bool();
	_useRed = rand_bool();
	_useGray = rand_bool();
	_useDiffBlue = rand_bool();
	_useDiffRed = rand_bool();
	_useDiffGray = rand_bool();
	_useReward = rand_bool();

	// Synapse weights
	_redWt = rand_double_uniform(-1.0, 1.0);
	_blueWt = rand_double_uniform(-1.0, 1.0);
	_grayWt = rand_double_uniform(-1.0, 1.0);
	_redDiffWt = rand_double_uniform(-1.0, 1.0);
	_blueDiffWt = rand_double_uniform(-1.0, 1.0);
	_grayDiffWt = rand_double_uniform(-1.0, 1.0);
	_rewardWt = 1.0; // Clamped

	// Action function parameters
	_actionM = rand_double_uniform(5.0, 45.0);
	_actionB = rand_double_uniform(0.0, 5.0);

	// Regular module learning rule
	_regLearnA = rand_double_uniform(-0.2, 0.2);
	_regLearnB = rand_double_uniform(-0.2, 0.2);
	_regLearnC = rand_double_uniform(-0.2, 0.2);
	_regLearnD = rand_double_uniform(-0.2, 0.2);

	// Differential module learning rule
	_diffLearnA = rand_double_uniform(-0.2, 0.2);
	_diffLearnB = rand_double_uniform(-0.2, 0.2);
	_diffLearnC = rand_double_uniform(-0.2, 0.2);
	_diffLearnD = rand_double_uniform(-0.2, 0.2);

	// Global learning rule
	_learnRate = rand_double_uniform(0.0, 1.0);

	// Modular learning dependencies
	_modDiffOnReg = rand_bool();
	_modRegOnDiff = rand_bool();
	_modDiffOnReward = rand_bool();
	_modRegOnReward = rand_bool();

	// State tracking
	_prevRed = 0.0;
	_prevBlue = 0.0;
	_prevGray = 0.0;
}

NivController::~NivController()
{
}

void NivController::ResetState()
{
	_prevRed = 0.0;
	_prevBlue = 0.0;
	_prevGray = 0.0;
}

double NivController::EvaluateNetwork(double red, double blue, double gray, double d_red, double d_blue, double d_gray, double reward)
{
	double out = reward;
	out += _redWt*red + _blueWt*blue + _grayWt*gray; // Regular
	out += _redDiffWt*d_red + _blueDiffWt*d_blue + _grayDiffWt*d_gray; // Differential
	return out;
}


double NivController::Train(double red, double blue, double gray, double d_red, double d_blue, double d_gray, double reward)
{
	double P = EvaluateNetwork(red, blue, gray, d_red, d_blue, d_gray, reward);
	
	// Update regular module weights
	if ( (!_modRegOnDiff || (d_red != 0.0)) && (!_modRegOnReward || (reward != 0.0)) )
		_redWt += _learnRate * (_regLearnA*red*P + _regLearnB*red + _regLearnC*P + _regLearnD);
	if ( (!_modRegOnDiff || (d_blue != 0.0)) && (!_modRegOnReward || (reward != 0.0)) )
		_blueWt += _learnRate * (_regLearnA*blue*P + _regLearnB*blue + _regLearnC*P + _regLearnD);
	if ( (!_modRegOnDiff || (d_gray != 0.0)) && (!_modRegOnReward || (reward != 0.0)) )
		_grayWt += _learnRate * (_regLearnA*gray*P + _regLearnB*gray + _regLearnC*P + _regLearnD);

	// Update differential module weights
	if ( (!_modDiffOnReg || (red != 0.0)) && (!_modDiffOnReward || (reward != 0.0)) )
		_redDiffWt += _learnRate * (_diffLearnA*d_red*P + _diffLearnB*d_red + _diffLearnC*P + _diffLearnD);

	if ( (!_modDiffOnReg || (blue != 0.0)) && (!_modDiffOnReward || (reward != 0.0)) )
		_blueDiffWt += _learnRate * (_diffLearnA*d_blue*P + _diffLearnB*d_blue + _diffLearnC*P + _diffLearnD);

	if ( (!_modDiffOnReg || (gray != 0.0)) && (!_modDiffOnReward || (reward != 0.0)) )
		_grayDiffWt += _learnRate * (_diffLearnA*d_gray*P + _diffLearnB*d_gray + _diffLearnC*P + _diffLearnD);

	return P;
}


bool NivController::ChangeDirection(double red_pct, double blue_pct, double gray_pct, double reward)
{
	double P = Train(red_pct, blue_pct, gray_pct, red_pct - _prevRed, blue_pct - _prevBlue, gray_pct - _prevGray, reward);
	_prevRed = red_pct;
	_prevBlue = blue_pct;
	_prevGray = gray_pct;

	double p_change = 1.0 / (1.0 + exp(_actionM*P + _actionB));
	double val = rand_double_uniform(0.0, 1.0);
	return (val < p_change);
}

/**********************************************************************/