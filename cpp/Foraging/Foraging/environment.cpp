#include <cmath>
#include <vector>
#include <cstdlib>
#include <ctime>

#include <environment.hpp>
#include <utils.hpp>

using namespace std;


// Implementation of the constant/deterministic target
ConstantTarget::ConstantTarget(double val)
{
	this->value = val;
}

double ConstantTarget::Payout()
{
	return this->value;
}
/***************************************************/


// Implementation of binary target
BinaryTarget::BinaryTarget(double p, double r)
{
	if (0.0 > p || p > 1.0)
		throw Exception("Probability outside [0,1]");

	this->prob = p;
	this->reward = r;
}

double BinaryTarget::Payout()
{
	srand((unsigned int)time(NULL));
	double value = rand() / (double(RAND_MAX)+1);
	if (value <= this->prob)
		return this->reward;
	else
		return 0.0;
}
/***************************************************************/


// Implementation of the discrete target
DiscreteTarget::DiscreteTarget(std::vector<double> probs, std::vector<double> rewards)
{
	if (probs.size() != rewards.size())
		throw Exception("Num rewards not equal to num probabilities");

	int N = probs.size();
	this->probs[0] = probs[0];
	for (int i=1; i<N; i++)
		this->probs[i] = this->probs[i-1] + probs[i];

	if (abs(1.0 - this->probs[N-1]) < 1e-6)    // Use epsilon threshold to check for equality, b/c dealing with floats
		this->probs[N-1] = 1.0;  // Adjust the probabilities so that they actually sum to 1.0
	else
		throw Exception("Probabilities don't sum to 1");
	
	this->rewards = rewards;
}

double DiscreteTarget::Payout()
{
	srand((unsigned int)time(NULL));
	double value = rand() / (double(RAND_MAX)+1);
	int index = 0;
	while (value < this->probs[index])
		index++;
	return this->rewards[index - 1];
}
/*************************************************************/