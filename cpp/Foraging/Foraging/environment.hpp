#include <vector>
#include <stdexcept>

using namespace std;


// A base class for a target of the foragers
class Target
{	
public:
	virtual double Payout() {return 0.0;};
};


// A constant source target, always providing the same reward
class ConstantTarget : public Target
{
private:
	double value;

public:
	ConstantTarget(double val);
	double Payout();
};
/*************************************/

// A binary target, producing either a reward or 0.0 according to some probability
class BinaryTarget : public Target
{
private:
	double prob;
	double reward;
public:
	BinaryTarget(double p, double r);
	double Payout();
};


// A discrete target, providing several possible rewards with different probabilities
class DiscreteTarget : public Target
{
private:
	std::vector<double> probs;
	std::vector<double> rewards;

public:
	DiscreteTarget(std::vector<double> probs, std::vector<double> rewards);
	double Payout();
};
/********************************************************************************/




// Simple target grid
class Field
{
private:
	Target* targets;
};
/**************************************************/



// A simple field with constant, binary, and neutral targets, as
// described in Niv, et. al. (2000)
class BasicField : public Field
{
public:
	BasicField(double certVal, double uncertVal, double uncertProb, double pctCert, double pctUncert);
};
/**********************************************************/