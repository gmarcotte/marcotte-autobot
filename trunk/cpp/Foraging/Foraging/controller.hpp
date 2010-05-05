#ifndef CONTROLLER_HPP
#define CONTROLLER_HPP

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>

class Controller
{
protected:

public:
	Controller() {};
	virtual ~Controller() {};
	virtual bool ChangeDirection(double red_pct, double blue_pct, double gray_pct, double reward) = 0;
	virtual void ResetState() = 0;
};


class RandomController : public Controller
{
private:
	double _prob;


public:
	RandomController(double prob);
	~RandomController();
	bool ChangeDirection();
};


class NivController : public Controller
{
protected:
	double _redWt, _blueWt, _grayWt;
	double _redDiffWt, _blueDiffWt, _grayDiffWt;
	double _rewardWt;
	bool _useRed, _useBlue, _useGray;
	bool _useDiffRed, _useDiffBlue, _useDiffGray;
	bool _useReward;

	double _actionM, _actionB;

	double _regLearnA, _regLearnB, _regLearnC, _regLearnD;
	double _diffLearnA, _diffLearnB, _diffLearnC, _diffLearnD;
	double _learnRate;

	bool _modRegOnReward, _modDiffOnReward;
	bool _modRegOnDiff, _modDiffOnReg;

	double _prevRed, _prevBlue, _prevGray;

public:
	NivController();
	~NivController();
	bool ChangeDirection(double red_pct, double blue_pct, double gray_pct, double reward);
	double Train(double red, double blue, double gray, double d_red, double d_blue, double d_gray, double reward);
	double EvaluateNetwork(double red, double blue, double gray, double d_red, double d_blue, double d_gray, double reward);
	void ResetState();
};





#endif //CONTROLLER_HPP