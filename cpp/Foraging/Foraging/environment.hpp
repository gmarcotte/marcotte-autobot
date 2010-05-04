
#ifndef ENVIRONMENT_HPP
#define ENVIRONMENT_HPP

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>

#include <windows.h>

#include <stdexcept>

using namespace std;

class Forager; // Declaration here because Field rendering needs Forager

const int NUM_COLORS = 3;

enum COLOR {
	RED,
	BLUE,
	GRAY
};

class BaseTarget
{
private:
	COLOR _color;

public:
	BaseTarget(COLOR clr = GRAY);
	virtual ~BaseTarget() {};
	char AsChar();
	void Render(HDC hdc, RECT box);
	virtual double Payout() = 0;
	COLOR GetColor();
};

// A constant source target, always providing the same reward
class ConstantTarget : public BaseTarget
{
private:
	double _value;

public:
	ConstantTarget(double val, COLOR clr = GRAY);
	double Payout() {return _value;};
};
/*************************************/

// A discrete target, providing several possible rewards with different probabilities
class DiscreteTarget : public BaseTarget
{
private:
	gsl_rng* _randGen;
	gsl_ran_discrete_t* _discretePdf;
	gsl_vector* _values;

public:
	DiscreteTarget(int N, double* probs, double* rewards, COLOR clr = GRAY);
	~DiscreteTarget();
	double Payout();
};
/********************************************************************************/


// A base field of targets
class Field
{
protected:
	int _height;
	int _width;	
	double _side;
	BaseTarget** _targets;

public:
	Field(int h, int w, double side);
	virtual ~Field();
	int GetIndex(int i, int j);
	int GetRowCoord(int index);
	int GetColCoord(int index);
	int GetRow(double y);
	int GetCol(double x);
	void AddTarget(int i, int j, BaseTarget* bt);
	void DeleteTarget(int i, int j);
	void ReplaceTarget(int i, int j, BaseTarget* bt);
	double Sample(int i, int j);
	double SampleCoord(double x, double y);
	COLOR GetColorByCoord(double x, double y);
	char* AsString();
	void Render(HDC hdc, RECT box);
	void Render(HDC hdc, int left, int top, int right, int bottom);
	void RenderForagerEllipse(HDC hdc, RECT box, Forager* fgr);
	void RenderForagerShadow(HDC hdc, RECT box, Forager* fgr);
};
/***************************/


// A square field of discrete targets and neutral areas, as described in Niv, et. al.
class BasicBinaryField : public Field
{
private:
	int _numNeutral;
	int _numRed;
	int _numBlue;
	
public:
	BasicBinaryField(int N, double nPct, double rPct, double rRwd, double rPrb, double bPct, double bRwd, double bPrb, double side = 10.0);
	int NumNeutral() {return _numNeutral;}
	int NumRed() {return _numRed;}
	int NumBlue() {return _numBlue;}
};



// A radial camera
class Camera
{
private:
	double _view_angle;
	double _focal_length;
	double _radius;
	double _deltaR;
	double _deltaTheta;
	int _numR;
	int _numTheta;
	COLOR* _visual_field;
	int* _color_counts;

public:
	Camera(double view_angle, double radius, double deltaR, double deltaTheta);
	~Camera();
	int GetNumR();
	int GetNumTheta();
	double GetR(int r_index, int th_index);
	double GetTheta(int r_index, int th_index);
	double GetFocalLength();
	double GetViewAngle();
	void SetVisualField(int r_index, int th_index, COLOR col);
	double GetPct(COLOR col);
	void Render(HDC hdc, POINT center, double rad);
	void PrintVisualField();
};
/********************************************************************************/


/* A base class for the forager */
class Forager
{
protected:
	double _reward;
	double _position[3];
	double _heading[3];
	double _speed;
	Camera* _camera;

public:
	Forager(double x, double y, double z,
		    double vx, double vy, double vz, double speed,
			double cone_angle, double proj_radius, double proj_dr, double proj_dth);
	~Forager();
	void GetPosition(double* pos);
	void GetHeader(double* head);
	double GetCameraAngle();
	void SetPosition(double x, double y, double z);
	void SetHeader(double vx, double vy, double vz);
	void SetSpeed(double speed);
	void GiveReward(double rwd);
	double GetReward();
	void Move(double dt);
	void UpdateVisualField(Field* field);
	void RenderCameraView(HDC hdc, RECT box);
	void PrintCameraView();
	double GetColorPct(COLOR clr);
};
/*****************************************************************************/

//class BasicBee : public Forager
//{
//};
/******************************************************************************/



#endif // ENVIRONMENT_HPP