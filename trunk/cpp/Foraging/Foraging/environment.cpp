#include <cmath>
#include <ctime>

#include <controller.hpp>
#include <environment.hpp>
#include <utils.hpp>

#include <iostream>
#include <windows.h>

using namespace std;

char COLOR_ABBR[3] = {'R', 'B', 'G'};
COLORREF COLOR_RGB[3] = {RGB(255, 0, 0),
						 RGB(0, 0, 255),
						 RGB(128, 128, 128)
};

// Implementation of the base target
BaseTarget::BaseTarget(COLOR clr)
{
	_color = clr;
}

char BaseTarget::AsChar()
{
	return COLOR_ABBR[_color];
}

COLOR BaseTarget::GetColor()
{
	return _color;
}

void BaseTarget::Render(HDC hdc, RECT box)
{
	HPEN NewPen = CreatePen(PS_SOLID, 1, COLOR_RGB[_color]);
	HPEN OldPen = (HPEN)SelectObject(hdc, NewPen);

	HBRUSH NewBrush = CreateSolidBrush(COLOR_RGB[_color]);
	HBRUSH OldBrush = (HBRUSH)SelectObject(hdc, NewBrush);

	Rectangle(hdc, box.left, box.top, box.right, box.bottom);

	SelectObject(hdc, OldBrush);
	SelectObject(hdc, OldPen);
	DeleteObject(NewBrush);
	DeleteObject(NewPen);
}
/**********************************************/

// Implementation of constant target
ConstantTarget::ConstantTarget(double val, COLOR clr)
: BaseTarget(clr)
{
	_value = val;
}
/**********************************************/

// Implementation of the discrete target
DiscreteTarget::DiscreteTarget(int N, double* probs, double* rewards, COLOR clr)
: BaseTarget(clr)
{
	for (int i=0; i<N; i++)
		if (probs[i] <= 0.0)
			throw Exception("Error: all probabilities must be positive");

	// Set up the discrete distribution
	_discretePdf = gsl_ran_discrete_preproc(N, probs);

	// Set up the rewards	
	_values = gsl_vector_alloc(N);
	for (int i=0; i<N; i++)
		gsl_vector_set(_values, i, rewards[i]);
}

DiscreteTarget::~DiscreteTarget()
{
	gsl_ran_discrete_free(_discretePdf);
	gsl_vector_free(_values);
}

double DiscreteTarget::Payout()
{
	size_t index = gsl_ran_discrete(get_rng(), _discretePdf);
	return gsl_vector_get(_values, index);
}
/*************************************************************/

// Implementation of Field
Field::Field(int h, int w, double side)
{
	_height = h;
	_width = w;
	_side = side;
	_targets = new BaseTarget*[h*w];
	for (int i=0; i<h*w; i++)
		_targets[i] = NULL;
}

Field::~Field()
{
	for (int i=0; i<_height; i++)
		for (int j=0; j<_width; j++)
			this->DeleteTarget(i, j);
	delete [] _targets;
}

int Field::GetIndex(int i, int j)
{
	if (i >= _height || i < 0)
		throw Exception("Error: row coord is not valid.");
	if (j >= _width || j < 0)
		throw Exception("Error: col coord is not valid.");
	return i*_width + j;
}

int Field::GetRowCoord(int index)
{
	if (index >= _width*_height || index < 0)
		throw Exception("Error: index is not valid.");
	return (index % _width);
}

int Field::GetColCoord(int index)
{
	int row = this->GetRowCoord(index);
	return index - row*_width;
}

int Field::GetRow(double y)
{
	int row;
	if (y >= 0 && y < _side*_width)
		row = (int)(y / _side);
	else
		row = -1;
	return row;
}

int Field::GetCol(double x)
{
	int col;
	if (x >= 0 && x < _side*_height)
		col = (int)(x / _side);
	else
		col = -1;
	return col;
}

void Field::AddTarget(int i, int j, BaseTarget *bt)
{
	int index = i*_width + j;
	if (_targets[index] != NULL)
		throw Exception("Error: adding target to non-null index.");
	else
		_targets[index] = bt;
}

void Field::DeleteTarget(int i, int j)
{
	int index = i*_width + j;
	if (_targets[index] == NULL)
		throw Exception("Error: trying to delete null index.");
	else
	{
		delete _targets[index];
		_targets[index] = NULL;
	}
}

void Field::ReplaceTarget(int i, int j, BaseTarget *bt)
{
	this->DeleteTarget(i, j);
	this->AddTarget(i, j, bt);
}

double Field::Sample(int i, int j)
{
	int index;
	try {
		index = this->GetIndex(i, j);
	}
	catch (Exception) {
		return 0.0;
	}
	if (_targets[index] == NULL)
		throw Exception("Error: Field target is uninitialized.");
	else
		return _targets[index]->Payout();
}

double Field::SampleCoord(double x, double y)
{
	int col = this->GetCol(x);
	int row = this->GetRow(y);
	return this->Sample(row, col);
}

COLOR Field::GetColorByCoord(double x, double y)
{
	int col = this->GetCol(x);
	int row = this->GetRow(y);
	int index;
	try {
		index = this->GetIndex(row, col);
	} catch (Exception) {
		return GRAY;
	}
	return _targets[index]->GetColor();
		
}

char* Field::AsString()
{
	char* out = new char[(_width+1)*_height];
	for (int i=0; i<_height; i++)
	{
		for (int j=0; j<_width; j++)
		{
			out[i*(_width+1)+j] = _targets[i*_width+j]->AsChar();
		}
		out[i*(_width+1)+_width] = '\n';
	}
	out[(_width+1)*_height - 1] = '\0';
	return out;
}

void Field::Render(HDC hdc, RECT box)
{
	double xStep = double((box.right - box.left)) / double(_width);
	double yStep = double((box.bottom - box.top)) / double(_height);

	RECT targetBox;
	for (int i=0; i<_height; i++)
	{
		for (int j=0; j<_width; j++)
		{
			int index = this->GetIndex(i, j);
			targetBox.left = (int)(box.left + j*xStep);
			targetBox.right = (int)(box.left + (j+1)*xStep);
			targetBox.top = (int)(box.top + i*yStep);
			targetBox.bottom = (int)(box.top + (i+1)*yStep);
			_targets[index]->Render(hdc, targetBox);
		}
	}
}

void Field::Render(HDC hdc, int left, int top, int right, int bottom)
{
	RECT box;
	box.left = left;
	box.top = top;
	box.right = right;
	box.bottom = bottom;
	this->Render(hdc, box);
}

void Field::RenderForagerShadow(HDC hdc, RECT box, Forager *fgr)
{
	double x_scale = (box.right - box.left) / (_side*_width);
	double y_scale = (box.bottom - box.top) / (_side*_height);

	HPEN GreenPen = CreatePen(PS_SOLID, 1, RGB(0, 255, 0));
	HBRUSH GreenBrush = CreateSolidBrush(RGB(0, 255, 0));
	HPEN OldPen = (HPEN)SelectObject(hdc, GreenPen);
	HBRUSH OldBrush = (HBRUSH)SelectObject(hdc, GreenBrush);

	double fgr_position[3];
	fgr->GetPosition(fgr_position);

	int rad = 2;
	RECT rect;
	rect.left = (int)(box.left + (fgr_position[0] - rad)*x_scale);
	rect.top = (int)(box.top + (fgr_position[1] - rad)*y_scale);
	rect.right = (int)(box.left + (fgr_position[0] + rad)*x_scale);
	rect.bottom = (int)(box.top + (fgr_position[1] + rad)*y_scale);
	Ellipse(hdc, rect.left, rect.top, rect.right, rect.bottom);

	SelectObject(hdc, OldPen);
	SelectObject(hdc, OldBrush);
	DeleteObject(GreenPen);
	DeleteObject(GreenBrush);
}

void Field::RenderForagerEllipse(HDC hdc, RECT box, Forager *fgr)
{
	double x_scale = (box.right - box.left) / (_side*_width);
	double y_scale = (box.bottom - box.top) / (_side*_height);

	HPEN YellowPen = CreatePen(PS_SOLID, 3, RGB(255, 255, 0));
	HPEN OldPen = (HPEN)SelectObject(hdc, YellowPen);

	// Get visual field of cone intersection
	// Define the cone
	double coneOrg[3];
	fgr->GetPosition(coneOrg);
	double coneDir[3];
	fgr->GetHeader(coneDir);
	double coneAngle = fgr->GetCameraAngle();
	//Ellipse(hdc, (int)coneOrg[0]-10, (int)coneOrg[1]-10, (int)coneOrg[0]+10, (int)coneOrg[1]+10);

	// Define the plane (xy)
	double planeOrg[3] = {0, 0, 0};
	double planeDir1[3] = {1, 0, 0};
	double planeDir2[3] = {0, 1, 0};
	
	double* ellipse = new double[5];
	conePlaneIntersection(coneOrg, coneDir, coneAngle, planeOrg, planeDir1, planeDir2, ellipse);
	
	double dx = ellipse[0];
	double dy = ellipse[1];
	double a = ellipse[2];
	double b = ellipse[3];
	double theta = ellipse[4];

	delete [] ellipse;

	// Create Ellipse
	RECT rect;
	rect.top = (int)(box.top + (dy - b)*y_scale);
	rect.left = (int)(box.left + (dx - a)*x_scale);
	rect.right = (int)(box.left + (dx + a)*x_scale);
	rect.bottom = (int)(box.top + (dy + b)*y_scale);
	POINT ellipsePts[13];
	EllipseToBezier(rect, ellipsePts);

	// Rotate
	POINT midPoint;
	midPoint.x = (int)((rect.left + rect.right)/2);
	midPoint.y = (int)((rect.top + rect.bottom)/2);
	Rotate(theta, midPoint, ellipsePts, 13);
	BeginPath(hdc);
	PolyBezier(hdc, ellipsePts, 13);
	EndPath(hdc);
	StrokePath(hdc);
	
	SelectObject(hdc, OldPen);
	DeleteObject(YellowPen);
}
/********************************************************************************/


// Implementation of BinaryBasicField
BasicBinaryField::BasicBinaryField(int N, double nPct, 
								   double rPct, double rRwd, double rPrb,
								   double bPct, double bRwd, double bPrb, double side)
: Field(N, N, side)
{
	// Initialize instance variables
	_numNeutral = 0;
	_numRed = 0;
	_numBlue = 0;
	
	// Initialize targets
	double probs[3] = {rPct, bPct, nPct};
	double types[3] = {0.0, 1.0, 2.0};
	DiscreteTarget* fieldPicker = new DiscreteTarget(3, probs, types);
	double nRwd = 0.0;
	double nPrb = 1.0;
	double rPrbs[2] = {rPrb, 1.0 - rPrb};
	double rRwds[2] = {rRwd, 0.0};
	double bPrbs[2] = {bPrb, 1.0 - bPrb};
	double bRwds[2] = {bRwd, 0.0};
	for (int i=0; i<_height; i++)
	{
		for (int j=0; j<_width; j++)
		{
			double type = fieldPicker->Payout();
			if (type == 0.0)
			{
				this->AddTarget(i, j, new DiscreteTarget(2, rPrbs, rRwds, RED));
				_numRed++;
			}
			else if (type == 1.0)
			{
				this->AddTarget(i, j, new DiscreteTarget(2, bPrbs, bRwds, BLUE));
				_numBlue++;
			}
			else if (type == 2.0)
			{
				this->AddTarget(i, j, new DiscreteTarget(1, &nPrb, &nRwd, GRAY));
				_numNeutral++;
			}
		}
	}
	delete fieldPicker;
}
/******************************************************************************************/



Camera::Camera(double view_angle, double radius, double deltaR, double deltaTheta)
{
	this->_view_angle = view_angle;
	this->_radius = radius;
	this->_focal_length = radius / tan(view_angle);
	this->_deltaR = deltaR;
	this->_deltaTheta = deltaTheta;
	this->_numR = (int)(radius / deltaR + 1);
	this->_numTheta = (int)(2*PI / deltaTheta + 1);
	this->_visual_field = new COLOR[_numR * _numTheta];
	this->_color_counts = new int[NUM_COLORS];
	for (int i=0; i<NUM_COLORS; i++)
		_color_counts[i] = 0;
	for (int i=0; i<_numR*_numTheta; i++)
		_visual_field[i] = RED;
	_color_counts[RED] = _numR*_numTheta;
}

Camera::~Camera()
{
	delete _visual_field;
	delete _color_counts;
}

double Camera::GetR(int r_index, int th_index)
{
	return r_index*_deltaR;	
}

double Camera::GetTheta(int r_index, int th_index)
{
	return th_index*_deltaTheta;
}

double Camera::GetFocalLength()
{
	return _focal_length;
}

int Camera::GetNumR()
{
	return _numR;
}

int Camera::GetNumTheta()
{
	return _numTheta;
}

double Camera::GetViewAngle()
{
	return _view_angle;
}

void Camera::SetVisualField(int r_index, int th_index, COLOR col)
{
	int index = r_index * _numTheta + th_index;
	_color_counts[_visual_field[index]]--;
	_visual_field[index] = col;
	_color_counts[_visual_field[index]]++;
}

void Camera::PrintVisualField()
{
	for (int i=0; i<_numR*_numTheta; i++)
	{
		cout << COLOR_ABBR[_visual_field[i]] << endl;	
	}
}

void Camera::Render(HDC hdc, POINT center, double rad)
{	
	double r, th;
	int x, y;
	HBRUSH NewBrush, OldBrush;
	HPEN NewPen, OldPen;
	for (int i=0; i<_numR; i++)
	{
		for (int j=0; j<_numTheta; j++)
		{
			r = rad * ((double)GetR(i, j)/ (double)_radius);
			th = -GetTheta(i, j);
			x = (int)(center.x + r*cos(th));
			y = (int)(center.y + r*sin(th));
			NewBrush = CreateSolidBrush(COLOR_RGB[_visual_field[i*_numTheta + j]]);
			NewPen = CreatePen(PS_SOLID, 1, COLOR_RGB[_visual_field[i*_numTheta + j]]);
			OldBrush = (HBRUSH)SelectObject(hdc, NewBrush);
			OldPen = (HPEN)SelectObject(hdc, NewPen);
			Ellipse(hdc, x-1, y-1, x+1, y+1);
			SelectObject(hdc, OldBrush);
			SelectObject(hdc, OldPen);
			DeleteObject(NewBrush);
			DeleteObject(NewPen);
		}
	}
}


double Camera::GetPct(COLOR col)
{
	return double(_color_counts[col]) / double(_numR * _numTheta);
}



// Implementation of the basic forager class
Forager::Forager(double x, double y, double z,
				 double vx, double vy, double vz,
				 double speed,
				 double cone_angle, double proj_radius, double proj_dr, double proj_dth)
{
	_reward = 0;
	this->SetPosition(x, y, z);
	this->SetHeader(vx, vy, vz);
	this->SetSpeed(speed);
	this->_camera = new Camera(cone_angle, proj_radius, proj_dr, proj_dth);
}

Forager::~Forager()
{
	delete _camera;
}

void Forager::SetController(Controller* ctrl)
{
	_controller = ctrl;
}

void Forager::GetPosition(double* pos)
{
	pos[0] = _position[0];
	pos[1] = _position[1];
	pos[2] = _position[2];
}

double Forager::GetCameraAngle()
{
	return _camera->GetViewAngle();
}

double Forager::GetColorPct(COLOR clr)
{
	return _camera->GetPct(clr);
}


void Forager::GetHeader(double* head)
{
	head[0] = _heading[0];
	head[1] = _heading[1];
	head[2] = _heading[2];
}

void Forager::SetSpeed(double speed)
{
	_speed = speed;
}

void Forager::SetPosition(double x, double y, double z)
{
	_position[0] = x;
	_position[1] = y;
	_position[2] = z;
}

void Forager::SetHeader(double vx, double vy, double vz)
{
	if (vx*vx + vy*vy + vz*vz < 1e-6)
		throw Exception("Error: Trying to set header to zero vector");
	_heading[0] = vx;
	_heading[1] = vy; 
	_heading[2] = vz;
	normalize(3, _heading, _heading);
}

void Forager::Move(double dt)
{
	_position[0] += _heading[0]*_speed*dt;
	_position[1] += _heading[1]*_speed*dt;
	_position[2] += _heading[2]*_speed*dt;
}

void Forager::ChangeToRandomHeading()
{
	double theta = gsl_rng_uniform(get_rng()) * (2*PI);
	double phi = gsl_rng_uniform(get_rng()) * (PI / 2.0);
	_heading[0] = cos(theta)*cos(phi);
	_heading[1] = sin(theta)*cos(phi);
	_heading[2] = -sin(phi);
	normalize(3, _heading, _heading); // Just to be sure it's a unit vector
}

void Forager::Update(Field* field)
{
	UpdateVisualField(field);
	if (_controller->ChangeDirection(_camera->GetPct(RED), _camera->GetPct(BLUE), _camera->GetPct(GRAY), 0.0))
		ChangeToRandomHeading();
}

void Forager::UpdateVisualField(Field* field)
{
	double center[3];
	center[0] = _position[0] - _camera->GetFocalLength()*_heading[0];
	center[1] = _position[1] - _camera->GetFocalLength()*_heading[1];
	center[2] = _position[2] - _camera->GetFocalLength()*_heading[2];

	double r1[3];
	double r2[3];
	double rho, theta;
	double P[3];

	double x_int, y_int;

	r1[0] = 1.0;
	r1[1] = 0.0;
	r1[2] = (_heading[0]*center[0] + _heading[1]*center[1] + _heading[2]*center[2] - _heading[0])/_heading[2];
	normalize(3, r1, r1);

	CrossProduct(_heading, r1, r2);
	normalize(3, r2, r2);
	for (int r_i=0; r_i<_camera->GetNumR(); r_i++)
	{
		for (int th_i=0; th_i<_camera->GetNumTheta(); th_i++)
		{
			rho = _camera->GetR(r_i, th_i);
			theta = _camera->GetTheta(r_i, th_i);
			P[0] = center[0] + rho*cos(theta)*r1[0] + rho*sin(theta)*r2[0];
			P[1] = center[1] + rho*cos(theta)*r1[1] + rho*sin(theta)*r2[1];
			P[2] = center[2] + rho*cos(theta)*r1[2] + rho*sin(theta)*r2[2];
			x_int = (P[0]*_position[2] - P[2]*_position[0])/(_position[2] - P[2]);
			y_int = (P[1]*_position[2] - P[2]*_position[1])/(_position[2] - P[2]);
			//cout << "Updating R=" << rho << ", theta=" << theta << " to " << COLOR_ABBR[field->GetColorByCoord(x_int, y_int)] << endl;
			_camera->SetVisualField(r_i, th_i, field->GetColorByCoord(x_int, y_int));
		}
	}
}

void Forager::RenderCameraView(HDC hdc, RECT box)
{
	POINT c;
	c.x = (box.left + box.right)/2;
	c.y = (box.top + box.bottom)/2;

	double rad;
	rad = min(box.right - box.left, box.bottom - box.top) / 2.0;
	_camera->Render(hdc, c, rad);
}

void Forager::PrintCameraView()
{
	_camera->PrintVisualField();
}

void Forager::GiveReward(double rwd)
{
	_reward += rwd;
	_controller->ChangeDirection(_camera->GetPct(RED), _camera->GetPct(BLUE), _camera->GetPct(GRAY), rwd);
	_controller->ResetState();
}

double Forager::GetReward()
{
	return _reward;
}