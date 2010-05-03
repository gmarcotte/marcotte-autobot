#include <utils.hpp>
#include <environment.hpp>
#include <forager.hpp>

#include <cmath>
using namespace std;

extern COLORREF COLOR_RGB[3];

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

void Camera::SetVisualField(int r_index, int th_index, COLOR col)
{
	int index = r_index * _numR + th_index;
	_visual_field[index] = col;
}


void Camera::Render(HDC hdc, POINT center, double rad)
{	
	double r, th;
	int x, y;
	HBRUSH NewBrush, OldBrush;
	HPEN NewPen, OldPen;
	for (int i=0; i<_numR; i++)
	{
		for (int j=0; j<_numR; j++)
		{
			r = rad * (GetR(i, j)/_radius);
			th = GetTheta(i, j);
			x = (int)(center.x + r*cos(th));
			y = (int)(center.y + r*sin(th));
			NewBrush = CreateSolidBrush(COLOR_RGB[_visual_field[i*_numR + j]]);
			NewPen = CreatePen(PS_SOLID, 1, COLOR_RGB[_visual_field[i*_numR + j]]);
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
	int ct = 0;
	for (int i=0; i<_numR*_numTheta; i++)
	{
		if (_visual_field[i] == col)
			ct++;
	}
	return double(ct) / double(_numR * _numTheta);
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

void Forager::UpdateVisualField()
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
			//_camera->SetVisualField(r_i, th_i, _task->GetFieldColor(x_int, y_int));
		}
	}
}

