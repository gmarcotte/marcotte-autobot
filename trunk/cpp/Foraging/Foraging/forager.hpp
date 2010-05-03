
#ifndef FORAGER_HPP
#define FORAGER_HPP

#include <windows.h>

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

public:
	Camera(double view_angle, double radius, double deltaR, double deltaTheta);
	int GetNumR();
	int GetNumTheta();
	double GetR(int r_index, int th_index);
	double GetTheta(int r_index, int th_index);
	double GetFocalLength();
	void SetVisualField(int r_index, int th_index, COLOR col);
	double GetPct(COLOR col);
	void Render(HDC hdc, POINT center, double rad);
};


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
	void SetPosition(double x, double y, double z);
	void SetHeader(double vx, double vy, double vz);
	void SetSpeed(double speed);
	void Move(double dt);
	void UpdateVisualField();
};


//class BasicBee : public Forager
//{
//};


#endif //FORAGER_HPP