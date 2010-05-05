
#ifndef TASKS_HPP
#define TASKS_HPP

#include <environment.hpp>

class ForagingTask
{
protected:
	Field* _field;
	Forager* _forager;
	double _xy_coord_max;
	double _z_coord_min, _z_coord_max;
	int _trials;
	int _completed;
	int _curr_its;

public:
	ForagingTask(Field* fld, Forager* frg, double init_pos_max, double init_height_min, double init_height_max, int trials);
	~ForagingTask();
	int CompletedTrials();
	void Update(double dt);
	void ResetPosition();
	void Run(double dt);

};


#endif // TASKS_HPP