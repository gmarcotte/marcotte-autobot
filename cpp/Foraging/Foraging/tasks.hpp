
#ifndef TASKS_HPP
#define TASKS_HPP

#include <environment.hpp>

class ForagingTask
{
protected:
	Field* _field;
	Forager* _forager;
	double _forager_init_pos[3];
	double _forager_init_dir[3];

public:
	ForagingTask(Field* fld, Forager* frg);
	~ForagingTask();
	void Update(double dt);

};


#endif // TASKS_HPP