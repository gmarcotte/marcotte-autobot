#include <tasks.hpp>
#include <environment.hpp>

#include <iostream>

using namespace std;

extern char COLOR_ABBR[3];

ForagingTask::ForagingTask(Field* fld, Forager* frg)
{
	_field = fld;
	_forager = frg;
	frg->GetPosition(_forager_init_pos);
	frg->GetHeader(_forager_init_dir);
}

ForagingTask::~ForagingTask()
{
}

void ForagingTask::Update(double dt)
{
	_forager->Move(dt);
	double curr_position[3];
	_forager->GetPosition(curr_position);
	if (curr_position[2] <= 0.0)
	{	
		double rwd = _field->SampleCoord(curr_position[0], curr_position[1]);
		cout << "Reward Given: " << rwd << " from " << COLOR_ABBR[_field->GetColorByCoord(curr_position[0], curr_position[1])] << " target." << endl;
		_forager->GiveReward(rwd);
		_forager->SetPosition(_forager_init_pos[0], _forager_init_pos[1], _forager_init_pos[2]);
		_forager->SetHeader(_forager_init_dir[0], _forager_init_dir[1], _forager_init_dir[2]);
	}
	_forager->UpdateVisualField(_field);
}

