#include <tasks.hpp>
#include <environment.hpp>
#include <utils.hpp>

#include <iostream>

using namespace std;

extern char COLOR_ABBR[3];

ForagingTask::ForagingTask(Field* fld, Forager* frg, double init_pos_max, double init_height_min, double init_height_max, int trials)
{
	_completed = 0;
	_trials = trials;
	_field = fld;
	_forager = frg;
	_xy_coord_max = init_pos_max;
	_z_coord_min = init_height_min;
	_z_coord_max = init_height_max;
	_curr_its = 0;
	ResetPosition();
}

ForagingTask::~ForagingTask()
{
}

int ForagingTask::CompletedTrials()
{
	return _completed;
}

void ForagingTask::ResetPosition()
{
	double new_pos[3];
	new_pos[0] = rand_double_uniform(0.0, _xy_coord_max);
	new_pos[1] = rand_double_uniform(0.0, _xy_coord_max);
	new_pos[2] = rand_double_uniform(_z_coord_min, _z_coord_max);
	_forager->SetPosition(new_pos[0], new_pos[1], new_pos[2]);
	_forager->ChangeToRandomHeading();
}

void ForagingTask::Update(double dt)
{
	_curr_its++;
	double curr_position[3];
	_forager->GetPosition(curr_position);
	if (curr_position[2] <= 0.0)
	{	
		_completed++;
		double rwd = _field->SampleCoord(curr_position[0], curr_position[1]);
		cout << _completed << "::" << rwd << "::" << COLOR_ABBR[_field->GetColorByCoord(curr_position[0], curr_position[1])] << "::" << _curr_its << endl;
		_forager->GiveReward(rwd);
		this->ResetPosition();
		_curr_its = 0;
	}
	else if (_curr_its > 50)
	{
		double curr_position[3];
		_forager->GetPosition(curr_position);
		this->ResetPosition();
		cout << "Max iterations reached -- Z=" << curr_position[2] << endl;
		_completed++;
		_curr_its = 0;
	}
	else
	{
		_forager->Update(_field);
		_forager->Move(dt);
	}
}


void ForagingTask::Run(double dt)
{
	//double time = 0.0;
	int switch_point = gsl_rng_uniform_int(get_rng(), _trials/2) + _trials/4;
	while (_completed < switch_point)
	{
		Update(dt);
	//	time += dt;
	}
	// Switch flower contingencies
	cout << "SWITCHING PROBABILITIES" << endl;
	_field->SwitchRedAndBlue();
	while (_completed < _trials)
	{
		Update(dt);
	//	time += dt;
	}
	cout <<"Total Reward: " << _forager->GetReward() << endl;
}
