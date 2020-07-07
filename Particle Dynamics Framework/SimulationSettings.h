#pragma once
#include "Vector3.h"

struct SimulationSettings
{
	int   types_num;
	int   units;
	float cutoff;
	float cutoff_inv;
	float dt;
	float half_dt;
	float half_dt2;
	float lambda_dt;
	float gamma;
	float sigma;

	int x_cells_num;
	int y_cells_num;
	int z_cells_num;
	int xy_cells_num;
	float half_size_squared;
	Vector3 step;
	Vector3 size;
};