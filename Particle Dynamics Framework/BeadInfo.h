#pragma once
#include "Vector3.h"

struct BeadInfo
{
	int			type;
	float		mass;
	float		mass_inv;
	float		lambda_dt_mass_inv;
	float		half_dt_mass_inv;
	float		half_dt2_mass_inv;
};