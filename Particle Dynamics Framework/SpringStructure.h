#pragma once
#include "Vector3.h"

struct SpringStructure
{	
	int						type;
	float					rest_length;
	float					stiffness;
	float					damping;
	int						p1;
	int						p2;
};