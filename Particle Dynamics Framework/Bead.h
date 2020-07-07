#pragma once
#include "Vector3.h"

struct Bead
{
	Vector3		position;
	Vector3		image_position;
	Vector3		velocity;
	Vector3		force;
	Vector3		velocity_old;
	Vector3		force_old;
};