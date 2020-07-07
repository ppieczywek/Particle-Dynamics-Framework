#pragma once
#include "stdafx.h"
#include <string>

struct BufferIndex
{
	static int bead_buffer;
	static int bead_info_buffer;
	static int spring_buffer;
	static int angle_buffer;
		
	static int grid_head_buffer;
	static int grid_tail_buffer;
	static int grid_cell_nh_buffer;

	static int random_buffer;
	static int random_buffer_stride;

	static int bead_nh_head_buffer;
	static int bead_nh_tail_buffer;


	static int bead_pair_coeff_buffer; //interactionType;
	static int bead_pair_table_buffer;// intercationCoeffs;
	//static int particleInfo;
	
};

