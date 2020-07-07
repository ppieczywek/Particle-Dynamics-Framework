// MSS_solver_gpu.cpp : Defines the entry point for the console application.
//
#pragma once
#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <tuple>
#include "CL/cl.hpp"
#include <math.h>
#include "SimulationManager.h"


int _tmain(int argc, _TCHAR* argv[])
{
	SimulationManager simulation_manager;
		
	if (simulation_manager.Initialize("input.txt"))
	{
		simulation_manager.Run();
		simulation_manager.Close();
	}
	
	//system("PAUSE");
	return 0;
}

