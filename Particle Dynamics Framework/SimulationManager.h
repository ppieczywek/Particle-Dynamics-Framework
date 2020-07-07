#pragma once
#include <typeindex> 
#include <iostream>
#include <fstream>
#include <memory>
#include <array>
#include <vector>
#include "openCLManager.h"
#include "Bead.h"
#include "SpringStructure.h"
#include "AngleBond.h"
#include "BufferChunk.h"
#include "BeadInfo.h"
#include "Type.h"
#include "Pair.h"
#include "SimulationSettings.h"
#include "GridParameters.h"
#include  <map>
#include  <set>
#include <algorithm>
#include <random>


class SimulationManager
{
	bool													init_status;
	int														simulation_steps;
	int														minimalization_steps;
	int														grid_update_freq;
	int														beads_num;
	int														springs_num;
	int														angles_num;
	//int														dump_counter;
	int														dump_frequency;
	std::string												dump_file_name;
	std::string												restart_file_name;
	std::string												vel_dump_file_name;
	SimulationSettings										settings;
	std::map<std::string,int>								type_id;
	std::vector<Type>										type_list;

	std::random_device					random_generator;
	std::mt19937						generator;//  (random_generator());
	std::normal_distribution<float>		toss; // (0.0, 1.0);

	std::vector<int>										export_bead_index;
	std::vector<std::string>								export_bead_name;


	std::vector<int>										collision_check_data;
	std::vector<InteractionNh>								random_check_data;
	std::vector<int>										grid_check_data;

	int														warning_msg_counter;

	bool													SetDumpFile(const std::string& value);
	bool													DumpData();
	bool													DumpRestartData();

		
	int														ScanText(std::vector<std::string>& text, std::string input, std::vector<std::string>& tokens, int start);

	bool													ReadInputFile(std::string inputFile);
	bool													ReadModelData(std::vector<std::string>& input_file_content);
	bool													ReadRestartFile(std::vector<std::string>& input_file_content);
	bool													ReadKernelData(std::string);
	bool													ReadTextTFileContent(std::string inputFile, std::vector<std::string>&	fileContents);
	bool													ReadBeadsData(std::vector<std::string>&	fileContents);
	bool													ReadSpringData(std::vector<std::string>&	fileContents);
	bool													ReadAngleData(std::vector<std::string>&	fileContents);

	bool													SetSimulationGrid();

	bool													RunSimulationLoop(int loop_size);
	bool													MinimizeSystemEnergy(int loop_size);
	bool													SimulationCheck();

	bool													GetSimulationSettings(std::vector<std::string>&	file_contents);
	bool													GetBeadInfo(std::vector<std::string>&	file_contents);
	bool													GetSimulationSteps(std::vector<std::string>& file_content);
	bool													GetOutputDataInfo(std::vector<std::string>&	file_contents);
	OpenClManager											openClManager;
public:

	
	bool													Initialize(std::string inputFile);
	bool													Run();
	void												    Close();

};