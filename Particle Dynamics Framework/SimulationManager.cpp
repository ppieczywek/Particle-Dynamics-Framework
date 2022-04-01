#pragma once
#include "stdafx.h"
#include "SimulationManager.h"
#include "HighResolutionTimer.h"
#include <map>

bool SimulationManager::Initialize(std::string input_file)
{
	std::cout << std::endl << "%------- Initializing Simulation Manager -------%" << std::endl << std::endl;
	simulation_steps		= -1;
	dump_frequency			= -1;
	beads_num				=  0;
	springs_num				=  0;
	grid_update_freq		=  8;
	minimalization_steps	=  1000;
	warning_msg_counter		=  0;

	generator = std::mt19937(random_generator());
	toss = std::normal_distribution<float>(0.0, 1.0);
	//toss = std::uniform_real_distribution<float>(0.0, 1.0);

	if (!openClManager.Initialize())
	{
		std::cout << "Error occured during initialization of openClManager!" << std::endl;
		std::cout << "Closing program!" << std::endl;
		return false;
	}

	if (!ReadInputFile(input_file)) return false;
	if (!ReadKernelData("kernel.cl")) return false;

	

	init_status = true;
	return true;
};


bool SimulationManager::Run()
{
	if (init_status)
	{
		float progress    = 0;
		int substep_num   = static_cast<int>(std::ceil(static_cast<float>(simulation_steps) / static_cast<float>(dump_frequency)));
		auto total_volume = (settings.size.x * settings.size.y * settings.size.z);
		auto density = static_cast<float>(beads_num) / (total_volume);
		HighResolutionTimer timer;

		std::cout << std::endl << std::endl << "%-------  Running Simulation  -------%" << std::endl << std::endl;
		
		std::cout << "\tBox size:\t" << settings.size.x << " "
			<< settings.size.y << " " << settings.size.z << std::endl;

		std::cout << "\tGrid cell size:\t" << settings.step.x << " x "
			<< settings.step.y << " x " << settings.step.z << std::endl;

		std::cout << "\tGrid cells:\t" << settings.x_cells_num << " x "
			<< settings.y_cells_num << " x " << settings.z_cells_num << std::endl;

		
		std::cout << "\tParticle density:\t" << density << std::endl;
		std::cout << "\tSimulation steps:\t" << simulation_steps << std::endl;
		std::cout << "\tData dump frequency:\t" << dump_frequency << std::endl;
		std::cout << "\tTime step:\t" << settings.dt << std::endl << std::endl << std::endl;
		
			
		//if (!DumpTopology()) return false;

		timer.Start();
		
		if (minimalization_steps > 0)
		{
			std::cout << "\tInitializing system energy minimalization." << std::endl; //zmieniæ to --
			if (!MinimizeSystemEnergy(minimalization_steps))
			{
				std::cout << "Error: An error occured during system energy minimalization." << std::endl;
				return false;
			}
			std::cout << "\tSystem energy minimalization finished." << std::endl;
		}


		for (int substep = 0; substep < substep_num; substep++)
		{
			if (!RunSimulationLoop(dump_frequency))
			{
				std::cout << "Error: Unable to start simulation main loop." << std::endl;
				return false;
			}
			if (!DumpData())
			{
				std::cout << "Error: An error occured during simulation data dump." << std::endl;
				return false;
			}
			progress = (static_cast<float>(substep) / static_cast<float>(substep_num))*100.0f;
			std::cout << "Simulation progress: " << progress << " %" << std::endl;
		}
		std::cout << "Simulation progress: 100.0%" << std::endl;

		DumpRestartData(); //dopracowaæ to

		double elapsed_seconds = timer.Stop();
		double elapsed_minutes = std::floor((elapsed_seconds / 60.0));
		elapsed_seconds = elapsed_seconds - elapsed_minutes * 60.0;
		std::cout << "Simulation finished after " << elapsed_minutes << " minutes and " << elapsed_seconds << " seconds." << std::endl;
		std::cout << "Closing program." << std::endl;
				


		return true;
	}
	return false;
};


bool	SimulationManager::MinimizeSystemEnergy(int loop_size)
{
	if (loop_size < 0)
	{
		std::cout << "Error: :MinimizeSystemEnergy function - loop size should be a positive number." << std::endl;
		return false;
	}

	int ResetBeadEnergy = openClManager.GetKernelId("ResetBeadEnergy");
	if (ResetBeadEnergy < 0)
	{
		std::cout << "Error: MinimizeSystemEnergy function - unable to get kernel id." << std::endl;
		return false;
	}
	
	if (!RunSimulationLoop(loop_size)) return false;
	
	std::vector<Bead> bead_data(openClManager.GetBufferLength(BufferIndex::bead_buffer), Bead());
	if (!openClManager.ReadBuffer(BufferIndex::bead_buffer, &bead_data.front())) return false;

	float average_velocity = 0.0f;
	for (int ii = 0; ii < bead_data.size(); ii++)
	{
		float length = bead_data[ii].velocity.x* bead_data[ii].velocity.x;
		length += bead_data[ii].velocity.y* bead_data[ii].velocity.y;
		length += bead_data[ii].velocity.z* bead_data[ii].velocity.z;
		
		if (length < 0.00001f)
		{
			length = 0.0f;
		}
		else
		{
			length = sqrt(length);
		}
		average_velocity += length;
	}
	average_velocity /= bead_data.size();

	
	
	for (int ii = 0; ii < bead_data.size(); ii++)
	{
		float magnitude = (toss(generator))*0.2f + average_velocity;
		Vector3 velocity;
		velocity.x = toss(generator);
		velocity.y = toss(generator);
		velocity.z = toss(generator);
		float length = (velocity.x*velocity.x + velocity.y*velocity.y + velocity.z*velocity.z);
		if (length < 0.00001f)
		{
			length = 0.0f;
		}
		else
		{
			length = 1.0f / sqrt(length);
		}

		velocity.x *= length*magnitude;
		velocity.y *= length*magnitude;
		velocity.z *= length*magnitude;
		bead_data[ii].velocity = velocity;
		bead_data[ii].velocity_old = velocity;
		bead_data[ii].image_position = bead_data[ii].position;

		bead_data[ii].force.x = 0.0f;
		bead_data[ii].force.y = 0.0f;
		bead_data[ii].force.z = 0.0f;

		bead_data[ii].force_old.x = 0.0f;
		bead_data[ii].force_old.y = 0.0f;
		bead_data[ii].force_old.z = 0.0f;
	}
	// nie zapisywac nowych predkosci
	//if (!openClManager.WriteBuffer<Bead>(BufferIndex::bead_buffer, bead_data)) return false;
	return true;
}


bool	SimulationManager::RunSimulationLoop(int loop_size)
{
	static int rnd_num_refresh_state = 0;
	static std::vector<InteractionNh>	random_coeff(openClManager.GetBufferLength(BufferIndex::random_buffer), InteractionNh());

	if (loop_size <= 0)
	{
		std::cout << "Error: RunSimulationLoop function - loop size should be a positive number." << std::endl;
		return false;
	}

	int		   grid_update_cnt			= 0;
	static int fill_grid_kr				= openClManager.GetKernelId("FillGrid");
	static int reset_grid_kr			= openClManager.GetKernelId("ResetGrid");
	static int resolve_springs_kr		= openClManager.GetKernelId("ResolveSprings");
	static int resolve_angles_kr		= openClManager.GetKernelId("ResolveAngles");
	static int resolve_contacts_kr		= openClManager.GetKernelId("ResolveContacts");
	static int advance_velocity_kr		= openClManager.GetKernelId("AdvanceVelocity");
	static int advance_position_kr		= openClManager.GetKernelId("AdvancePosition");
	static int init_contact_list_kr		= openClManager.GetKernelId("InitContactList");
	static int reset_contact_list_kr	= openClManager.GetKernelId("ResetContactList");
	static int build_contact_list_kr	= openClManager.GetKernelId("BuildContactList");
	
	for (auto step = 0; step < loop_size; step++)
	{
		if (grid_update_cnt == grid_update_freq) grid_update_cnt = 0;

		if (!openClManager.executeKernel(advance_position_kr)) return false;
		//if (resolve_springs_kr > 0)
		if (springs_num > 0)
		{
			if (!openClManager.executeKernel(resolve_springs_kr)) return false;
		}

		//if (resolve_angles_kr > 0)
		if (angles_num > 0)
		{
			if (!openClManager.executeKernel(resolve_angles_kr)) return false;
		}

		if (grid_update_cnt == 0)
		{
			if (!openClManager.executeKernel(reset_grid_kr)) return false;
			if (!openClManager.executeKernel(reset_contact_list_kr)) return false;
			if (!openClManager.executeKernel(fill_grid_kr)) return false;
			if (!openClManager.executeKernel(build_contact_list_kr)) return false;
			if (!openClManager.executeKernel(init_contact_list_kr)) return false;
		}

		if (!openClManager.executeKernel(resolve_contacts_kr)) return false;
		if (!openClManager.executeKernel(advance_velocity_kr)) return false;
		if (!openClManager.runKernels()) return false;

		grid_update_cnt++;
		
		if (rnd_num_refresh_state == 128)
		{
			float avg = 0.0f;
			float	factor = 1.0f / sqrt(settings.dt);

			for (int ii = 0; ii < random_coeff.size(); ii++)
			{
				for (int jj = 0; jj < 256; jj++)
				{
					random_coeff[ii].n[jj] = factor * toss(generator);
					avg = avg + random_coeff[ii].n[jj];
				}

			}

			avg /= random_coeff.size() * 256;
			for (int ii = 0; ii < random_coeff.size(); ii++)
			{
				for (int jj = 0; jj < 256; jj++)
				{
					random_coeff[ii].n[jj] -= avg;
				}
			}

			if (!openClManager.WriteBuffer<InteractionNh>(BufferIndex::random_buffer, random_coeff)) return false;
			rnd_num_refresh_state = 0;
		}
		rnd_num_refresh_state += 1;
		
		
		//SimulationCheck();
		
	}
	return true;
};



int SimulationManager::ScanText(std::vector<std::string>& text, std::string input, std::vector<std::string>& tokens, int startLine)
{
	std::string	 str;
	tokens.clear();
	if (startLine >= 0)
	{
		if (!text.empty())
		{
			for (int ii = startLine; ii < static_cast<int>(text.size()); ii++)
			{
				str = text[ii];
				std::istringstream iss(str);
				std::copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(), std::back_inserter(tokens));
				if (!tokens.empty())
				{
					if (tokens[0] == input)
					{
						return ii+1;
					}
				}
				tokens.clear();
			}
		}
	}
	return -1;
};


bool SimulationManager::ReadTextTFileContent(std::string inputFile, std::vector<std::string>&	fileContents)
{
	std::ifstream	file(inputFile);

	if (file.is_open())
	{
		fileContents.clear();
		std::string str;
		while (std::getline(file, str))
		{
			fileContents.push_back(str);
		}
		file.close();
		return true;
	}
	else
	{
		return false;
	}
}


bool SimulationManager::GetSimulationSettings(std::vector<std::string>&	file_contents)
{
	std::vector<std::string>	tokens;

	if (ScanText(file_contents, "step", tokens, 0))
	{
		if (tokens.size() > 1)
		{
			settings.dt = std::stof(tokens[1]);
			if (settings.dt > 0.0f)
			{
				settings.half_dt = settings.dt * 0.5f;
				settings.half_dt2 = settings.dt * settings.dt * 0.5f;
			}
			else
			{
				std::cout << "Error: incorrect value for time step!" << std::endl; //linia
				return false;
			}
		}
		else
		{
			std::cout << "Error: time step was not specified." << std::endl;
			return false;
		}
	}

	if (ScanText(file_contents, "lambda", tokens, 0))
	{
		if (tokens.size() > 1)
		{
			settings.lambda_dt = std::stof(tokens[1]);
			if (settings.lambda_dt > 0.0f)
			{
				settings.lambda_dt *= settings.dt;
			}
			else
			{
				std::cout << "Error: incorrect value for lambda!" << std::endl; //linia
				return false;
			}
		}
		else
		{
			std::cout << "Info: Lambda was not specified in the input file. using deafult value 0.5" << std::endl;
			return false;
		}
	}


	if (ScanText(file_contents, "sigma", tokens, 0))
	{
		if (tokens.size() > 1)
		{
			settings.sigma = std::stof(tokens[1]);
			if (settings.sigma < 0.0f)
			{
				std::cout << "Error: incorrect value for sigma!" << std::endl; //linia
				return false;
			}
		}
		else
		{
			std::cout << "Error: sigma was not specified in the input file!" << std::endl;
			return false;
		}
	}


	if (ScanText(file_contents, "gamma", tokens, 0))
	{
		if (tokens.size() > 1)
		{
			settings.gamma = std::stof(tokens[1]);
			if (settings.gamma < 0.0f)
			{
				std::cout << "Error: incorrect value for gamma!" << std::endl; //linia
				return false;
			}
		}
		else
		{
			std::cout << "Error: gamma was not specified in the input file!" << std::endl;
			return false;
		}
	}

	if (ScanText(file_contents, "units", tokens, 0))
	{
		if (tokens.size() > 1)
		{
			if (tokens[1] == "native")
			{
				settings.units = 1;
			}
			else if (tokens[1] == "reduced")
			{
				settings.units = 0;
			}
			else
			{
				std::cout << "Error: parameter 'units' takes two values - native or reduced!" << std::endl;
				return false;
			}
		}
		else
		{	
			settings.units = 0; //by default //std::cout << "Error: wrong definition of 'units' parameter!" << std::endl;
			//return false;
		}
	}
	

	if (ScanText(file_contents, "cutoff", tokens, 0))
	{
		if (tokens.size() > 1)
		{
			settings.cutoff = std::stof(tokens[1]);
			if (settings.cutoff > 0.0f)
			{
				settings.cutoff_inv = 1.0f / settings.cutoff;
			}
			else
			{
				std::cout << "Error: incorrect value for cutoff radius!" << std::endl;
				return false;
			}
		}
		else
		{
			std::cout << "Error: cutoff radius was not specified in the input file!" << std::endl;
			return false;
		}
	}

	if (ScanText(file_contents, "box", tokens, 0))
	{
		if (tokens.size() == 4)
		{
			settings.size.x = std::stof(tokens[1]);
			settings.size.y = std::stof(tokens[2]);
			settings.size.z = std::stof(tokens[3]);

			if (settings.size.x < 0.0f || settings.size.y < 0.0f || settings.size.z < 0.0f)
			{
				std::cout << "Error: box dimensions must have positive values!" << std::endl;
				return false;
			}

			if ((settings.size.x / settings.cutoff) < 3.0f ||
				(settings.size.y / settings.cutoff) < 3.0f ||
				(settings.size.z / settings.cutoff) < 3.0f)
			{
				std::cout << "Error: box dimensions must be at least three times bigger than cutoff radius!" << std::endl;
				return false;
			}
		}
		else
		{
			std::cout << "Error: box dimensions were not specified in the input file!" << std::endl;
			return false;
		}
	}
		

	int	types_num = 0;
	int my_start = 0;
	while ((my_start = ScanText(file_contents, "type", tokens, my_start)) != -1)
	{
		if (tokens.size() == 3)
		{
			types_num++;
		}
	}

	if (types_num <= 0)
	{
		std::cout << "Error: input file requires bead type definition!" << std::endl;
		return false;
	}

	settings.types_num = types_num;
	
	if (!SetSimulationGrid()) return false;

	return true;
};


bool SimulationManager::GetBeadInfo(std::vector<std::string>&	file_contents)
{
	std::vector<std::string>	tokens;

	int my_start = 0;
	int types_num = 0;
	while ((my_start = ScanText(file_contents, "type", tokens, my_start)) != -1)
	{
		if (tokens.size() == 3)
		{
			auto id = tokens[1];
			auto mass = std::stof(tokens[2]);

			auto it = type_id.find(id);
			if (it == type_id.end())
			{
				Type type;
				type.id = id;
				type.mass = mass;
				type_list.push_back(type);
				type_id[id] = types_num;
				types_num++;
			}
			else
			{
				std::cout << "Error: in input file line -  " << my_start << "\t" << " - bead type already defined!" << std::endl;
				return false;
			}
		}
		else
		{
			std::cout << "Error: in input file line -  " << my_start << "\t" << " - incorrect type definition!" << std::endl;
			return false;
		}
	}

	if (types_num <= 0)
	{
		std::cout << "Error: input file requires bead type definition!" << std::endl;
		return false;
	}

	
	int pairs_num = ((types_num*(types_num - 1)) / 2) + types_num;
	
	std::vector<Pair> pair_list(pairs_num, Pair());
	for (int ii = 0; ii < pair_list.size(); ii++)
	{
		pair_list[ii].p1 = -1;
		pair_list[ii].p2 = -1;
	}
	
	my_start = 0;
	while ((my_start = ScanText(file_contents, "pair", tokens, my_start)) != -1)
	{
		if (tokens.size() == 5)
		{
			std::string type_name_1 = tokens[1];
			std::string type_name_2 = tokens[2];

			if (type_id.find(type_name_1) == type_id.end())
			{
				std::cout << "Error: input file line - " << my_start << "\t" << " - first bead type id deos not match any of the defined types.!" << std::endl;
				return false;
			}
			if (type_id.find(type_name_2) == type_id.end())
			{
				std::cout << "Error: input file line - " << my_start << "\t" << " - second bead type id deos not match any of the defined types.!" << std::endl;
				return false;
			}

			int type_1 = type_id.find(type_name_1)->second;
			int type_2 = type_id.find(type_name_2)->second;

			if (type_1 > type_2) std::swap(type_1, type_2);

			int pair_list_index = (types_num * type_1) + type_2 - ((type_1*(type_1 + 1)) / 2);

			if (pair_list_index >= 0 && pair_list_index < pair_list.size())
			{
				Pair interaction_pair;
				interaction_pair.p1 = type_1;
				interaction_pair.p2 = type_2;
				interaction_pair.c1 = std::stof(tokens[3]);
				interaction_pair.c2 = -1.0f / (std::stof(tokens[4]) / settings.cutoff);
				pair_list[pair_list_index] = interaction_pair;
			}
		}
	}

	for (int ii = 0; ii < pair_list.size(); ii++)
	{
		if (pair_list[ii].p1 == -1 && pair_list[ii].p2 == -1)
		{
			std::cout << "Error: pair list contains undefined interactions!" << std::endl;
			return false;
		}
	}

	if (pairs_num > 0)
	{
		BufferIndex::bead_pair_table_buffer = openClManager.CreateBuffer<Pair>(pair_list, CL_MEM_READ_ONLY);
	}
	return true;
};

bool SimulationManager::GetSimulationSteps(std::vector<std::string>& file_contents)
{
	std::vector<std::string>	tokens;

	if (ScanText(file_contents, "run", tokens, 0) < 0)
	{
		std::cout << "Error:  number of simalation steps was not specified." << std::endl;
		return false;
	}
	else
	{
		if (tokens.size() != 2)
		{
			std::cout << "Error: incorrect value for data dump freqeancy!" << std::endl; //linia
			return false;
		}
		simulation_steps = std::stoi(tokens[1]);
	}
	

	if (ScanText(file_contents, "minimize", tokens, 0) < 0)
	{
		std::cout << "Error:  number of simalation steps was not specified." << std::endl;
		return false;
	}
	else
	{
		if (tokens.size() != 2)
		{
			std::cout << "Error:  number of simalation steps was not specified." << std::endl;
			return false;
		}
		minimalization_steps = std::stoi(tokens[1]);
	}

	return true;
};

bool SimulationManager::GetOutputDataInfo(std::vector<std::string>&	file_contents)
{

	std::vector<std::string>	tokens;

	if (ScanText(file_contents, "output", tokens, 0))
	{
		if (tokens.size() > 1)
		{
			dump_file_name = tokens[1] + "_coord.xyz";
			vel_dump_file_name = tokens[1] + "_vel.xyz";
			restart_file_name = tokens[1] + "_restart.res";
		}
		else
		{
			std::cout << "Error: dump file was not specified." << std::endl;
			return false;
		}
	}

	if (!SetDumpFile(dump_file_name)) return false;

	if (!SetDumpFile(vel_dump_file_name)) return false;

	if (!SetDumpFile(restart_file_name)) return false;


	if (ScanText(file_contents, "output_freq", tokens, 0))
	{
		if (tokens.size() > 1)
		{
			dump_frequency = std::stoi(tokens[1]);
			if (dump_frequency <= 0)
			{
				std::cout << "Error: incorrect value for data dump freqeancy!" << std::endl; //linia
				return false;
			}
		}
		else
		{
			std::cout << "Error:  data dump freqeancy was not specified." << std::endl;
			return false;
		}
	}

	std::vector<int> export_type_list;
	int my_start = 0;
	while ((my_start = ScanText(file_contents, "export", tokens, my_start)) != -1)
	{
		if (tokens.size() == 3)
		{
			if (tokens[1] == "type")
			{
				std::string type = tokens[2];
				if (type_id.find(type) == type_id.end())
				{
					std::cout << "Error: input file line - " << my_start << "\t" << " - export bead type id deos not match any of the defined types.!" << std::endl;
					return false;
				}

				int export_type = type_id.find(type)->second;
				export_type_list.push_back(export_type);
			}
		}
		else
		{
			std::cout << "Error: cutoff radius was not specified in the input file!" << std::endl;
			return false;
		}
	}

	if (!openClManager.GetStatus()) return false;
	if (BufferIndex::bead_info_buffer < 0) return false;

	std::vector<BeadInfo> results_particle_info(openClManager.GetBufferLength(BufferIndex::bead_info_buffer), BeadInfo());
	if (results_particle_info.size() == 0)
	{
		std::cout << "Error: DumpTopology function - unable to load buffer data!" << std::endl;
		return false;
	}

	if (!openClManager.ReadBuffer(BufferIndex::bead_info_buffer, &results_particle_info.front()))
	{
		std::cout << "Error: DumpTopology function - unable to load buffer data!" << std::endl;
		return false;
	}
	
	if (export_type_list.empty()) return false;
	for (auto jj = 0; jj < openClManager.GetBufferElements(BufferIndex::bead_info_buffer); jj++)
	{
		if (std::find(export_type_list.begin(), export_type_list.end(), results_particle_info[jj].type) != export_type_list.end())
		{
			auto type = results_particle_info[jj].type;
			if (type < 0 || type > type_list.size()) return false;

			export_bead_index.push_back(jj);
			export_bead_name.push_back(type_list[type].id);
		}
	}
	
	return true;
};

bool	SimulationManager::ReadInputFile(std::string input_file)
{	
	std::vector<std::string>	file_contents;
	std::cout << std::endl << std::endl << "%------ Reading simulation input file ------%" << std::endl << std::endl;
		
	if (!ReadTextTFileContent(input_file, file_contents))
	{
		std::cout << "Error: unable to open input file:" << input_file << std::endl;
		std::cout << "Closing program!" << std::endl;
		return false;
	}
			
	if (!GetSimulationSettings(file_contents)) return false;
	if (!GetSimulationSteps(file_contents)) return false;
	if (!GetBeadInfo(file_contents)) return false;
	if (!ReadModelData(file_contents)) return false;
	if (!ReadRestartFile(file_contents)) return false;
	if (!GetOutputDataInfo(file_contents)) return false;
				
	return true;
};


bool SimulationManager::ReadKernelData(std::string input_file)
{
	std::ifstream file(input_file);
	if (!file.is_open())	return false;
	if (!file.good())		return false;

	std::string	 kernel_code((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());

	std::vector<std::string> kernel_names;
	kernel_names.push_back("AdvancePosition");
	kernel_names.push_back("ResetGrid");
	kernel_names.push_back("ResetContactList");
	kernel_names.push_back("FillGrid");
	kernel_names.push_back("BuildContactList");
	kernel_names.push_back("InitContactList");
	kernel_names.push_back("ResolveContacts");
	kernel_names.push_back("AdvanceVelocity");
	kernel_names.push_back("ResetBeadEnergy");
	kernel_names.push_back("ResolveSprings");
	kernel_names.push_back("ResolveAngles");

	if (!openClManager.registerKernel(kernel_names, kernel_code))
	{
		return false;
	}
	auto bead_buffer = openClManager.GetBuffer(BufferIndex::bead_buffer);
	auto angle_buffer = openClManager.GetBuffer(BufferIndex::angle_buffer);
	auto random_buffer = openClManager.GetBuffer(BufferIndex::random_buffer);
	auto random_buffer_stride = openClManager.GetBuffer(BufferIndex::random_buffer_stride);
	auto spring_buffer = openClManager.GetBuffer(BufferIndex::spring_buffer);
	auto bead_info_buffer = openClManager.GetBuffer(BufferIndex::bead_info_buffer);
	auto grid_head_buffer = openClManager.GetBuffer(BufferIndex::grid_head_buffer);
	auto grid_tail_buffer = openClManager.GetBuffer(BufferIndex::grid_tail_buffer);
	auto grid_cell_nh_buffer = openClManager.GetBuffer(BufferIndex::grid_cell_nh_buffer);
	auto bead_nh_head_buffer = openClManager.GetBuffer(BufferIndex::bead_nh_head_buffer);
	auto bead_nh_tail_buffer = openClManager.GetBuffer(BufferIndex::bead_nh_tail_buffer);
	auto bead_pair_coeff_buffer = openClManager.GetBuffer(BufferIndex::bead_pair_coeff_buffer);
	auto bead_pair_table_buffer = openClManager.GetBuffer(BufferIndex::bead_pair_table_buffer);

	if (bead_buffer == nullptr) return false;
	if (random_buffer == nullptr) return false;
	if (grid_head_buffer == nullptr) return false;
	if (grid_tail_buffer == nullptr) return false;
	if (grid_cell_nh_buffer == nullptr) return false;
	if (bead_nh_head_buffer == nullptr) return false;
	if (bead_nh_tail_buffer == nullptr) return false;
	if (bead_pair_coeff_buffer == nullptr) return false;
	if (bead_pair_table_buffer == nullptr) return false;

	int beads_num = openClManager.GetBufferElements(BufferIndex::bead_buffer);
	int grid_size = openClManager.GetBufferElements(BufferIndex::grid_head_buffer);
	int list_size = openClManager.GetBufferElements(BufferIndex::bead_nh_head_buffer);

	if (!openClManager.SetKernelArg("AdvancePosition", 0, *bead_buffer)) return false;
	if (!openClManager.SetKernelArg("AdvancePosition", 1, *bead_info_buffer)) return false;
	if (!openClManager.SetKernelArg("AdvancePosition", 2, settings)) return false;
	if (!openClManager.SetKernelArg("AdvancePosition", 3, beads_num)) return false;
	openClManager.SetKernelGlobalWorkRange("AdvancePosition", openClManager.GetBufferLength(BufferIndex::bead_buffer));

	if (spring_buffer != nullptr)
	{
		int springs_num = openClManager.GetBufferElements(BufferIndex::spring_buffer);
		if (!openClManager.SetKernelArg("ResolveSprings", 0, *bead_buffer)) return false;
		if (!openClManager.SetKernelArg("ResolveSprings", 1, *spring_buffer)) return false;
		if (!openClManager.SetKernelArg("ResolveSprings", 2, springs_num)) return false;
		openClManager.SetKernelGlobalWorkRange("ResolveSprings", openClManager.GetBufferLength(BufferIndex::spring_buffer));
	}

	if (angle_buffer != nullptr)
	{
		//int springs_num = openClManager.GetBufferElements(BufferIndex::spring_buffer);
		if (!openClManager.SetKernelArg("ResolveAngles", 0, *angle_buffer)) return false;
		if (!openClManager.SetKernelArg("ResolveAngles", 1, *bead_buffer)) return false;
		if (!openClManager.SetKernelArg("ResolveAngles", 2, angles_num)) return false;
		openClManager.SetKernelGlobalWorkRange("ResolveAngles", openClManager.GetBufferLength(BufferIndex::angle_buffer));

	}

	if (!openClManager.SetKernelArg("ResetGrid", 0, *grid_head_buffer)) return false;
	if (!openClManager.SetKernelArg("ResetGrid", 1, grid_size)) return false;
	openClManager.SetKernelGlobalWorkRange("ResetGrid", openClManager.GetBufferLength(BufferIndex::grid_head_buffer));

	if (!openClManager.SetKernelArg("ResetContactList", 0, *bead_nh_head_buffer)) return false;
	if (!openClManager.SetKernelArg("ResetContactList", 1, list_size)) return false;
	openClManager.SetKernelGlobalWorkRange("ResetContactList", openClManager.GetBufferLength(BufferIndex::bead_nh_head_buffer));

	if (!openClManager.SetKernelArg("FillGrid", 0, *bead_buffer)) return false;
	if (!openClManager.SetKernelArg("FillGrid", 1, *grid_head_buffer)) return false;
	if (!openClManager.SetKernelArg("FillGrid", 2, *grid_tail_buffer)) return false;
	if (!openClManager.SetKernelArg("FillGrid", 3, settings)) return false;   //tu tez zmiana
	if (!openClManager.SetKernelArg("FillGrid", 4, beads_num)) return false;
	openClManager.SetKernelGlobalWorkRange("FillGrid", openClManager.GetBufferLength(BufferIndex::bead_buffer));

	if (!openClManager.SetKernelArg("BuildContactList", 0, *grid_head_buffer)) return false;
	if (!openClManager.SetKernelArg("BuildContactList", 1, *grid_tail_buffer)) return false;
	if (!openClManager.SetKernelArg("BuildContactList", 2, *grid_cell_nh_buffer)) return false;
	if (!openClManager.SetKernelArg("BuildContactList", 3, *bead_nh_head_buffer)) return false;
	if (!openClManager.SetKernelArg("BuildContactList", 4, *bead_nh_tail_buffer)) return false;
	if (!openClManager.SetKernelArg("BuildContactList", 5, *bead_buffer)) return false;
	if (!openClManager.SetKernelArg("BuildContactList", 6, settings)) return false;
	if (!openClManager.SetKernelArg("BuildContactList", 7, grid_size)) return false;
	openClManager.SetKernelGlobalWorkRange("BuildContactList", openClManager.GetBufferLength(BufferIndex::grid_head_buffer));

	if (!openClManager.SetKernelArg("InitContactList", 0, *bead_nh_head_buffer)) return false;
	if (!openClManager.SetKernelArg("InitContactList", 1, *bead_nh_tail_buffer)) return false;
	if (!openClManager.SetKernelArg("InitContactList", 2, *bead_info_buffer)) return false;
	if (!openClManager.SetKernelArg("InitContactList", 3, *bead_pair_coeff_buffer)) return false;
	if (!openClManager.SetKernelArg("InitContactList", 4, *bead_pair_table_buffer)) return false;
	if (!openClManager.SetKernelArg("InitContactList", 5, settings.types_num)) return false;
	if (!openClManager.SetKernelArg("InitContactList", 6, list_size)) return false;
	openClManager.SetKernelGlobalWorkRange("InitContactList", openClManager.GetBufferLength(BufferIndex::bead_nh_head_buffer));

	if (!openClManager.SetKernelArg("ResolveContacts", 0, *bead_nh_head_buffer)) return false;
	if (!openClManager.SetKernelArg("ResolveContacts", 1, *bead_nh_tail_buffer)) return false;
	if (!openClManager.SetKernelArg("ResolveContacts", 2, *bead_buffer)) return false;
	if (!openClManager.SetKernelArg("ResolveContacts", 3, *random_buffer)) return false;
	if (!openClManager.SetKernelArg("ResolveContacts", 4, *bead_pair_coeff_buffer)) return false;
	if (!openClManager.SetKernelArg("ResolveContacts", 5, *bead_pair_table_buffer)) return false;
	if (!openClManager.SetKernelArg("ResolveContacts", 6, settings)) return false;
	if (!openClManager.SetKernelArg("ResolveContacts", 7, list_size)) return false;
	openClManager.SetKernelGlobalWorkRange("ResolveContacts", openClManager.GetBufferLength(BufferIndex::bead_nh_head_buffer));

	if (!openClManager.SetKernelArg("AdvanceVelocity", 0, *bead_buffer)) return false;
	if (!openClManager.SetKernelArg("AdvanceVelocity", 1, *bead_info_buffer)) return false;
	if (!openClManager.SetKernelArg("AdvanceVelocity", 2, *random_buffer)) return false;
	if (!openClManager.SetKernelArg("AdvanceVelocity", 3, beads_num)) return false;
	if (!openClManager.SetKernelArg("AdvanceVelocity", 4, openClManager.GetBufferLength(BufferIndex::bead_buffer))) return false;
	openClManager.SetKernelGlobalWorkRange("AdvanceVelocity", openClManager.GetBufferLength(BufferIndex::bead_buffer));
	//openClManager.SetKernelGlobalWorkRange("AdvanceVelocity", );

	if (!openClManager.SetKernelArg("ResetBeadEnergy", 0, *bead_buffer)) return false;
	if (!openClManager.SetKernelArg("ResetBeadEnergy", 1, beads_num)) return false;
	openClManager.SetKernelGlobalWorkRange("ResetBeadEnergy", openClManager.GetBufferLength(BufferIndex::bead_buffer));

	return true;
};



bool  SimulationManager::ReadBeadsData(std::vector<std::string>&	file_contents)
{	
	std::vector<Bead>					bead_buffer;
	std::vector<BeadInfo>				bead_info_buffer;
	std::vector<std::string>		    tokens;
	int									my_start = 0;
	
	bead_buffer.reserve(20000);
	bead_info_buffer.reserve(20000);

	while ((my_start = ScanText(file_contents, "BEAD", tokens, my_start)) != -1)
	{
		if (tokens[0] == "BEAD")
		{
			if (tokens.size() == 5 || tokens.size() == 8) //iloæ pól
			{
				BeadInfo		bead_info;
				Bead			bead;
				
				bead.position.x = std::stof(tokens[2]);
				bead.position.y = std::stof(tokens[3]);
				bead.position.z = std::stof(tokens[4]);
						
				if (settings.units == 0)
				{
					bead.position.x *= settings.cutoff_inv;
					bead.position.y *= settings.cutoff_inv;
					bead.position.z *= settings.cutoff_inv;
				}

				if (tokens.size() == 5)
				{
					bead.image_position = bead.position;
				}

				if (tokens.size() == 8)
				{					
					bead.image_position.x = std::stof(tokens[5]);
					bead.image_position.y = std::stof(tokens[6]);
					bead.image_position.z = std::stof(tokens[7]);
					
					if (settings.units == 0)
					{
						bead.image_position.x *= settings.cutoff_inv;
						bead.image_position.y *= settings.cutoff_inv;
						bead.image_position.z *= settings.cutoff_inv;
					}
				}
				


				if (bead.position.x > settings.size.x ||
					bead.position.y > settings.size.y ||
					bead.position.z > settings.size.z ||
					bead.position.x < 0.0f ||
					bead.position.y < 0.0f ||
					bead.position.z < 0.0f)
				{
					std::cout << "Error: model geometry does not fit the simulation box size." << std::endl;
					return false;
				}

				bead.velocity.x = 0.0f;
				bead.velocity.y = 0.0f;
				bead.velocity.z = 0.0f;
				
				bead.velocity_old.x = 0.0f;
				bead.velocity_old.y = 0.0f;
				bead.velocity_old.z = 0.0f;

				bead.force.x = 0.0f;
				bead.force.y = 0.0f;
				bead.force.z = 0.0f;

				bead.force_old.x = 0.0f;
				bead.force_old.y = 0.0f;
				bead.force_old.z = 0.0f;

				auto it = type_id.find(tokens[1]);
				if (it == type_id.end())
				{
					std::cout << "Error: model file line - " << my_start << "\t" 
						      << " - bead type does not match any of defined types." << std::endl;
					return false; 
				}
				
				auto index = it->second;

				bead_info.type = index;
				bead_info.mass = type_list[index].mass;
				
				if (bead_info.mass > 0.0f)
					bead_info.mass_inv = 1.0f / bead_info.mass;
				else
					bead_info.mass_inv = 0.0f;

				bead_info.lambda_dt_mass_inv = settings.lambda_dt * bead_info.mass_inv;
				bead_info.half_dt_mass_inv = settings.half_dt * bead_info.mass_inv;
				bead_info.half_dt2_mass_inv = settings.half_dt2 * bead_info.mass_inv;

				bead_buffer.push_back(bead);
				bead_info_buffer.push_back(bead_info);
			}
			else
			{
				std::cout << "Erro: invalid Particle data definition at line - "
					      << my_start << " - solver terminates." << std::endl;
				//system("PAUSE");
				return false;
			}
		}
	}

	beads_num = static_cast<int>(bead_buffer.size());

	if (beads_num > 0)
	{		
		//std::random_device					random_generator;
		//std::mt19937						generator(random_generator());
		//std::normal_distribution<float>		toss(0.0, 1.0);


		if ((BufferIndex::bead_buffer = openClManager.CreateBuffer<Bead>(bead_buffer, CL_MEM_READ_WRITE)) < 0)
		{
			std::cout << "Error: SetSimulationGrid function - unable to load BufferIndex::particles !" << std::endl;
			return false;
		}

		float avg = 0.0f;
		float	factor = 1.0f / sqrt(settings.dt);
		std::vector<InteractionNh>	random_coeff(openClManager.GetBufferLength(BufferIndex::bead_buffer), InteractionNh());

		for (int ii = 0; ii < random_coeff.size(); ii++)
		{
			for (int jj = 0; jj < 256; jj++)
			{
				random_coeff[ii].n[jj] = factor * toss(generator);
				avg = avg + random_coeff[ii].n[jj];
			}
			
		}

		avg /= random_coeff.size()*256; 
		for (int ii = 0; ii < random_coeff.size(); ii++)
		{
			for (int jj = 0; jj < 256; jj++)
			{
				random_coeff[ii].n[jj] -= avg;
			}
		}



		if ((BufferIndex::bead_info_buffer = openClManager.CreateBuffer<BeadInfo>(bead_info_buffer, CL_MEM_READ_ONLY)) < 0)
		{
			std::cout << "Error: SetSimulationGrid function - unable to load BufferIndex::particles_info !" << std::endl;
			return false;
		}


		//zmiana
		if ((BufferIndex::random_buffer = openClManager.CreateBuffer<InteractionNh>(random_coeff, CL_MEM_READ_WRITE)) < 0)
		{
			std::cout << "Error: SetSimulationGrid function - unable to load BufferIndex::randomBuffer !" << std::endl;
			return false;
		}
		else
		{
			random_check_data.resize(openClManager.GetBufferLength(BufferIndex::random_buffer));
		}
		if ((BufferIndex::bead_nh_head_buffer = openClManager.CreateBuffer<int>(beads_num, CL_MEM_READ_WRITE)) < 0)
		{
			std::cout << "Error: SetSimulationGrid function - unable to load BufferIndex::listHead !" << std::endl;
			return false;
		}
		else
		{
			collision_check_data.resize(openClManager.GetBufferLength(BufferIndex::bead_nh_head_buffer));
		}
		
		if ((BufferIndex::bead_nh_tail_buffer = openClManager.CreateBuffer<ParticleNh>(beads_num, CL_MEM_READ_WRITE)) < 0)
		{
			std::cout << "Error: SetSimulationGrid function - unable to load BufferIndex::listTail  !" << std::endl;
			return false;
		}

		if ((BufferIndex::bead_pair_coeff_buffer = openClManager.CreateBuffer<ParticleNh>(beads_num, CL_MEM_READ_WRITE)) < 0) //zmienione z InteractionNH na ParticleNh
		{
			std::cout << "Error: SetSimulationGrid function - unable to load BufferIndex::intercationCoeffs !" << std::endl;
			return false;
		}
		
		std::cout << "\tBead data uploaded succesfully." << std::endl;
		return true;
	}
	else
	{
		std::cout << "Error: model data requires bead definition." << std::endl;
		return false;
	}
}

bool  SimulationManager::ReadSpringData(std::vector<std::string>&	file_contents)
{
	std::vector<SpringStructure>	spring_data;
	std::vector<std::string>		tokens;
	int								my_start = 0;

	spring_data.reserve(2000);
	while ((my_start = ScanText(file_contents, "SPRING", tokens, my_start)) != -1)
	{
		if (tokens[0] == "SPRING")
		{
			if (tokens.size() == 5)
			{
				SpringStructure spring;
				spring.p1 = std::stoi(tokens[1]);
				if (spring.p1 < 0 && spring.p1 >= beads_num)
				{
					std::cout << "Error: model file line - " << my_start << "\t"
						<< " - wrong bead index reference." << std::endl;
					return false;
				}
				
				spring.p2 = std::stoi(tokens[2]);
				if (spring.p2 < 0 && spring.p2 >= beads_num)
				{
					std::cout << "Error: model file line - " << my_start << "\t"
						<< " - wrong bead index reference." << std::endl;
					return false;
				}

				if (settings.units == 0)
				{
					spring.rest_length = std::stof(tokens[3]) * settings.cutoff_inv;
				}
				else
				{
					spring.rest_length = std::stof(tokens[3]);
				}
				
				if (spring.rest_length < 0.0f)
				{
					std::cout << "Error: model file line - " << my_start << "\t"
						<< " - wrong rest length value." << std::endl;
					return false;
				}

				spring.stiffness = std::stof(tokens[4]);
				if (spring.stiffness < 0.0f)
				{
					std::cout << "Error: model file line - " << my_start << "\t"
						<< " - wrong stiffness value." << std::endl;
					return false;
				}

				spring_data.push_back(spring);
			}
			else
			{
				std::cout << "Error: invalid spring data definition at line - "
					<< my_start << " - solver terminates." << std::endl;
				//system("PAUSE");
				return false;
			}
		}
	}
	if (spring_data.size() > 0)
	{
		springs_num = static_cast<int>(spring_data.size());
		BufferIndex::spring_buffer = openClManager.CreateBuffer<SpringStructure>(spring_data, CL_MEM_READ_ONLY);
		if (BufferIndex::spring_buffer >= 0)
		{
			std::cout << "\tSpring bonds data uploaded succesfully." << std::endl;
			return true;
		}
		else
		{
			std::cout << "Error: unable to upload spring bonds data." << std::endl;
			return false;
		}
	}
	else
	{
		return true;
	}
};

bool  SimulationManager::ReadAngleData(std::vector<std::string>&	file_contents)
{
	std::vector<AngleBond>			angle_data;
	std::vector<std::string>		tokens;
	int								my_start = 0;

	angle_data.reserve(2000);
	while ((my_start = ScanText(file_contents, "ANGLE", tokens, my_start)) != -1)
	{
		if (tokens[0] == "ANGLE")
		{
			if (tokens.size() == 7)
			{
				AngleBond angle_bond;
				
				angle_bond.type = std::stoi(tokens[1]);

				angle_bond.b1 = std::stoi(tokens[2]);
				if (angle_bond.b1 < 0 && angle_bond.b1 >= beads_num)
				{
					std::cout << "Error: model file line - " << my_start << "\t"
						<< " - wrong bead index reference." << std::endl;
					return false;
				}

				angle_bond.b2 = std::stoi(tokens[3]);
				if (angle_bond.b2 < 0 && angle_bond.b2 >= beads_num)
				{
					std::cout << "Error: model file line - " << my_start << "\t"
						<< " - wrong bead index reference." << std::endl;
					return false;
				}
				
				angle_bond.b3 = std::stoi(tokens[4]);
				if (angle_bond.b3 < 0 && angle_bond.b3 >= beads_num)
				{
					std::cout << "Error: model file line - " << my_start << "\t"
						<< " - wrong bead index reference." << std::endl;
					return false;
				}
				
				angle_bond.c1 = std::stof(tokens[5]);
				angle_bond.angle = std::stof(tokens[6]);
				

				angle_data.push_back(angle_bond);
			}
			else
			{
				std::cout << "Error: invalid spring data definition at line - "
					<< my_start << " - solver terminates." << std::endl;
				//system("PAUSE");
				return false;
			}
		}
	}
	if (angle_data.size() > 0)
	{
		angles_num = static_cast<int>(angle_data.size());
		BufferIndex::angle_buffer = openClManager.CreateBuffer<AngleBond>(angle_data, CL_MEM_READ_ONLY);
		if (BufferIndex::angle_buffer >= 0)
		{
			std::cout << "\tAngle bonds data uploaded succesfully." << std::endl;
			return true;
		}
		else
		{
			std::cout << "Error: unable to upload angle bonds data." << std::endl;
			return false;
		}
	}
	else
	{
		return true;
	}
};



bool SimulationManager::ReadRestartFile(std::vector<std::string>& input_file_content)
{
	std::string					restart_file_name;
	std::vector<std::string>	restart_file_content;
	std::vector<std::string>	tokens;

	if (!openClManager.GetStatus())
	{
		return false;
	}

	if (ScanText(input_file_content, "restart", tokens, 0))
	{
		if (tokens.size() > 1)
		{
			restart_file_name = tokens[1];
			tokens.clear();
		}
		else
		{
			std::cout << "Restart data file was not specified." << std::endl;
			return true;
		}
	}

	std::cout << "\tReading restart data from file:  " << restart_file_name << std::endl;
	if (!ReadTextTFileContent(restart_file_name, restart_file_content))
	{
		std::cout << "Error: unable to open model file: " << restart_file_name << std::endl;
		system("PAUSE");
		return false;
	}


	if (!restart_file_content.empty())
	{
		std::vector<Bead> bead_data(openClManager.GetBufferLength(BufferIndex::bead_buffer), Bead());
		if (!openClManager.ReadBuffer(BufferIndex::bead_buffer, &bead_data.front())) return false;
		int bead_index = 0;

		for (int ii = 0; ii < static_cast<int>(restart_file_content.size()); ii++)
		{
			auto str = restart_file_content[ii];
			std::istringstream iss(str);
			std::copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(), std::back_inserter(tokens));
			if (!tokens.empty())
			{
				if (tokens.size() == 12)
				{
					if (bead_index < bead_data.size())
					{
						bead_data[bead_index].position.x = std::stof(tokens[0]);
						bead_data[bead_index].position.y = std::stof(tokens[1]);
						bead_data[bead_index].position.z = std::stof(tokens[2]);

						bead_data[bead_index].image_position.x = std::stof(tokens[3]);
						bead_data[bead_index].image_position.y = std::stof(tokens[4]);
						bead_data[bead_index].image_position.z = std::stof(tokens[5]);

						bead_data[bead_index].velocity.x = std::stof(tokens[6]);
						bead_data[bead_index].velocity.y = std::stof(tokens[7]);
						bead_data[bead_index].velocity.z = std::stof(tokens[8]);

						bead_data[bead_index].force.x = std::stof(tokens[9]);
						bead_data[bead_index].force.y = std::stof(tokens[10]);
						bead_data[bead_index].force.z = std::stof(tokens[11]);
					}
					bead_index++;
				}
			}
			tokens.clear();
		}

		if (!openClManager.WriteBuffer<Bead>(BufferIndex::bead_buffer, bead_data)) return false;

		return true;
	}
	else
	{
		return false;
	}
};

bool SimulationManager::ReadModelData(std::vector<std::string>& input_file_content)
{
	std::string					model_file_name;
	std::vector<std::string>	model_file_content;
	std::vector<std::string>	tokens;

	if (!openClManager.GetStatus())
	{
		return false;
	}

	if (ScanText(input_file_content, "model", tokens, 0))
	{
		if (tokens.size() > 1)
		{
			model_file_name = tokens[1];
		}
		else
		{
			std::cout << "Error: model data file was not specified in the input file!" << std::endl;
			return false;
		}
	}
		
	std::cout << "\tReading model data from file:  " << model_file_name << std::endl;
	if (!ReadTextTFileContent(model_file_name, model_file_content))
	{
		std::cout << "Error: unable to open model file: " << model_file_name << std::endl;
		system("PAUSE");
		return false;
	}
	
	if (model_file_content.empty())
	{
		std::cout << "Model file is empty. Solver terminates." << std::endl;
		system("PAUSE");
		return false;
	}
	
	if (!ReadBeadsData(model_file_content)) return false;
	if (!ReadSpringData(model_file_content)) return false;
	if (!ReadAngleData(model_file_content)) return false;

	std::cout << "\tModel data succesfully uploaded." << std::endl;
	return true;	
};


bool SimulationManager::SetDumpFile(const std::string& value)
{
	std::cout << "\tOutput file: " << value << std::endl;
	if (std::ifstream(value))
	{
		std::cout << "\tWarning - output file already exists - overwriting existing file" << std::endl;
	}
	std::fstream dump_file;
		
	dump_file.open(value, std::ios::out);
	if (dump_file.is_open())
	{
		if (dump_file.good())
		{
			dump_file.close();
			//dump_file_name = value;
			return true;
		}
	}
	std::cout << "Error: Unable to set output file!" << std::endl;
	return false;
};


bool SimulationManager::SimulationCheck()
{
	if (openClManager.GetStatus())
	{
		if (BufferIndex::bead_nh_head_buffer != -1)
		{			
			openClManager.ReadBuffer(BufferIndex::random_buffer, &random_check_data.front());

			if (openClManager.ReadBuffer(BufferIndex::bead_nh_head_buffer, &collision_check_data.front()))
			{

				return true;
			}
			else
			{
				return false;
			}

			if (openClManager.ReadBuffer(BufferIndex::grid_head_buffer, &grid_check_data.front()))
			{
				int size = openClManager.GetBufferElements(BufferIndex::grid_head_buffer);
				for (int ii = 0; ii < size; ii++)
				{
					if (grid_check_data[ii] >= 30)	warning_msg_counter++;
				}
				return true;
			}
			else
			{
				return false;
			}

		}
		else
		{
			return false;
		}
	}
	else
	{
		return false;
	}
}
bool	SimulationManager::DumpRestartData()
{
	if (!openClManager.GetStatus())
	{
		std::cout << "Error: DumpRestartData() - dump flag is seto to false or openClManager is not initialized!" << std::endl;
		return false;
	}
	if (BufferIndex::bead_buffer < 0) return false;

	std::vector<Bead> bead_data(openClManager.GetBufferLength(BufferIndex::bead_buffer), Bead());
	std::fstream dump_file(restart_file_name, std::ios_base::app);
	
	if (!openClManager.ReadBuffer(BufferIndex::bead_buffer, &bead_data.front())) return false;
	if (!dump_file.is_open()) return false;
	if (!dump_file.good()) return false;

	dump_file << bead_data.size() << std::endl;
	dump_file << "FRAME" << std::endl;
	for (auto ii = 0; ii < beads_num; ii++) /// tyle ile rzeczywicie jest cz¹stek
	{
		dump_file << "\t" << bead_data[ii].position.x;
		dump_file << "\t" << bead_data[ii].position.y;
		dump_file << "\t" << bead_data[ii].position.z;

		dump_file << "\t" << bead_data[ii].image_position.x;
		dump_file << "\t" << bead_data[ii].image_position.y;
		dump_file << "\t" << bead_data[ii].image_position.z;
			
		dump_file << "\t" << bead_data[ii].velocity.x;
		dump_file << "\t" << bead_data[ii].velocity.y;
		dump_file << "\t" << bead_data[ii].velocity.z;

		dump_file << "\t" << bead_data[ii].force.x;
		dump_file << "\t" << bead_data[ii].force.y;
		dump_file << "\t" << bead_data[ii].force.z << std::endl;
	}
	return true;
}



bool SimulationManager::DumpData()
{
	if (!openClManager.GetStatus())
	{
		std::cout << "Error: DumpTopology function - dump flag is seto to false or openClManager is not initialized!" << std::endl;
		return false;
	}
	if (BufferIndex::bead_buffer < 0) return false;
	if (export_bead_index.empty()) return false;

	static int dump_counter = 1;
	static std::vector<Bead> bead_data(openClManager.GetBufferLength(BufferIndex::bead_buffer), Bead());
	static std::fstream dump_file(dump_file_name, std::ios_base::app);
	static std::fstream vel_dump_file(vel_dump_file_name, std::ios_base::app);
	
	if (!openClManager.ReadBuffer(BufferIndex::bead_buffer, &bead_data.front())) return false;
	if (!dump_file.is_open() && !vel_dump_file.is_open()) return false;
	if (!dump_file.good() && !vel_dump_file.good()) return false;

	dump_file << export_bead_index.size() << std::endl;
	dump_file << "FRAME" << "\t" << dump_counter*dump_frequency << "\t" << settings.cutoff << "\t" << settings.sigma << "\t"
		<< settings.gamma << "\t" << settings.dt << "\t" << settings.size.x << "\t" << settings.size.y << "\t" << settings.size.z << std::endl;

	for (auto ii = 0; ii < export_bead_index.size(); ii++)
	{
		dump_file << export_bead_name[ii];
		//dump_file << "\t" << bead_data[export_bead_index[ii]].image_position.x;
		//dump_file << "\t" << bead_data[export_bead_index[ii]].image_position.y;
		//dump_file << "\t" << bead_data[export_bead_index[ii]].image_position.z << std::endl;
		dump_file << "\t" << bead_data[export_bead_index[ii]].position.x;
		dump_file << "\t" << bead_data[export_bead_index[ii]].position.y;
		dump_file << "\t" << bead_data[export_bead_index[ii]].position.z << std::endl;
	}
	
	vel_dump_file << export_bead_index.size() << std::endl;
	vel_dump_file << "FRAME" << std::endl;
	for (auto ii = 0; ii < export_bead_index.size(); ii++)
	{
		vel_dump_file << export_bead_name[ii];
		vel_dump_file << "\t" << bead_data[export_bead_index[ii]].velocity.x;
		vel_dump_file << "\t" << bead_data[export_bead_index[ii]].velocity.y;
		vel_dump_file << "\t" << bead_data[export_bead_index[ii]].velocity.z << std::endl;
	}

	dump_counter++;
	return true;
};

bool SimulationManager::SetSimulationGrid()
{
	//ustalanie iloci komórek w siatce na podstawie max promienia cz¹steczki

	if (settings.size.x <= 0.0f || settings.size.y <= 0.0f || settings.size.z <= 0.0f)
	{
		std::cout << "Error: SetSimulationGrid function - invalid grid size data!" << std::endl;
		return false;
	}

	std::vector<GridCellNh>			grid_cell_nhood;

	if (settings.units == 0)
	{
		settings.size.x *= settings.cutoff_inv;
		settings.size.y *= settings.cutoff_inv;
		settings.size.z *= settings.cutoff_inv;

		settings.x_cells_num = static_cast<int>(floor(settings.size.x / 1.1f));
		settings.y_cells_num = static_cast<int>(floor(settings.size.y / 1.1f));
		settings.z_cells_num = static_cast<int>(floor(settings.size.z / 1.1f));

		settings.xy_cells_num = settings.x_cells_num * settings.y_cells_num;

		settings.step.x = settings.size.x / static_cast<float>(settings.x_cells_num);
		settings.step.y = settings.size.y / static_cast<float>(settings.y_cells_num);
		settings.step.z = settings.size.z / static_cast<float>(settings.z_cells_num);
	}
	else
	{
		settings.x_cells_num = static_cast<int>(floor(settings.size.x / (settings.cutoff*1.1f)));
		settings.y_cells_num = static_cast<int>(floor(settings.size.y / (settings.cutoff*1.1f)));
		settings.z_cells_num = static_cast<int>(floor(settings.size.z / (settings.cutoff*1.1f)));

		settings.xy_cells_num = settings.x_cells_num * settings.y_cells_num;

		settings.step.x = settings.size.x / static_cast<float>(settings.x_cells_num);
		settings.step.y = settings.size.y / static_cast<float>(settings.y_cells_num);
		settings.step.z = settings.size.z / static_cast<float>(settings.z_cells_num);
	}

	settings.half_size_squared = std::min<float>(std::min<float>(settings.size.x, settings.size.y), settings.size.z)*0.5f;
	settings.half_size_squared *= settings.half_size_squared;

	//wielkoæ siatki
	int grid_size = settings.x_cells_num * settings.y_cells_num * settings.z_cells_num;

	grid_cell_nhood.resize(grid_size);

	//ustalanie adresów komórek s¹siaduj¹cych 
	for (int zz = 0; zz < settings.z_cells_num; zz++)
	{
		for (int yy = 0; yy < settings.y_cells_num; yy++)
		{
			for (int xx = 0; xx < settings.x_cells_num; xx++)
			{
				int targetCell = xx + yy * settings.x_cells_num + zz * settings.x_cells_num * settings.y_cells_num;
				int cnt = 0;
				for (int dx = -1; dx < 2; dx++)
				{
					for (int dy = -1; dy < 2; dy++)
					{
						int nx = xx + dx;
						int ny = yy + dy;
						int nz = zz - 1;

						if (nx < 0) nx = settings.x_cells_num - 1;
						if (nx >= settings.x_cells_num) nx = 0;

						if (ny < 0) ny = settings.y_cells_num - 1;
						if (ny >= settings.y_cells_num) ny = 0;

						if (nz < 0) nz = settings.z_cells_num - 1;
						if (nz >= settings.z_cells_num) nz = 0;

						int nh = nx + ny * settings.x_cells_num + nz * settings.x_cells_num * settings.y_cells_num;
						grid_cell_nhood[targetCell].n[cnt] = nh;
						cnt++;
					}
				}

				for (int dy = -1; dy < 2; dy++)
				{
					int nx = xx + 1;
					int ny = yy + dy;
					int nz = zz;

					if (nx < 0) nx = settings.x_cells_num - 1;
					if (nx >= settings.x_cells_num) nx = 0;

					if (ny < 0) ny = settings.y_cells_num - 1;
					if (ny >= settings.y_cells_num) ny = 0;

					if (nz < 0) nz = settings.z_cells_num - 1;
					if (nz >= settings.z_cells_num) nz = 0;

					int nh = nx + ny * settings.x_cells_num + nz * settings.x_cells_num * settings.y_cells_num;
					grid_cell_nhood[targetCell].n[cnt] = nh;
					cnt++;
				}

				int nx = xx;
				int ny = yy - 1;
				int nz = zz;

				if (nx < 0) nx = settings.x_cells_num - 1;
				if (nx >= settings.x_cells_num) nx = 0;

				if (ny < 0) ny = settings.y_cells_num - 1;
				if (ny >= settings.y_cells_num) ny = 0;

				if (nz < 0) nz = settings.z_cells_num - 1;
				if (nz >= settings.z_cells_num) nz = 0;

				int nh = nx + ny * settings.x_cells_num + nz * settings.x_cells_num * settings.y_cells_num;
				grid_cell_nhood[targetCell].n[cnt] = nh;
				cnt++;
			}
		}
	}


	BufferIndex::grid_head_buffer = openClManager.CreateBuffer<int>(grid_size, CL_MEM_READ_WRITE);
	if (BufferIndex::grid_head_buffer < 0)
	{
		std::cout << "Error: SetSimulationGrid function - unable to load BufferIndex::gridHead!" << std::endl;
		return false;
	}
	else
	{
		grid_check_data.resize(openClManager.GetBufferLength(BufferIndex::grid_head_buffer));
	}

	BufferIndex::grid_tail_buffer = openClManager.CreateBuffer<GridCell>(grid_size, CL_MEM_READ_WRITE);
	if (BufferIndex::grid_tail_buffer < 0)
	{
		std::cout << "Error: SetSimulationGrid function - unable to load BufferIndex::gridTail!" << std::endl;
		return false;
	}

	BufferIndex::grid_cell_nh_buffer = openClManager.CreateBuffer<GridCellNh>(grid_cell_nhood, CL_MEM_READ_ONLY);
	if (BufferIndex::grid_cell_nh_buffer < 0)
	{
		std::cout << "Error: SetSimulationGrid function - unable to load BufferIndex::nHoodBuffer!" << std::endl;
		return false;
	}

	return true;
};

void SimulationManager::Close()
{
	//if (dump_file.is_open()) dump_file.close();
};
