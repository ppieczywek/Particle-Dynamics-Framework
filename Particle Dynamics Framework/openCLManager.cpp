#pragma once
#include "stdafx.h"
#include "openCLManager.h"




bool OpenClManager::Initialize()
{
	localWorkGroupSize = 64;

	std::vector<cl::Platform> all_platforms;
	cl::Platform::get(&all_platforms);

	if (all_platforms.size() == 0)
	{
		std::cout << "\tOpenClManager:" << "\t\t" << "No platforms found. Check OpenCL installation!\n";
		return false;
	}

	for (int ii = 0; ii < all_platforms.size(); ii++)
	{
		std::cout << "\t[" << ii + 1 << "]" << "\t" << all_platforms[ii].getInfo<CL_PLATFORM_NAME>() << std::endl;
	}
	std::cout << "\t[" << all_platforms.size()+1 << "]" << "\t" << "terminate program" << std::endl;
		
	int input = 0;
	while (1)
	{
		//std::cin >> input;
		input = 1;

		if ((all_platforms.size() + 1) == input)
		{
			return false;
		}
		if ((input-1) >= 0 && (input-1) < all_platforms.size())
		{
			default_platform = all_platforms[static_cast<size_t>(input-1)];
			std::cout << "\tSelected platform:\t" << default_platform.getInfo<CL_PLATFORM_NAME>() << "\r";
			break;
		}
		
	}
	

	

	std::vector<cl::Device> all_devices;
	default_platform.getDevices(CL_DEVICE_TYPE_GPU, &all_devices);
	if (all_devices.size() == 0)
	{
		std::cout << "\tOpenClManager:" << "\t\t" << " No GPU devices found !\n";
		default_platform.getDevices(CL_DEVICE_TYPE_CPU, &all_devices);
		if (all_devices.size() == 0)
		{
			return false;
		}
	}

	default_device = all_devices[0];
	std::cout << "\tSelected device:\t" << default_device.getInfo<CL_DEVICE_NAME>() << std::endl;

	

	global_mem_avalilable =  default_device.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>();
	global_mem_avalilable /= (1024 * 1024);
	std::cout << "\tDevice memory:\t\t" << global_mem_avalilable << " MB" << std::endl << std::endl << std::endl;

//	cl_int err = default_device.getInfo(CL_DEVICE_GLOBAL_MEM_SIZE, &size);
//	if (err == CL_SUCCESS)
//	{
//		global_mem_avalilable = size;
//	
//		global_mem_allocated = 0;
//	}
	

	
	

//#define CL_DEVICE_GLOBAL_MEM_SIZE                   0x101F
//#define CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE          0x1020
//#define CL_DEVICE_MAX_CONSTANT_ARGS                 0x1021
//#define CL_DEVICE_LOCAL_MEM_TYPE                    0x1022
//#define CL_DEVICE_LOCAL_MEM_SIZE                    0x1023


	context = cl::Context({ default_device });
	queue = cl::CommandQueue(context, default_device);

	managerStatus = true;
	return true;
};

bool OpenClManager::registerKernel(const std::vector<std::string> &kernelNames, const std::string &kernelCode)
{
	if (!managerStatus)
	{
		std::cout << "OpenClManager:	manager not initialized" << std::endl;
	}
	
	cl::Program::Sources sources;
	sources.push_back({ kernelCode.c_str(), kernelCode.length() });
	std::string options = "-cl-mad-enable";
	
	cl::Program program(context, sources);
	if (program.build({ default_device },options.c_str(), NULL,NULL) != CL_SUCCESS)
	{
		std::cout << " Error building: " << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(default_device) << "\n";
		return false;
	}

	for (auto ii = 0; ii < kernelNames.size(); ii++)
	{
		cl_int err;
		cl::Kernel my_kernel(program, kernelNames[ii].c_str(), &err);
		if (err == CL_SUCCESS)
		{
			//std::cout << "Kernel registered - " << kernelNames[ii] << std::endl;
			kernels.push_back(my_kernel);
			kernelsNames.push_back(kernelNames[ii]);
			kernelGlobalWorkRange.push_back(0);
		}
		else
		{
			checkErr(err);
			return false;
		}
	}
	return true;
};

cl::Buffer* OpenClManager::GetBuffer(const int& bufferId)
{
	if (!clBuffer.empty())
	{
		if (bufferId < clBuffer.size())
		{
			size_t size = 0;
			clBuffer[bufferId].getInfo(CL_MEM_SIZE, &size);
			if (size != 0)
			{
				return &clBuffer[bufferId];
			}
		}
	}
	return nullptr;
};

int OpenClManager::GetLocalWorkGroupSize()
{
	return localWorkGroupSize;
}

int OpenClManager::GetBufferSize(const int& bufferId)
{
	if (!clBuffer.empty())
	{
		if (bufferId < clBuffer.size())
		{
			size_t size = 0;
			clBuffer[bufferId].getInfo(CL_MEM_SIZE, &size);
			if (size != 0)
			{
				return clBufferSize[bufferId];
			}
		}
	}
	return -1;
};

int OpenClManager::GetBufferElements(const int& bufferId)
{
	if (!clBuffer.empty())
	{
		if (bufferId < clBuffer.size())
		{
			size_t size = 0;
			clBuffer[bufferId].getInfo(CL_MEM_SIZE, &size);
			if (size != 0)
			{
				return clBufferElements[bufferId];
			}
		}
	}
	return -1;
};


int OpenClManager::GetBufferLength(const int& bufferId)
{
	if (!clBuffer.empty())
	{
		if (bufferId < clBuffer.size())
		{
			size_t size = 0;
			clBuffer[bufferId].getInfo(CL_MEM_SIZE, &size);
			if (size != 0)
			{
				return clBufferLength[bufferId];
			}
		}
	}
	return -1;
};


int OpenClManager::GetBufferStructureSize(const int& bufferId)
{
	if (!clBuffer.empty())
	{
		if (bufferId < clBuffer.size())
		{
			size_t size = 0;
			clBuffer[bufferId].getInfo(CL_MEM_SIZE, &size);
			if (size != 0)
			{
				return clBufferStructureSize[bufferId];
			}
		}
	}
	return -1;
};



int OpenClManager::LoadBuffer(const void* ptr, unsigned structureSize, int bufferSize, int bufferElements, cl_mem_flags flag)
{
	cl::Buffer buffer(context, flag, structureSize *  bufferSize);
	auto err = queue.enqueueWriteBuffer(buffer, CL_TRUE, 0, structureSize * bufferSize, ptr);
	if (err == CL_SUCCESS)
	{
		clBuffer.push_back(buffer);
		clBufferSize.push_back(bufferElements*structureSize);
		clBufferElements.push_back(bufferElements);
		clBufferLength.push_back(bufferSize);
		clBufferStructureSize.push_back(structureSize);

		return static_cast<int>(clBuffer.size() - 1);
	}
	else
	{
		std::cout << "OpenClManager:	error while loading buffer data" << std::endl;
		checkErr(err);
		return -1;
	}
};


bool OpenClManager::ReadBuffer(int bufferId, void* ptr)
{
	if (!clBuffer.empty())
	{
		if (bufferId < clBuffer.size())
		{
			size_t size = 0;
			clBuffer[bufferId].getInfo(CL_MEM_SIZE, &size);
			if (size != 0)
			{
				cl_int err = queue.enqueueReadBuffer(clBuffer[bufferId], 
													 CL_TRUE,
													 0,
													 clBufferLength[bufferId]*clBufferStructureSize[bufferId],
													 ptr);	
				if (err == CL_SUCCESS)
				{
					return true;
				}
				else
				{
					std::cout << "OpenClManager:	error while reading buffer data" << std::endl;
					checkErr(err);
					return false;
				}
			}
		}
	}
	return false;
};

bool OpenClManager::runKernels()
{
	if (managerStatus)
	{
		cl_int err = queue.finish();
		if (err == CL_SUCCESS)
		{
			return true;
		}
		else
		{
			checkErr(err);
			return false;
		}
	}
	else
	{
		std::cout << "OpenClManager:	manager not initialized" << std::endl;
		return false;
	}
}

bool OpenClManager::executeKernel(int kernelId)
{
	if (managerStatus)
	{
		if (kernelId >= 0 && kernelId < kernels.size())
		{
			if (kernelGlobalWorkRange[kernelId] > 0)
			{
				cl_int err = queue.enqueueNDRangeKernel(kernels[kernelId],
														cl::NullRange,
														cl::NDRange(kernelGlobalWorkRange[kernelId]),
														cl::NDRange(localWorkGroupSize),
														NULL,
														NULL);
				if (err == CL_SUCCESS)
				{
					return true;
				}
				else
				{
					std::cout << "OpenClManager:	error during execution of kernel:" << kernelsNames[kernelId] << std::endl;
					checkErr(err);
					return false;
				}
			}
			else
			{
				std::cout << "OpenClManager:	incorrect kernel work range" << std::endl;
				return false;
			}
		}
		else
		{
			std::cout << "OpenClManager:	kernel structure is empty" << std::endl;
			return false;
		}
	}
	else
	{
		std::cout << "OpenClManager:	manager not initialized" << std::endl;
		return false;
	}
};

bool OpenClManager::GetStatus()
{
	return managerStatus;
};


bool	OpenClManager::SetKernelArg(std::string name, cl_uint argId, size_t size, const void *argPtr)
{
	if (managerStatus)
	{
		int kernelId = -1;
		
		for (int ii = 0; ii < kernelsNames.size(); ii++)
		{
			if (name == kernelsNames[ii])
			{
				kernelId = ii;
			}
		}


		if (kernelId < kernels.size() && kernelId >= 0)
		{
			cl_int err = kernels[kernelId].setArg(argId, size, argPtr);
			if (err != CL_SUCCESS)
			{
				std::cout << "\t" << "unable to set kernel argument" << std::endl;
				std::cout << "\t" << "Kernel id: " << kernelId << "\t" << "Argument id: " << argId << std::endl;
				return false;
			}
			else
			{
				return true;
			}
		}
		else
		{
			std::cout <<  "\t" << "unable to set kernel argument - system not initialized" << std::endl;
			return false;
		}
	}
	else
	{
		return false;
	}
};

int	OpenClManager::GetNumberOfKernels()
{
	if (!kernels.empty())
	{
		return static_cast<int>(kernels.size());
	}
	else
	{
		return -1;
	}
}

void	OpenClManager::SetKernelGlobalWorkRange(std::string name, int workRange)
{
	if (managerStatus)
	{
		int kernelId = -1;

		for (int ii = 0; ii < kernelsNames.size(); ii++)
		{
			if (name == kernelsNames[ii])
			{
				kernelId = ii;
			}
		}

		if (!kernels.empty())
		{
			if (kernelId != -1 && kernelId < kernels.size())
			{
				kernelGlobalWorkRange[kernelId] = workRange;
			}
			else
			{
				std::cout << "Wrong kernel ID invoked!!" << std::endl;
			}
		}
		else
		{
			std::cout << "There are no kernels present in this system!!" << std::endl;
		}
	}
};

void OpenClManager::Close()
{
	std::cout << "OpenClManager:" << "\t" << "manager cleanup...." << std::endl;
};


int	OpenClManager::GetKernelId(std::string name)
{
	if (managerStatus)
	{
		for (int ii = 0; ii < kernelsNames.size(); ii++)
		{
			if (name == kernelsNames[ii])
			{
				return ii;
			}
		}
		return -1;
	}
	else
	{
		return -1;
	}
}

void OpenClManager::checkErr(cl_int err)
{
	switch (err)
	{
		case CL_INVALID_COMMAND_QUEUE:
			std::cout << "command_queue is not a valid command-queue" << std::endl;
			break;
			
		case CL_INVALID_CONTEXT:
			std::cout << "context associated with command_queue and buffer are not the	same or if the context associated with command_queue and events in	event_wait_list are not the same" << std::endl;
			break;

		case CL_INVALID_MEM_OBJECT:
			std::cout << "buffer is not a valid buffer object" << std::endl;
			break;

		case CL_INVALID_VALUE:
			std::cout << "the region being read specified by (offset, cb) is out of bounds or if ptr is a NULL value" << std::endl;
			break;
		
		case CL_INVALID_EVENT_WAIT_LIST:
			std::cout << "event_wait_list is NULL and num_events_in_wait_list greater than 0, or event_wait_list is not NULL and num_events_in_wait_list is 0, or if event objects in event_wait_list are not valid events" << std::endl;
			break;

		case CL_MEM_OBJECT_ALLOCATION_FAILURE:
			std::cout << "there is a failure to allocate memory for data store associated with buffer" << std::endl;
			break;

		case CL_OUT_OF_HOST_MEMORY:
			std::cout << "there is a failure to allocate resources required by the OpenCL implementation on the host" << std::endl;
			break;

		case CL_OUT_OF_RESOURCES:
			std::cout << "out of resources" << std::endl;
			break;

		case CL_INVALID_PROGRAM:
			std::cout << "OpenClManager:	program is not a valid program object!" << std::endl;
			break;
	
		case CL_INVALID_PROGRAM_EXECUTABLE:
			std::cout << "OpenClManager:	there is no successfully built executable for program!" << std::endl;
			break;
	
		case CL_INVALID_KERNEL_NAME:
			std::cout << "OpenClManager:	 kernel_name is not found in program!" << std::endl;
			break;
				
		case CL_INVALID_KERNEL:
			std::cout << "OpenClManager:	kernel is not a valid kernel object" << std::endl;
			break;

		case CL_INVALID_KERNEL_ARGS:
			std::cout << "OpenClManager:	the kernel argument values have not been specified" << std::endl;
			break;
		
		case CL_INVALID_WORK_DIMENSION:
				std::cout << "OpenClManager:	work_dim is not a valid value(i.e.a value between 1 and 3)" << std::endl;
			break;
		
		case CL_INVALID_WORK_GROUP_SIZE:
			std::cout << "OpenClManager:	local_work_size is specified and number of work - items specified by global_work_size is not evenly divisable by size of work - group given by local_work_size or does not match the work - group size specified for kernel using the __attribute__((reqd_work_group_size(X, Y, Z))) qualifier in program source" << std::endl;
			break;

		case CL_INVALID_WORK_ITEM_SIZE:
			std::cout << "OpenClManager:	the number of work - items specified in any of local_work_size[0], ... local_work_size[work_dim - 1] is greater than the corresponding values specified by CL_DEVICE_MAX_WORK_ITEM_SIZES[0], ....CL_DEVICE_MAX_WORK_ITEM_SIZES[work_dim - 1]" << std::endl;
			break;

		case CL_INVALID_GLOBAL_OFFSET:
			std::cout << "OpenClManager:	global_work_offset is not NULL" << std::endl;
			break;

		default:
			std::cout << "OpenClManager:	unknown error code" << std::endl;
			break;
	}
}