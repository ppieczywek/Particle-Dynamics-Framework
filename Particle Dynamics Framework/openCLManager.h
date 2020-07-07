#pragma once
#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <tuple>
#include <conio.h>
#include "CL/cl.hpp"
//#include "TypeCodeInitializer.h"
#include "BufferIndex.h"

class OpenClManager
{
	bool												managerStatus;
	int													localWorkGroupSize;
	unsigned long long									global_mem_avalilable;
	unsigned long long									global_mem_allocated;

	std::vector<std::string>							kernelsNames;
	std::vector<cl::Kernel>								kernels;
	std::vector<int>									kernelGlobalWorkRange;

	std::vector<cl::Buffer>								clBuffer;
	std::vector<int>									clBufferSize; 
	std::vector<int>									clBufferElements;
	std::vector<int>									clBufferLength;
	std::vector<int>									clBufferStructureSize;

	cl::Device											default_device;
	cl::Platform										default_platform;
	cl::Context											context;
	cl::CommandQueue									queue;

	void												checkErr(cl_int err);

public:

	bool												GetStatus();
	bool												registerKernel(const std::vector<std::string> &kernelNames, const std::string &kernelCode);
	bool												executeKernel(int kernelId);
	bool												runKernels();

	bool												Initialize();
	cl::Buffer*											GetBuffer(const int&);
	int													LoadBuffer(const void* ptr, unsigned, int, int, cl_mem_flags flag);
	bool												ReadBuffer(int, void*);
	int													GetBufferSize(const int&);
	int													GetBufferElements(const int&);
	int													GetBufferLength(const int&);
	int													GetBufferStructureSize(const int&);
	int													GetLocalWorkGroupSize();
	int													GetNumberOfKernels();

	template<class systemType> bool						SetKernelArg(std::string name, cl_uint argId, const systemType &argPtr);

	template<class type>       int	                    CreateBuffer(std::vector<type>& data, cl_mem_flags flag);
	template<class type>       int						CreateBuffer(int size, cl_mem_flags flag);
	template<class type>	   int				        WriteBuffer(int bufferId, std::vector<type>& data);

	bool												SetKernelArg(std::string name, cl_uint argId, size_t size, const void *argPtr);
	void												SetKernelGlobalWorkRange(std::string name, int workRange);
	int													GetKernelId(std::string name);
	void												Close();
	
};


template<class systemType> bool	OpenClManager::SetKernelArg(std::string name, cl_uint argId, const systemType &argPtr)
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
			cl_int err = kernels[kernelId].setArg<systemType>(argId, argPtr);
			if (err != CL_SUCCESS)
			{
				//-38
				checkErr(err);
				std::cout << "\t" << "unable to set kernel argument" << std::endl;
				std::cout << "\t" << "Kernel id: " << kernelsNames[kernelId] << "\t" << "Argument id: " << argId << std::endl;
				return false;
			}
			else
			{
				return true;
			}
		}
		else
		{
			std::cout << "\t" << "unable to set kernel argument - system not initialized" << std::endl;
			return false;
		}
	}
	else
	{
		return false;
	}
};


template<class type> int OpenClManager::WriteBuffer(int bufferId, std::vector<type>& data)
{
	if (data.size() > 0)
	{
		unsigned structureSize = sizeof(type);
		int      bufferLength = static_cast<int>(data.size());
		//int      adjustedSize = static_cast<int>(localWorkGroupSize * std::ceil((float)bufferElements / localWorkGroupSize));

		//data.resize(adjustedSize);

		//if (clBufferSize[bufferId] != (bufferLength*structureSize)) return false;
		//if (clBufferElements[bufferId] != bufferElements) return false;
		if (clBufferLength[bufferId] != bufferLength ) return false;
		if (clBufferStructureSize[bufferId] != structureSize) return false;

		auto err = queue.enqueueWriteBuffer(clBuffer[bufferId], CL_TRUE, 0, structureSize * bufferLength, &data.front());
		if (err == CL_SUCCESS)
		{			
			return static_cast<int>(1);
		}
		else
		{
			std::cout << "OpenClManager:	error while loading buffer data" << std::endl;
			checkErr(err);
			return -1;
		}
	}
	else
	{
		std::cout << "OpenClManager:	error while loading buffer data - buffer is empty" << std::endl;
		return -2;
	}
};

template<class type> int OpenClManager::CreateBuffer(std::vector<type>& data, cl_mem_flags flag)
{
	if (data.size() > 0)
	{
		unsigned structureSize = sizeof(type);
		int      bufferElements = static_cast<int>(data.size());
		int      adjustedSize = static_cast<int>(localWorkGroupSize * std::ceil((float)bufferElements / localWorkGroupSize));

		data.resize(adjustedSize);

		cl::Buffer buffer(context, flag, structureSize *  adjustedSize);
		auto err = queue.enqueueWriteBuffer(buffer, CL_TRUE, 0, structureSize * adjustedSize, &data.front());
		if (err == CL_SUCCESS)
		{
			clBuffer.push_back(buffer);
			clBufferSize.push_back(bufferElements*structureSize);
			clBufferElements.push_back(bufferElements);
			clBufferLength.push_back(adjustedSize);
			clBufferStructureSize.push_back(structureSize);
			return static_cast<int>(clBuffer.size() - 1);
		}
		else
		{
			std::cout << "OpenClManager:	error while loading buffer data" << std::endl;
			checkErr(err);
			return -1;
		}
	}
	else
	{
		std::cout << "OpenClManager:	error while loading buffer data - buffer is empty" << std::endl;
		return -2;
	}
};


template<class type> int OpenClManager::CreateBuffer(int size, cl_mem_flags flag)
{
	if (size > 0)
	{
		unsigned			structureSize  = sizeof(type);
		int					bufferElements = size;
		int					adjustedSize   = static_cast<int>(localWorkGroupSize * std::ceil((float)bufferElements / localWorkGroupSize));
		std::vector<type>	empty_buffer;
				
		empty_buffer.resize(adjustedSize);

		cl::Buffer buffer(context, flag, structureSize *  adjustedSize);
		auto err = queue.enqueueWriteBuffer(buffer, CL_TRUE, 0, structureSize * adjustedSize, &empty_buffer.front());
		if (err == CL_SUCCESS)
		{
			clBuffer.push_back(buffer);
			clBufferSize.push_back(bufferElements*structureSize);
			clBufferElements.push_back(bufferElements);
			clBufferLength.push_back(adjustedSize);
			clBufferStructureSize.push_back(structureSize);
			return static_cast<int>(clBuffer.size() - 1);
		}
		else
		{
			std::cout << "OpenClManager:	error while loading buffer data - buffer is empty" << std::endl;
			return -1;
		}
	}
	else
	{
		std::cout << "OpenClManager:	error while loading buffer data - buffer is empty" << std::endl;
		return -1;
	}
};