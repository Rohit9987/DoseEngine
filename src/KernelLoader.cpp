#include "dose_engine/kernel/KernelLoader.h"

#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <ios>
#include <regex>
#include <stdexcept>

namespace doseengine::kernel
{
	
	namespace
	{
		template <typename T>
		void readValue(std::ifstream& in, T& value)
		{
			in.read(reinterpret_cast<char*>(&value), sizeof(T));
			if(!in)
				throw std::runtime_error("KernelLoader: failed to read binary value");
		}


		template <typename T, std::size_t N>
		void readArray(std::ifstream& in, std::array<T, N>& values)
		{
			in.read(reinterpret_cast<char*>(values.data()), sizeof(T) * N);
			if(!in)
				throw std::runtime_error("KernelLoader: failed to read binary array");
		}

		std::uint64_t voxelCount(const std::array<std::uint64_t, 3>& shape)
		{
			return shape[0] * shape[1] * shape[2];
		}

		double parseEnergyMeVFromFilename(const std::string& path)
		{
			const std::string name = std::filesystem::path(path).filename().string();
			
			const std::regex rgx(R"(^([0-9]+(?:\.[0-9]+)?)MeV)");
			std::smatch match;

			if(!std::regex_search(name, match, rgx))
				return 0.0;	// allow unknown energy for now

			return std::stod(match[1].str());
		}

	}	// anonymous namespace

	Kernel3D KernelLoader::loadKernelBin(
			const std::string& path,
			bool normalizeToUnitSum)
	{
		std::ifstream in(path, std::ios::binary);
		if(!in)
			throw std::runtime_error("KernelLoader: count not open file: " + path);

		std::uint64_t size = 0;
		std::array<std::uint64_t, 3> shape{};
		std::array<double, 3> half{};
		std::array<double, 3> offset{};

		readValue(in, size);
		readValue(in, shape);
		readValue(in, half);
		readValue(in, offset);

		const std::uint64_t expected = voxelCount(shape);

		if(size != expected)
			throw std::runtime_error("KernelLoader: header size does no match shape product");

		const double dx_mm = 2.0 * half[0];
		const double dy_mm = 2.0 * half[1];
		const double dz_mm = 2.0 * half[2];

		core::Grid3D<float> grid(
				static_cast<std::size_t>(shape[0]),
				static_cast<std::size_t>(shape[1]),
				static_cast<std::size_t>(shape[2]),
				dx_mm, dy_mm, dz_mm,
				offset[0], offset[1], offset[2]
				);

		grid.data().resize(static_cast<std::size_t>(size));

		in.read(
				reinterpret_cast<char*>(grid.data().data()),
				static_cast<std::streamsize>(size * sizeof(float))
			   );

		if(!in)
			throw std::runtime_error("KernelLoader: failed to read kernel voxel data");

		const double energy_MeV = parseEnergyMeVFromFilename(path);
		Kernel3D  kernel(std::move(grid), energy_MeV);

		if(normalizeToUnitSum)
			kernel.normalizeToUnitSum();

		return kernel;
	}
}
