#pragma once

#include <bits/types/cookie_io_functions_t.h>
#include <vector>
#include <stdexcept>
#include <cstddef>

namespace doseengine::core
{

	template<typename T>
	class Grid3D
	{
	public:
		Grid3D() = default;

		Grid3D(
				std::size_t nx,
				std::size_t ny,
				std::size_t nz,
				double dx_mm,
				double dy_mm,
				double dz_mm,
				double x0_mm = 0.0,
				double y0_mm = 0.0,
				double z0_mm = 0.0
			  )
			    : nx_(nx), ny_(ny), nz_(nz),
          dx_mm_(dx_mm), dy_mm_(dy_mm), dz_mm_(dz_mm),
          x0_mm_(x0_mm), y0_mm_(y0_mm), z0_mm_(z0_mm),
          data_(nx * ny * nz)
		{
				if (nx == 0 || ny == 0 || nz == 0)
					throw std::runtime_error("Grid3D: dimensions must be non-zero");

				if (dx_mm <= 0.0 || dy_mm <= 0.0 || dz_mm <= 0.0)
					throw std::runtime_error("Grid3D: spacing must be positive");
		}

		T& operator()(std::size_t i, std::size_t j, std::size_t k)
		{
			return data_.at(index(i, j, k));
		}

		const T& operator()(std::size_t i, std::size_t j, std::size_t k) const
		{
			return data_.at(index(i, j, k));
		}
		
		std::size_t index(std::size_t i, std::size_t j, std::size_t k) const
		{
			return k * nx_ * ny +
				   j * nx_ +
				   i;
		}

		double x(std::size_t i) const { return x0_mm + static_cast<double>(i) * dx_mm_; }
		double y(std::size_t j) const { return y0_mm + static_cast<double>(j) * dy_mm_; }
		double z(std::size_t k) const { return z0_mm + static_cast<double>(k) * dz_mm_; }

		void fill(const T& value)
		{
			std::fill(data_.begin(), data_.end(), value);
		}

		std::size_t nx() const { return nx_; }
		std::size_t ny() const { return ny_; }
		std::size_t nz() const { return nz_; }

		double dx_mm() const { return dx_mm_; }
		double dy_mm() const { return dy_mm_; }
		double dz_mm() const { return dz_mm_; }

		double x0_mm() const { return x0_mm_; }
		double y0_mm() const { return y0_mm_; }
		double z0_mm() const { return z0_mm_; }


		std::vector<T>& data() { return data_; }
		const std::vector<T>& data() const { return data_; }

	private:
		std::size_t nx_ = 0;
		std::size_t ny_ = 0;
		std::size_t nz_ = 0;

		double dx_mm_ = 1.0;
		double dy_mm_ = 1.0;
		double dz_mm_ = 1.0;

		double x0_mm_ = 0.0;
		double y0_mm_ = 0.0;
		double z0_mm_ = 0.0;

		std::vector<T> data_;
	};

}
