#pragma once

#include "dose_engine/core/Grid3D.h"

#include <stdexcept>
#include <limits>
#include <cmath>

namespace doseengine::kernel {

	class Kernel3D
	{
	public:
		Kernel3D() = default;

		Kernel3D(
				core::Grid3D<float> values,
				double energy_MeV
				)
			: values_(std::move(values)),
			  energy_MeV_(energy_MeV)
		{}

		const core::Grid3D<float>& values() const { return values_;}
		core::Grid3D<float>& values() { return values_;}

		double energyMeV() const { return energy_MeV_;}

		double sum() const
		{
			double s = 0.0;
			for(float v: values_.data())
				s += static_cast<double>(v);
			return s;
		}

		void normalizeToUnitSum()
		{
			const double s = sum();

			if (!(s > 0.0))
				throw std::runtime_error("Kernel3D: cannot normalize kernel with non-positive sum");

			for (auto& v: values_.data())
				v = static_cast<float>(static_cast<double>(v)/s);
		}

		std::size_t nx() const { return values_.nx(); }
		std::size_t ny() const { return values_.ny(); }
		std::size_t nz() const { return values_.nz(); }

		double dx_mm() const { return values_.dx_mm(); }
		double dy_mm() const { return values_.dy_mm(); }
		double dz_mm() const { return values_.dz_mm(); }

		double x0_mm() const { return values_.x0_mm(); }
		double y0_mm() const { return values_.y0_mm(); }
		double z0_mm() const { return values_.z0_mm(); }

	private:
		core::Grid3D<float> values_;
		double energy_MeV_ = 0.0;
	};
}
