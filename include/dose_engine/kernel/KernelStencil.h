#pragma once

#include "dose_engine/kernel/Kernel3D.h"

#include <vector>
#include <cstddef>

namespace doseengine::kernel{

	struct KernelTap
	{
		int dx = 0,
			dy = 0,
			dz = 0;
		float weight = 0.0f;
	};

	class KernelStencil
	{

	public:
		static Kernel3D binToCoarseKernel(
				const Kernel3D& fineKernel,
				double coarseSpacing_mm,
				std::size_t nx,
				std::size_t ny,
				std::size_t nz,
				double x0_mm,
				double y0_mm,
				double z0_mm);

		static KernelStencil fromKernel(
				const Kernel3D& coarseKernel,
				double keepFraction = 0.99,
				bool renormalizeKeptWeights = true);

		const std::vector<KernelTap>& taps() const {return taps_;}

		double spacingMm() const { return spacing_mm_;}

		std::size_t size() const { return taps_.size(); }

		double retainedFraction() const { return retainedFraction_; }

	private:
		std::vector<KernelTap> taps_;
		double spacing_mm_ = 1.0;
		double retainedFraction_ = 1.0;
	};
}
