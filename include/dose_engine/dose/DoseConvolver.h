#pragma once

#include "dose_engine/core/Grid3D.h"
#include "dose_engine/kernel/Kernel3D.h"
#include <cstddef>


namespace doseengine::dose {

	class DoseConvolver
	{
	public:
		static core::Grid3D<float> convolveNaive(
				const core::Grid3D<float>& terma,
				const kernel::Kernel3D& kernel
				);

	private:
		static std::size_t originIndex(double offset_mm, double spacing_mm);
	};
}
