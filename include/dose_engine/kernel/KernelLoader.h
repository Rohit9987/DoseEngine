#pragma once

#include "dose_engine/kernel/Kernel3D.h"

#include <string>

namespace doseengine::kernel {
	class KernelLoader
	{
	public:
		static Kernel3D loadKernelBin(
				const std::string& path,
				bool normalizeToUnitSum = true
				);
	};
}
