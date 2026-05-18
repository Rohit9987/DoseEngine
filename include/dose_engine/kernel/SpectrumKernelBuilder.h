#pragma once

#include "dose_engine/kernel/Kernel3D.h"
#include "dose_engine/kernel/KernelLoader.h"
#include "dose_engine/physics/BeamSpectrum.h"
#include "dose_engine/physics/PhotonAttenuationTable.h"


#include <string>
#include <vector>

namespace doseengine::kernel {

	class SpectrumKernelBuilder
	{
	public:
		static Kernel3D buildWeightedKernel(
				const std::vector<SpectrumKernelFile>& kernelFiles,
				const physics::BeamSpectrum& spectrum,
				const physics::PhotonAttenuationTable& attenuation,
				double referenceDepth_mm,
				bool normalizeEachKernel = true
				);
	};
}
