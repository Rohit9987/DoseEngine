#include "dose_engine/kernel/SpectrumKernelBuilder.h"
#include "dose_engine/kernel/Kernel3D.h"
#include "dose_engine/physics/BeamSpectrum.h"
#include "dose_engine/physics/PhotonAttenuationTable.h"

#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <iostream>

namespace doseengine::kernel {

	namespace {
		const physics::SpectrumPoint* findSpectrumPoint(
				const physics::BeamSpectrum& spectrum,
				double energy_MeV)
		{
			constexpr double tol = 1.0e-6;

			for(const auto& p: spectrum.points())
			{
				if(std::abs(p.energy_MeV - energy_MeV) < tol)
					return &p;
			}

			return nullptr;
		}


		void assetSameGeometry(
				const Kernel3D& a,
				const Kernel3D& b)
		{
			if (a.nx() != b.nx() || a.ny() != b.ny() || a.nz() != b.nz())
				throw std::runtime_error("SpectrumKernelBuilder: kernel shape mismatch");

			if (std::abs(a.dx_mm() - b.dx_mm()) > 1e-9 ||
				std::abs(a.dy_mm() - b.dy_mm()) > 1e-9 ||
				std::abs(a.dz_mm() - b.dz_mm()) > 1e-9)
				throw std::runtime_error("SpectrumKernelBuilder: kernel spacing mismatch");

			if (std::abs(a.x0_mm() - b.x0_mm()) > 1e-9 ||
				std::abs(a.y0_mm() - b.y0_mm()) > 1e-9 ||
				std::abs(a.z0_mm() - b.z0_mm()) > 1e-9)
				throw std::runtime_error("SpectrumKernelBuilder: kernel offset mismatch");
		}
	}	// anon namespace

	Kernel3D SpectrumKernelBuilder::buildWeightedKernel(
			const std::vector<SpectrumKernelFile>& kernelFiles,
			const physics::BeamSpectrum& spectrum,
			const physics::PhotonAttenuationTable& attenuation,
			double referenceDepth_Mm,
			bool normalizeEachKernel)
	{
		if(kernelFiles.empty())
			throw std::runtime_error("SpectrumKernelBuilder: no kernel files provided");

		std::vector<double> weights;
		weights.reserve(kernelFiles.size());

		double weightSum = 0.0;


		for(const auto& f: kernelFiles)
		{
			const auto* sp = findSpectrumPoint(spectrum, f.energy_MeV);

			if(!sp)
			{
				std::cout << "skipping kernel energy not in spectrum: "
						  << f.energy_MeV << " MeV\n";
				weights.push_back(0.0);
				continue;
			}

			const double mu = attenuation.muPerMm(f.energy_MeV, 1.0);
			const double muen = attenuation.muenPerMm(f.energy_MeV, 1.0);

			const double w = sp->relativeWeight
							* f.energy_MeV
							* std::exp(-mu * referenceDepth_Mm)
							* muen;

			weights.push_back(w);
			weightSum += w;
		}

		if(!weightSum > 0.0)
			throw std::runtime_error("SpectrumKernelBuilder: total spectral kernel weight is zero");

		for(auto& w: weights)
			w/= weightSum;

		bool initialized = false;

		Kernel3D weightedKernel;

		for(std::size_t n = 0; n < kernelFiles.size(); ++n)
		{
			const double w = weights[n];

			if(w == 0.0)
				continue;

			auto kernel = KernelLoader::loadKernelBin(
					kernelFiles[n].path,
					normalizeEachKernel);


			if(!initialized)
			{
				weightedKernel = kernel;
				weightedKernel.values().fill(0.0f);
				initialized = true;
			}
			else
				assetSameGeometry(weightedKernel, kernel);

			const auto& src = kernel.values().data();
			auto& dst = weightedKernel.values().data();

			for(std::size_t i = 0; i < dst.size(); ++i)
			{
				dst[i] += static_cast<float>(w * static_cast<double>(src[i]));
			}

			std::cout << "Added kernel " << kernelFiles[n].energy_MeV << " MeV with weight " << w << "\n";
		}

		if(!initialized)
			throw std::runtime_error("SpectrumKernelBuilder: no matching kernels were loaded.");

		weightedKernel.normalizeToUnitSum();


		std::cout << "Spectrum-weighted kernel built\n";
		std::cout << "  Sum = " << weightedKernel.sum() << "\n";
		std::cout << "  Reference depth = "
				  << referenceDepth_Mm << " mm\n";

		return weightedKernel;
	}
}
