#pragma once

#include "dose_engine/physics/BeamSpectrum.h"
#include "dose_engine/physics/PhotonAttenuationTable.h"
#include <memory>

namespace doseengine::physics {

	class EffectiveAttenuation
	{
	public:
		static double weightedMuPerMm(
				const BeamSpectrum& spectrum,
				const PhotonAttenuationTable& material)
		{
			double sum = 0.0;

			for(const auto& p: spectrum.points())
			{
				sum += p.relativeWeight *
					   material.muPerMm(p.energy_MeV);
			}

			return sum;
		}

		static double weightedMuenPerMm(
				const BeamSpectrum& spectrum,
				const PhotonAttenuationTable& material)
		{
			double sum = 0.0;

			for(const auto& p: spectrum.points())
			{
				sum += p.relativeWeight *
					   material.muenPerMm(p.energy_MeV);
			}

			return sum;
		}
	};
}
