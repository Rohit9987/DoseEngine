#pragma once

#include "dose_engine/core/Grid3D.h"
#include "dose_engine/core/Volume.h"
#include "dose_engine/geometry/BeamGeometry.h"

#include "dose_engine/physics/BeamSpectrum.h"
#include "dose_engine/physics/PhotonAttenuationTable.h"

#include "headmodel/grid/Grid2D.h"

namespace doseengine::dose
{

	class TermaCalculator
	{
	public:
		struct Params		// TODO: can be deleted now, waiting for confirm
		{
			double mu_per_mm = 0.005;	// simple first-pass attenuation coefficient
			double mu_en_per_mm = 0.005;
			bool useInverseSquare = false;
		};

		static core::Grid3D<float> computeWaterTerma(
				const core::Volume& volume,
				const headmodel::grid::Grid2D<float>& fluence,
				const geometry::BeamGeometry& beam,
				const physics::BeamSpectrum& spectrum,
				const physics::PhotonAttenuationTable& attenuation,
				double density_g_per_cm3 = 1.0,
				bool useInverseSquare = false);

	private:
		static float sampleFluenceNearest(
				const headmodel::grid::Grid2D<float>& fluence,
				double x_mm,
				double y_mm);
	};
}
