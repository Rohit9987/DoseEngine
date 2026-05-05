#pragma once

#include "dose_engine/core/Grid3D.h"
#include "dose_engine/core/Volume.h"
#include "dose_engine/geometry/BeamGeometry.h"

#include "headmodel/grid/Grid2D.h"

namespace doseengine::dose
{

	class TermaCalculator
	{
	public:
		struct Params
		{
			double mu_per_mm = 0.005;	// simple first-pass attenuation coefficient
			double mu_en_per_mm = 0.005;
			bool useInverseSquare = false;
		};

		static core::Grid3D<float> computeWaterTerma(
				const core::Volume& volume,
				const headmodel::grid::Grid2D<float>& fluence,
				const geometry::BeamGeometry& beam,
				const Params& Params);

	private:
		static float sampleFluenceNearest(
				const headmodel::grid::Grid2D<float>& fluence,
				double x_mm,
				double y_mm);
	};
}
