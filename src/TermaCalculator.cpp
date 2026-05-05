#include "dose_engine/dose/TermaCalculator.h"

#include <cmath>
#include <algorithm>
#include <ctime>
#include <limits>
#include <stdexcept>
#include <sys/types.h>

namespace doseengine::dose
{

	core::Grid3D<float> TermaCalculator::computeWaterTerma(
			const core::Volume& volume,
			const headmodel::grid::Grid2D<float>& fluence,
			const geometry::BeamGeometry& beam,
			const Params& params)
	{
		if(params.mu_per_mm < 0.0 || params.mu_en_per_mm < 0.0)
			throw std::runtime_error("TermaCalculator: attenuation coefficients must be non-negative");

		const auto& density = volume.density();

		core::Grid3D<float> terma(
				density.nx(), density.ny(), density.nz(),
				density.dx_mm(), density.dy_mm(), density.dz_mm(),
				density.x0_mm(), density.y0_mm(), density.z0_mm()
				);

		terma.fill(0.0f);

		const double sourceToIso_mm = beam.sad_mm;


		for(std::size_t k = 0; k < density.nz(); ++k)
		{
			const double z_mm = density.z(k);
			const double attenuation = std::exp(-params.mu_per_mm * z_mm);

			double invSq = 1.0;
			if(params.useInverseSquare)
			{
				const double sourceToPoint_mm = sourceToIso_mm + z_mm;
				invSq = (sourceToIso_mm * sourceToIso_mm)/
							(sourceToPoint_mm * sourceToPoint_mm);
			}

			for(std::size_t j = 0; j < density.ny(); ++j)
			{
				const double y_mm = density.y(j);

				for(std::size_t i = 0; i < density.nx(); ++i)
				{
					const double x_mm = density.x(i);

					const float phi = sampleFluenceNearest(fluence, x_mm, y_mm);

					const double rho  = density(i, j, k);

					const double value = 
						static_cast<double>(phi)
						* attenuation
						* invSq
						* params.mu_en_per_mm
						* rho;

					terma(i, j, k) = static_cast<float>(value);
				}
			}
		}

		return terma;
	}

	// TODO: double x0() const { return m_x0; }
	//		 double y0() const { return m_y0; }
	//		 add these functions to Grid2D.h in  headmodel


	float TermaCalculator::sampleFluenceNearest(
		const headmodel::grid::Grid2D<float>& fluence,
		double x_mm,
		double y_mm
	)
	{
		int bestI = -1;
		int bestJ = -1;

		double bestDx = 1.0e30;
		double bestDy = 1.0e30;

		for (int i = 0; i < fluence.nx(); ++i)
		{
			const double d = std::abs(x_mm - fluence.xCenter(i));
			if (d < bestDx)
			{
				bestDx = d;
				bestI = i;
			}
		}

		for (int j = 0; j < fluence.ny(); ++j)
		{
			const double d = std::abs(y_mm - fluence.yCenter(j));
			if (d < bestDy)
			{
				bestDy = d;
				bestJ = j;
			}
		}

		if (bestI < 0 || bestJ < 0)
			return 0.0f;

		if (bestDx > 0.5 * fluence.dx() || bestDy > 0.5 * fluence.dy())
			return 0.0f;

		return fluence(bestI, bestJ);
	}
}
