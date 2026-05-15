#include "dose_engine/dose/StencilConvolver.h"

#include <cmath>
#include <stdexcept>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace doseengine::dose
{

core::Grid3D<float> StencilConvolver::convolveGather(
    const core::Grid3D<float>& terma,
    const kernel::KernelStencil& stencil
)
{
    if (std::abs(stencil.spacingMm() - terma.dx_mm()) > 1e-9 ||
        std::abs(stencil.spacingMm() - terma.dy_mm()) > 1e-9 ||
        std::abs(stencil.spacingMm() - terma.dz_mm()) > 1e-9)
    {
        throw std::runtime_error("StencilConvolver: stencil spacing must match TERMA grid spacing");
    }

    core::Grid3D<float> dose(
        terma.nx(), terma.ny(), terma.nz(),
        terma.dx_mm(), terma.dy_mm(), terma.dz_mm(),
        terma.x0_mm(), terma.y0_mm(), terma.z0_mm()
    );

    dose.fill(0.0f);

    const auto& taps = stencil.taps();


	// Convolution
	#pragma omp parallel for collapse(3) schedule(dynamic)
	for (long tz = 0; tz < static_cast<long>(dose.nz()); ++tz)
	{
		for (long ty = 0; ty < static_cast<long>(dose.ny()); ++ty)
		{
			for (long tx = 0; tx < static_cast<long>(dose.nx()); ++tx)
			{
				double sum = 0.0;

				for (const auto& tap : taps)
				{
					const long sx = tx - tap.dx;
					const long sy = ty - tap.dy;
					const long sz = tz - tap.dz;

					if (sx < 0 || sy < 0 || sz < 0 ||
						sx >= static_cast<long>(terma.nx()) ||
						sy >= static_cast<long>(terma.ny()) ||
						sz >= static_cast<long>(terma.nz()))
					{
						continue;
					}

					sum += static_cast<double>(
						terma(
							static_cast<std::size_t>(sx),
							static_cast<std::size_t>(sy),
							static_cast<std::size_t>(sz)
						)
					) * static_cast<double>(tap.weight);
				}

				dose(
					static_cast<std::size_t>(tx),
					static_cast<std::size_t>(ty),
					static_cast<std::size_t>(tz)
				) = static_cast<float>(sum);
			}
		}
	}
    return dose;
}

}
