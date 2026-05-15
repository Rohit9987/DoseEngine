#include "dose_engine/dose/StencilConvolver.h"

#include <cmath>
#include <stdexcept>
#include <iostream>

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

    for (std::size_t tz = 0; tz < dose.nz(); ++tz)
    {
		std::cout << tz << "/" << dose.nz() << '\n';
        for (std::size_t ty = 0; ty < dose.ny(); ++ty)
        {
            for (std::size_t tx = 0; tx < dose.nx(); ++tx)
            {
                double sum = 0.0;

                for (const auto& tap : taps)
                {
                    const long sx = static_cast<long>(tx) - tap.dx;
                    const long sy = static_cast<long>(ty) - tap.dy;
                    const long sz = static_cast<long>(tz) - tap.dz;

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

                dose(tx, ty, tz) = static_cast<float>(sum);
            }
        }
    }

    return dose;
}

}
