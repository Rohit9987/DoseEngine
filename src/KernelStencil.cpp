#include "dose_engine/kernel/KernelStencil.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <iostream>

namespace doseengine::kernel
{

Kernel3D KernelStencil::binToCoarseKernel(
    const Kernel3D& fineKernel,
    double coarseSpacing_mm,
    std::size_t nx,
    std::size_t ny,
    std::size_t nz,
    double x0_mm,
    double y0_mm,
    double z0_mm
)
{
    if (coarseSpacing_mm <= 0.0)
        throw std::runtime_error("KernelStencil: invalid coarse spacing");

    core::Grid3D<float> coarse(
        nx, ny, nz,
        coarseSpacing_mm, coarseSpacing_mm, coarseSpacing_mm,
        x0_mm, y0_mm, z0_mm
    );

    coarse.fill(0.0f);

    const auto& fine = fineKernel.values();

    for (std::size_t kz = 0; kz < fine.nz(); ++kz)
    {
        const double z = fine.z(kz);

        const long cz = static_cast<long>(
            std::llround((z - z0_mm) / coarseSpacing_mm)
        );

        if (cz < 0 || cz >= static_cast<long>(nz))
            continue;

        for (std::size_t ky = 0; ky < fine.ny(); ++ky)
        {
            const double y = fine.y(ky);

            const long cy = static_cast<long>(
                std::llround((y - y0_mm) / coarseSpacing_mm)
            );

            if (cy < 0 || cy >= static_cast<long>(ny))
                continue;

            for (std::size_t kx = 0; kx < fine.nx(); ++kx)
            {
                const double x = fine.x(kx);

                const long cx = static_cast<long>(
                    std::llround((x - x0_mm) / coarseSpacing_mm)
                );

                if (cx < 0 || cx >= static_cast<long>(nx))
                    continue;

                coarse(
                    static_cast<std::size_t>(cx),
                    static_cast<std::size_t>(cy),
                    static_cast<std::size_t>(cz)
                ) += fine(kx, ky, kz);
            }
        }
    }

    Kernel3D out(std::move(coarse), fineKernel.energyMeV());

    const double before = out.sum();
    if (before <= 0.0)
        throw std::runtime_error("KernelStencil: coarse kernel has zero sum");

    out.normalizeToUnitSum();

    std::cout << "Coarse kernel created\n";
    std::cout << "  Coarse shape = " << nx << " x " << ny << " x " << nz << "\n";
    std::cout << "  Sum before renormalisation = " << before << "\n";
    std::cout << "  Sum after renormalisation = " << out.sum() << "\n";

    return out;
}

KernelStencil KernelStencil::fromKernel(
    const Kernel3D& coarseKernel,
    double keepFraction,
    bool renormalizeKeptWeights
)
{
    if (keepFraction <= 0.0 || keepFraction > 1.0)
        throw std::runtime_error("KernelStencil: keepFraction must be in (0, 1]");

    if (std::abs(coarseKernel.dx_mm() - coarseKernel.dy_mm()) > 1e-9 ||
        std::abs(coarseKernel.dx_mm() - coarseKernel.dz_mm()) > 1e-9)
    {
        throw std::runtime_error("KernelStencil: only isotropic kernels supported for stencil");
    }

    const auto& K = coarseKernel.values();

    const int ox = static_cast<int>(std::llround(-coarseKernel.x0_mm() / coarseKernel.dx_mm()));
    const int oy = static_cast<int>(std::llround(-coarseKernel.y0_mm() / coarseKernel.dy_mm()));
    const int oz = static_cast<int>(std::llround(-coarseKernel.z0_mm() / coarseKernel.dz_mm()));

    struct Candidate
    {
        KernelTap tap;
        double absWeight;
    };

    std::vector<Candidate> candidates;
    candidates.reserve(K.nx() * K.ny() * K.nz());

    double totalAbs = 0.0;

    for (std::size_t kz = 0; kz < K.nz(); ++kz)
    {
        for (std::size_t ky = 0; ky < K.ny(); ++ky)
        {
            for (std::size_t kx = 0; kx < K.nx(); ++kx)
            {
                const float w = K(kx, ky, kz);

                if (w == 0.0f)
                    continue;

                KernelTap tap;
                tap.dx = static_cast<int>(kx) - ox;
                tap.dy = static_cast<int>(ky) - oy;
                tap.dz = static_cast<int>(kz) - oz;
                tap.weight = w;

                const double aw = std::abs(static_cast<double>(w));
                totalAbs += aw;

                candidates.push_back({tap, aw});
            }
        }
    }

    std::sort(
        candidates.begin(),
        candidates.end(),
        [](const Candidate& a, const Candidate& b)
        {
            return a.absWeight > b.absWeight;
        }
    );

    KernelStencil stencil;
    stencil.spacing_mm_ = coarseKernel.dx_mm();

    double runningAbs = 0.0;
    double keptSignedSum = 0.0;

    for (const auto& c : candidates)
    {
        if (runningAbs / totalAbs >= keepFraction)
            break;

        stencil.taps_.push_back(c.tap);
        runningAbs += c.absWeight;
        keptSignedSum += c.tap.weight;
    }

    stencil.retainedFraction_ = runningAbs / totalAbs;

    if (renormalizeKeptWeights)
    {
        if (std::abs(keptSignedSum) < 1e-30)
            throw std::runtime_error("KernelStencil: kept stencil sum is zero");

        for (auto& tap : stencil.taps_)
            tap.weight = static_cast<float>(tap.weight / keptSignedSum);
    }

    std::cout << "Kernel stencil created\n";
    std::cout << "  Taps kept = " << stencil.taps_.size() << "\n";
    std::cout << "  Retained fraction = " << stencil.retainedFraction_ << "\n";

    return stencil;
}

}
