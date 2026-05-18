#pragma once

#include <vector>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <string>

namespace doseengine::kernel
{

struct SpectrumKernelFile
{
    double energy_MeV = 0.0;
    std::string path;
};

}

namespace doseengine::physics
{

struct SpectrumPoint
{
    double energy_MeV;
    double relativeWeight;
};

class BeamSpectrum
{
public:
    explicit BeamSpectrum(std::vector<SpectrumPoint> points)
        : points_(std::move(points))
    {
        if (points_.empty())
            throw std::runtime_error("BeamSpectrum: empty spectrum");

        normalize();
    }

    const std::vector<SpectrumPoint>& points() const
    {
        return points_;
    }

    static BeamSpectrum sixMV()
    {
        return BeamSpectrum({
            {0.10, 0.055},
            {0.20, 0.098},
            {0.30, 0.136},
            {0.40, 0.168},
            {0.50, 0.196},
            {0.60, 0.220},
            {0.80, 0.256},
            {1.00, 0.279},
            {1.25, 0.292},
            {1.50, 0.292},
            {2.00, 0.266},
            {3.00, 0.172},
            {4.00, 0.086},
            {5.00, 0.025},
            {6.00, 0.003}
        });
    }

    static std::vector<kernel::SpectrumKernelFile> sixMVKernelFiles(
        const std::string& kernelDir
    )
    {
        return {
            {0.10, kernelDir + "0.1MeV_kernel.bin"},
            {0.20, kernelDir + "0.2MeV_kernel.bin"},
            {0.30, kernelDir + "0.3MeV_kernel.bin"},
            {0.40, kernelDir + "0.4MeV_kernel.bin"},
            {0.50, kernelDir + "0.5MeV_kernel.bin"},
            {0.60, kernelDir + "0.6MeV_kernel.bin"},
            {0.80, kernelDir + "0.8MeV_kernel.bin"},
            {1.00, kernelDir + "1.0MeV_kernel.bin"},
            {1.25, kernelDir + "1.25MeV_kernel.bin"},
            {1.50, kernelDir + "1.50MeV_kernel.bin"},
            {2.00, kernelDir + "2.0MeV_kernel.bin"},
            {3.00, kernelDir + "3.0MeV_kernel.bin"},
            {4.00, kernelDir + "4.0MeV_kernel.bin"},
            {5.00, kernelDir + "5.0MeV_kernel.bin"},
            {6.00, kernelDir + "6.0MeV_kernel.bin"}
        };
    }

private:
    void normalize()
    {
        double sum = 0.0;

        for (const auto& p : points_)
            sum += p.relativeWeight;

        if (sum <= 0.0)
            throw std::runtime_error("BeamSpectrum: invalid spectrum sum");

        for (auto& p : points_)
            p.relativeWeight /= sum;
    }

    std::vector<SpectrumPoint> points_;
};

}
