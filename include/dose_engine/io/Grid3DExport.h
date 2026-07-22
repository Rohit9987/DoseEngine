#pragma once

#include "dose_engine/core/Grid3D.h"

#include <fstream>
#include <string>
#include <stdexcept>
#include <cmath>

namespace doseengine::io
{

inline void exportCaxProfileCsv(
    const std::string& filename,
    const core::Grid3D<float>& grid,
    const std::string& valueName
)
{
    std::ofstream out(filename);

    if (!out)
        throw std::runtime_error("exportCaxProfileCsv: could not open file: " + filename);

    out << "k,z_mm," << valueName << "\n";

    const std::size_t cx = grid.nx() / 2;
    const std::size_t cy = grid.ny() / 2;

    for (std::size_t k = 0; k < grid.nz(); ++k)
    {
        out << k << ","
            << grid.z(k) << ","
            << grid(cx, cy, k) << "\n";
    }
}

inline void exportCrosslineProfileCsv(
    const std::string& filename,
    const core::Grid3D<float>& grid,
    const std::string& valueName,
    std::size_t zIndex
)
{
    if (zIndex >= grid.nz())
        throw std::runtime_error("exportCrosslineProfileCsv: zIndex out of range");

    std::ofstream out(filename);

    if (!out)
        throw std::runtime_error("exportCrosslineProfileCsv: could not open file: " + filename);

    out << "i,x_mm," << valueName << "\n";

    const std::size_t cy = grid.ny() / 2;

    for (std::size_t i = 0; i < grid.nx(); ++i)
    {
        out << i << ","
            << grid.x(i) << ","
            << grid(i, cy, zIndex) << "\n";
    }
}

inline void exportCrosslineProfileAtDepthCsv(
    const std::string& filename,
    const core::Grid3D<float>& grid,
    const std::string& valueName,
    double depth_mm
)
{
    const long zIndex = static_cast<long>(
        std::llround((depth_mm - grid.z0_mm()) / grid.dz_mm())
    );

    if (zIndex < 0 || zIndex >= static_cast<long>(grid.nz()))
        throw std::runtime_error("exportCrosslineProfileAtDepthCsv: depth outside grid");

    exportCrosslineProfileCsv(
        filename,
        grid,
        valueName,
        static_cast<std::size_t>(zIndex)
    );
}

}
