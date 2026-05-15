#pragma once

#include "dose_engine/core/Grid3D.h"
#include "dose_engine/kernel/KernelStencil.h"

namespace doseengine::dose
{

class StencilConvolver
{
public:
    static core::Grid3D<float> convolveGather(
        const core::Grid3D<float>& terma,
        const kernel::KernelStencil& stencil
    );
};

}
