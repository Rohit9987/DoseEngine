#pragma once

#include "dose_engine/phantom/Phantom.h"

namespace doseengine::phantom {

	class WaterPhantom: public Phantom
	{
	public:
		WaterPhantom(
				std::size_t nx,
				std::size_t ny,
				std::size_t nz,
				double spacing_mm
				)
				{
				core::Grid3D<float> density(
						nx, ny, nz,
						spacing_mm, spacing_mm, spacing_mm,
						-0.5 * static_cast<double>(nx) * spacing_mm,
						-0.5 * static_cast<double>(ny) * spacing_mm,
						0.0
						);

				density.fill(1.0f);
				volume_ = core::Volume(std::move(density));
				}

		const core::Volume& volume() const override { return volume_;}
		core::Volume& volume() override { return volume_;}

	private:
		core::Volume volume_;
	};
}
