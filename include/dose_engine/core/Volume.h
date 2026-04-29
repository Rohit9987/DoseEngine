#pragma once

#include "dose_engine/core/Grid3D.h"

namespace doseengine::core {

	class Volume
	{
	public:
		Volume() = default;

		explicit Volume(Grid3D<float> density_g_per_cm3)
			: density_(std::move(density_g_per_cm3))
		{}

		const Grid3D<float>& density() const {return density_;}
		Grid3D<float>& density() { return density_;}


		std::size_t nx() const { return density_.nx(); }
		std::size_t ny() const { return density_.ny(); }
		std::size_t nz() const { return density_.nz(); }

	private:
		Grid3D<float> density_;

	};
}
