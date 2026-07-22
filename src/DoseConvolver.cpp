#include "dose_engine/dose/DoseConvolver.h"
#include "dose_engine/core/Grid3D.h"

#include <cmath>
#include <stdexcept>


namespace doseengine::dose {
	std::size_t DoseConvolver::originIndex(double offset_mm, double spacing_mm)
	{
		if(spacing_mm <= 0.0)
			throw std::runtime_error("DoseConvolver: invalid kernel spacing");

		const double idx = -offset_mm / spacing_mm;

		if (idx < 0.0)
			throw std::runtime_error("DoseConvolver: kernel origin lies outside grid");

		return static_cast<std::size_t>(std::llround(idx));
	}

	core::Grid3D<float> DoseConvolver::convolveNaive(
			const core::Grid3D<float>& terma,
			const kernel::Kernel3D& kernel)
	{
		core::Grid3D<float> dose(
			terma.nx(), terma.ny(), terma.nz(),
			terma.dx_mm(), terma.dy_mm(), terma.dz_mm(),
			terma.x0_mm(), terma.y0_mm(), terma.z0_mm()
		);

		dose.fill(0.0f);

		const auto& K = kernel.values();

		const std::size_t kx0 = originIndex(kernel.x0_mm(), kernel.dx_mm());
		const std::size_t ky0 = originIndex(kernel.y0_mm(), kernel.dy_mm());
		const std::size_t kz0 = originIndex(kernel.z0_mm(), kernel.dz_mm());



		for(std::size_t sz = 0; sz < terma.nz(); ++sz)
		{
			for(std::size_t sy = 0; sy < terma.ny(); ++sy)
			{
				for(std::size_t sx = 0; sx < terma.nx(); ++sx)
				{
					const float t = terma(sx, sy, sz);

					if(t == 0.0f)
						continue;

					for(std::size_t kz = 0; kz < K.nz(); ++kz)
					{
						const long dz = static_cast<long>(kz) - static_cast<long>(kz0);
						const long tz = static_cast<long>(sz) + dz;
						if(tz < 0 || tz >= static_cast<long>(dose.nz()))
							continue;

						for(std::size_t ky = 0; ky < K.ny(); ++ky)
						{
							const long dy = static_cast<long>(ky) - static_cast<long>(ky0);
							const long ty = static_cast<long>(sy) + dy;
							if(ty < 0 || ty >= static_cast<long>(dose.ny()))
								continue;

							for(std::size_t kx = 0; kx < K.nx(); ++kx)
							{
								const long dx = static_cast<long>(kx) - static_cast<long>(kx0);
								const long tx = static_cast<long>(sx) + dx;
								if(tx < 0 || tx >= static_cast<long>(dose.nx()))
									continue;

								if(K(kx, ky, kz) < 1e-8f)
									continue;

								dose(
									static_cast<std::size_t>(tx),
									static_cast<std::size_t>(ty),
									static_cast<std::size_t>(tz)
									)
									+= t * K(kx, ky, kz);
							}
						}
					}
				}
			}
		}

		return dose;
	}
}

