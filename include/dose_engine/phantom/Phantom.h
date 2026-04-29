#pragma once

#include "dose_engine/core/Volume.h"

namespace doseengine::phantom
{
	class Phantom
	{
	public:
		virtual ~Phantom() = default;

		virtual const core::Volume& volume() const = 0;
		virtual core::Volume& volume() = 0;

	};
}

