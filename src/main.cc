#include <iostream>

#include "dose_engine/phantom/WaterPhantom.h"
#include "dose_engine/geometry/BeamGeometry.h"

int main()
{
	doseengine::phantom::WaterPhantom phantom(
			150, 150, 150,
			2.0);

	doseengine::geometry::BeamGeometry beam;
	beam.fieldX_mm = 100.0;
	beam.fieldY_mm = 100.0;

	std::cout << "Water phantom created\n";
    std::cout << "Grid: "
              << phantom.volume().nx() << " x "
              << phantom.volume().ny() << " x "
              << phantom.volume().nz() << "\n";

    return 0;
}

