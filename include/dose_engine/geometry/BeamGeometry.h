#pragma one

#include <cmath>

namespace doseengine::geometry
{

struct Vec3
{
	double x_mm = 0.0;
	double y_mm = 0.0;
	double z_mm = 0.0;
};

struct BeamGeometry
{
	Vec3 sourcePosition_mm {0.0, 0.0, -1000.0};

	double sad_mm = 1000.0;
	double ssd_mm = 1000.0;

	double gantry_deg = 0.0;
	double collimator_deg = 0.0;
	double couch_deg = 0.0;

	double fieldX_mm = 100.0;
	double fieldY_mm = 100.0;

	Vec3 isocentre_mm {0.0, 0.0, 0.0};
};
}



