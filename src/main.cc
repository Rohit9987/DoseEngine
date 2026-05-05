#include <iostream>
#include <chrono>

#include "headmodel/HeadModel.h"
#include "headmodel/HeadModelConfig.h"

#include "dose_engine/phantom/WaterPhantom.h"
#include "dose_engine/geometry/BeamGeometry.h"
#include "dose_engine/dose/TermaCalculator.h"


#include "dose_engine/physics/BeamSpectrum.h"
#include "dose_engine/physics/PhotonAttenuationTable.h"
#include "dose_engine/physics/EffectiveAttenuation.h"

void test_photonAttenuation()
{
auto spectrum = doseengine::physics::BeamSpectrum::sixMV();
auto waterAtt = doseengine::physics::PhotonAttenuationTable::water();

double sumWeights = 0.0;
double meanEnergy = 0.0;

for (const auto& p : spectrum.points())
{
    sumWeights += p.relativeWeight;
    meanEnergy += p.relativeWeight * p.energy_MeV;
}

double mu_eff =
    doseengine::physics::EffectiveAttenuation::weightedMuPerMm(
        spectrum, waterAtt
    );

double muen_eff =
    doseengine::physics::EffectiveAttenuation::weightedMuenPerMm(
        spectrum, waterAtt
    );

std::cout << "Spectrum sanity check\n";
std::cout << "  Sum weights = " << sumWeights << "\n";
std::cout << "  Mean energy = " << meanEnergy << " MeV\n";
std::cout << "  Effective mu = " << mu_eff << " mm^-1\n";
std::cout << "  Effective mu_en = " << muen_eff << " mm^-1\n";

}

int main()
{
	auto start = std::chrono::high_resolution_clock::now();

    headmodel::HeadModelConfig headConfig;
    headConfig.nx = 201;
    headConfig.ny = 201;
    headConfig.dx_mm = 1.0;
    headConfig.dy_mm = 1.0;

    headmodel::HeadModel head(headConfig);

    auto fluence = head.computeOpenField(100.0); // 10 x 10 cm
												 
	//time
	auto end = std::chrono::high_resolution_clock::now();

	auto duration_ms =
		std::chrono::duration_cast<std::chrono::seconds>(end - start);

	std::cout << "HEAD MODEL calculation time: "
			  << duration_ms.count()
			  << " s\n";

    doseengine::phantom::WaterPhantom phantom(
        150, 150, 150,
        2.0
    );

    doseengine::geometry::BeamGeometry beam;
    beam.sad_mm = 1000.0;
    beam.fieldX_mm = 100.0;
    beam.fieldY_mm = 100.0;

    doseengine::dose::TermaCalculator::Params termaParams;
    termaParams.mu_per_mm = 0.005;
    termaParams.mu_en_per_mm = 0.005;
    termaParams.useInverseSquare = false;

	start = std::chrono::high_resolution_clock::now();

	auto terma = doseengine::dose::TermaCalculator::computeWaterTerma(
		phantom.volume(),
		fluence.total,
		beam,
		termaParams
	);

	end = std::chrono::high_resolution_clock::now();

	duration_ms =
		std::chrono::duration_cast<std::chrono::seconds>(end - start);

	std::cout << "TERMA calculation time: "
			  << duration_ms.count()
			  << " s\n";

	test_photonAttenuation();
    return 0;
}
