#include <iostream>
#include <fstream>
#include <chrono>

#include "headmodel/HeadModel.h"
#include "headmodel/HeadModelConfig.h"

#include "dose_engine/phantom/WaterPhantom.h"

#include "dose_engine/geometry/BeamGeometry.h"

#include "dose_engine/physics/BeamSpectrum.h"
#include "dose_engine/physics/PhotonAttenuationTable.h"
#include "dose_engine/physics/EffectiveAttenuation.h"

#include "dose_engine/kernel/KernelLoader.h"

#include "dose_engine/dose/TermaCalculator.h"
#include "dose_engine/dose/DoseConvolver.h"

#include "dose_engine/kernel/KernelStencil.h"
#include "dose_engine/dose/StencilConvolver.h"
#include "dose_engine/kernel/SpectrumKernelBuilder.h"

#include "dose_engine/io/Grid3DExport.h"

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
        100, 100, 120,
        2.0
    );

    doseengine::geometry::BeamGeometry beam;
    beam.sad_mm = 1000.0;
    beam.fieldX_mm = 100.0;
    beam.fieldY_mm = 100.0;

	auto spectrum = doseengine::physics::BeamSpectrum::sixMV();
	auto waterAtt = doseengine::physics::PhotonAttenuationTable::water();

	start = std::chrono::high_resolution_clock::now();

	auto terma = doseengine::dose::TermaCalculator::computeWaterTerma(
		phantom.volume(),
		fluence.total,
		beam,
		spectrum,
		waterAtt,
		1.0,
		true
	);

	end = std::chrono::high_resolution_clock::now();
	duration_ms =
		std::chrono::duration_cast<std::chrono::seconds>(end - start);

	std::cout << "TERMA calculation time: "
			  << duration_ms.count()
			  << " s\n";
	std::cout << "TERMA size: (" 
			  << terma.nx() << ','
			  << terma.ny() << ','
			  << terma.nz() << ")\n";

	auto t10 = terma(terma.nx()/2, terma.ny()/2, 50);
	auto t20 = terma(terma.nx()/2, terma.ny()/2, 100);

	std::cout << "terma at 10 cm : " << t10 << "\n"
			  << "terma at 20 cm : " << t20 << "\n"
			  << "terma_TPR20/10 : " << t20/t10 << '\n';

	start = std::chrono::high_resolution_clock::now();

	const std::string kernelDir = "../../Kernel/kernel-downsampled/";

	auto kernelFiles =
		doseengine::physics::BeamSpectrum::sixMVKernelFiles(kernelDir);

	auto weightedKernel =
		doseengine::kernel::SpectrumKernelBuilder::buildWeightedKernel(
			kernelFiles,
			spectrum,
			waterAtt,
			100.0,  // reference depth = 10 cm
			true    // normalize each monoenergetic kernel first
		);

	auto coarseKernel = doseengine::kernel::KernelStencil::binToCoarseKernel(
		weightedKernel,
		2.0,        // match TERMA voxel spacing
		51, 51, 51,
		-50.0, -50.0, -20.0
	);

	auto stencil = doseengine::kernel::KernelStencil::fromKernel(
		coarseKernel,
		0.95,   // keep 99.5% of kernel magnitude
		true     // renormalize kept taps to unit sum
	);

	auto dose = doseengine::dose::StencilConvolver::convolveGather(
		terma,
		stencil
	);

	std::cout << "Stencil dose calculated\n";
	std::cout << "Stencil taps = " << stencil.size() << "\n";


	end = std::chrono::high_resolution_clock::now();
	duration_ms =
		std::chrono::duration_cast<std::chrono::seconds>(end - start);

	std::cout << "DOSE calculation time: "
			  << duration_ms.count()
			  << " s\n";

	auto d10 = dose(dose.nx()/2, dose.ny()/2, 50);
	auto d20 = dose(dose.nx()/2, dose.ny()/2, 100);

	std::cout << "Dose at 10 cm : " << d10 << "\n"
			  << "Dose at 20 cm : " << d20 << "\n"
			  << "TPR20/10 : " << d20/d10 << '\n';


	doseengine::io::exportCaxProfileCsv(
		"terma_cax.csv",
		terma,
		"terma"
	);

	doseengine::io::exportCrosslineProfileAtDepthCsv(
		"terma_crossline_100mm.csv",
		terma,
		"terma",
		100.0
	);

	doseengine::io::exportCaxProfileCsv(
		"dose_cax.csv",
		dose,
		"dose"
	);

	doseengine::io::exportCrosslineProfileAtDepthCsv(
		"dose_crossline_100mm.csv",
		dose,
		"dose",
		100.0
	);
	
	

    return 0;
}
