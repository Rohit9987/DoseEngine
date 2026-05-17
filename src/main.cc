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


	start = std::chrono::high_resolution_clock::now();

	std::vector<doseengine::kernel::SpectrumKernelFile> kernelFiles = {
		{0.10, "../../Kernel/kernel-downsampled/0.1MeV_kernel.bin"},
		{0.20, "../../Kernel/kernel-downsampled/0.2MeV_kernel.bin"},
		{0.30, "../../Kernel/kernel-downsampled/0.3MeV_kernel.bin"},
		{0.40, "../../Kernel/kernel-downsampled/0.4MeV_kernel.bin"},
		{0.50, "../../Kernel/kernel-downsampled/0.5MeV_kernel.bin"},
		{0.60, "../../Kernel/kernel-downsampled/0.6MeV_kernel.bin"},
		{0.80, "../../Kernel/kernel-downsampled/0.8MeV_kernel.bin"},
		{1.00, "../../Kernel/kernel-downsampled/1.0MeV_kernel.bin"},
		{1.25, "../../Kernel/kernel-downsampled/1.25MeV_kernel.bin"},
		{1.50, "../../Kernel/kernel-downsampled/1.50MeV_kernel.bin"},
		{2.00, "../../Kernel/kernel-downsampled/2.0MeV_kernel.bin"},
		{3.00, "../../Kernel/kernel-downsampled/3.0MeV_kernel.bin"},
		{4.00, "../../Kernel/kernel-downsampled/4.0MeV_kernel.bin"},
		{5.00, "../../Kernel/kernel-downsampled/5.0MeV_kernel.bin"},
		{6.00, "../../Kernel/kernel-downsampled/6.0MeV_kernel.bin"}
	};

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
			  << "Dose at 10 cm : " << d20 << "\n"
			  << "TPR20/10 : " << d20/d10 << '\n';



	
	std::ofstream out("dose_cax.csv");
	out << "k,z_mm,dose\n";

	const std::size_t cx = dose.nx() / 2;
	const std::size_t cy = dose.ny() / 2;
	for (std::size_t k = 0; k < dose.nz(); ++k)
	{
		out << k << "," << dose.z(k) << ","
			<< dose(cx, cy, k) << "\n";
	}
	out.close();

	
	std::ofstream out2("dose.csv");
	out << "k,x_mm,terma\n";

	const std::size_t cz = 50;
	for (std::size_t k = 0; k < dose.nx(); ++k)
	{
		out2 << k << "," << dose.x(k) << ","
			<< dose(k, cy, cz) << "\n";
	}
	out2.close();
	

    return 0;
}
