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



/*
	auto print_terma = [terma](const int k)-> double
	{
		double t = terma(75, 75, k);
		std::cout << "CAX depth "<< terma.z(k) << " mm TERMA = "
              << t << "\n";
		return t;
	};
	double t = print_terma(100);
	t /= print_terma(50);
	std::cout << "TPR20/10 is approximately: " << t << '\n';
	


	// csv export
	std::ofstream out("terma_cax.csv");
	out << "k,z_mm,terma\n";

	const std::size_t cx = terma.nx() / 2;
	const std::size_t cy = terma.ny() / 2;
	for (std::size_t k = 0; k < terma.nz(); ++k)
	{
		out << k << "," << terma.z(k) << ","
			<< terma(cx, cy, k) << "\n";
	}
	out.close();

	
	std::ofstream out2("terma_crossline.csv");
	out << "k,x_mm,terma\n";

	const std::size_t cz = 50;
	for (std::size_t k = 0; k < terma.nx(); ++k)
	{
		out2 << k << "," << terma.x(k) << ","
			<< terma(k, cy, cz) << "\n";
	}
	out2.close();
	

	auto kernel1 = []()
	{
		auto kernel = doseengine::kernel::KernelLoader::loadKernelBin(
				"../../Kernel/kernel-downsampled/1.0MeV_kernel.bin",
			true
		);

		std::cout << "Kernel loaded\n";
		std::cout << "Energy: " << kernel.energyMeV() << " MeV\n";
		std::cout << "Shape: "
				  << kernel.nx() << " x "
				  << kernel.ny() << " x "
				  << kernel.nz() << "\n";

		std::cout << "Spacing: "
				  << kernel.dx_mm() << ", "
				  << kernel.dy_mm() << ", "
				  << kernel.dz_mm() << " mm\n";

		std::cout << "Offset: "
				  << kernel.x0_mm() << ", "
				  << kernel.y0_mm() << ", "
				  << kernel.z0_mm() << " mm\n";

		std::cout << "Kernel sum: " << kernel.sum() << "\n";

		// save csv
		std::ofstream out("kernel_cax.csv");
		out << "k,z,kernel\n";

		const std::size_t cx = kernel.nx() / 2;
		const std::size_t cy = kernel.ny() / 2;
		for (std::size_t k = 0; k < kernel.nz(); ++k)
		{
			out << k << "," <<  k*kernel.dz_mm() + kernel.z0_mm() << ","
				<< kernel.values()(cx, cy, k) << "\n";
		}
		out.close();

		
		std::ofstream out2("kernel_crossline.csv");
		out2 << "k,x_mm,kernel\n";

		const std::size_t cz = 40;
		for (std::size_t k = 0; k < kernel.nx(); ++k)
		{
			out2 << k << "," << k * kernel.dx_mm() + kernel.x0_mm() << ","
				<< kernel.values()(k, cy, cz) << "\n";
		}
		out2.close();
			
	};

	kernel1();


*/




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
		false
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

	auto kernel = doseengine::kernel::KernelLoader::loadKernelBin(
			"../../Kernel/kernel-downsampled/1.0MeV_kernel.bin",
			true
		);

	auto coarseKernel = doseengine::kernel::KernelStencil::binToCoarseKernel(
		kernel,
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

	std::cout << "Dose at 10 cm : " << dose(dose.nx()/2, dose.ny()/2, 50) << "\n"
			  << "Dose at 10 cm : " << dose(dose.nx()/2, dose.ny()/2, 100) << "\n";

    return 0;
}
