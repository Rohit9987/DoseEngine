#pragma once
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <cmath>

namespace doseengine::physics {
	
	struct AttenuationPoint
	{
		double energy_MeV;
		double mu_over_rho;
		double muen_over_rho;
	};

	class PhotonAttenuationTable
	{
	public:
		explicit PhotonAttenuationTable(std::vector<AttenuationPoint> points)
			: points_(std::move(points))
		{
			if(points_.size() < 2)
				throw std::runtime_error("PhotonAttenuationTable: need at least two points");

			std::sort(points_.begin(), points_.end(),
					[](const auto& a, const auto& b)
						{
							return a.energy_MeV < b.energy_MeV;
						});
		}

		double muOverRho(double energy_MeV) const
		{
			return interpolateLogLog(
					energy_MeV,
					[](const AttenuationPoint& p)
					{ return p.mu_over_rho;}
					);
		}

		double muenOverRho(double energy_MeV) const
		{
			return interpolateLogLog(
					energy_MeV,
					[](const AttenuationPoint& p)
					{ return p.muen_over_rho;}
					);
		}

		double muPerMm(double energy_MeV, double density_g_per_cm3 = 1.0) const
		{
			// μ = (μ/ρ) ρ
			// cm^-1 -> mm^-1 by dividing by 10
			
			return muOverRho(energy_MeV) * density_g_per_cm3/10.0;
		}


		double muenPerMm(double energy_MeV, double density_g_per_cm3 = 1.0) const
		{
			return muenOverRho(energy_MeV) * density_g_per_cm3/10.0;
		}

		static PhotonAttenuationTable water()
		{
			return PhotonAttenuationTable({
				{0.001, 4.078E+03, 4.065E+03},
				{0.0015, 1.376E+03, 1.372E+03},
				{0.002, 6.173E+02, 6.152E+02},
				{0.003, 1.929E+02, 1.917E+02},
				{0.004, 8.278E+01, 8.191E+01},
				{0.005, 4.258E+01, 4.188E+01},
				{0.006, 2.464E+01, 2.405E+01},
				{0.008, 1.037E+01, 9.915E+00},
				{0.010, 5.329E+00, 4.944E+00},
				{0.015, 1.673E+00, 1.374E+00},
				{0.020, 8.096E-01, 5.503E-01},
				{0.030, 3.756E-01, 1.557E-01},
				{0.040, 2.683E-01, 6.947E-02},
				{0.050, 2.269E-01, 4.223E-02},
				{0.060, 2.059E-01, 3.190E-02},
				{0.080, 1.837E-01, 2.597E-02},
				{0.100, 1.707E-01, 2.546E-02},
				{0.150, 1.505E-01, 2.764E-02},
				{0.200, 1.370E-01, 2.967E-02},
				{0.300, 1.186E-01, 3.192E-02},
				{0.400, 1.061E-01, 3.279E-02},
				{0.500, 9.687E-02, 3.299E-02},
				{0.600, 8.956E-02, 3.284E-02},
				{0.800, 7.865E-02, 3.206E-02},
				{1.000, 7.072E-02, 3.103E-02},
				{1.250, 6.323E-02, 2.965E-02},
				{1.500, 5.754E-02, 2.833E-02},
				{2.000, 4.942E-02, 2.608E-02},
				{3.000, 3.969E-02, 2.281E-02},
				{4.000, 3.403E-02, 2.066E-02},
				{5.000, 3.031E-02, 1.915E-02},
				{6.000, 2.770E-02, 1.806E-02},
				{8.000, 2.429E-02, 1.658E-02},
				{10.000, 2.219E-02, 1.566E-02},
				{15.000, 1.941E-02, 1.441E-02},
				{20.000, 1.813E-02, 1.382E-02}
			});
		}

	private:
		template<typename Getter>
		double interpolateLogLog(double energy_MeV, Getter getter) const
		{
			if(energy_MeV <= points_.front().energy_MeV)
				return getter(points_.front());
			if(energy_MeV >= points_.back().energy_MeV)
				return getter(points_.back());

			for(std::size_t i = 0; i + 1 < points_.size(); i++)
			{
				const auto& a = points_[i];
				const auto& b = points_[i+1];

				if(energy_MeV >= a.energy_MeV && energy_MeV <= b.energy_MeV)
				{
					const double logE = std::log(energy_MeV);
					const double logE1 = std::log(a.energy_MeV);
					const double logE2 = std::log(b.energy_MeV);

					const double logY1 = std::log(getter(a));
					const double logY2 = std::log(getter(b));

					const double t = (logE - logE1) / (logE2 - logE1);
					return std::exp(logY1 + t * (logY2 - logY1));
				}
			}
			throw std::runtime_error("PhotonAttenuationTable: interpolation failed");
		}

		std::vector<AttenuationPoint> points_;
	};
}
