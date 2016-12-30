#include <cmath>
#include "melfilter.hpp"

namespace parole {
	
	MELFilter::MELFilter(const int Nf, const double fs) : coordinates(0.0, Nf+2), fs(fs) {
		const double W = f2F(fs/2)/(Nf + 1); //(half) width of a filter in the MEL frequency domain
		for(int i = 1; i < Nf+2; ++i)
			coordinates[i] = coordinates[i-1] + W;
		coordinates = F2f(coordinates);
	}

	double MELFilter::f2F(double f) {
		return 2595*std::log10(1.0 + f/700.0);
	}

	std::valarray<double> MELFilter::F2f(const std::valarray<double>& F) {
		return 700*(std::pow(10.0, F/2959.0) - 1);
	}

	std::valarray<double> MELFilter::operator()(const std::valarray<double>& spectrum) const {
		const int Nf = coordinates.size() - 2;
		const double step = fs/(2*spectrum.size());
		std::valarray<double> result(0.0, Nf);
		for(int i = 0; i < Nf; ++i) {
			const double& fmin = coordinates[i];
			const double& fmid = coordinates[i+1];
			const double& fmax = coordinates[i+2];
			const double al = 1/(fmid - fmin);
			const double bl = -al*fmin;
			const double ar = 1/(fmid - fmax);
			const double br = -ar*fmax;
			const int jmid = std::floor(fmid/step);
			const int jmax = std::floor(fmax/step);
			int j;
			for(j = std::floor(fmin/step); j <= jmid; ++j)
				result[i] += al*spectrum[j] + bl;
			for(; j <= jmax; ++j)
				result[i] += ar*spectrum[j] + br;
		}
		return result;
	}
}
