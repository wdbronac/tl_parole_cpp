#pragma once
#include <valarray>

namespace parole {
	class MELFilter {
		private:
			std::valarray<double> coordinates;
			const double fs;
			static double f2F(double f); //linear frequency domain to MEL frequency domain conversion
			static std::valarray<double> F2f(const std::valarray<double>& F); //MEL frequency domain to linear frequency domain conversion
			//int Nf = coordinates.size() - 2;

		public:
			std::valarray<double> operator()(const std::valarray<double>& spectrum) const; //The spectrum is assumed to have the same sampling frequency
			MELFilter(const int Nf, const double fs);
	};
}
