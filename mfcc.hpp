//mfcc.hpp

#pragma once
#include <valarray>
#include <ostream>

namespace parole {

	class MFCC;
	std::ostream& operator<<(std::ostream& os, const MFCC& mfcc);

	class MFCC {
		private:
			std::valarray<std::valarray<double>> m_samples;
			const int Nf;
			friend std::ostream& operator<<(std::ostream& os, const MFCC& mfcc);

		public:
			MFCC& operator=(MFCC&& mfcc);
			MFCC& operator=(const MFCC& mfcc);
			void process(const std::valarray<double>& meloutput, const double& energy);
			MFCC(int Nf, int size);
			MFCC(const MFCC& mfcc);
			MFCC(MFCC&& mfcc);
			MFCC();
			const std::valarray<std::valarray<double>>& samples() const;
	};
	
	std::valarray<double> linear(const int N);
	std::valarray<double> idct(const std::valarray<double>& signal);
	MFCC make_mfcc(double* samples, const int size, const int sample_rate);
}
