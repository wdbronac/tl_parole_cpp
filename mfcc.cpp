#include <cmath>
#include <fftw3.h>
#include "mfcc.hpp"
#include "melfilter.hpp"
#include "parameters.hpp"

namespace parole {
	//Helper function to get a valarray from 0 to N
	std::valarray<double> linear(const int N) {
		std::valarray<double> result(0.0,N);
		for(int i = 1; i < N; ++i)
			result[i] = i;
		return result;
	}

	//inverse discrete cosine transform
	std::valarray<double> idct(const std::valarray<double>& signal) {
		const int N = signal.size();
		const auto k = linear(N);
		std::valarray<double> result(0.0, N);
		double n = 0.0;
		std::valarray<double> w(std::sqrt(std::valarray<double>(2.0/N, N)));
		w[0] = std::sqrt(1.0/N);
		for(auto& x : result) {
			x = (w*signal*std::cos((n+0.5)*k*PI/N)).sum();
			++n;
		}
		return result;
	}

	MFCC::MFCC(int Nf, int size) : m_samples(std::valarray<double>(3*(Nf+1)), size), Nf(Nf) {
	}

	MFCC::MFCC() : m_samples(), Nf(20) {
	}

	const std::valarray<std::valarray<double>>& MFCC::samples() const {
		return m_samples;
	}

	void MFCC::process(const std::valarray<double>& meloutput, const double& energy) {
		auto coeffs = std::slice(0,Nf,1);
		auto d_coeffs = std::slice(Nf,Nf,1);
		auto dd_coeffs = std::slice(2*Nf,Nf,1);

		m_samples = m_samples.cshift(1);
		m_samples[m_samples.size()-1][coeffs] = idct(std::log(std::abs(meloutput)));
		m_samples[m_samples.size()-1][3*Nf] = energy;
		m_samples[m_samples.size()-2][3*Nf+1] = energy - m_samples[m_samples.size()-2][3*Nf];
		m_samples[m_samples.size()-3][3*Nf+2] = m_samples[m_samples.size()-2][3*Nf+1] - m_samples[m_samples.size()-3][3*Nf+1];
		m_samples[m_samples.size()-2][d_coeffs] = static_cast<std::valarray<double>>(m_samples[m_samples.size()-1][coeffs]) -
			static_cast<std::valarray<double>>(m_samples[m_samples.size()-2][coeffs]);
		m_samples[m_samples.size()-3][dd_coeffs] = static_cast<std::valarray<double>>(m_samples[m_samples.size()-2][d_coeffs]) -
			static_cast<std::valarray<double>>(m_samples[m_samples.size()-3][d_coeffs]);
	}

	MFCC& MFCC::operator=(const MFCC& mfcc) {
		if(Nf != mfcc.Nf) throw "Wrong value of Nf";
		m_samples = mfcc.m_samples;
		return *this;
	}

	MFCC& MFCC::operator=(MFCC&& mfcc) {
		if(Nf != mfcc.Nf) throw "Wrong value of Nf";
		m_samples = std::move(mfcc.m_samples);
		return *this;
	}

	MFCC::MFCC(const MFCC& mfcc) : m_samples(mfcc.m_samples), Nf(mfcc.Nf) {
	}

	MFCC::MFCC(MFCC&& mfcc) : m_samples(std::move(mfcc.m_samples)), Nf(mfcc.Nf) {
	}

	void subsampling(double* samples, const int size, const int sample_rate) {
		if(sample_rate == FREQUENCE) return;
		const int r = sample_rate % int(FREQUENCE);
		const int q = sample_rate / FREQUENCE;
		double subsamples[size/q+1];
		for(int i = 0; i < size/q; ++i)
			for(int j = 0; j < q; ++j) subsamples[i] += samples[j+i]/q;
		for(int j = 0; j < r; ++j) subsamples[size/q] += samples[size-r+j]/r;
		for(int i = 0; i < size/q; ++i) samples[i] = subsamples[i];
	}

	MFCC make_mfcc(double* samples, const int size, const int sample_rate) {
		if(sample_rate > FREQUENCE)
			subsampling(samples, size, sample_rate);
		else if(sample_rate < FREQUENCE)
			throw	"The reference signal's sample rate is too low,\n\
				please reduce the incoming signal's sample rate";
		MFCC mfcc(NBFILTRES, size/DECALAGE_FEN - NSEGMENT);
		MELFilter melfilter(NBFILTRES, FREQUENCE);
		fftw_plan p;
		fftw_complex cfft[TAILLE_FEN/2 + 1];
		std::valarray<double> spectrum(TAILLE_FEN/2);
		for(int i = 0; i < size-TAILLE_FEN; i += DECALAGE_FEN) {
			p = fftw_plan_dft_r2c_1d(TAILLE_FEN, &samples[i], cfft, FFTW_ESTIMATE);
			fftw_execute(p);
			for(int j = 0; j < TAILLE_FEN/2; ++j)
				spectrum[j] = std::sqrt(cfft[j][0]*cfft[j][0] + cfft[j][1]*cfft[j][1]);

			mfcc.process(melfilter(spectrum), (spectrum*spectrum).sum());
			fftw_destroy_plan(p);
		}
		return mfcc;
	}

	std::ostream& operator<<(std::ostream& os, const MFCC& mfcc) {
		for(const auto& x : mfcc.m_samples) {
			for(const auto& y : x) os << y << ',';
			os << std::endl;
		}
		return os;
	}
}
