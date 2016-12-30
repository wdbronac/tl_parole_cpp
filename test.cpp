#include <iostream>
#include <cmath>
#include <fftw3.h>
#include "mfcc.hpp"
#include "melfilter.hpp"
#include "parameters.hpp"

namespace test {
	void valarray() {
		//Tests valarray deep copying
		std::valarray<std::valarray<double>> M(std::valarray<double>(0.0,10),10);
		std::valarray<std::valarray<double>> X(M);
		M[0][0] = 1.0;
		for(const auto& x : X) {
			for(const auto& y : x)
				std::cout << y << ", ";
			std::cout << std::endl;
		}
	}

	void mfcc() {
		constexpr int size = TAILLE_FEN*NSAMPLES;
		double signal[size];
		for(int i = 0; i < size; ++i) signal[i] = std::sin(400*i*PI/FREQUENCE);
		auto mfcc = parole::make_mfcc(signal, size, FREQUENCE);
		std::cout << mfcc << std::endl;
	}

	void idct() {
		std::valarray<double> signal(TAILLE_BUFFER*2);
		for(int i = 0; i < TAILLE_BUFFER; ++i) signal[i] = i;
		for(int i = 0; i < TAILLE_BUFFER; ++i) signal[TAILLE_BUFFER+i] = TAILLE_BUFFER-i;
		for(int i = 0; i < TAILLE_BUFFER*2; ++i) std::cout << signal[i] << ", ";
		std::cout << std::endl;

		auto idct = parole::idct(signal);

		for(const auto& x : idct) std::cout << x << ", ";
		std::cout << std::endl;
	}

	void fft() {
		double signal[TAILLE_FEN];
		fftw_plan p;
		constexpr int s = TAILLE_FEN/2 + 1; 
		fftw_complex cfft[s];

		p = fftw_plan_dft_r2c_1d(TAILLE_FEN, signal, cfft, 0);
		for(int i = 0; i < TAILLE_FEN; ++i) {
			signal[i] = std::sin(400*i*PI/FREQUENCE);
		}

		fftw_execute(p);

		for(int j = 0; j < s-1; ++j)
			std::cout << std::sqrt(cfft[j][0]*cfft[j][0] + cfft[j][1]*cfft[j][1]) << ", ";
		std::cout << std::endl;
		fftw_destroy_plan(p);
	}

	void mel() {
		double signal[TAILLE_FEN];
		fftw_plan p;
		fftw_complex cfft[TAILLE_FEN/2 + 1];
		p = fftw_plan_dft_r2c_1d(TAILLE_FEN, signal, cfft, 0);
		for(int i = 0; i < TAILLE_FEN; ++i) signal[i] = std::sin(400*i*PI/FREQUENCE);

		fftw_execute(p);
		std::valarray<double> spectrum(TAILLE_FEN/2);
		for(int j = 0; j < TAILLE_FEN/2; ++j)
			spectrum[j] = std::sqrt(cfft[j][0]*cfft[j][0] + cfft[j][1]*cfft[j][1]);
		parole::MELFilter melfilter(NBFILTRES, FREQUENCE);

		for(const auto& x : melfilter(spectrum)) std::cout << x << ", ";
		fftw_destroy_plan(p);
		std::cout << std::endl;
	}

	void dtw() {
		std::cerr << "\e[1;31mNot implemented\e[0m" << std::endl;
	}
}

int main(int argc, char** argv) {
	std::string error_message("Please give one of the following arguments: mfcc, idct, fft, mel, dtw, valarray");
	if(argc == 2) {
		std::string arg(argv[1]);
		if(arg == "mfcc")	test::mfcc();
		else if(arg == "idct") test::idct();
		else if(arg == "fft") test::fft();
		else if(arg == "mel") test::mel();
		else if(arg == "dtw") test::dtw();
		else if(arg == "valarray") test::valarray();
		else {
			std::cerr <<  error_message << std::endl; 
			return 1;
		}
		return 0;
	} else {
		std::cerr << error_message << std::endl;
		return 1;
	}
}
