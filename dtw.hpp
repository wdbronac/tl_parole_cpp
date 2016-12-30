//dtw.hpp

#pragma once
#include <valarray>

namespace parole {
	double dtw(const std::valarray<std::valarray<double> >& t1, const std::valarray<std::valarray<double> >& t2);
	std::valarray<std::valarray<double> > dtw_distance(const std::valarray<std::valarray<double> >& X, const std::valarray<std::valarray<double> >& Y);
}

