#include <cmath>
#include "dtw.hpp"
#include "parameters.hpp"

namespace parole {

	std::valarray<std::valarray<double> > dtw_distance(const std::valarray<std::valarray<double> >& X,
			const std::valarray<std::valarray<double> >& Y) {
		const int n_col_X = X.size()-2;
		const int n_col_Y = Y.size()-2;
		std::valarray<std::valarray<double> > dist(std::valarray<double>(0.0, n_col_Y), n_col_X);
		int ligne;

    #pragma omp parallel for private(ligne)
		for(ligne = 0; ligne < n_col_Y; ++ligne) {
			for(int colonne = 0; colonne < n_col_X; ++colonne) {
				auto epsilon = X[colonne] - Y[ligne];
				dist[colonne][ligne] = std::sqrt((epsilon*epsilon).sum());
			}
		}
		return dist;
	}

	double dtw(const std::valarray<std::valarray<double> >& X,
			const std::valarray<std::valarray<double> >& Y) {
		const int n_col = X.size()-2;
		const int n_lign = Y.size()-2;
		auto dist = dtw_distance(X,Y);

		// create cost matrix
		std::valarray<std::valarray<double> > cost(std::valarray<double>(0.0, n_lign + 1), n_col + 1);

		cost[0][0] = 0;
		// calculate first row
		for(int i = 1; i <= n_col; i++)
			cost[i][0] = infinity;
		// calculate first column
		for(int j = 1; j <= n_lign; j++)
			cost[0][j] = infinity;

		// fill matrix
		for(int i = 1; i <= n_col; i++) {
			for(int j = 1; j <= n_lign; j++) {
				cost[i][j] = std::min({cost[i-1][j], cost[i][j-1], cost[i-1][j-1]}) + dist[i-1][j-1];
			}
		}
		return cost[n_col][n_lign];
	}
}
