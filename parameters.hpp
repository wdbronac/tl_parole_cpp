#include <boost/math/constants/constants.hpp>
#include <limits>
#include <chrono>

constexpr int TAILLE_BUFFER = 128;			// Nombre d'echantillons dans un buffer
constexpr int NSAMPLES = 128*2;						// Nombre de mfccs en mémoire 
constexpr int CANAUXIN = 1;							// Nombre de canaux en entree (1)
constexpr int CANAUXOUT = 2;						// Nombre de canaux en sortie (2)
constexpr int NSEGMENT = 4;							// Nombre de segements de la fenetre
constexpr int NBFILTRES = 20;						// Nombre de filtres de MEL
constexpr int TAILLE_FEN = TAILLE_BUFFER*NSEGMENT;
constexpr int DECALAGE_FEN = TAILLE_BUFFER;

constexpr double PI = boost::math::constants::pi<double>();
constexpr double FREQUENCE = 16000;			// Frequence d'echantillonnage.
constexpr double THRESHOLD = 1000;				//TODO : trouver le bon threshold
//constexpr double ENERGY_THRESHOLD = 4000;	//TODO : trouver le bon threshold
constexpr double infinity = std::numeric_limits<double>::infinity();
constexpr std::chrono::duration<double, std::milli> WAIT_DURATION(200);
