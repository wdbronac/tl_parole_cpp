//C++ std libs
#include <string>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <valarray>
#include <thread>
#include <mutex>
#include <chrono>

//ROS
#include <ros/ros.h>
#include <geometry_msgs/Twist.h>

//Boost
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>

//C libs
#include <sndfile.h>
#include <fftw3.h>
#include <supelec-audio2.h>

//local includes
#include "dtw.hpp"
#include "melfilter.hpp"
#include "mfcc.hpp"
#include "parameters.hpp"

// le dernier segment est mis a jour a chaque iteration
const std::string PATHREF("./audio");
std::mutex run, mfcc_lock;
bool mfcc_running;

void mfcc_processing(parole::MFCC& mfcc) {
  int i,j;
  int bufsize = TAILLE_BUFFER;
  int frequence = FREQUENCE;

  // Ouverture audio :
  if (!(frequence = initialiseSonCanaux(frequence, CANAUXIN, CANAUXOUT, TAILLE_BUFFER))) {
    std::cerr << "Problème d'initialisation du système audio." << std::endl;
    mfcc_running = false;
    return;
  }

  std::cout << "Fréquence : " << frequence << ", taille entrée : " << recupereTailleBufferEntree() << std::endl
  echantillon* buffer;
  double signal[NSEGMENT*TAILLE_BUFFER*CANAUXIN];
  fftw_plan p;
  fftw_complex cfft[NSEGMENT*TAILLE_BUFFER*CANAUXIN/2 + 1];
  p = fftw_plan_dft_r2c_1d(NSEGMENT*TAILLE_BUFFER*CANAUXIN, signal, cfft, 0);
  std::valarray<double> spectrum(TAILLE_BUFFER*CANAUXIN*NSEGMENT/2);
  double energy(0);
  parole::MELFilter melfilter(NBFILTRES, frequence);

  for(i = 0; i < NSEGMENT-1; ++i) {
    const int jmax = (i+1)*bufsize*CANAUXIN;
    buffer = recupereBufferEntreeEtTaille(&bufsize);
    for(j = i*bufsize*CANAUXIN; j < jmax; ++j)
      signal[j] = static_cast<double>(buffer[j % bufsize]);
    libereBufferEntree();
  }

  for(; j < NSEGMENT*TAILLE_BUFFER*CANAUXIN; ++j)
    signal[j] = 0.0; //zero padding in case the real buffer is smaller than intended

  while(run.try_lock()) {
    i = (NSEGMENT-1)*bufsize*CANAUXIN;
    for(buffer = recupereBufferEntreeEtTaille(&bufsize); i < NSEGMENT*bufsize*CANAUXIN; ++i)
      signal[i] = static_cast<double>(buffer[i % bufsize]);
    libereBufferEntree();

    fftw_execute(p);

    for(int j = 0; j < NSEGMENT*bufsize*CANAUXIN/2; ++j)
      spectrum[j] = std::sqrt(cfft[j][0]*cfft[j][0] + cfft[j][1]*cfft[j][1]);

    std::move(signal + bufsize*CANAUXIN, signal + NSEGMENT*bufsize*CANAUXIN, signal);
    energy = (spectrum*spectrum).sum();
    mfcc_lock.lock();
    mfcc.process(melfilter(spectrum), energy);
    mfcc_lock.unlock();
    run.unlock();
  }
  fftw_destroy_plan(p);
}


class RobotDriver {
  private:
    //! The node handle we'll be using
    ros::NodeHandle nh_;
    //! We will be publishing to the "/base_controller/command" topic to issue commands
    ros::Publisher cmd_vel_pub_;

  public:
    //! ROS node initialization
    RobotDriver(ros::NodeHandle &nh) {
      nh_ = nh;
      //set up the publisher for the cmd_vel topic
      cmd_vel_pub_ = nh_.advertise<geometry_msgs::Twist>("/turtle1/cmd_vel", 1000);
    }

		void do_for(const geometry_msgs::Twist& action, const std::chrono::duration<double, std::milli>& duration) {
			auto start = std::chrono::system_clock::now();
		while(std::chrono::system_clock::now() - start < duration) {
				cmd_vel_pub_.publish(action);
				std::this_thread::sleep_for(WAIT_DURATION);
			}
		}

    bool driveParole() {
      parole::MFCC mfcc(NBFILTRES, NSAMPLES);
      mfcc_running = true;

      geometry_msgs::Twist base_cmd;
      //On met les MFCC de référence dans la mémoire pour pouvoir comparer
      //ouvre tous les fichiers un a un
      boost::filesystem::path path(PATHREF);
      std::unordered_map<std::string, parole::MFCC> mfccs;
      if(!boost::filesystem::exists(path)) {
        std::cerr << "Invalid path!" << std::endl;
        return 0;
      } else if(boost::filesystem::is_directory(path)) {
				std::ofstream output;
        std::cout << "Found the following files in " << path << ":" << std::endl;
        for(boost::filesystem::directory_iterator it(path), end; it != end; ++it) {
          if(boost::filesystem::is_regular_file(it->status())) {
            std::string filename(it->path().filename().string());
            std::cout << filename << std::endl;
            SNDFILE *infile = NULL;
            SF_INFO sfinfo;
            infile = sf_open(it->path().string().c_str(), SFM_READ, &sfinfo);
            if(sfinfo.channels != 1) throw "Too many channels";

            int num_items = sfinfo.channels*sfinfo.frames;
            int sample_rate = sfinfo.samplerate;
            int N = num_items;
            double samples[N];
            sf_read_double(infile, samples, num_items);
						std::size_t pos = filename.find(".wav");
						std::string name;
						if(pos != std::string::npos) name = filename.substr(0,pos);
						std::cout << "Writing MFCC in " << name << ".csv" << std::endl;
            mfccs[name] = parole::make_mfcc(samples, N, sample_rate);
						output.open(name + ".csv");
						output << mfccs[name] << std::endl;
						output.close();
          }
        }
      }
			parole::MFCC mfcc_copy;
      std::thread mfcc_thread(mfcc_processing, std::ref(mfcc));
			int step = 6;
			while(nh_.ok() && step) {
				double best_dtw(infinity);
				std::string best_order;
				std::string last_order = "stop";
				if(mfcc_running) {
					mfcc_lock.lock();
					mfcc_copy = mfcc;
					mfcc_lock.unlock();
				} else {
					switch(step) {
						case 6:
							std::cerr << "MFCC processing stopped unexpectedly. Launching test routine..." << std::endl;
							std::cout << "EN AVANT !" << std::endl;
							mfcc_copy = mfccs["enavant"];
							break;
						case 5:
							std::cout << "A DROITE !" << std::endl;
							mfcc_copy = mfccs["adroite"];
							break;
						case 4:
							std::cout << "EN AVANT !" << std::endl;
							mfcc_copy = mfccs["enavant"];
							break;
						case 3:
							std::cout << "A GAUCHE !" << std::endl;
							mfcc_copy = mfccs["agauche"];
							break;
						case 2:
							std::cout << "EN AVANT !" << std::endl;
							mfcc_copy = mfccs["enavant"];
							break;
						case 1:
							std::cout << "STOP !" << std::endl;
							mfcc_copy = mfccs["stop"];
							break;
					}
					--step;
				}

				for(const auto& item : mfccs) {
					const auto& ref = item.second;
					const auto& filename = item.first;
					if(last_order == "stop" || last_order == "enavant")
						if(boost::starts_with(filename, last_order)) continue;
					if(last_order == "adroite" || last_order == "agauche")
						if(boost::starts_with(filename, "enavant")) continue;
					double current_dtw = parole::dtw(mfcc_copy.samples(), ref.samples());
					if(current_dtw < best_dtw) {
						best_dtw = current_dtw;
						best_order = filename;
					}
				}

				if(best_dtw < THRESHOLD) {
					std::cout << best_order << ": " << best_dtw << std::endl;
					if(boost::starts_with(best_order, "stop")) {
						last_order = "stop";
						base_cmd.linear.x = 0.0;
						base_cmd.angular.z = 0.0;
						do_for(base_cmd, WAIT_DURATION*25);
					} else if(boost::starts_with(best_order, "agauche")) {
						last_order = "agauche";
						base_cmd.angular.z = 1;
						base_cmd.linear.x = 0.25;
						do_for(base_cmd, WAIT_DURATION*5);
						base_cmd.angular.z = 0.0;
						cmd_vel_pub_.publish(base_cmd);
					} else if(boost::starts_with(best_order, "enavant")) {
						last_order = "enavant";
						base_cmd.linear.x = 0.25;
						base_cmd.angular.z = 0.0;
						do_for(base_cmd, WAIT_DURATION*25);
					} else if(boost::starts_with(best_order, "adroite")) {
						last_order = "adroite";
						base_cmd.angular.z = -1;
						base_cmd.linear.x = 0.25;
						do_for(base_cmd, WAIT_DURATION*5);
						base_cmd.angular.z = 0.0;
						cmd_vel_pub_.publish(base_cmd);
					}
				}
			}

			base_cmd.linear.x = 0.0;
			base_cmd.angular.z = 0.0;
			cmd_vel_pub_.publish(base_cmd);

      run.lock();
      mfcc_thread.join();
      run.unlock();
      return true;
    }
};

int main(int argc, char** argv) {
  //init the ROS node
  ros::init(argc, argv, "robot_driver");
  ros::NodeHandle nh;
  RobotDriver driver(nh);
  return !driver.driveParole();
}
