#include "StatisticLogger.h"
#include <cmath>

using namespace std;

StatisticLogger::StatisticLogger(int nx, int ny, string fileName) : nx_(nx), ny_(ny) {
    open(fileName);
}

void StatisticLogger::open(string fileName) {
    isOpened_ = true;
    log_.open(fileName.c_str());
}

void StatisticLogger::close() {
    if(isOpened_) {
        log_.close();
        isOpened_ = false;
    }
}

void StatisticLogger::addHeadder() {
    log_ << "#time "
         << "potential_energy "
         << "kinetic_energy "
         << "mean_hight "
         << "skewness "
         << "kurtosis"
         << endl;
}

void StatisticLogger::addLog(double time, const double* eta, double potentialEnergy, double kineticEnergy) {
    double eta2 = 0;
    double eta3 = 0;
    double eta4 = 0;
    for(int i=0; i<nx_*ny_; i++) {
        eta2 += pow(eta[i], 2);
        eta3 += pow(eta[i], 3);
        eta4 += pow(eta[i], 4);
    }
    double rms = sqrt(eta2/nx_/ny_);
    double skewness = eta3/pow(rms, 3)/nx_/ny_;
    double kurtosis = eta4/pow(rms, 4)/nx_/ny_;
    log_ << time << " "
         << potentialEnergy << " "
         << kineticEnergy << " "
         << rms << " "
         << skewness << " "
         << kurtosis
         << endl;

}
