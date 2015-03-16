#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <complex>
#include "HosmCalculator.h"
#include <boost/timer.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include "Property.h"
#include "Exception.h"
#include "util.h"
#include "Logger.h"
#include "Define.h"
//#include <fftw3.h>
//#include <mpi.h>
//#include <fftw3-mpi.h>
#include "SpectrumCreater.h"
#include "WaveBinData.h"
#include "StatisticLogger.h"


using namespace boost;
using namespace std;


int main(int argc, char **argv) {

    bool printlog = false;
    if(argc > 1 && string(argv[1]) == "--print") printlog = true;
//    MPI::Init(argc, argv);
//    fftw_mpi_init();
//    int myrank = MPI::COMM_WORLD.Get_rank();


    boost::timer timer;

    /**          load settings     START        */
    string inputFile = "input.in";
    flow::util::Property property;
    try {
        property.read( inputFile );
    } catch( flow::exception::IO& e) {
        cerr << e.info() << endl;
        cerr << "file open failed." << endl;
        return 0;
    }

    double lx;
    double ly;
    int nx;
    int ny;
    double dt;
    double endTime;
    string spectrum;
    double jonswap_E;
    double jonswap_wp;
    double gravity;
    int outputStep;

    int level;
    HosmCalculator::TimeIntegration timeIntegration;

    try {
        lx = property.get<double>("Lx");
        ly = property.get<double>("Ly");
        nx = property.get<int>("nx");
        ny = property.get<int>("ny");
        level = property.get<int>("level");
        dt      = property.get<double>("dt");
        endTime = property.get<double>("endTime");
        string ti = property.get<string>("timeIntegration");
        spectrum = property.get<string>("spectrum", "simple");
        outputStep = property.get<int>("outputStep", 1);
        if(spectrum == "jonswap" || spectrum == "JONSWAP") {
            jonswap_E = property.get<double>("jonswap_E");
            jonswap_wp = property.get<double>("jonswap_wp");
        }
        gravity = property.get<double>("gravity", 9.8);
        cout << "Time integration : " << ti << endl;
        if(ti == "Euler") {
            timeIntegration = HosmCalculator::Euler;
        } else if (ti == "RK4") {
            timeIntegration = HosmCalculator::RK4;
        } else if (ti == "RK3") {
            timeIntegration = HosmCalculator::RK3;
        } else if (ti == "RK2") {
            timeIntegration = HosmCalculator::RK2;
        } else {
            cerr << "no such time integration" << endl;
            throw;
        }
    } catch (flow::exception::Parse& e) {
        cerr << e.info() << endl;
    }
    double dx = lx/nx;
    double dy = ly/ny;

    /**          load settings    END           */


//    ptrdiff_t local_nx, local_x_start;
//    fftw_mpi_local_size_2d(nx, ny, MPI::COMM_WORLD, &local_nx, &local_x_start);

    ubvector_d eta(nx*ny, 0);
    ubvector_d phi(nx*ny, 0);
    double time = 0;
    int step = 0;

    HosmCalculator calculator(nx, ny, dx, dy, level, timeIntegration, gravity);
    const SpectrumProperty& specProp = calculator.specProp();

    /**          create initial wave field       */
    if(spectrum == "jonswap" || spectrum == "JONSWAP") {
        createJONSWAPSpectrum(&eta[0], &phi[0], nx, ny, specProp, jonswap_E, jonswap_wp, gravity);
    } else {
        createSpectrum("spectrum.in", &eta[0], &phi[0], nx, ny, specProp, gravity);
    }
    cout << "Initial Spectrum created." << endl;

    // output
    writeGrid(&(eta[0]),"surf_ini.out", nx, ny);
    writeGrid(&(phi[0]),"phi_ini.out", nx, ny);
    writeSpectrum(&(eta[0]), "spec_ini.out", nx, ny, specProp);


    /** trace each targetting waves */
    Logger logger;
    logger.open("logging.in");
    logger.set(nx, ny);
    logger.addLog(time, eta, phi, specProp);

    WaveBinData wavebin("wave.bin", nx, ny, specProp);
    wavebin.addHeader(level, lx, ly, dt, endTime, outputStep, gravity);
    wavebin.addData(&eta[0],  &phi[0], time);

    StatisticLogger statisticLogger(nx, ny, "statistic.log");
    statisticLogger.addHeadder();


    cout << "Start main loop" << endl;
    while(time <= endTime) {
        calculator.next(eta, phi, dt);

        step++;
        time += dt;

        logger.addLog(time, eta, phi, specProp);
        statisticLogger.addLog(time, &eta[0], calculator.getPotentialEnergy(), calculator.getKineticEnergy());

        if(step % outputStep == 0) {
            wavebin.addData(&eta[0],  &phi[0], time);
        }

        if(printlog) {
            cout << step << " "
                 << time << " "
                 << calculator.getPotentialEnergy() + calculator.getKineticEnergy()
                 << endl;
        }

        if(isnan( calculator.getPotentialEnergy() ) ) {
            cerr << "NaN value is detected." << endl;
            break;
        }

    }

    statisticLogger.close();
    wavebin.close();
    logger.close();


    // output
    writeGrid(&(eta[0]),"eta.out", nx, ny);
    writeGrid(&(phi[0]),"phi.out", nx, ny);
    writeSpectrum(&(eta[0]), "spec.out", nx, ny, specProp);


    cout << "Total " << step << " steps" << endl;
    cout << timer.elapsed() << "[sec]" << endl;
    cout << timer.elapsed()/step << "[sec] for 1 step" << endl;


    return 0;
}
