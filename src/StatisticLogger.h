#ifndef STATISTIC_LOGGER__
#define STATISTIC_LOGGER__

#include <fstream>
#include <string>

class StatisticLogger
{
public:
    explicit StatisticLogger(int nx, int ny, std::string fileName);
    virtual ~StatisticLogger() { close(); }

    void addHeadder();
    void addLog(double time, const double* eta, double potentialEnergy, double kineticEnergy);
    void open(std::string fileName);
    void close();
private:
    int nx_, ny_;
    std::ofstream log_;
    bool isOpened_;

};

#endif /**  STATISTIC_LOGGER__ */
