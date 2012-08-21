#ifndef MYTIME_CC
#define MYTIME_CC


#include <iostream>
#include <time.h>
using namespace std;

void convertToTime(double etc) {
    int dd, hh, mm, ss;
    if (etc < 60) {
        std::cout << int(etc) << " secs";
    } else {
        ss = int(etc) % 60;
        mm = int(etc / 60);
        etc = etc / 60;
        if (etc < 60) {
            std::cout << mm << " mins " << ss << " secs";
        } else {
            hh = int(mm / 60);
            mm = mm % 60;
            etc = etc / 60;
            if (etc < 24) {
                std::cout << hh << " hrs " << mm << " mins " << ss << " secs";
            }
            else{
                dd = int(hh / 24);
                hh = hh % 24;
                std::cout << dd << " days " << hh << " hrs " << mm << " mins ";
            }
        }
    }
}

template <class T, class S>
void calculateETA(time_t _startTime, T calculationsCompleted, S totalCalculations) {
    std::cout << int((calculationsCompleted*100)/totalCalculations) << "% Complete";

    time_t end;
    time(&end);
    double etc = difftime(end, _startTime);
    std::cout << "\tRunning: ";
    convertToTime(etc);
    etc = (etc * (totalCalculations - calculationsCompleted)) / calculationsCompleted;

    std::cout << "\tRemaining: ";
    convertToTime(etc);
    std::cout << std::endl;
}

#endif
