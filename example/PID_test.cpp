//---------------------------------------------------------------------------------------------------------------------
//  Vertical Engineering Solutions
//---------------------------------------------------------------------------------------------------------------------
// 
//  Copyright 2020 Vertical Engineering Solutions  - All Rights Reserved
// 
//  Unauthorized copying of this file, via any medium is strictly prohibited Proprietary and confidential.
// 
//  All information contained herein is, and remains the property of Vertical Engineering Solutions.  The 
//  intellectual and technical concepts contained herein are proprietary to Vertical Engineering Solutions 
//  and its suppliers and may be covered by UE and Foreign Patents, patents in process, and are protected 
//  by trade secret or copyright law. Dissemination of this information or reproduction of this material is 
//  strictly forbidden unless prior written permission is obtained from Adobe Systems Incorporated.
//
//---------------------------------------------------------------------------------------------------------------------
//
//  Maintainer: pramon@vengineerings.com
//
//---------------------------------------------------------------------------------------------------------------------

#include <QApplication>
#include <QVector>
#include "qcustomplot.h"

#include <pidpp/PID.h>

#include<iostream>

#include <boost/numeric/odeint.hpp>
using namespace boost::numeric::odeint;

typedef std::vector<double> state_type;
/* The rhs of x' = f(x) */
class FirstOrderSystem {
    public:
        FirstOrderSystem( double _c ) : c_(_c) { }

        void operator() ( const state_type &x , state_type &dxdt, const double _t ) {
            dxdt[0] = -1/c_ * x[0] + x[1];
        }

    private:
        double c_;
};

class SecondOrderSystem {
    public:
        SecondOrderSystem( double _a, double _b ) : a_(_a), b_(_b) { }

        void setU(double _u){
            u_ = _u;
        }

        void operator() ( const state_type &x , state_type &dxdt, const double _t ) {
            dxdt[0] = x[1];
            dxdt[1] = u_ - b_*x[0] -a_* x[1];
        }

    private:
        double a_, b_;
        double u_;
};


int main(int _argc, char** _argv){
    QApplication a(_argc, _argv);

    //-----------------------------------------------------------------------------------------------------------------
    // Simulation
    //-----------------------------------------------------------------------------------------------------------------
    
    
    float kp = 4.0;
    float ki = 5.0;
    float kd = 0.0;

    //-----------------------------------------------------------------------------------------------------------------
    // Free run
    //-----------------------------------------------------------------------------------------------------------------
    std::vector<double> stateFree = {0, 0};
    std::vector<double> valuesFree;
    std::vector<double> stepFree;
    std::vector<double> timeFree;
    
    {
        SecondOrderSystem system(1,1);
        runge_kutta4< state_type > stepper;
        system.setU(0);
        integrate_const( stepper ,  system, stateFree , 0.0 , 10.0 , 0.01 );

        const double dt = 0.01;
        for( double t=0.0 ; t<10.0 ; t+= dt ){
            if(t >1){
                system.setU(1);
                stepFree.push_back(1);
            }else{
                stepFree.push_back(0);
            }

            stepper.do_step( system , stateFree , t , dt );
            timeFree.push_back(t);
            valuesFree.push_back(stateFree[0]);
        }
    }

    //-----------------------------------------------------------------------------------------------------------------
    // PID no windup
    //-----------------------------------------------------------------------------------------------------------------
    std::vector<double> stateNoWindup = {0, 0};
    std::vector<double> valuesNoWindup;
    std::vector<double> stepNoWindup;
    std::vector<double> timeNoWindup;
    
    {
        SecondOrderSystem system(1,1);
        runge_kutta4< state_type > stepper;
        system.setU(0);
        integrate_const( stepper ,  system, stateNoWindup , 0.0 , 10.0 , 0.01 );
        pidpp::PID pid(kp, ki, kd,-2,2);
        pid.reference(1);

        const double dt = 0.01;
        for( double t=0.0 ; t<10.0 ; t+= dt ){
            if(t >1){
                stepNoWindup.push_back(1);
                system.setU(pid.update(stateNoWindup[0], dt));
            }else{
                stepNoWindup.push_back(0);
            }

            stepper.do_step( system , stateNoWindup , t , dt );
            timeNoWindup.push_back(t);
            valuesNoWindup.push_back(stateNoWindup[0]);
        }
    }

    //-----------------------------------------------------------------------------------------------------------------
    // PID windup sat
    //-----------------------------------------------------------------------------------------------------------------
    std::vector<double> stateWindupSat = {0, 0};
    std::vector<double> valuesWindupSat;
    std::vector<double> stepWindupSat;
    std::vector<double> timeWindupSat;
    
    {
        SecondOrderSystem system(1,1);
        runge_kutta4< state_type > stepper;
        system.setU(0);
        integrate_const( stepper ,  system, stateWindupSat , 0.0 , 10.0 , 0.01 );
        pidpp::PID pid(kp, ki, kd,-2,2);
        pid.reference(1);
        pid.setAntiWindup(pidpp::PID::AntiWindupMethod::Saturation, {-0.2, 0.2});

        const double dt = 0.01;
        for( double t=0.0 ; t<10.0 ; t+= dt ){
            if(t >1){
                stepWindupSat.push_back(1);
                system.setU(pid.update(stateWindupSat[0], dt));
            }else{
                stepWindupSat.push_back(0);
            }

            stepper.do_step( system , stateWindupSat , t , dt );
            timeWindupSat.push_back(t);
            valuesWindupSat.push_back(stateWindupSat[0]);
        }
    }
    
    //-----------------------------------------------------------------------------------------------------------------
    // PID windup Back calculation
    //-----------------------------------------------------------------------------------------------------------------
    std::vector<double> stateBackCalc = {0, 0};
    std::vector<double> valuesBackCalc;
    std::vector<double> stepBackCalc;
    std::vector<double> timeBackCalc;
    
    {
        SecondOrderSystem system(1,1);
        runge_kutta4< state_type > stepper;
        system.setU(0);
        integrate_const( stepper ,  system, stateBackCalc , 0.0 , 10.0 , 0.01 );
        pidpp::PID pid(kp, ki, kd,-2,2);
        pid.reference(1);
        pid.setAntiWindup(pidpp::PID::AntiWindupMethod::BackCalculation, {1.2});

        const double dt = 0.01;
        for( double t=0.0 ; t<10.0 ; t+= dt ){
            if(t >1){
                stepBackCalc.push_back(1);
                system.setU(pid.update(stateBackCalc[0], dt));
            }else{
                stepBackCalc.push_back(0);
            }

            stepper.do_step( system , stateBackCalc , t , dt );
            timeBackCalc.push_back(t);
            valuesBackCalc.push_back(stateBackCalc[0]);
        }
    }

    //-----------------------------------------------------------------------------------------------------------------
    // PID windup clamp
    //-----------------------------------------------------------------------------------------------------------------
    std::vector<double> stateClamp = {0, 0};
    std::vector<double> valuesClamp;
    std::vector<double> stepClamp;
    std::vector<double> timeClamp;
    
    {
        SecondOrderSystem system(1,1);
        runge_kutta4< state_type > stepper;
        system.setU(0);
        integrate_const( stepper ,  system, stateClamp , 0.0 , 10.0 , 0.01 );
        pidpp::PID pid(kp, ki, kd,-2,2);
        pid.reference(1);
        pid.setAntiWindup(pidpp::PID::AntiWindupMethod::Clamping, {1});

        const double dt = 0.01;
        for( double t=0.0 ; t<10.0 ; t+= dt ){
            if(t >1){
                stepClamp.push_back(1);
                system.setU(pid.update(stateClamp[0], dt));
            }else{
                stepClamp.push_back(0);
            }

            stepper.do_step( system , stateClamp , t , dt );
            timeClamp.push_back(t);
            valuesClamp.push_back(stateClamp[0]);
        }
    }

    //-----------------------------------------------------------------------------------------------------------------
    // PLOT
    //-----------------------------------------------------------------------------------------------------------------
    QCustomPlot plot;
    plot.setInteraction(QCP::iSelectPlottables, true);
    plot.legend->setVisible(true);
    plot.legend->setBrush(QBrush(QColor(255,255,255,230)));
    plot.setMinimumWidth(800);
    plot.setMinimumHeight(500);
    // Step
    QPen graphPen;
    graphPen.setColor(QColor(255,0,0));
    graphPen.setWidthF(4);
    plot.addGraph()->setPen(graphPen);
    plot.graph(0)->setData(QVector<double>::fromStdVector(timeFree), QVector<double>::fromStdVector(stepFree));
    plot.graph(0)->setName("Target Step");
    // Open step response
    graphPen.setColor(QColor(0,255,0));
    graphPen.setWidthF(2);
    plot.addGraph()->setPen(graphPen);
    plot.graph(1)->setData(QVector<double>::fromStdVector(timeFree), QVector<double>::fromStdVector(valuesFree));
    plot.graph(1)->setName("Free form");
    // No wind up PID
    graphPen.setColor(QColor(0,0,255));
    graphPen.setWidthF(2);
    plot.addGraph()->setPen(graphPen);
    plot.graph(2)->setData(QVector<double>::fromStdVector(timeNoWindup), QVector<double>::fromStdVector(valuesNoWindup));
    plot.graph(2)->setName("PID no antiwindup");
    
    // wind up sat PID
    graphPen.setColor(QColor(0,255,255));
    graphPen.setWidthF(2);
    plot.addGraph()->setPen(graphPen);
    plot.graph(3)->setData(QVector<double>::fromStdVector(timeWindupSat), QVector<double>::fromStdVector(valuesWindupSat));
    plot.graph(3)->setName("with saturation");

    // wind up back calculation
    graphPen.setColor(QColor(255,0,255));
    graphPen.setWidthF(2);
    plot.addGraph()->setPen(graphPen);
    plot.graph(4)->setData(QVector<double>::fromStdVector(timeBackCalc), QVector<double>::fromStdVector(valuesBackCalc));
    plot.graph(4)->setName("Back Calculation");

    // wind up Clamp
    graphPen.setColor(QColor(255,255,0));
    graphPen.setWidthF(2);
    plot.addGraph()->setPen(graphPen);
    plot.graph(5)->setData(QVector<double>::fromStdVector(timeClamp), QVector<double>::fromStdVector(valuesClamp));
    plot.graph(5)->setName("Clamping");

    // give the axes some labels:
    plot.xAxis->setLabel("t (s)");
    plot.yAxis->setLabel("value (u)");
    // set axes ranges, so we see all data:
    plot.xAxis->setRange(0, 10);
    plot.yAxis->setRange(   -0.5, 2);
    plot.replot();
    plot.show();




    return a.exec();

}