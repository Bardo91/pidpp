//---------------------------------------------------------------------------------------------------------------------
//  PID++
//---------------------------------------------------------------------------------------------------------------------
//  Copyright 2020 Pablo Ramon Soria (a.k.a. bardo91) pabramsor@gmail.com
//---------------------------------------------------------------------------------------------------------------------
//  Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
//  and associated documentation files (the "Software"), to deal in the Software without restriction, 
//  including without limitation the rights to use, copy, modify, merge, publish, distribute, 
//  sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is 
//  furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included in all copies or substantial 
//  portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING 
//  BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
//  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES 
//  OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
//  CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//---------------------------------------------------------------------------------------------------------------------

#include <pidpp/PID.h>
#include <algorithm>
#include <thread>
#include <chrono>
#include <iostream>
#include <cmath>
#include <Eigen/Eigen>

#include <cmath>

namespace pidpp{

    float PID::EuclideanDistance(float _a, float _b) { return _b - _a; };
    float PID::AngularDistance(float _a, float _b){
        Eigen::Vector3f v1 = {cos(_a), sin(_a), 0};
        Eigen::Vector3f v2 = {cos(_b), sin(_b), 0};
        auto v3 = v1.cross(v2);
        float magn = acos(v1.dot(v2));
        float angle = v3[2];
        return angle;
    };


    PID::PID(float _kp, float _ki, float _kd, float _minSat, float _maxSat) {
        kp_ = _kp;
        ki_ = _ki;
        kd_ = _kd;
        minSat_ = _minSat;
        maxSat_ = _maxSat;
        reference_ = 0;
        updateFn_ = std::bind(&PID::updateAWU_None, this, std::placeholders::_1, std::placeholders::_2);
	}

    void PID::setAntiWindup(AntiWindupMethod _antiWindup, std::vector<float> _params){
        switch (_antiWindup) {
        case AntiWindupMethod::None:
            updateFn_ = std::bind(&PID::updateAWU_None,             this, std::placeholders::_1, std::placeholders::_2);
            break;
        case AntiWindupMethod::Saturation:
            assert(_params.size() == 2);
            updateFn_ = std::bind(&PID::updateAWU_Saturation,       this, std::placeholders::_1, std::placeholders::_2);
            minWindup_ = _params[0];
            maxWindup_ = _params[1];
            break;
        case AntiWindupMethod::BackCalculation:
            assert(_params.size() == 1);
            updateFn_ = std::bind(&PID::updateAWU_BackCalculation,  this, std::placeholders::_1, std::placeholders::_2);
            backCalculationCte_ = _params[0];
            break;
        case AntiWindupMethod::Clamping:
            updateFn_ = std::bind(&PID::updateAWU_Clamping,         this, std::placeholders::_1, std::placeholders::_2);
            break;
        
        }
    }


    void PID::reference(float _ref, bool _reset) { 
        reference_ = _ref; 
        if(_reset){
            accumErr_ = 0; 
            lastError_ = 0; 
            lastResult_ = 0; 
            bouncingFactor_ = 0.1;
        }
    }


    float PID::update(float _val, float _incT) {
        return updateFn_(_val, _incT);
    }

    //----------------------------------------------------
    float PID::updateAWU_None(float _val, float _incT){
        float dt = _incT; // 666 input arg?

        float err = distanceFn_(_val, reference_);//reference_ - _val;
        accumErr_ += err * dt;


        // Compute PID
        float unsaturated =     kp_ * err + 
                                ki_ * accumErr_ + 
                                kd_ * (err - lastError_) / dt;
        lastError_ = err;

        // Saturate signal
        float saturated = std::min(std::max(unsaturated, minSat_), maxSat_);
        lastResult_ = saturated * bouncingFactor_;
        bouncingFactor_ *= 2.0;
        bouncingFactor_ = bouncingFactor_ > 1.0 ? 1.0 : bouncingFactor_;
        return lastResult_;
    }

    float PID::updateAWU_Saturation(float _val, float _incT){
        float dt = _incT; // 666 input arg?

        float err = distanceFn_(_val, reference_);//reference_ - _val;
        accumErr_ += err * dt;
        accumErr_ = std::min(std::max(accumErr_, minWindup_), maxWindup_);

        // Compute PID
        float unsaturated =     kp_ * err + 
                                ki_ * accumErr_ + 
                                kd_ * (err - lastError_) / dt;

        lastError_ = err;

        // Saturate signal
        float saturated = std::min(std::max(unsaturated, minSat_), maxSat_);
        lastResult_ = saturated * bouncingFactor_;
        bouncingFactor_ *= 2.0;
        bouncingFactor_ = bouncingFactor_ > 1.0 ? 1.0 : bouncingFactor_;
        return lastResult_;
    }

    float PID::updateAWU_BackCalculation(float _val, float _incT){
        float dt = _incT; // 666 input arg?

        float err = distanceFn_(_val, reference_);//reference_ - _val;
        accumErr_ += err * dt +  backCalculationCte_*lastBackCalculation_*dt;

        // Compute PID
        float unsaturated =     kp_ * err + 
                                ki_ * accumErr_ + backCalculationCte_*lastBackCalculation_ + 
                                kd_ * (err - lastError_) / dt;
        lastError_ = err;

        // Saturate signal
        float saturated = std::min(std::max(unsaturated, minSat_), maxSat_);
        lastBackCalculation_ = saturated - unsaturated;
        lastResult_ = saturated * bouncingFactor_;
        bouncingFactor_ *= 2.0;
        bouncingFactor_ = bouncingFactor_ > 1.0 ? 1.0 : bouncingFactor_;
        return lastResult_;

    }

    float PID::updateAWU_Clamping(float _val, float _incT){
        float dt = _incT; // 666 input arg?

        float err = distanceFn_(_val, reference_);//reference_ - _val;
        accumErr_ += err * dt;

        // Clamping
        bool isSaturated = (lastResult_ == minSat_ || lastResult_ == maxSat_);
        bool sameSign = std::signbit(err) == std::signbit(lastResult_);
        clampFactor_ = (isSaturated && sameSign);

        // Compute PID
        float unsaturated =     kp_ * err + 
                                clampFactor_* ki_ * accumErr_ + 
                                kd_ * (err - lastError_) / dt;
        lastError_ = err;

        // Saturate signal
        float saturated = std::min(std::max(unsaturated, minSat_), maxSat_);
        lastResult_ = saturated * bouncingFactor_;
        bouncingFactor_ *= 2.0;
        bouncingFactor_ = std::min(bouncingFactor_, 1.0);
        return lastResult_;
    }


}

