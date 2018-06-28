/*
 * dt_var.h
 *
 *
 *  Created by Joseph Tibbs on 6/27/18.
 *  
 *
 */

//=====================================
//include guard
#ifndef __DT_VAR_H_INCLUDED__
#define __DT_VAR_H_INCLUDED__

//=====================================
//included dependences
#include "globals.h"
#include "vector"


//variable timestep operator class

class dt_var
{
    public:
        dt_var::update_dt_var(double& t, );

    public:
        double tcurr, stable_thresh;
        int slow_param, net_status, myosins_status, crosslks_status, stable_checks,
         slow_down, slowed_down;
    
};

#endif