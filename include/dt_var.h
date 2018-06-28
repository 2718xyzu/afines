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
        dt_var(double final_time, int num_msgs, int chk_steps,  ostream& account_file);
        int update_dt_var(double& t, double& dt, int net_status, int myosins_status, int crosslks_status);

    public:
        double tcurr, stable_thresh, tfinal;
        int slow_param, net_status, myosins_status, crosslks_status, stable_checks,
         slow_down, slowed_down, nmsgs, check_steps;
        ostream& file_counts;
    
};

#endif