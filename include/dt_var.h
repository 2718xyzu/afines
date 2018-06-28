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
        dt_var(double final_time, int num_msgs, int chk_steps);

        int update_dt_var(double& t, double& dt, int& count, int net_status, int myosins_status, int crosslks_status, 
            ostream& account_file);

        void clear_all(vector<string> &time_past, vector<string> &count_diff, vector<string> &count_past,
            vector<string> &actins_past, vector<string> &links_past, vector<string> &motors_past, vector<string> &crosslks_past,
            vector<string> &thermo_past, vector<double> &stretching_energy_past, vector<double> &bending_energy_past, 
            vector<double> &potential_energy_motors_past,vector<double> &potential_energy_crosslks_past);

    public:
        double tcurr, stable_thresh, tfinal;
        int slow_param, net_status, myosins_status, crosslks_status, stable_checks,
         slow_down, slowed_down, nmsgs, check_steps;
        
    
};

#endif