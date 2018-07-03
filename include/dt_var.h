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
        dt_var(double final_time, int num_msgs, int chk_steps, double stable_threshold);

        int update_dt_var(double& t, double& dt, int& count, int net_status, int myosins_status, int crosslks_status, 
            ostream& account_file);

        void clear_all(vector<double> &time_past, vector<double> &count_past,
            vector<string> &actins_past, vector<string> &links_past, vector<string> &motors_past, vector<string> &crosslks_past,
            vector<string> &thermo_past, vector<double> &stretching_energy_past, vector<double> &bending_energy_past, 
            vector<double> &potential_energy_motors_past,vector<double> &potential_energy_crosslks_past);

        void set_test(int test_param);

    public:
        double tcurr, dtcurr, stable_thresh, tfinal, slow_amount;
        int slow_param, slow_threshold, countcurr, net_status, myosins_status, crosslks_status, stable_checks,
         slow_down, slowed_down, nmsgs, check_steps, test_check;
        
    
};

#endif