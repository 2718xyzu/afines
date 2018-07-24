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
#include "filament_ensemble.h"
#include "motor_ensemble.h"


//variable timestep operator class

class dt_var
{
    public:
        dt_var(int method, double final_time, int num_msgs, int chk_steps, double stable_threshold, double initDt, int numRetry,
            array<double, 3> thresholds);

        int update_dt_var(double& t, double& dt, int& count, ostream& file_counts);

        void clear_all(vector<double> &time_past, vector<double> &count_past,
            vector<string> &actins_past, vector<string> &links_past, vector<string> &motors_past, vector<string> &crosslks_past,
            vector<string> &thermo_past, vector<double> &stretching_energy_past, vector<double> &bending_energy_past, 
            vector<double> &potential_energy_motors_past,vector<double> &potential_energy_crosslks_past);

        void set_test(int test_param);

        void erase1(vector<double> &time_past, vector<double> &count_past,
            vector<string> &actins_past, vector<string> &links_past, vector<string> &motors_past, vector<string> &crosslks_past,
            vector<string> &thermo_past, vector<double> &stretching_energy_past, vector<double> &bending_energy_past, 
            vector<double> &potential_energy_motors_past,vector<double> &potential_energy_crosslks_past);

        void update_thresholds();

        void check_energies(filament_ensemble * network, motor_ensemble * myosins, motor_ensemble * crosslks);

    public:
        double tcurr, dtcurr, stable_thresh, tfinal, slow_amount, minDt, net_thresh, myosins_thresh, crosslks_thresh;
        int slow_param, slow_threshold, countcurr, net_status, myosins_status, crosslks_status, stable_checks,
         slow_down, slowed_down, nmsgs, check_steps, test_check, retries, backed_up, check_count, var_dt_meth;

        array<int, 10> n_s, m_s, c_s;
        array<int, 3> obj_statuses;
        array<double, 3> obj_thresholds;
        
    
};

#endif