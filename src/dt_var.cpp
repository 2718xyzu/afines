/*-------------------------------------
    dt_var is a class which contains methods
    relating to implementing a variable timestep within
    a simulation

    Created by Joseph Tibbs, 06/27/18
-------------------------------------*/


#include "globals.h"
#include "dt_var.h"
#include "vector"

dt_var::dt_var(int method, double final_time, int num_msgs, int chk_steps, double stable_threshold, double initDt, int numRetry){ 

    tfinal = final_time;
    nmsgs = num_msgs;
    check_steps = chk_steps;
    stable_checks = 0;
    stable_thresh = stable_threshold;
    test_check = 0;
    if (method == 1){ //conservative method
        slow_threshold = 4;
        slow_amount = 1.5;
    }else if(method ==2){ //aggressive method
        slow_threshold = 6;
        slow_amount = 1.95;
    }
    minDt = initDt;
    retries = numRetry;
    backed_up = 2;

}


int dt_var::update_dt_var(double& t, double& dt, int& count, int net_status, int myosins_status, 
        int crosslks_status, ostream& file_counts){

    int returned_int;
    tcurr = t;
    dtcurr = dt;
    countcurr = count;
    slow_param = 0;
    if ((net_status >= 2 || myosins_status >= 2 || crosslks_status >= 2) && backed_up<0) {
        if (slowed_down<(retries+1)){
        t -= check_steps * dt;
        count -= check_steps;
        stable_checks = floor(stable_checks/2);
        if (net_status+myosins_status+crosslks_status < slow_threshold && slow_down == 0){ 
            slow_down = 0;
        }else{
            slow_down = 1;
        }
        if (slow_down && dt>=minDt) {
            dt /= 2;
            slowed_down = 1;
            file_counts<<"\nt = " <<tcurr<<"\tSlow Down, status:\tn_s = "<<net_status<<"\tm_s = "<<myosins_status<<"\tc_s = "<<crosslks_status;
        } else if(slow_down){
            dt /= 1.5;
            slowed_down++;
            file_counts<<"\nt = " <<tcurr<<"\tSlow Down (floor), status:\tn_s = "<<net_status<<"\tm_s = "<<myosins_status<<"\tc_s = "<<crosslks_status;
        } else {
            file_counts<<"\nt = " <<tcurr<<"\tEnergy exceeded, status:\tn_s = "<<net_status<<"\tm_s = "<<myosins_status<<"\tc_s = "<<crosslks_status;
        }
        slow_down = 1;
        //something has blown up, so we go back
        cout<<"t = "<<tcurr<<" back up to "<<t<<endl;
        returned_int = 1;
        }else{
            cout<<"t = "<<tcurr<<" back up farther to ";
            returned_int = 5;
        }
    } else {
        returned_int = 2;
        backed_up--;
        if (slowed_down==1) {
            dt *= slow_amount;
            slowed_down = 0;
        }else if(slowed_down > 1){
            dt *= 1.5;
            slowed_down--;
        }
        if (net_status+myosins_status+crosslks_status < 2) {
            if (stable_checks>stable_thresh && dt<(tfinal/double(nmsgs*2)) && dt<get_upper_dt(minDt, tfinal, t)){
                dt *= 1.1;
                stable_checks = 0;
                file_counts<<"\nt = "<<tcurr<<"\tAll energies are low, speeding up now";
            }
            stable_checks += 10;
            file_counts<<"\nAll energies are low (no change)";
        } else {
            file_counts<<"\nNo change in time step, non-relaxed state";
            stable_checks += (5-(net_status+myosins_status+crosslks_status));
        }
        slow_down = 0;
                
            } 
        
    return returned_int;
}

void dt_var::clear_all(vector<double> &time_past, vector<double> &count_past,
    vector<string> &actins_past, vector<string> &links_past, vector<string> &motors_past, vector<string> &crosslks_past,
    vector<string> &thermo_past, vector<double> &stretching_energy_past, vector<double> &bending_energy_past, 
    vector<double> &potential_energy_motors_past,vector<double> &potential_energy_crosslks_past){

    time_past.clear();
    // count_diff.clear();
    count_past.clear();
    actins_past.clear();
    links_past.clear();
    motors_past.clear();
    crosslks_past.clear();
    thermo_past.clear();
    stretching_energy_past.clear();
    bending_energy_past.clear();
    potential_energy_motors_past.clear();
    potential_energy_crosslks_past.clear();


}

void dt_var::set_test(int test_param){
    test_check = test_param;
}

void dt_var::erase1(vector<double> &time_past, vector<double> &count_past,
    vector<string> &actins_past, vector<string> &links_past, vector<string> &motors_past, vector<string> &crosslks_past,
    vector<string> &thermo_past, vector<double> &stretching_energy_past, vector<double> &bending_energy_past, 
    vector<double> &potential_energy_motors_past,vector<double> &potential_energy_crosslks_past){

    time_past.erase(time_past.begin());
    count_past.erase(count_past.begin());
    actins_past.erase(actins_past.begin());
    links_past.erase(links_past.begin());
    motors_past.erase(motors_past.begin());
    crosslks_past.erase(crosslks_past.begin());
    thermo_past.erase(thermo_past.begin());
    stretching_energy_past.erase(stretching_energy_past.begin());
    bending_energy_past.erase(bending_energy_past.begin());
    potential_energy_motors_past.erase(potential_energy_motors_past.begin());
    potential_energy_crosslks_past.erase(potential_energy_crosslks_past.begin());


}

void dt_var::update_thresholds(array<int, 3> statuses){
    count = (count+1)%10;
    n_s[count] = statuses[0];
    m_s[count] = statuses[1];
    c_s[count] = statuses[2]; 

    int max_n_s = *max_element(n_s, n_s+10);
    int max_m_s = *max_element(m_s, m_s+10);
    int max_c_s = *max_element(c_s, c_s+10);

}