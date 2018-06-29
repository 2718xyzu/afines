/*-------------------------------------
    dt_var is a class which contains methods
    relating to implementing a variable timestep within
    a simulation

    Created by Joseph Tibbs, 06/27/18
-------------------------------------*/


#include "globals.h"
#include "dt_var.h"
#include "vector"

dt_var::dt_var(double final_time, int num_msgs, int chk_steps, double stable_threshold){ 

    tfinal = final_time;
    nmsgs = num_msgs;
    check_steps = chk_steps;
    stable_checks = 0;
    stable_thresh = stable_threshold;
    test_check = 0;
}


int dt_var::update_dt_var(double& t, double& dt, int& count, int net_status, int myosins_status, 
    int crosslks_status, ostream& file_counts){
    int returned_int;
    tcurr = t;
    dtcurr = dt;
    countcurr = count;
    slow_param = 0;
    if ((net_status == 2 || myosins_status == 2 || crosslks_status == 2) && slowed_down<2 ) {
        t -= check_steps * dt;
        count -= check_steps;
        stable_checks = 0;
        if (net_status+myosins_status+crosslks_status<3 && slow_down == 0){ 
            slow_down = 0;
        }else{
            slow_down = 1;
        }
        if (slow_down && dt>6E-5) {
            dt /= 2;
            slowed_down = 1;
            file_counts<<"\nt = " <<tcurr<<"\tSlow Down, status:\tn_s = "<<net_status<<"\tm_s = "<<myosins_status<<"\tc_s = "<<crosslks_status;
        } else if(slow_down){
            dt /= 2;
            slowed_down = 2;
            file_counts<<"\nt = " <<tcurr<<"\tSlow Down (floor), status:\tn_s = "<<net_status<<"\tm_s = "<<myosins_status<<"\tc_s = "<<crosslks_status;
        } else {
            file_counts<<"\nt = " <<tcurr<<"\tEnergy exceeded, status:\tn_s = "<<net_status<<"\tm_s = "<<myosins_status<<"\tc_s = "<<crosslks_status;
        }
        slow_down = 1;
        //something has blown up, so we back up
        cout<<"t = "<<tcurr<<" back up to "<<t<<endl;
        returned_int = 1;
    } else {
        returned_int = 2;
        if (slowed_down==1) {
            dt *= 1.5;
            slowed_down = 0;
        }else if(slowed_down==2){
            dt *= 2;
            slowed_down = 0;
        }
        if (net_status == 0 && myosins_status == 0 && crosslks_status == 0) {
            if (stable_checks>stable_thresh && dt<(tfinal/double(nmsgs*2)) && dt<get_upper_dt(tfinal, t)){
                dt *= 1.1;
                stable_checks = 0;
                file_counts<<"\nt = "<<tcurr<<"\tAll energies are low, speeding up now";
            }
            stable_checks += 5;
            file_counts<<"\nAll energies are low (no change)";
        } else {
            file_counts<<"\nNo change in time step, non-relaxed state";
            stable_checks += 1;
        }
        slow_down = 0;
                
            } 
    if (abs(test_check)==1){
        dt = dtcurr;
        count = countcurr;
        t = tcurr;
        returned_int = 2;
        if (count%140000 == 0 && test_check == 1){
            t -= check_steps * dt;
            count -= check_steps;
            test_check = -1;
        }
        else test_check = 1;



    }else if(abs(test_check) ==2){
        dt = dtcurr;
        count = countcurr;
        t = tcurr;
        returned_int = 2;
        if (count%140000 == 0 && test_check == 2){
            t -= check_steps * dt;
            count -= check_steps;
            dt /= 2;
            slowed_down = 2;
            test_check = -2;
        }
        else test_check = 2;

    }else if(test_check == 3){
        dt = dtcurr;

    }else if(test_check == 4){
        returned_int = 2;
        dt = dtcurr;
        count = countcurr;
        t = tcurr;
        if (count%140000 == 0){
            returned_int = 10;
        }
    }else if(abs(test_check)==5){
        dt = dtcurr;
        count = countcurr;
        t = tcurr;
        returned_int = 2;
        if (count%140000 == 0 && test_check == 5){
            t -= check_steps * dt;
            count -= check_steps;
            test_check = -5;
            returned_int = 1;
        }
        else test_check = 5;



    return returned_int;
}

// void dt_var::say_my_name(){
//     returned_int = var_dt.update_dt_var(&t, &dt, net_status, myosins_status, crosslks_status);

//     if (returned_int == 1){
//         delete net;
//         net = new filament_ensemble(stored_actin_pos_vec, {xrange, yrange}, {xgrid, ygrid}, dt,
//             temperature, viscosity, link_length,
//             link_stretching_stiffness, fene_pct, link_bending_stiffness,
//             link_fracture_force, bnd_cnd);
//         delete myosins;
//         myosins = new motor_ensemble( stored_a_motor_pos_vec, {xrange, yrange}, dt, temperature,
//             a_motor_length, net, a_motor_v, a_motor_stiffness, fene_pct, a_m_kon, a_m_koff,
//             a_m_kend, a_m_stall, a_m_cut, viscosity, a_motor_lcatch, a_m_fracture_force, bnd_cnd, light_param);
//         delete crosslks;
//         crosslks = new motor_ensemble( stored_p_motor_pos_vec, {xrange, yrange}, dt, temperature,
//             p_motor_length, net, p_motor_v, p_motor_stiffness, fene_pct, p_m_kon, p_m_koff,
//             p_m_kend, p_m_stall, p_m_cut, viscosity, p_motor_lcatch, p_m_fracture_force, bnd_cnd, {0,0});

//         int temp_size = time_past.size();
//         for (int i = 0; i<temp_size; i++){
//             print_times.push_back(time_past.back());
//             time_past.pop_back(); 
//         }

//     clear_all(&time_past, &count_diff, &count_past, &actins_past, &links_past, &motors_past, &crosslks_past,
//        &thermo_past, &stretching_energy_past, &bending_energy_past, &potential_energy_motors_past, &potential_energy_crosslks_past);

//     }else if(returned_int == 2){
//         net->set_dt(dt);
//         myosins->set_dt(dt);
//         crosslks->set_dt(dt);

//         stored_actin_pos_vec = net->get_vecvec();
//         stored_a_motor_pos_vec = myosins->get_vecvec();
//         stored_p_motor_pos_vec = crosslks->get_vecvec();
//         int flush = 0;
//         for (unsigned int i = 0; i<time_past.size(); i++){

//             cout<<"Wrote out t = "<<to_string(time_past[i])<<endl;
//             if (time_past[i]>tinit) time_str ="\n";
//             time_str += "t = "+to_string(time_past[i]);

//             file_a << time_str<<"\tN = "<<to_string(net->get_nactins());
//             file_a << actins_past[i];

//             file_l << time_str<<"\tN = "<<to_string(net->get_nlinks());
//             file_l << links_past[i];

//             file_am << time_str<<"\tN = "<<to_string(myosins->get_nmotors());
//             file_am << motors_past[i];

//             file_pm << time_str<<"\tN = "<<to_string(crosslks->get_nmotors());
//             file_pm << crosslks_past[i];

//             file_th << time_str<<"\tN = "<<to_string(net->get_nfilaments());
//             file_th << thermo_past[i];

//             file_pe << to_string(stretching_energy_past[i]) + "\t" + to_string(bending_energy_past[i]) + "\t" +
//                 to_string(potential_energy_motors_past[i]) + "\t" + to_string(potential_energy_crosslks_past[i]) << endl;
            
//             flush = 1;

//         }
//         if (flush) {
//             var_dt.clear_all(&time_past, &count_diff, &count_past, &actins_past, &links_past, &motors_past, &crosslks_past,
//                 &thermo_past, &stretching_energy_past, &bending_energy_past, &potential_energy_motors_past, &potential_energy_crosslks_past);

//             file_a<<std::flush;
//             file_l<<std::flush;
//             file_am<<std::flush;
//             file_pm<<std::flush;
//             file_th<<std::flush;
//             file_pe<<std::flush;
//             file_time<<std::flush;
//             file_counts<<std::flush;
//             flush = 0;
//         }

//     }
//     file_time<<t<<"\t"<<dt<<endl;
// }

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

