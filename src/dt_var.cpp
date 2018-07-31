/*-------------------------------------
    dt_var is a class which contains methods
    relating to implementing a variable timestep within
    a simulation

    Created by Joseph Tibbs, 06/27/18
-------------------------------------*/


#include "globals.h"
#include "dt_var.h"
#include "vector"
// #include "filament_ensemble.h"
// #include "motor_ensemble.h"

dt_var::dt_var(int method, double final_time, int num_msgs, int chk_steps, double stable_threshold, double initDt, int numRetry,
        array<double, 3> thresholds){ 

    tfinal = final_time;
    nmsgs = num_msgs;
    check_steps = chk_steps;
    stable_checks = 0;
    stable_thresh = stable_threshold;
    test_check = 0;
    if (method == 1){ //conservative method
        slow_threshold = 4;
        slow_amount = 1.5; //slow_amount is the reciprocal of the net proportional change after a single reversal (assuming time step changes)
    }else if(method ==2){ //aggressive method
        slow_threshold = 6;
        slow_amount = 1.95;
    }
    minDt = initDt;
    retries = numRetry;
    backed_up = 2;
    obj_thresholds = thresholds;
    check_count = 0;
    var_dt_meth = method;

}


int dt_var::update_dt_var(double& t, double& dt, int& count, ostream& file_counts){

    int returned_int;
    tcurr = t;
    dtcurr = dt;
    countcurr = count;
    slow_param = 0;
    if ((obj_statuses[0]>=2 || obj_statuses[1]>=2  || obj_statuses[2]>=2) && !(backed_up==1 && slowed_down>=(retries+1))) {
        if (slowed_down<(retries+1)){
        t -= check_steps * dt;
        count -= check_steps;
        stable_checks = floor(stable_checks/2);
        if (obj_statuses[0]+obj_statuses[1]+obj_statuses[2] < slow_threshold && slow_down == 0){ 
            slow_down = 0;
        }else{
            slow_down = 1;
        }
        if (slow_down && dt>=minDt) {
            dt /= 2;
            slowed_down = 1;
            file_counts<<"\nt = " <<tcurr<<"\tSlow Down, status:\tn_s = "<<obj_statuses[0]<<"\tm_s = "<<obj_statuses[1]<<"\tc_s = "<<obj_statuses[2];
        } else if(slow_down){
            dt /= 1.5;
            slowed_down++;
            file_counts<<"\nt = " <<tcurr<<"\tSlow Down (floor), status:\tn_s = "<<obj_statuses[0]<<"\tm_s = "<<obj_statuses[1]<<"\tc_s = "<<obj_statuses[2];
        } else {
            file_counts<<"\nt = " <<tcurr<<"\tEnergy exceeded, status:\tn_s = "<<obj_statuses[0]<<"\tm_s = "<<obj_statuses[1]<<"\tc_s = "<<obj_statuses[2];
            slow_down = 1;
        }
        
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
        if (obj_statuses[0]+obj_statuses[1]+obj_statuses[2] < 2) {
            if (stable_checks>stable_thresh && dt<(tfinal/double(nmsgs*2)) && dt<get_upper_dt(minDt, tfinal, t)){
                dt *= 1.1;
                stable_checks = 0;
                file_counts<<"\nt = "<<tcurr<<"\tAll energies are low, speeding up now";
            }
            stable_checks += 10;
            file_counts<<"\nAll energies are low (no change)";
        } else {
            file_counts<<"\nNo change in time step, non-relaxed state";
            stable_checks += (5-(obj_statuses[0]+obj_statuses[1]+obj_statuses[2]));
        }
        slow_down = 0;
                
            } 
        
    return returned_int;
}

// int dt_var::update_objects(){

//             if (returned_int == 1){
//                 file_counts<<"\n dt is now "<<dt<<endl;
            
//                 delete net;
//                 net = new filament_ensemble(*net2);
                
//                 delete myosins;
//                 myosins = new motor_ensemble(*myosins2);

//                 delete crosslks;
//                 crosslks = new motor_ensemble(*crosslks2);

//                 myosins->set_fil_ens(net);
//                 crosslks->set_fil_ens(net);
//                 net->set_dt(dt);
//                 myosins->set_dt(dt);
//                 crosslks->set_dt(dt);

//                 for (unsigned int i = 0; i<time_past.size(); i++){
//                     print_times.push_back(time_past.back());
//                     time_past.pop_back(); 
//                 }

//             var_dt.clear_all(time_past, count_past, actins_past, links_past, motors_past, crosslks_past,
//                 thermo_past, stretching_energy_past, bending_energy_past, potential_energy_motors_past, potential_energy_crosslks_past);

//             }else if(returned_int == 2){

//                 net->set_dt(dt);
//                 myosins->set_dt(dt);
//                 crosslks->set_dt(dt);

//                 delete backupNet1;
//                 delete backupMyosins1;
//                 delete backupCrosslks1;
//                 backupNet1 = net2;
//                 backupMyosins1 = myosins2;
//                 backupCrosslks1 = crosslks2;

//                 net2 = new filament_ensemble(*net);
//                 myosins2 = new motor_ensemble(*myosins);
//                 crosslks2 = new motor_ensemble(*crosslks);
//                 myosins2->set_fil_ens(net2);
//                 crosslks2->set_fil_ens(net2);

//                 int flush = 0;
//                 int lastPrint = (int) time_past.size();


//                 for (unsigned int i = 0; i<time_past.size(); i++){
//                     if(time_past[i]>backupNet1->t) continue;
//                     lastPrint = (int) i;
//                     cout<<"Wrote out t = "<<to_string(time_past[i])<<endl;
//                     if (time_past[i]>tinit) time_str ="\n";
//                     time_str += "t = "+to_string(time_past[i]);

//                     file_a << time_str<<"\tN = "<<to_string(net->get_nactins());
//                     file_a << actins_past[i];

//                     file_l << time_str<<"\tN = "<<to_string(net->get_nlinks());
//                     file_l << links_past[i];

//                     file_am << time_str<<"\tN = "<<to_string(myosins->get_nmotors());
//                     file_am << motors_past[i];

//                     file_pm << time_str<<"\tN = "<<to_string(crosslks->get_nmotors());
//                     file_pm << crosslks_past[i];

//                     file_th << time_str<<"\tN = "<<to_string(net->get_nfilaments());
//                     file_th << thermo_past[i];

//                     file_pe << to_string(stretching_energy_past[i]) + "\t" + to_string(bending_energy_past[i]) + "\t" +
//                         to_string(potential_energy_motors_past[i]) + "\t" + to_string(potential_energy_crosslks_past[i]) << endl;
                    
//                     flush = 1;
//                 }
//                 if (lastPrint<(((int) time_past.size())-1)){
//                     cout<<"erasing first elements"<<endl;
//                     for (int i = 0; i<lastPrint-1; i++){
//                         var_dt.erase1(time_past, count_past, actins_past, links_past, motors_past, crosslks_past,
//                         thermo_past, stretching_energy_past, bending_energy_past, potential_energy_motors_past, potential_energy_crosslks_past);
//                     }
//                 }else if(lastPrint == (((int) time_past.size())-1)) {
//                         var_dt.clear_all(time_past, count_past, actins_past, links_past, motors_past, crosslks_past,
//                             thermo_past, stretching_energy_past, bending_energy_past, potential_energy_motors_past, potential_energy_crosslks_past);
//                 }
//                 if (flush) {
                    
//                     file_a<<std::flush;
//                     file_l<<std::flush;
//                     file_am<<std::flush;
//                     file_pm<<std::flush;
//                     file_th<<std::flush;
//                     file_pe<<std::flush;
//                     file_time<<std::flush;
//                     file_counts<<std::flush;
//                     flush = 0;
//                 }

//             } else if(returned_int == 5){
//                 cout<<(backupNet1->t)<<endl;
//                 file_counts<<"\n dt is now "<<dt<<endl;
            
//                 delete net;
//                 net = new filament_ensemble(*backupNet1);
//                 delete net2;
//                 net2 = new filament_ensemble(*backupNet1);
                
//                 delete myosins;
//                 myosins = new motor_ensemble(*backupMyosins1);
//                 delete myosins2;
//                 myosins2 = new motor_ensemble(*backupMyosins1);

//                 delete crosslks;
//                 crosslks = new motor_ensemble(*backupCrosslks1);
//                 delete crosslks2;
//                 crosslks2 = new motor_ensemble(*backupCrosslks1);


//                 myosins->set_fil_ens(net);
//                 crosslks->set_fil_ens(net);
//                 net->set_dt(dt);
//                 myosins->set_dt(dt);
//                 crosslks->set_dt(dt);               
//                 myosins2->set_fil_ens(net2);
//                 crosslks2->set_fil_ens(net2);
//                 net2->set_dt(dt);
//                 myosins2->set_dt(dt);
//                 crosslks2->set_dt(dt);
                
//                 t = net->t;

//                 var_dt.backed_up = 2;

//                 for (unsigned int i = 0; i<time_past.size(); i++){
//                     print_times.push_back(time_past.back());
//                     time_past.pop_back(); 
//                 }

//             var_dt.clear_all(time_past, count_past, actins_past, links_past, motors_past, crosslks_past,
//                 thermo_past, stretching_energy_past, bending_energy_past, potential_energy_motors_past, potential_energy_crosslks_past);

//             }


//             file_time<<t<<"\t"<<dt<<endl;
//             obj_statuses = {0,0,0};
// }


void dt_var::clear_all(vector<double> &time_past, vector<double> &count_past, //in the interests of keeping network.cpp cleaner
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


void dt_var::check_energies(filament_ensemble * network, motor_ensemble * myosins, motor_ensemble * crosslks){
        obj_statuses[0] = max(obj_statuses[0], network->check_link_energies(var_dt_meth, obj_thresholds[0]));
        obj_statuses[1] = max(obj_statuses[1], myosins->check_energies(0, obj_thresholds[1]));
        obj_statuses[2] = max(obj_statuses[2],crosslks->check_energies(0, obj_thresholds[2]));
}

void dt_var::check_blowups(filament_ensemble * network, motor_ensemble * myosins, motor_ensemble * crosslks, double timestamp){
        obj_statuses[0] = max(obj_statuses[0], network->check_blowups(timestamp));
        obj_statuses[1] = max(obj_statuses[1], myosins->check_blowups(timestamp));
        obj_statuses[2] = max(obj_statuses[2],crosslks->check_blowups(timestamp));
}


string dt_var::update_thresholds(){
    string out = "";
    check_count = (check_count+1)%10;
    n_s[check_count] = obj_statuses[0];
    m_s[check_count] = obj_statuses[1];
    c_s[check_count] = obj_statuses[2]; 

    if(check_count == 0){

    int max_n_s = *max_element(n_s.begin(), n_s.begin()+10);
    int max_m_s = *max_element(m_s.begin(), m_s.begin()+10);
    int max_c_s = *max_element(c_s.begin(), c_s.begin()+10);

    int sum_n_s = accumulate(n_s.begin(), n_s.begin()+10, 0);
    int sum_m_s = accumulate(m_s.begin(), m_s.begin()+10, 0);
    int sum_c_s = accumulate(c_s.begin(), c_s.begin()+10, 0);

    if(sum_n_s < 4){
        obj_thresholds[0] *= .95;
    }else if(max_n_s >= 3){
        int increase_threshold = 0;
        for (int i = 0; i<10; i++){
            if(n_s[i]>=3) increase_threshold++;
        }
        if (increase_threshold>1) obj_thresholds[0] *= 1.05;
    }
    if(sum_m_s < 4){
        obj_thresholds[1] *= .95;
    }else if(max_m_s >= 3){
        int increase_threshold = 0;
        for (int i = 0; i<10; i++){
            if(m_s[i]>=3) increase_threshold++;
        }
        if (increase_threshold>1) obj_thresholds[1] *= 1.05;
    }

    if(sum_c_s < 4){
        obj_thresholds[2] *= .95;
    }else if(max_c_s >= 3){
        int increase_threshold = 0;
        for (int i = 0; i<10; i++){
            if(c_s[i]>=3) increase_threshold++;
        }
        if (increase_threshold>1) obj_thresholds[2] *= 1.05;
    }
    out = "\t"+to_string(obj_thresholds[0])+"\t"+to_string(obj_thresholds[1])+"\t"+to_string(obj_thresholds[2]);
    }

    return out;
}