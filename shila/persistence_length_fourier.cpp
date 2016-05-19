#include "actin_myosin_flexible.cpp"


#define xrange 25.0
#define yrange 25.0
#define xgrid 100.0
#define ygrid 100.0

#define tinit 0.0
#define tfinal 100
// #define dt 0.001 -- defined previously
#define print_dt 1000

int main(int argc, char* argv[]){
    
    int seed=time(NULL);
    srand(seed);
    
    std::string output_file="output.txt";
    std::string actin_output="actin_final.txt";
    std::string myosin_output="myosin_final.txt";
    std::string persistence_length_output = "fourier_modes.txt";
    
    // CONTROLS :
    double motor_length=0.5;
    double motor_density=0;
    double motor_stiffness=50.0; //pN / um
    double vmotor=1.0;
    double m_kon=90.0;          
    double m_kend=5.0;
    double m_koff=1.0; 
    double viscosity=0.5;
    

    // VARIABLES :
    double npolymer = 1, nmonomer = 100;
    double link_stiffness = motor_stiffness/2;
    
    std::string link_color = "1"; //"blue";
    
    if (argc>1) {
        nmonomer = std::atof(argv[1]);
        link_stiffness = std::atof(argv[2]);
    }
    
    // DERIVED QUANTITIES :
    double actin_length=10/nmonomer; //length of a monomer
    double actin_density= npolymer*nmonomer/(xrange*yrange);//0.65;
    double link_length = actin_length/10; 
    
    int n_modes = 20;
    std::map<int, std::vector<double> > fm;

    std::ofstream a_final, m_final, p_final;
    a_final.open(actin_output.c_str());
    m_final.open(myosin_output.c_str());
    p_final.open(persistence_length_output.c_str());
    std::ofstream o_file;
    o_file.open(output_file.c_str());
    o_file << " FILE: " << output_file <<"\n"<< "Actin Density: " << actin_density  << ", Actin Mean Length: " << actin_length << "\n";
    o_file << " Motor Density: " << motor_density << ", Motor Rest Length: " << motor_length << ", Motor Stiffness: " << motor_stiffness<<"\n";
    o_file << ", Motor unloaded speed: " << vmotor << ", Motor binding rate: " << m_kon <<"\n"<<"Motor unbinding rate: " << m_koff << ", Motor end detachment rate: " << m_kend<<", Viscosity: "<<viscosity<<"\n";
    o_file << " Link Rest Length: "<< link_length <<", Link Stiffness: " << link_stiffness <<"\n";
    o_file << " Simulation time: " << tfinal - tinit << ", dt: " << dt <<", Number of time steps between output files: "<< print_dt<<"\n";
    o_file.close();
    
    
    int count=0;
    int timesteps=ceil(tfinal/dt);
	char afile[timesteps+2];
	char mfile[timesteps+2];
    char tfile[timesteps+2];
    
    
    std::cout<<"\nCreating actin network..";
	actin_ensemble net=actin_ensemble(actin_density,xrange,yrange,xgrid,ygrid,actin_length,viscosity,nmonomer,link_length);
    std::cout<<"\nAdding motors..";
    motor_ensemble myosins=motor_ensemble(motor_density,xrange,yrange,motor_length,&net,vmotor,motor_stiffness,m_kon,m_koff,m_kend,actin_length,viscosity);
    std::cout<<"\nAdding links to connect actin filament monomers...";
    net.connect_polymers( &myosins, link_length, link_stiffness, link_color );

    double t=tinit;
    std::cout<<"\nUpdating motors, filaments and crosslinks in the network..";
    
    while (t<=tfinal) {
        //print time count
		if (count%1000==0) {
			std::cout<<"\nTime counts: "<<count;
		}
        
        //update network
        net.update();
        net.quad_update();
        
        //print to file
	    
        if (count%print_dt==0) {
		
            sprintf (afile, "afile%d.txt", count/print_dt);
			sprintf (mfile, "mfile%d.txt", count/print_dt);
            sprintf (tfile, "tfile%d.txt", count/print_dt);
			std::ofstream file_a, file_m, file_t;
			file_a.open(afile);
			file_m.open(mfile);
            file_t.open(tfile);
			net.write(file_a);
			myosins.motor_write(file_m);
            myosins.motor_tension(file_t);
			file_a.close();
			file_m.close();
            file_t.close();
            
		}
        
        for (int n = 1; n <= n_modes; n++){
            //assuming only one polymer
            fm[n].push_back(net.get_fourier_mode(n, 0));
        
        }
        
        //myosins.reshape();
        myosins.motor_walk();
        //
        t+=dt;
		count++;
    }
    net.write(a_final);
    myosins.motor_write(m_final);
    
    // write the fourier mode file 
    for (int n = 1; n<= n_modes; n++){
        p_final<< ( 1/(n*n*pi*pi) )<<"\t"<<mode_var(fm[n],0)<<"\n"; 
    }
    
    a_final.close();
    m_final.close();
    p_final.close(); 

    std::cout<<"\nTime counts: "<<count;
	std::cout<<"\nExecuted";
	std::cout<<"\n Done\n";
    return 0;
}
