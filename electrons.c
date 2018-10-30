//Program to calculate the electron positions if the electrons are released in a circular disk
//Invocation: ./<execname>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gmp.h>
#include "electrons.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

//Global variables
electron* state;
electron* newstate;

//System constants;
mpf_t  force_multiplier; 
mpf_t  field_multiplier;
mpf_t  electric_permittivity;
mpf_t  electron_charge;
double step_size;
double d = 1E-9;
char isVideo = 0;

//Function declarations
void system_step(unsigned long int noelectrons, double radius);
void calc_pot(	 unsigned long int noelectrons, double plate_distance, unsigned long int nostepspotential,  double x, double y, 			 mpf_t* potential);
void ef_z(		 unsigned long int noelectrons, double x, 			   double y, 							double z, double plate_distance, mpf_t* ef_z);
void save_state( const char* filename, 			unsigned long int noelectrons);

//Main function
int main(int argc, char** argv){
	//setup variables
	double 				radius;
	unsigned long int 	noelectrons, noiterations;
	char 				filename[101];
	char 				buffer[101];
	char 				cmd[101];
	char 				outfile[101];
	unsigned long int 	h, i, j, k;
	electron 			el;
	double 				angle, angle_per_electron;
	
	unsigned long int 	nostepspotential;
	double 				plate_distance; 
	double 				total_charge;
	double 				capacity;
	
	//SETUP GMP
	mpf_set_default_prec(128);
	
	//SETUP constants
	mpf_init_set_str(electric_permittivity, "8.854187817@-12", 10);
	mpf_init_set_str(electron_charge, "-1.60217657@-19", 10);
	
	mpf_t denum;
	mpf_init_set_ui(denum, 4);
	
	mpf_t mpfpi;
	mpf_init_set_str(mpfpi, "3.14159265359", 10);
	
	mpf_mul(denum, denum, mpfpi);
	mpf_mul(denum, denum, electric_permittivity);

	mpf_t negone;
	mpf_init_set_d(negone, -1.0);
				
	mpf_init(force_multiplier);
	mpf_init(field_multiplier);
	mpf_div(force_multiplier, negone, denum);
	
	mpf_clear(mpfpi);
	
	mpf_mul(force_multiplier, force_multiplier, electron_charge);
	mpf_set(field_multiplier, force_multiplier);
	
	mpf_mul(field_multiplier, field_multiplier, negone);
	mpf_mul(force_multiplier, force_multiplier, electron_charge);
	
	mpf_clear(negone);
				
	//Gather configuration information	
	if(argc < 9){
		printf("You can also call this executable: \n ./<exec> radius platedistance noelectrons noiterations nostepspotential <v/i> filebuffer outfile\n");
		//Prompt the user for the configuration
		printf("Enter the plate radius: ");
		scanf("%lf", &radius);
	
		printf("Enter the plate distance: ");
		scanf("%lf", &plate_distance);
	
		printf("Enter the number of electrons simulated: ");
		scanf("%li", &noelectrons);
	
		printf("Enter the number of iterations for the simulation: ");
		scanf("%li", &noiterations);
	
		printf("Enter the number of steps for calculating the potential: ");
		scanf("%li", &nostepspotential);
		
		char input;
		printf("Do you want to generate a video or an image? [y/n]: ");
		scanf("%c", &input);
		if(input == 'v') isVideo = 1;
	
		printf("Enter the file name to save to: ");
		scanf("%s", filename);
		
		printf("Enter the image file name: ");
		scanf("%s", outfile);
	}else{
		radius 				= atof( argv[1]);
		plate_distance		= atof( argv[2]);
		noelectrons 		= atoll(argv[3]);
		noiterations 		= atoll(argv[4]);
		nostepspotential 	= atoll(argv[5]);
		char videoimage		= argv[6][0];
		if(videoimage == 'v'){
			isVideo = 1;
		}
		strcpy(filename, argv[7]);
		strcpy(outfile, argv[8]);
	}
	
	//setup step size
	step_size = radius / noiterations;
	
	//Setup the buffers!
	state 		= (electron*) malloc(sizeof(electron) * noelectrons);
	newstate 	= (electron*) malloc(sizeof(electron) * noelectrons);
	
	//Move the electrons to their initial positions.
	angle_per_electron = (2 * M_PI) / noelectrons;
	for(h = 0; h < noelectrons; h++){
		angle = angle_per_electron * h;
		state[h].x = radius * d * cos(angle);
		state[h].y = radius * d * sin(angle);
	}
	
	char databuffer[201];
	char vidbuffer[201];
	char folderbuffer[101];
	
	if(isVideo){
		sprintf(folderbuffer, "VID_%s", outfile);
		mkdir(folderbuffer, 0777);
		sprintf(folderbuffer, "VID_%s/DATA", outfile);
		mkdir(folderbuffer, 0777);
		sprintf(folderbuffer, "VID_%s/IMG", outfile);
		mkdir(folderbuffer, 0777);
		
		sprintf(databuffer, "VID_%s/DATA/0.txt", outfile);
		save_state(databuffer, noelectrons);

		sprintf(vidbuffer, "VID_%s/IMG/000000.png", outfile);
		sprintf(cmd, "python plot.py %lf %s %s %i", radius, databuffer, vidbuffer, 0);
		system(cmd);
	}
	
	//Run the simulation
	for(i = 0; i < noiterations; i++){
		printf("\rIteration: %li", i+1); fflush(stdout);
		system_step(noelectrons, radius);
	
		//Copy newstate to state.
		for(j = 0; j < noelectrons; j++){
			state[j] = newstate[j];
		}
		
		if(isVideo){
			sprintf(databuffer, "VID_%s/DATA/%i.txt", outfile, i+1);
			save_state(databuffer, noelectrons);
		
			sprintf(vidbuffer, "VID_%s/IMG/%06i.png", outfile, i+1);
			sprintf(cmd, "python plot.py %lf %s %s %i", radius, databuffer, vidbuffer, i+1);
			system(cmd);
		}
	}
	printf("\n");
	
	//We use the results to calculate the potential.
	mpf_t potential;
	mpf_init(potential);
	calc_pot(noelectrons, plate_distance, nostepspotential, 0, 0, &potential);
	gmp_printf("Found potential: %.*Ff Volts\n", 20, potential);
	mpf_t mp_tc;
	mpf_init(mp_tc);
	mpf_mul_ui(mp_tc, electron_charge, noelectrons);
	mpf_t cap;
	mpf_init(cap);
	mpf_div(cap, mp_tc, potential);
	
	//Convert to PicoFarads
	mpf_mul_ui(cap, cap, 10E12);
	
	gmp_printf("Found capacity: %.*Ff pF\n", 20, cap);
	
	mpf_clear(potential);
	mpf_clear(cap);
	mpf_clear(mp_tc);

	if(!isVideo){
		save_state(filename, noelectrons);
	
		printf("Plotting results...\n");
		//We call on the plotting script
		sprintf(cmd, "python plot.py %lf %s %s", radius, filename, outfile);
		system(cmd);
	}else{
		//we call the rendering script
		sprintf(cmd, "ffmpeg -i /VID_%s/IMG/\%06i.png %s.mpg", outfile, outfile);
		system(cmd); 
	}
	
	free(state);
	free(newstate);
	
	return 0;
}

void save_state(const char* filename, unsigned long int noelectrons){
	//We need to output state to filename
	FILE* output = fopen(filename, "w+");
	char buffer[101];
	unsigned long int k;
	electron el;
	
	for(k = 0; k < noelectrons; k++){
		el = state[k];
		sprintf(buffer, "%lf %lf\n", el.x, el.y); 
		fwrite(buffer, sizeof(char), strlen(buffer), output);
	}
	fclose(output);
}

void system_step(unsigned long int noelectrons, double radius){
	//Setup the required buffers
	double x1, y1, x2, y2, x3, y3;
	double xdist, ydist, dist2, dist;
	electron el1, el2, el3;
	unsigned long int i, j, k;
	double t1, t2;
		
	mpf_t force_mpz;
	mpf_t mpfxdist, mpfydist, mpfdist2, mpfdist, mpfxdist2, mpfydist2;
	mpf_t force_x, force_y, total_force_x, total_force_y, force_len, force_len2, force;
	mpf_t distratx, distraty;
	mpf_t force_x2, force_y2;
	mpf_t step_x, step_y, mpfstep_size;
	mpf_t xforcerat, yforcerat;
		
	mpf_init_set_str(force_mpz, "1", 10);
	mpf_init(mpfxdist);
	mpf_init(mpfxdist2);
	mpf_init(mpfydist);
	mpf_init(mpfydist2);
	mpf_init(mpfdist2);
	mpf_init(mpfdist);
	mpf_inits(force_x, force_y, total_force_x, total_force_y, force_len, force_len2, force, NULL);
	mpf_inits(distratx, distraty, NULL);
	
	mpf_inits(force_x2, force_y2, NULL);
	mpf_inits(step_x, step_y, mpfstep_size, NULL);
	mpf_inits(xforcerat, yforcerat, NULL);

	mpf_set_d(mpfstep_size, step_size);
	
	//We calculate nett force for each electron, then place that information in newstate.
	for(i = 0; i < noelectrons; i++){
		el1 = state[i];
		x1 = el1.x;
		y1 = el1.y;
		
		mpf_set_d(total_force_x, 0.0);
		mpf_set_d(total_force_y, 0.0);
		
		for(j = 0; j < noelectrons; j++){
			if( i != j){
				//We calculate the force from this electron onto our current electron
				el2 = state[j];
				x2 = el2.x;
				y2 = el2.y;
				
				xdist = x2 - x1;
				ydist = y2 - y1;
				
				mpf_set_d(mpfxdist, xdist);
				mpf_set_d(mpfydist, ydist);
				
				mpf_mul(mpfxdist2, mpfxdist, mpfxdist);
				mpf_mul(mpfydist2, mpfydist, mpfydist);
								
				mpf_add(mpfdist2, mpfxdist2, mpfydist2);
				mpf_sqrt(mpfdist, mpfdist2);
				
				//gmp_printf("o: %.*Ff\n", 20, mpfdist2);
								
				mpf_div(force_mpz, force_multiplier, mpfdist2);
				
				//gmp_printf("f: %.*Ff\n", 20, force_mpz);

				mpf_div(distratx, mpfxdist, mpfdist);
				mpf_div(distraty, mpfydist, mpfdist);
				
				mpf_mul(force_x, force_mpz, distratx);
				mpf_mul(force_y, force_mpz, distraty);
				
				//gmp_printf("x: %.*Ff\n", 20, force_x);
				
				mpf_add(total_force_x, total_force_x, force_x);
				mpf_add(total_force_y, total_force_y, force_y);
				
				//gmp_printf("t: %.*Ff\n", 20, total_force_x);
			}
		}			
		mpf_mul(force_x2, total_force_x, total_force_x);
		mpf_mul(force_y2, total_force_y, total_force_y);
		
		mpf_add(force_len2, force_x2, force_y2);
		mpf_sqrt(force_len, force_len2);	
		
		mpf_div(xforcerat, total_force_x, force_len);
		mpf_div(yforcerat, total_force_y, force_len);
		
		mpf_mul(step_x, mpfstep_size, xforcerat);
		mpf_mul(step_y, mpfstep_size, yforcerat);
				
		//Now that we have the total force on our electron, we move it a bit.
		t1 = mpf_get_d(step_x); t2 = mpf_get_d(step_y);
		//printf("%f\n", t1); fflush(stdout);
		x1 = x1 + t1;
		y1 = y1 + t2;
		
		//We need to make sure that the electron is not leaving the disk.
		dist2 = x1 * x1 + y1 * y1;
		if(dist2 > (radius * radius)){
			//We need to move it back in! (to the nearest edge)
			dist = sqrt(dist2);
			x1 = x1 / dist * radius;
			y1 = y1 / dist * radius; 
		}
		
		//We need to prevent collisions. If there is an electron already at this position, we need to move it slightly
		RESTART_CHECK: ; //no-op 
		for(k = 0; k < i; k++){
			if( k != i){
				el3 = newstate[k];
				x3 = el3.x;
				y3 = el3.y;
				
				if(x3 == x1 && y3 == y1){
					//We move new_x and new_y slightly away from the other electron
					if(x1 >= 0){
						x1 -= d; 
					}else{
						x1 += d;
					}
						
					if(y1 >= 0){
						y1 -= d;
					}else{
						y1 += d;
					}
					
					//We restart the check.
					goto RESTART_CHECK;
				}
			}
		}
		
		//We now need to store the new state.
		newstate[i].x = x1;
		newstate[i].y = y1;
	}
	
	mpf_clear(force_mpz);
	mpf_clear(mpfdist2);
	mpf_clear(mpfxdist);
	mpf_clear(mpfxdist2);
	mpf_clear(mpfydist);
	mpf_clear(mpfydist2);
	mpf_clear(mpfdist);
	mpf_clears(force_x, force_y, total_force_x, total_force_y, force_len, force_len2, force, NULL);
	mpf_clears(distratx, distraty, NULL);
	mpf_clears(force_x2, force_y2, NULL);
	mpf_clears(step_x, step_y, mpfstep_size, NULL);
	mpf_clears(xforcerat, yforcerat, NULL);
}

void calc_pot(unsigned long int noelectrons, double plate_distance, unsigned long int nostepspotential, double x, double y, mpf_t* potential){	
	double z = 0;
	double pot_step_size;
		
	mpf_t addition;
	mpf_t addpot;
	mpf_t mpfpot_step_size;

	mpf_init(addition);
	mpf_init(addpot);
	mpf_init(mpfpot_step_size);

	pot_step_size = plate_distance / nostepspotential;
	mpf_set_d(mpfpot_step_size, pot_step_size);

	//Calculate the capacity
	mpf_set_ui(*potential, 0);
	
	for(z = d; z <= plate_distance; z+= pot_step_size){
		ef_z(noelectrons, x, y, z, plate_distance, &addpot);
		mpf_mul(addition, addpot, mpfpot_step_size);
		mpf_add(*potential, *potential, addition);
	}
	
	mpf_clear(mpfpot_step_size);
	mpf_clear(addpot);
	mpf_clear(addition);
}

void ef_z(unsigned long int noelectrons, double x, double y, double z, double plate_distance, mpf_t* ef_z){
	unsigned long int i;
	electron el;
	double elx, ely;
	double x_dist, y_dist, z_dist, dist2, dist, distrat;
	
	mpf_t el_field;
	mpf_t mpfdist2;
	mpf_t addition;
	
	mpf_init(addition);
	mpf_init(mpfdist2);
	mpf_init(el_field);
	
	mpf_set_d(*ef_z, 0.0);
	for(i = 0; i < noelectrons; i++){
		el  = state[i];
		elx = el.x;
		ely = el.y;
				
		x_dist = elx - x;
		y_dist = ely - y;
		z_dist = z;
		
		dist2 		= z_dist * z_dist + x_dist * x_dist + y_dist * y_dist;
		dist		= sqrt(dist2);
		
		mpf_set_d(mpfdist2, dist2);
		mpf_div(el_field, field_multiplier, mpfdist2);
	
		distrat = z_dist / dist;
		mpf_set_d(addition, distrat);
								
		mpf_mul(addition, addition, el_field);
		mpf_add(*ef_z, *ef_z, addition);
	}
	
	mpf_clear(mpfdist2);
	mpf_clear(addition);
	mpf_clear(el_field);
}
