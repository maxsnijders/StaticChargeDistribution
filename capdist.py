from __future__ import division
import numpy as np
import math

#ALL UNITS IMPLICITLY SI DEFAULT WITHOUT PREFIX UNLESS NOTED OTHERWISE
#ALL COORDINATES CARTESIAN

#CONFIG SECTION
plate_radius 	= 0.1;

#SYSTEM SETTINGS
no_electrons    = 5000;
electron_charge	= -1.60217657E-19;
plate_charge 	= no_electrons * electron_charge;
d				= 1E-9;
el_step_size	= plate_radius / 10;
rel_el_perm		= 1;
coulomb_const	= 1/(4 * math.pi * 8.854187817E-12 * rel_el_perm);
sys_max_steps 	= 50;
bin_side		= plate_radius / 10;

def run():
	#SYSTEM STATE
	electrons = [0] * no_electrons;

	#Set up the initial positions for all electrons
	for i in np.arange(0, no_electrons, 1):
		r = plate_radius * d;
		a = i * (2 * math.pi / no_electrons);
		x = r * np.cos(a);
		y = r * np.sin(a);
		electrons[i] = [x, y];

	#Run the simulation
	for s in np.arange(0, sys_max_steps, 1):
		print "Iteration: ", s+1;
		electrons = electrons_step(electrons);
	
	return electrons;
#===============FUNCTIONS:

#Calculate a step for the current electron state.
def electrons_step(electron_state):
	new_electron_state = electron_state[:];
	
	for i in np.arange(0, no_electrons, 1):
		electron = electron_state[i];
		x1 = electron[0];
		y1 = electron[1];
		
		net_force = np.zeros(2);
		for j in np.arange(0, no_electrons, 1):	
			if( j == i):
				continue;
			electron2 = electron_state[j];
			x2		= electron2[0];
			y2		= electron2[1];
			x_dist = x1 - x2;
			y_dist = y1 - y2;
			dist 	= math.sqrt(x_dist**2 + y_dist**2);
			force   = coulomb_const * electron_charge**2 / (dist**2);
			force_x = force * (x_dist / dist);
			force_y = force * (y_dist / dist);
			net_force[0] += force_x;
			net_force[1] += force_y;
			
		net_force_len = math.sqrt(net_force[0]**2 + net_force[1]**2);
		
		step_x = net_force[0] / net_force_len * el_step_size;
		step_y = net_force[1] / net_force_len * el_step_size;
		
		new_x = x1 + step_x;
		new_y = y1 + step_y;
		
		#Prevent running off of the plate
		dist_from_origin = math.sqrt(new_x**2 + new_y**2);
		if( dist_from_origin > plate_radius ):
			new_x = new_x / dist_from_origin * plate_radius;
			new_y = new_y / dist_from_origin * plate_radius;
	
		#Prevent two electrons at the exact same position.
		should_restart = True;
		while(should_restart):
			should_restart = False;	
			for j in np.arange(0, no_electrons, 1):
				electron3 = new_electron_state[j];
				x3 = electron3[0];
				y3 = electron3[1];
				if(x3 == new_x and y3 == new_y):
					if(new_x >= 0):
						new_x -= d;
					if(new_y >= 0):
						new_y -= d;
					if(new_x < 0):
						new_x += d;
					if(new_y < 0):
						new_y += d;
						
					should_restart = True;
					break;
		
		new_electron_state[i] = [new_x, new_y];
	return new_electron_state;
	
electrons = run();
np.savetxt("electrons.txt", np.array(electrons));
