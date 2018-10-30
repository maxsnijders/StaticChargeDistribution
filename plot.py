from __future__ import division
import numpy as np
import math
import sys

import matplotlib
from mpl_toolkits.mplot3d import axes3d
from matplotlib import pyplot as plt

def run():
	if(len(sys.argv) < 4):
		print "You need to specify 3 arguments: plateradius filenameload outputfilename"
		return;
	
	plate_radius = float(sys.argv[1]);
	electrons = np.loadtxt(sys.argv[2]);
	no_electrons = electrons.shape[0];
	outfile = sys.argv[3];
	
	#We need to output the results		
	fig2 = plt.figure();
	ax2 = fig2.add_subplot(111);
	el_x = [electrons[i][0] for i in range(0,no_electrons)]
	el_y = [electrons[i][1] for i in range(0,no_electrons)]
	
	circle1=plt.Circle((0,0),plate_radius,color='r',fill=False)
	fig2.gca().add_artist(circle1)
	
	
	ax2.scatter(el_x, el_y);
	if(len(sys.argv) >= 5):
		ax2.set_title("Electron Positions: {N} electrons, t={T}".format(N = no_electrons, T=sys.argv[4]));
	else:	
		ax2.set_title("Electron Positions: {N} electrons".format(N = no_electrons));
	ax2.set_xlabel("X pos (m)");
	ax2.set_ylabel("Y pos (m)", rotation="horizontal", labelpad=30);
	ax2.set_aspect('equal');
	ax2.set_xlim([-plate_radius * 1.5, plate_radius * 1.5]);
	ax2.set_ylim([-plate_radius * 1.5, plate_radius * 1.5]);
	

	
	fig2.savefig(outfile, dpi=300);
		
run();
