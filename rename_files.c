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

int main(int argc, char** argv){
	char oldname[200];
	char newname[200];
	char padded[5];
	
	for(int i = 0; i <= 1000; i++){
		sprintf(oldname, "VID_1000-1000/IMG/%i.png", i);
		sprintf(newname, "VID_1000-1000/IMG/%04i.png", padded);
		rename(oldname, newname);
	}
}
