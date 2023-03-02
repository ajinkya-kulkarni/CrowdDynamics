
// Code for Wet Granular Flow and Self Proplled particles for MONODISPERSE particles. 

// This is the universal code for elliptical particles using the multi-sphere model, and self propelled particles //

// Developed at Multiphase Flows Lab, Department of Applied Mechanics, IIT Madras.

// Last modified on Feb 5, 2016 by Ajinkya Kulkarni. Email: kulkajinkya@gmail.com

#include<stdio.h>
#include<math.h>
#include<time.h>
#include<stdlib.h> // This library is used for drand48 function which generated random numbers later in the force_calc function.

/*******************************************************************************************
Domain size is the span of the domain. Diameter of the circular cavity or Length of the square cavity.

Modification in the code: Viscous drag relation updated to include mass averaged relative velocity instead of number averaged, since the heavier particles' velocity have a greater say in the fluid velocity value.

Additional monitor: Mass mixing parameter added in addition to the existing number mixing parameter. Total energy directly available that can be plotted with mixing parameter.

Polydisperse material property with different stiffness constants as inputs and evaluating equivalent stiffness coefficient in order to calculate the force on the particles.

!!! Don't change the parameters 0.008 from the force_calc function and the parameter 0.01 from the post processing function. Ever. !!!

Also contains the self propelled force parameter 'beta' 

*******************************************************************************************/

#define tf 20 // Flow time
#define dt 0.00005  // time step

#define pt 7289 // TOTAL number of particles (wall plus free particles)
int core=6120; // free particles

int multisphere_model=0; // if 1, activates the multisphere model for elliptical bacteria, currently comprising of 3 spheres

int modified_beta=0; // if 1, activates the modified beta model for elliptical bacteria, which modifies the self-propelling force into the direction of the orientation of the chain. Only to be used in conjunction with the multi-sphere model.

int random_orientation=0; //if 1, activates the multisphere model for elliptical bacteria with random initial conditions. If 0, the chains are arranged in a rotary sense. Only to be used in conjunction with the multi-sphere model.

float limit_overlap=0; // lower limit for the overlap between 2 particles, usually 1e-7. Eliminates the numerical error associated with determining the particle overlap.

float rotating_wall=1; // switch for the rotating wall condition. ON if 1, else OFF
/*************************************************************************************************************************************************************************/


// We start the simulation in the following way.

//Taking random initial velocities (whose net vector sum is zero) as the input, but no spring forces for the first iteration. From second iteration onwards, the spring forces are calculated.


/*************************** part pertaining to multi-sphere model ********************/

float ang_acc[pt];
float ang_acc_prev[pt];

float omega_chain[pt]; // calculate the omega of the chain
float omega_chain_prev[pt];

float net_moment[pt];

float theta_initial[pt];
float theta_chain[pt];
float theta_original[pt];

float mominertia = 1e-5; // moment of inertia for the chain

float rand_or[pt]; // random array for random initial locations AND orientations

int sign=1; // signum function used in initial co-ordinates section.

/*************************************************************************************/

#define anifrq 1000 //frequency for writing data for simulation
#define cellfrq 1 //frequency for cell information updates

#define mw 1169  // Sum total of all the wall particles + baffle particles
#define mw1 284 // 1st wall particles
#define mw2 290 // 2nd wall particles
#define mw3 295 // 3rd wall particles
#define mw4 300 // 4th wall particles

/******************************************* part pertaining to Rotating walls and baffle(s) if any ******************************************/

float theta_wall[pt]; // angle measurement for the wall particles. To move the wall.
float r_wall=0;

#define mbaff 29 // particles in the baffle(s) if present

/*********************************************************************************************************************************************/

#define ga 0 // Gravitational accelaration <--- switched off if 0, on if >0
#define beta 1 // Self-propelling parameter <--- switched off if 0, on if 1
float gamma_limit = 1;

float dm = 0.00500; // particle diameter
float mass = 0.004189; // mass of each particle

float cv; // Hydrodynamic coordination coefficient which is mu/rho.

//float limit; // the bounding factor for the random langevin forces. The value is taken from the limit.txt file

float video_time; // seconds of video generated in realtime.

#define del 0.000625 // del square which is equal to (5*(particle dia))^2. It is the radius of influence

#define rho 100 // rho equals difference of densities between the fluids.

/********************************************** circular wall related stuff ************************************************/

float omega=0; // Omega in radians/sec of the rotating drum.
float pi=3.14159265358979323846264338327; // pi
float x_c=0.4; // x-coord of the centre of the cavity
float y_c=0.4; // y-coord of the centre of the cavity

float r_cavity1=0.2257; // radius for the 1st moving wall only if the domain is circular!

float r_cavity2; //radius for the 2nd moving wall which is radius of 1st+sqrt(3)*(radius of the particle).

float r_cavity3; // radius for the 3rd moving wall which is radius of 2nd+sqrt(3)*(radius of the particle).

float r_cavity4;; // radius for the 4th moving wall which is radius of 3rd+sqrt(3)*(radius of the particle).

/***************************************************************************************************************************/


/************************************* domain generalisation *****************************/

#define length 0.8 // total width of domain. 

float domain_length;

int nfc=length/0.008;  // used in the force calc function which equals width/0.008. Formerly used to be 52.
//nfc=length/0.008;

int npp=length/0.01; // used in the post processing function which equals width/0.01. Formerly used to be 40.
//npp=length/0.01;

int nfcsq; // square of the nfc value, formerly known as 2704.
int nppsq; // square of the npp value formerly known as 1600.

int cell_pt[10000][200];  // !!! Change this value for every different domain. It is the square of the nfc value, formerly known as 2704.

//float grid[62500][3]; // !!! Change this value for every different domain. It is the square of the nfc value, formerly known as 2704.

/***************************************************************/

#define sigma 213.34

float V_surr = del;

int timestep;

void force_calc_initial(float fx[],float fy[],float x[],float y[],float vx[],float vy[], float k_stiff[], float cv);

void force_calc(float fx[],float fy[],float x[],float y[],float vx[],float vy[], float k_stiff[], float cv, float gamma[], float video_time);

int cj=0,cf=0;
int tempg1,tempg2,tempg3;

float exf[801];

float ds; // distance between two particles

/***************************************************************************************************************************/


int main()

{



  FILE *fda,*fdv;//, *force_particles;
  
 // FILE *fdmeanv;

  FILE *cordin,*velin;
  //FILE *lin;
  FILE *cvin;
  
  FILE *gammain;

  float vx[pt],vy[pt]; //velocity

  float v_tot[pt]; // squares of the velocities
  float v_tot_final[pt]; // total velocity.

  float r_tot[pt]; // squares of the distances
  float r_tot_final[pt]; // total distance.
  
  float x[pt],y[pt]; // coordinates

  float xi[pt],yi[pt];//initial coordinates

  float vxi[pt],vyi[pt];//initial velocity

  float k_stiff[pt]; // individual stiffness coefficient (spring constant of the particles)

  float gamma[pt];

  nppsq=npp*npp;

  float time1,dx,dy;

  float fxi[pt],fyi[pt];///previous force
  float fx[pt],fy[pt]; //current force

  int i=0,timestep,k,a,ts;
  
  float d=0;
  float di=0.00005;
  
  for(i=0;i<801;i++)
    {
      exf[i]=exp(-(d*d)/del);
      d=d+di;
    }  


/************************ STRUCTURE IMPORT ********************************/ 

  cordin=fopen("init_positions.txt","r");
  velin=fopen("init_velocities.txt","r"); 
  //lin=fopen("limit.txt","r");
  cvin=fopen("cv.txt","r");
   

  gammain=fopen("gamma.txt","r");


  for(i=0;i<pt;i++)
    {
      if(fscanf(cordin,"%d %f %f %f\n",&i,&x[i],&y[i],&k_stiff[i])){};
      if(fscanf(velin,"%d %f %f\n",&i,&vx[i],&vy[i])){};  

	if(fscanf(gammain,"%d %f\n",&i,&gamma[i])){};

  
    }
  //if(fscanf(lin,"%f",&limit)){};

  if(fscanf(cvin,"%f",&cv)){};
    
  fclose(cordin);
  fclose(velin); 
  //fclose(lin); 
  fclose(cvin);

  fclose(gammain);


/****************** hardcoding the locations of the core particles, so as to initially keep them as a chain ********************/ 

if (multisphere_model == 1)
	{ 
		srand48(time(NULL));

		for(i=0;i<core;i++)
   		{					
		rand_or[i] = drand48();
   		}
 
		for(i=0;i<(core-1);i++)
   			{					
				theta_initial[i] = atan2((y[i]-y_c),(x[i]-x_c));
   			}

		for(i=0;i<core;i++)
      		{

      		  if(i%3 == 0)
	    

	    		{

			if(random_orientation == 1)

				{ 
	    				sign = sign*(-1);
	    				
	    				x[i+1] = x[i] + (sign)*dm*cos((sign)*rand_or[i]);   // hardcoded the co-ordinates for the random orientation case.
	    				y[i+1] = y[i] + (sign)*dm*sin((sign)*rand_or[i]); 
	     		
	     				x[i+2] = x[i] - (sign)*dm*cos((sign)*rand_or[i] + pi); 
	    				y[i+2] = y[i] - (sign)*dm*sin((sign)*rand_or[i] + pi); 
				}

			else
				{
	    				sign = sign*(1);
	    				
	    				x[i+1] = x[i] + (sign)*dm*cos((sign)*theta_initial[i]+pi/2);  // hardcoded the co-ordinates for the rotary sense case.
	    				y[i+1] = y[i] + (sign)*dm*sin((sign)*theta_initial[i]+pi/2);
	    		
	     				x[i+2] = x[i] - (sign)*dm*cos((sign)*theta_initial[i]+pi/2);
	    				y[i+2] = y[i] - (sign)*dm*sin((sign)*theta_initial[i]+pi/2);
				}					

      			}


		}

	for(i=0;i<(core-1);i++)
   	{
					
	theta_original[i] = atan2((y[i+1]-y[i]),(x[i+1]-x[i])); // preserving the original orientations of the chains, so that the simulations start from this orientation and NOT from zero angle.

   	}


}
/*********************************************************************************************************************/


dx=length/npp;  // !!!
dy=length/npp;  // !!!
  
nfcsq=nfc*nfc;


force_calc_initial(fx,fy,x,y,vx,vy,k_stiff,cv); // first call of the function. Starts the simulation.


/************************************************************ Opening the text files to be written ********************************************************/
fda=fopen("data.txt","w");
fdv=fopen("velocity.txt","w");
//fdmeanv=fopen("meanvel.txt","w");
/*********************************************************************************************************************/


ts = tf/dt + 1;

if (multisphere_model == 1)
	{

	for(i=0;i<pt;i++) // initialising all the elements to be zero. 

		{

 			ang_acc[i] = 0;
 			ang_acc_prev[i] = 0;

 			omega_chain[i] = 0; // calculate the omega of the chain
 			omega_chain_prev[i] = 0;

 			net_moment[i] = 0;
 			theta_chain[i] = 0;

		}
	}


//***************************************** cicular moving wall(s) section (velocity updation section) **************************************//

for(i=0;i<pt;i++) // initializing the array to zero. Necessary evil.
{
theta_wall[i]=0;
}

if(omega==0)

{
	for(i=core;i<core+mw;i++) // since omega of the circular wall is zero, the wall velocities are zero
		{
			vx[i] = 0;
			vy[i] = 0;		
		}
}

else

{

float alpha[pt]; // Used in calculating the (moving) wall velocities.

for(i=0;i<pt;i++) // initializing all alphas to be zero.
	{
		alpha[i]=0;
	}

for(i=core;i<core+mw1;i++) // moving wall no 1
    
    {
    

      alpha[i] = atan2((y[i]-y_c),(x[i]-x_c)); //x_c, y_c is the centre of the cavity
   
	  vx[i] = -omega*r_cavity1*sin(alpha[i]);
	  
	  vy[i] = omega*r_cavity1*cos(alpha[i]);
                               	
	}


for(i=core+mw1;i<core+mw1+mw2;i++) // moving wall no 2
    
    {
    
    alpha[i] = atan2((y[i]-y_c),(x[i]-x_c)); //x_c, y_c is the centre of the cavity
   

          r_cavity2 = r_cavity1 + sqrt(3)*(dm/2); 

	  vx[i] = -omega*r_cavity2*sin(alpha[i]);
	  
	  vy[i] = omega*r_cavity2*cos(alpha[i]);
	                                   	
	}


for(i=core+mw1+mw2;i<core+mw1+mw2+mw3;i++) // moving wall no 3
{

alpha[i] = atan2((y[i]-y_c),(x[i]-x_c)); //x_c, y_c is the centre of the cavity
   
	  
          r_cavity3 = r_cavity2 + sqrt(3)*(dm/2);

	  vx[i] = -omega*r_cavity3*sin(alpha[i]);
	  
	  vy[i] = omega*r_cavity3*cos(alpha[i]);
	  
                                	
	}


for(i=core+mw1+mw2+mw3;i<core+mw1+mw2+mw3+mw4;i++) // moving wall no 4
{

alpha[i] = atan2((y[i]-y_c),(x[i]-x_c)); //x_c, y_c is the centre of the cavity
   
	  
          r_cavity4 = r_cavity3 + sqrt(3)*(dm/2);

	  vx[i] = -omega*r_cavity4*sin(alpha[i]);
	  
	  vy[i] = omega*r_cavity4*cos(alpha[i]);
                                 	
	}
}


/************************************************************** end of moving wall(s) section ***********************************************************/


// Start of the time loop //

  for(timestep=0;timestep<ts;timestep++)
    {
      time1=(timestep+1)*dt;
      cf=timestep/cellfrq;
      cj=timestep;
      a=timestep/anifrq;
      
      if(timestep==a*anifrq)

	{
   	  
	  video_time = timestep*dt;

          if (omega==0)

          	{
          	
          		 if (multisphere_model == 1)
					{
						
						printf("Self Propelled mode is ON for elliptical particles\n");
					if (random_orientation == 1)

						{
						printf("Intitial orientations are random\n");	
						}

					else
						{
			
						printf("Intitial orientations are ordered in a rotary sense\n");

						}
					if (modified_beta == 1)
						{
						printf("Thrust vector is in the direction of orientation\n");
						}	
					else
						{
						printf("Thrust vector is in the direction of velocity\n");		
						}

						printf("Moment of Inertia is %.2e\n",mominertia);
					}
				else
					{
             					printf("Self Propelled mode is ON for circular particles\n");
            				}
            }

          else

          	{
            	 printf("Rotating Drum mode is ON, with an omega of %0.2f radians/sec\n",omega); 
			if(rotating_wall==1)
				{
             				printf("Walls are rotating\n");					
				}
			else
				{
				        printf("Walls are not rotating\n");					
				}
		}
			
          printf("Current CV is %0.2f\n",cv);
          
          if (omega==0)
          	{
		
          	printf("Current Beta is %d\n",beta);
		
		printf("Nucleation till %0.2f secs\n",gamma_limit);
          	
		}
		
          	
          printf("Number of particles simulated: %d\n",core);

          printf("Total Flow time is %d secs\n",tf);

          domain_length=length;

          printf("Diameter of the domain is %0.1fm\n",domain_length);
                             

          printf("Seconds of video generated till now: %0.2f secs\n",video_time);

          printf("\n");

	  for(k=0;k<pt;k++)
	    {
		

	      v_tot[k]=(vx[k]*vx[k]+vy[k]*vy[k]);
 	      v_tot_final[k]=sqrt(v_tot[k]); // total velocity of the particles.

	      r_tot[k]=((x[k]-x_c)*(x[k]-x_c)+(y[k]-y_c)*(y[k]-y_c));
 	      r_tot_final[k]=sqrt(r_tot[k]); // 
 	      
 	      fprintf(fda,"%f\t%f\t%f\t%f\t%f\n",time1-dt,x[k],y[k],dm/2,r_tot_final[k]);

	      fprintf(fdv,"%f\t%f\t%f\t%f\n",time1-dt,vx[k],vy[k],v_tot_final[k]);

	    }


	 // for(k=0;k<core;k++)
	   // {
		  
 	  //    fprintf(fdmeanv,"%f\t%f\t%f\n",time1-dt,vx[k],vy[k]);

	  //  }

	  fprintf(fda,"\n");
	  fprintf(fdv,"\n");
	//fprintf(fdmeanv,"\n");

	}



/***************************************** Velocity VERLET ALGORITHM **************************************/

	
for(i=0;i<core;i++)
{
	
     	 x[i]=x[i]+(vx[i]*dt)+(((dt*dt)*(fx[i]/mass))/2);
      	 y[i]=y[i]+(vy[i]*dt)+(((dt*dt)*(fy[i]/mass))/2);

// Pseudo v(t + del_t) calculation
      
	  vxi[i]=vx[i];
	  vyi[i]=vy[i];
	  	  
	  vx[i]=vx[i]+((fx[i]/mass)*dt);
	  vy[i]=vy[i]+((fy[i]/mass)*dt);
		  	
	  fxi[i]=fx[i];
	  fyi[i]=fy[i];
			
			
if (multisphere_model == 1)
	{
 	  
	  	if(i%3 == 0)
	  	
			{
	 		
	  		
	  				theta_chain[i] = theta_chain[i] + omega_chain[i]*dt + ang_acc[i]*dt*dt*0.5;
	  				
					x[i+1] = x[i] + dm*cos(theta_chain[i] + theta_original[i]); // rotation relative to x[i]
					y[i+1] = y[i] + dm*sin(theta_chain[i] + theta_original[i]);
						    
					x[i+2] = x[i] - dm*cos(theta_chain[i] + theta_original[i]);
					y[i+2] = y[i] - dm*sin(theta_chain[i] + theta_original[i]);
							
					
	  				omega_chain_prev[i] = omega_chain[i];
	  
	  				omega_chain[i] = omega_chain[i] + (ang_acc[i]*dt);
	  		
	  				ang_acc_prev[i] = ang_acc[i];
	  		}
	  		
	  }	
	
}



/*************************************************************************************************************/


force_calc(fx,fy,x,y,vx,vy,k_stiff,cv,gamma,video_time);


/********************************* part pertaining to rotating wall and rotating baffle(s) **********************************/

if ((rotating_wall==1) && (omega!=0))
{

for(i=core+mw1;i<pt;i++)

{
	   						
		
//if(i<core+mw1)
//{
//r_wall=r_cavity1;
//theta_wall[i] = atan2((y[i]-y_c),(x[i]-x_c));
//theta_wall[i] = theta_wall[i] + 2*asin(0.5*dm/r_cavity1) + omega*dt;
//}

if(i<core+mw1+mw2) // moving the 2nd wall
{
r_wall = r_cavity1 + sqrt(3)*(dm/2);
theta_wall[i] = atan2((y[i]-y_c),(x[i]-x_c));
theta_wall[i] = theta_wall[i] + 2*asin(0.5*dm/r_cavity2) + omega*dt;
}

if((i>=core+mw1+mw2) && (i<core+mw1+mw2+mw3)) // moving the 3rd wall
{
r_wall = r_cavity1 + 2*sqrt(3)*(dm/2);
theta_wall[i] = atan2((y[i]-y_c),(x[i]-x_c));
theta_wall[i] = theta_wall[i] + 2*asin(0.5*dm/r_cavity3) + omega*dt;
}

if((i>=core+mw1+mw2+mw3) && (i<core+mw1+mw2+mw3+mw4)) //moving the 4th wall
{
r_wall = r_cavity1 + 3*sqrt(3)*(dm/2);
theta_wall[i] = atan2((y[i]-y_c),(x[i]-x_c));
theta_wall[i] = theta_wall[i] + 2*asin(0.5*dm/r_cavity4) + omega*dt;
}


		x[i] = x_c + r_wall*cos(theta_wall[i]); // rotation relative to x_c
		y[i] = y_c + r_wall*sin(theta_wall[i]);  
	

//for(i=0;i<pt;i++)
//{
//theta_wall[i]=0;
//}

//for(i=core+mw1+mw2+mw3+mw4;i<core+mw1+mw2+mw3+mw4+mbaff;i++) // rotate the particles in the baffles
  // 	{					
	//   theta_wall[i] = atan2((y[i]-y_c),(x[i]-x_c));
  //	}

//for(i=core+mw1+mw2+mw3+mw4;i<core+mw1+mw2+mw3+mw4+mbaff;i++) // rotate the particles from the 1st wall

//	{	
//		r_baff[i]=

//		theta_wall[i] = theta_wall[i] + omega*dt; //theta_wall + 1*asin(0.5*dm/r_cavity1);

//		x[i] = x_c + r_baff*cos(theta_wall[i]); // rotation relative to x[i]
//		y[i] = y_c + r_baff*sin(theta_wall[i]);  
//	}

}						

}	


/**************************************************** End of rotating walls and baffle(s) section *********************************************************/


	
for(i=0;i<core;i++)
    
	{
		
	  vx[i]=vxi[i]+((fxi[i]/mass)+(fx[i]/mass))*dt*0.5; // Bug fixed for the first term on the RHS, vx[i] replaced with vxi[i] 
	  
	  vy[i]=vyi[i]+((fyi[i]/mass)+(fy[i]/mass))*dt*0.5; // Bug fixed for the first term on the RHS, vy[i] replaced with vyi[i]
	  
	if (multisphere_model == 1)
		{
		if(i%3 == 0)
	
			{

				omega_chain[i] = omega_chain_prev[i] + (ang_acc_prev[i] + ang_acc[i])*dt*0.5;
						      				 		
				vx[i+1] = vx[i] - dm*omega_chain[i]*sin(theta_chain[i] + theta_original[i]); //reassigned velocities using omega and angular accelaration
	    			vy[i+1] = vy[i] + dm*omega_chain[i]*cos(theta_chain[i] + theta_original[i]);
	    					    				
				vx[i+2] = vx[i] + dm*omega_chain[i]*sin(theta_chain[i] + theta_original[i]);
	    			vy[i+2] = vy[i] - dm*omega_chain[i]*cos(theta_chain[i] + theta_original[i]);
	    			
					
			}
		}
	
	



	 }


/***************************************************************************************************************************/

  }

  fclose(fda);
  fclose(fdv);
//fclose(fdmeanv);

}


/************************************************* INITIAL force calculation function *********************************************/

void force_calc_initial(float fx[],float fy[],float x[],float y[],float vx[],float vy[],float k_stiff[],float cv)

{
  
  int i,j,k,id;
  float dx,dy,d,f;
  float k_eq; // equivalent stiffness coefficient
  float vel_mag, v_i = 0.0, v_j = 0.0; 
  float fvx,fvy;
  int bin[9],cal[500];
  int pointer;

  float exv,ex,eyv,ey,tmpexp;
  float m_part, m_fluid;
  int temp1,temp2,temp3,extemp;


/******************** Random number generation for random forces (to be activated only if there is an existence of a random langevin force)********************/

	 	   
//printf("f_rand is %0.20f\n",f_rand);
//f_rand_x = ((float)rand()/(float)(RAND_MAX/limit));// - (limit/2);
//f_rand_y = ((float)rand()/(float)(RAND_MAX/limit));//((float)rand()/(float)(RAND_MAX/limit1));// - (limit1/2);//f_rand_x;//(float)rand()/(float)(RAND_MAX/limit1);

//f_rand_x = r[1];
//f_rand_y = r[2];

//float a = 1.0;
   // for (int i=0;i<2;i++)
   //  rand_nos = ((float)rand()/(float)(RAND_MAX)) * a;

//printf("x random force is %0.20f\n",f_rand_x);
//printf("y random force is %0.20f\n",f_rand_y);

//**********************************************************************************//


if(cf==cj*cellfrq)
    {
      for(i=0;i<nfcsq;i++)  

/*************************************************** Grid generation ***************************************************/
	{	
	  cell_pt[i][1]=0;

	  for(k=2;k<npp;k++)
	    {
	      cell_pt[i][k]=-1;
	    } 
	}   

   
      for(i=0;i<pt;i++) //Extration of cell ID
	{      
	  tempg1=x[i]/0.008;
	  tempg2=y[i]/0.008;
	  
	  tempg3=((tempg2+1)*nfc)+(tempg1+1);  
	  cell_pt[tempg3][1]=cell_pt[tempg3][1]+1;
	  

	  tempg1=cell_pt[tempg3][1];
	  cell_pt[tempg3][tempg1+1]=i;
	   
	}
	 
    }
   
/*******************************************************************************************************************************/


  for(i=0;i<core+mw;i++) //start of for loop
    {
    	
      fx[i]=0;
      fy[i]=-mass*ga;	
      temp1=x[i]/0.008;
      temp2=y[i]/0.008;      
      id=((temp2+1)*nfc)+(temp1+1); //cell ID

/******************************************* Binning the particles ****************************************************/
     
      bin[0]=id+nfc+1;    
      bin[1]=id+nfc;
      bin[2]=id+nfc-1;
      bin[3]=id+1;
      bin[4]=id;
      bin[5]=id-1;
      bin[6]=id-(nfc-1);
      bin[7]=id-nfc;
      bin[8]=id-(nfc+1);
      
/*******************************************************************************************************************************/
      
      pointer=0;
   
      for(k=0;k<9;k++)
	{
          if(cell_pt[bin[k]][1]>0)
	    {
	      for(j=1;j<=cell_pt[bin[k]][1];j++)
		{                             
		  cal[j-1+pointer]=cell_pt[bin[k]][j+1];            
		}
	      pointer=pointer+cell_pt[bin[k]][1];
	    }
	}

      exv=0;
      ex=0;
      eyv=0;
      ey=0;
      m_part=0;

for(j=0;j<pointer;j++) //inner for loop
{     
   	   if(i!=cal[j])

	 { 
	      dx=x[cal[j]]-x[i];
	      dy=y[cal[j]]-y[i];
	      		   
	      d=(sqrt((dx*dx)+(dy*dy)));
	  
/***************************************** Fluid drag force part ****************************************/	
	

extemp=((d/del)*100)/8;
	  
tmpexp=exf[extemp];

ex=ex + tmpexp* mass; // mass weighted instead of simply number weighted - bug fixed w.r.t. v3
	      
m_part = m_part + mass;

exv=exv + tmpexp * (mass*vx[cal[j]]); // mass weighted instead of simply number weighted - bug fixed w.r.t. v3
 	      
eyv=eyv + tmpexp* (mass*vy[cal[j]]); // mass weighted instead of simply number weighted - bug fixed w.r.t. v3
	      
	 }
}
	
	
m_fluid = rho*(V_surr - m_part/sigma);
    

	
      if(ex==0)

	{
	  fvx=cv*dm*(0 - vx[i]);
	  fvy=cv*dm*(0 - vy[i]);
	}
	

	
      if(ex!=0)

	{
	  fvx=cv*dm*((m_part/(m_part + m_fluid))*(exv/ex) - vx[i]);
	  fvy=cv*dm*((m_part/(m_part + m_fluid))*(eyv/ex) - vy[i]);


	}

	
/******************************************************* Self Propelled forces part **************************************************/

if(modified_beta == 1)

		{
				
      				fx[i]=fx[i]+ fvx + beta*mass*cos(theta_original[i]); // self propelling force is in the direction of orientation vector. (i+1)st particle is considered leading always

      				fy[i]=fy[i]+ fvy + beta*mass*sin(theta_original[i]);

		}

else
		{
				
				vel_mag = sqrt(vx[i]*vx[i] + vy[i]*vy[i]);

       				 if (vel_mag > 1e-30) 
          			  {
            				v_i = vx[i]/sqrt(vx[i]*vx[i] + vy[i]*vy[i]);
            				v_j = vy[i]/sqrt(vx[i]*vx[i] + vy[i]*vy[i]);
            			  }				


				fx[i]=fx[i]+ fvx + beta*mass*v_i; //self propelling force is in the direction of velocity vector.
      				fy[i]=fy[i]+ fvy + beta*mass*v_j;

			
	
		}


			
/*************************************************** moment calculation for the chain ***********************************************/

if (multisphere_model == 1)

{

if((i%3 == 0) && (i<core))

		{


			
			net_moment[i] = ((((x[i+2]-x[i])*fy[i+2]) -  ((y[i+2]-y[i])*fx[i+2])) + (((x[i+1]-x[i])*fy[i+1]) - ((y[i+1]-y[i])*fx[i+1])));

			
			ang_acc[i] = (net_moment[i]/mominertia);


	    		fx[i] = fx[i] + fx[i+1] + fx[i+2];
	    		fy[i] = fy[i] + fy[i+1] + fy[i+2];


	     	}

      		
      	}

/***************************************************************************************************************************************/

     }


}



/************************************************* force calculation function *********************************************/

void force_calc(float fx[],float fy[],float x[],float y[],float vx[],float vy[],float k_stiff[],float cv, float gamma[], float video_time)
{
  
  int i,j,k,id;
  float dx,dy,d,f;
  float k_eq; // equivalent stiffness coefficient
  float vel_mag, v_i = 0.0, v_j = 0.0; 
  float fvx,fvy;
  int bin[9],cal[500];
  int pointer;

  float exv,ex,eyv,ey,tmpexp;
  float m_part, m_fluid;
  int temp1,temp2,temp3,extemp;



/******************** Random number generation for random forces (to be activated only if there is an existence of a random langevin force)********************/

	 	   
//printf("f_rand is %0.20f\n",f_rand);
//f_rand_x = ((float)rand()/(float)(RAND_MAX/limit));// - (limit/2);
//f_rand_y = ((float)rand()/(float)(RAND_MAX/limit));//((float)rand()/(float)(RAND_MAX/limit1));// - (limit1/2);//f_rand_x;//(float)rand()/(float)(RAND_MAX/limit1);

//f_rand_x = r[1];
//f_rand_y = r[2];

//float a = 1.0;
   // for (int i=0;i<2;i++)
   //  rand_nos = ((float)rand()/(float)(RAND_MAX)) * a;

//printf("x random force is %0.20f\n",f_rand_x);
//printf("y random force is %0.20f\n",f_rand_y);

//**********************************************************************************//


if(cf==cj*cellfrq)
    {
      for(i=0;i<nfcsq;i++)  

/*************************************************** Grid generation ***************************************************/
	{	
	  cell_pt[i][1]=0;

	  for(k=2;k<npp;k++)
	    {
	      cell_pt[i][k]=-1;
	    } 
	}   

   
      for(i=0;i<pt;i++) //Extration of cell ID
	{      
	  tempg1=x[i]/0.008;
	  tempg2=y[i]/0.008;
	  
	  tempg3=((tempg2+1)*nfc)+(tempg1+1);  
	  cell_pt[tempg3][1]=cell_pt[tempg3][1]+1;
	  

	  tempg1=cell_pt[tempg3][1];
	  cell_pt[tempg3][tempg1+1]=i;
	   
	}
	 
    }
   
/*******************************************************************************************************************************/


  for(i=0;i<core+mw;i++) //start of for loop
    {
    	
      fx[i]=0;
      fy[i]=-mass*ga;	
      temp1=x[i]/0.008;
      temp2=y[i]/0.008;      
      id=((temp2+1)*nfc)+(temp1+1); //cell ID

/******************************************* Binning the particles ****************************************************/
     
      bin[0]=id+nfc+1;    
      bin[1]=id+nfc;
      bin[2]=id+nfc-1;
      bin[3]=id+1;
      bin[4]=id;
      bin[5]=id-1;
      bin[6]=id-(nfc-1);
      bin[7]=id-nfc;
      bin[8]=id-(nfc+1);
      
/*******************************************************************************************************************************/
      
      pointer=0;
   
      for(k=0;k<9;k++)
	{
          if(cell_pt[bin[k]][1]>0)
	    {
	      for(j=1;j<=cell_pt[bin[k]][1];j++)
		{                             
		  cal[j-1+pointer]=cell_pt[bin[k]][j+1];            
		}
	      pointer=pointer+cell_pt[bin[k]][1];
	    }
	}

      exv=0;
      ex=0;
      eyv=0;
      ey=0;
      m_part=0;

      for(j=0;j<pointer;j++) //inner for loop
	{     
   	   if(i!=cal[j])
	    { 
	      dx=x[cal[j]]-x[i];
	      dy=y[cal[j]]-y[i];
	      		   
	     d=(sqrt((dx*dx)+(dy*dy)));
	  
	     ds = dm;//((dm[i]+dm[cal[j]])/2); // !!!!! Major bug fixed!!!!! "dm[j]" replaced by "dm[cal[j]]"
	      


float overlap = (ds-d);	// compressive overlap between the particles	


 if(d<0.2*ds) // preventing excess compression. Manually assigning the value to 0.1*ds;
	{
	  d=0.1*ds; // decide upon the limit for this. why 0.1?
	}

 if(overlap>limit_overlap)

	{
	
//if((j==i+1) || (j==i-1))
//{
//printf("i is %d\n",i);
//printf("j is %d\n",j);
//printf("overlap is %0.20f\n",overlap);
//}
       	  k_eq = (k_stiff[i]*k_stiff[cal[j]])/(k_stiff[i] + k_stiff[cal[j]]);

	  f=(-k_eq*(ds-d));
		  
	 fx[i]=fx[i]+(f*(dx/d));// + (f_rand_x*(dx/d)); // f_rand_x are the pairwise random forces.

	 fy[i]=fy[i]+(f*(dy/d));// + (f_rand_y*(dy/d));




	}


/***************************************** Fluid drag force part ****************************************/	
	

extemp=((d/del)*100)/8;
	  
tmpexp=exf[extemp];

ex=ex + tmpexp* mass; // mass weighted instead of simply number weighted - bug fixed w.r.t. v3
	      
m_part = m_part + mass;

exv=exv + tmpexp * (mass*vx[cal[j]]); // mass weighted instead of simply number weighted - bug fixed w.r.t. v3
 	      
eyv=eyv + tmpexp* (mass*vy[cal[j]]); // mass weighted instead of simply number weighted - bug fixed w.r.t. v3
	      
	 }
}
	
	
m_fluid = rho*(V_surr - m_part/sigma);
    

	
      if(ex==0)

	{
	  fvx=cv*dm*(0 - vx[i]);
	  fvy=cv*dm*(0 - vy[i]);
	}
	

	
      if(ex!=0)

	{
	  fvx=cv*dm*((m_part/(m_part + m_fluid))*(exv/ex) - vx[i]);
	  fvy=cv*dm*((m_part/(m_part + m_fluid))*(eyv/ex) - vy[i]);


	}

	
/******************************************************* Self Propelled forces part **************************************************/

if(modified_beta == 1)

		{
				
      				fx[i]=fx[i]+ fvx + beta*mass*cos(theta_original[i]); // self propelling force is in the direction of orientation vector. (i+1)st particle is considered leading always
      				fy[i]=fy[i]+ fvy + beta*mass*sin(theta_original[i]);

		}

else
		{
				
				vel_mag = sqrt(vx[i]*vx[i] + vy[i]*vy[i]);

       				 if (vel_mag > 1e-30) 
          			  {
            				v_i = vx[i]/sqrt(vx[i]*vx[i] + vy[i]*vy[i]);
            				v_j = vy[i]/sqrt(vx[i]*vx[i] + vy[i]*vy[i]);
            			  }				


				//fx[i]=fx[i]+ fvx + beta*mass*v_i; //self propelling force is in the direction of velocity vector.
      				//fy[i]=fy[i]+ fvy + beta*mass*v_j;



//float gamma=sqrt(1 - beta*beta); // Self-propelling MOTIVE parameter from 0 to 1


				//float theta_gamma[pt]; float vradial_single[pt]; float vtheta_single[pt];
				float theta_gamma = 0;

				/*for(i=0;i<pt;i++)
					{
						theta_gamma[i] = 0;
						vradial_single[i] = 0;
						vtheta_single[i] = 0;
					}

 				*/

				theta_gamma = atan2((y[i]-y_c),(x[i]-x_c));
				//vradial_single[i] = vy[i]*sin(theta_gamma[i]) + vx[i]*cos(theta_gamma[i]);
				//vtheta_single[i] = vy[i]*cos(theta_gamma[i]) - vx[i]*sin(theta_gamma[i]);


if (video_time<gamma_limit)
{
fx[i]=fx[i]+ fvx + beta*mass*v_i - mass*gamma[i]*(sin(theta_gamma));

fy[i]=fy[i]+ fvy + beta*mass*v_j + mass*gamma[i]*(cos(theta_gamma));

//printf("i is %d and gamma is %f\n",i,gamma);

}

else

{
fx[i]=fx[i]+ fvx + beta*mass*v_i; 

fy[i]=fy[i]+ fvy + beta*mass*v_j;
}


//printf("i is %d and gamma is %0.2f\n",i,gamma[i]);


			
	
		}


			
/*************************************************** moment calculation for the chain ***********************************************/

if (multisphere_model == 1)

{

if((i%3 == 0) && (i<core))

		{


			
			net_moment[i] = ((((x[i+2]-x[i])*fy[i+2]) -  ((y[i+2]-y[i])*fx[i+2])) + (((x[i+1]-x[i])*fy[i+1]) - ((y[i+1]-y[i])*fx[i+1])));

			
			ang_acc[i] = (net_moment[i]/mominertia);


	    		fx[i] = fx[i] + fx[i+1] + fx[i+2];
	    		fy[i] = fy[i] + fy[i+1] + fy[i+2];


	     	}

      		
      	}

/***************************************************************************************************************************************/

     }


}


