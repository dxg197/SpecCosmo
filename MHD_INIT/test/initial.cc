
/************************************************************************************************************
* By John Hamilton and Cindi Ballard Fall 2005 University Houston Clear Lake                                 *
* Initial conditions for GRMHD under perturbations from gravitational waves.                                *
************************************************************************************************************/

#include <math.h>
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#ifndef PI             //This is included in case this code is tested on a Microsoft compiler
  #define PI 3.14159265358979323846
#endif

/************************************************************************************************************
* List assumptions first                                                                                    *
*	Use Geometrized Units:                                                                                  *
*		Gravitational Constant               G = 1                                                          *
*		Speed of Light                       c = 1                                                          *
*		Boltzmann Constant                   k = 1                                                          *
*		Coulomb Force Constant  epislon_naught = 1/(4 pi)      //D J Griffiths "Intro to E&M" 3rd Ed pg 558 *
*	from en.wikipedia.org                                                                                   *
************************************************************************************************************/

/************************************************************************************************************
List constants for conversions
	Uncomment if required; not current used in this version of code

#ifndef Grav_Const
  #define Grav_Const 6.6727-11      //units are meter^3 kilogram^-1 second^-2
#endif

#ifndef Speed_Light
  #define Speed_Light 2.9979e9      //units are meter second^-1
#endif

#ifndef Boltz_Const
  #define Boltz_Const 1.3807e-23    //units are kilogram meter^2 second^-2 Kelvin^-1
#endif

#ifndef Eps_Naught
  #define Eps_Naught 8.8542e-12     //units are Amps^2 second^4 kilogram^-1 meter^-3
#endif

#ifndef Solar_Mass
  #define Solar_Mass 1.9989e30       //units are kilogram
#endif
************************************************************************************************************/

/************************************************************************************************************
* SI Variables in terms of Solar Mass; by J. Hamilton and C. Ballard                                        *
*	for more information on actual variable manipulation                                                    * 
*	see C. Ballard Fall 2005 Research Paper                                                                 *
*	 (Calculated by Mathematica)                                                                            *
************************************************************************************************************/

#ifndef Kilogram
  #define Kilogram 5.003e-31             //units are Solar_Mass
#endif

#ifndef Second
  #define Second 2.020e5                 //units are Solar_Mass
#endif

#ifndef Meter
  #define Meter 6.738e-4                 //units are Solar_Mass
#endif

#ifndef Kelvin
  #define Kelvin 7.685e-71               //units are Solar_Mass
#endif

#ifndef Density
  #define Density 1.635e-21              //units are Solar_Mass^-2
#endif

#ifndef SIE                              //Specific Internal Energy
  #define SIE 2.172e53                   //units are Solar_Mass^-1  //Must use Temperature in Solar Mass units
#endif

#ifndef Pressure
  #define Pressure 1.819e-38             //units are Solar_Mass^-2
#endif

#ifndef Amp
  #define Amp 2.874e-26                  //unitless in Geometrized Units
#endif

#ifndef Tesla
  #define Tesla 4.265e-12                //units are Solar_Mass^-1
#endif

#ifndef E_Field
  #define E_Field 1.423e-24              //units are Solar_Mass^-1
#endif

/************************************************************************************************************
* Solar Mass Variables in terms of SI units; by C. Ballard and J.Hamilton                                   *
*	for more information on actual variable manipulation                                                    * 
*	see C. Ballard Fall 2005 Research Paper                                                                 *
*	 (Calculated by Mathematica)                                                                            *
************************************************************************************************************/

#ifndef SM_Kilogram
  #define SM_Kilogram 1.999e30           //units are kilogram
#endif

#ifndef SM_Second
  #define SM_Second 4.950e-6             //units are seconds
#endif

#ifndef SM_Meter
  #define SM_Meter 1.484e3               //units are meters
#endif

#ifndef SM_Kelvin
  #define SM_Kelvin 1.301e70             //units are Kelvin
#endif

#ifndef SM_Density
  #define SM_Density 6.116e20            //units are kilogram meter^-3
#endif

#ifndef SM_SIE                           //Specific Internal Energy
  #define SM_SIE 4.605e-54               //units are meter^2 second ^-2 Kelvin^-1
#endif

#ifndef SM_Pressure
  #define SM_Pressure 5.497e37           //units are kilogram meter^-1 seconds^-2
#endif

#ifndef SM_Amp
  #define SM_Amp 3.479e25                //units are sqrt(kilogram meter^3 second^-4)
#endif

#ifndef SM_Tesla
  #define SM_Tesla 2.344e11              //units are kilogram Amp^-1 second^-2
#endif

#ifndef SM_E_Field
  #define SM_E_Field 7.028e23            //units are kilogram meter Amp^-1 second^-3
#endif


#ifndef CONSTANTS          
  #define CONSTANTS            //currently not reserved for any thing
  #define cc 1                 //speed of light is 299792458 m/s
  #define GG 1                 //Newton's Gravitational Constant 6.67e-11 kg m^3/s^2
  #define kB 1                 //Boltzmann Constant 1.38e-23 (kg m^2)/(s^2 K)
/*********************************************************************************************************************	*Coulomb's Force Constant 1/4*pi*epsilon_naught                                                                       *
*	{epsilon_naught 8.85e-12 (s^4*A^2)/(kg*m^3)}                                                                 *
*********************************************************************************************************************/
  #define cfC 1/4*PI           //Coloumb Force Constant
    //FIX in a future release  #define qq 0                 //value of the exponent in line element scaling factor function see Carroll (2.54) unitless
  #define i_rho pow(10.,8)     //init rho_naught of the universe Mathematical Cosmology third frame page 130 units kg/m^3
  #define i_temp pow(10.,10)   //init temp of the universe Mathematical Cosmology third frame page 130 units K
  #define i_time 1.1           //time at start of early universe from Mathematical Cosmology by Islam
  #define i_vel0 1             //init spatial t vel=1 we are comoving; units are m/s  
  #define i_vel1 0             //init spatial x vel=0 we are comoving; units are m/s
  #define i_vel2 0             //init spatial y vel=0 we are comoving; units are m/s
  #define i_vel3 0             //init spatial z vel=0 we are comoving; units are m/s
  #define i_bo0 0              //init Bo field t in Tesla=kg/s^2*A
  #define i_bo1 0              //init Bo field x in Tesla=kg/s^2*A
  #define i_bo2 8.26e-5        //init Bo field y in Tesla=kg/s^2*A (Duez 420)
  #define i_bo3 0              //init Bo field z in Tesla=kg/s^2*A
  #define i_eo0 0              //init Eo field t in E_Field=kg m/s^3*A  (co-moving frame)
  #define i_eo1 0              //init Eo field x in E_Field=kg m/s^3*A
  #define i_eo2 0              //init Eo field y in E_Field=kg m/s^3*A
  #define i_eo3 0              //init Eo field z in E_Field=kg m/s^3*A
  #define i_pressure 9e24      //init Pressure in units of kg/(m*s^2) derived from rho_naught
  #define i_eps 3/2*1.38e-23*i_temp //internal energy units kg*m^2/s^2 or 3/2*k*T    //End of FIX in a future release
#endif


/************************************************************************************************************
* Line Element similar to [Carroll "An Intro to GR" eq. (2.53)]                                             *
*             ds^2=-(c*dx0)^2+[a*(dx1)^2+b*(dx2)^2+d*(dx3)^2]Y(t)^2                                         *
*                                                                                                           *
* Scaling Function in the line element [Carroll "An Intro to GR" eq. (2.54)]                                *
*             Y(t)= t^q  q={0-1}                                                                            *
*                                                                                                           *
* Gravity wave: currently it is linear [Duez arXiv:astro-ph/0503420 v2  8 Jul 2005]                         *
*             g_wave_plus = A+ = Ao*cos(kt)sin(kz)          eq. (66)                                         *
*             g_wave_cross = Ax = Bo*cos(kt)sin(kz)         eq. (67)                                        *
*                                                                                                           *
* Perturbed Metric from [Duez arXiv:astro-ph/0503421 v2  8 Jul 2005 eq. (24)]                               *
*       g   = n   +  h                                                                                      *
*        uv    uv     uv                                                                                    *
*                                                                                                           *
*       /-1  0  0  0 \              / 0     0        0     0 \                                              *
*      |  0  1  0  0  |            |  0  A(t,z)   B(t,z)   0  |                                             *
* n  = |  0  0  1  0  |       h  = |  0  B(t,z)  -A(t,z)   0  |                                             *
*  uv   \ 0  0  0  1 /         uv   \ 0     0        0     0 /                                              *
*                                                                                                           *
*                                                                                                           *
************************************************************************************************************/ 

/************************************************************************************************************
* Stress Energy Tensor [Duez arXiv:astro-ph/0503420 v2  8 Jul 2005 (33)]                                    *
*  munu        mu    mu nu     mu        munu  mu nu                                                        *
* T    =(p *h+b b )*U *U +[P+(b *b )/2]*g    -b  b                                                          *
*         o      mu               mu                                                                        *
* h =1+eps+P/p                                               //Specific Enthalpy                            *
*  o          o                                                                                             *
*                                                                                                           *
*  mu  mu    /                                                                       mu     mu              *
* b  =B     / sqrt(4PI)                        //Conversion of 4vector b   to B   as a function of 4velocity*
*      (UU)/                                                                                                *
************************************************************************************************************/
/*******************************************************************************************************************
* The following variables are defined as                                                                           *
* rho = density                                                                                                    *
* BB_mu = four vector magnetic field                                                                               *
* UU = four velocity                                                                                               *
* PP = pressure                                                                                                    *
* eps = specific internal energy (SIE)                                                                             *
*******************************************************************************************************************/
//4velocity of an observer in units of Solar_Mass^0
CCTK_REAL UU_mu[4]={i_vel0*Meter/Second,i_vel1*Meter/Second,i_vel2*Meter/Second,i_vel3*Meter/Second};
CCTK_REAL BB_mu[4]={i_bo0*Tesla,i_bo1*Tesla,i_bo2*Tesla,i_bo3*Tesla};    //Init 4magnetic field in units of Solar_Mass^-1
CCTK_REAL rho=i_rho*Density;                                  //Init Density of rho in Solar_Mass
CCTK_REAL PP=i_pressure*Pressure;                             //Init Pressure in units of Solar_Mass^-2
CCTK_REAL eps=i_eps*SIE;                                      //specific internal energy density
CCTK_REAL h=1+eps+PP/rho;                                     //add description
CCTK_REAL BB0=(UU_mu[1]*BB_mu[1]+UU_mu[2]*BB_mu[2]+UU_mu[3]*BB_mu[3]);//intermediate step in calculation of BB_ofU[mu]
CCTK_REAL lapse();
CCTK_REAL ULapse=(1/(UU_mu[0]*lapse()));                      //intermediate step in calculation of BB_ofU[mu]
CCTK_REAL BB_ofU[4]={                                         /*See Duez 420 eq 23,24*/
                    BB0/lapse(),                              /*because of cartesian coordinates there is no difference*/
                    ULapse*(BB_mu[1]+BB0*UU_mu[1]),           /*between contra/cov tensers so */
                    ULapse*(BB_mu[2]+BB0*UU_mu[2]),           /*b_cov_mu * b_con_mu = b^2*/
                    ULapse*(BB_mu[3]+BB0*UU_mu[3])            /**/
                    }; 
//CCTK_REAL find_Tmunu(TT_contra_munu);                         //this will find the value of 4X4 stress energy tensor
/*
CCTK_REAL time is in seconds\ 
CCTK_REAL z_dir in in meters |--with units of solar masses
CCTK_REAL plusA is in meters/
CCTK_REAL kk is the wave number and is unitless
*/

extern "C" void Calc_Init(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Set up shorthands

  CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  CCTK_REAL dz = CCTK_DELTA_SPACE(2);
  CCTK_REAL dt = CCTK_DELTA_TIME;

  CCTK_REAL dx2 = dx*dx;
  CCTK_REAL dy2 = dy*dy;
  CCTK_REAL dz2 = dz*dz;
  CCTK_REAL dt2 = dt*dt;

  CCTK_REAL dx2i = 1.0/dx2;
  CCTK_REAL dy2i = 1.0/dy2;
  CCTK_REAL dz2i = 1.0/dz2;

  CCTK_REAL factor = 2*(1 - (dt2)*(dx2i + dy2i + dz2i));

  int istart = 1;
  int jstart = 1;
  int kstart = 1;

  int iend = cctk_lsh[0]-1;
  int jend = cctk_lsh[1]-1;
  int kend = cctk_lsh[2]-1;

  for (int k=kstart; k<kend; k++)
  {
    for (int j=jstart; j<jend; j++)
    {
      for (int i=istart; i<iend; i++)
      {
      
      int vindex =  CCTK_GFINDEX3D(cctkGH,i,j,k);
      
      gxx[vindex] = GG_contra_munu[1][1]; 
      gxy[vindex] = GG_contra_munu[1][2];
      gxz[vindex] = GG_contra_munu[1][3];
      gyy[vindex] = GG_contra_munu[2][2];
      gyz[vindex] = GG_contra_munu[2][3];
      gzz[vindex] = GG_contra_munu[3][3];
 
 	/*We made the assumption that our shift will always be a constant.  
 	Thus the Lie derivative wrt the shift will produce a zero value.  
 	Future expansion can be done to eliminate this assumption in the extrinsic curvature.*/
 	
    CCTK_REAL wave_term = kk(time,z_dir)*sin(kk(time,z_dir)*t)*sin(kk(time,z_dir)*z_dir)/(2*lapse());
    
    //intermediate step      
    
      kxx[vindex] = -plusA(time,z_dir)*wave_term;        		// Duez 420 (2) Extrinsic Curvature
      kxy[vindex] = -crossB(time,z_dir)*wave_term;       		//               ||
      kxz[vindex] = 0;                                   		//               ||
      kyy[vindex] = plusA(time,z_dir)*wave_term;         		//               ||
      kyz[vindex] = 0;                                   		//               ||
      kzz[vindex] = 0;                                   		//               ||
    
      alp[vindex] = GG_contra_munu[0][0];   					//This is the lapse()
      
      betax[vindex] = GG_contra_munu[0][1]; 					//This is the shift()
      betay[vindex] = GG_contra_munu[0][2];
      betaz[vindex] = GG_contra_munu[0][3]; 					/* Next three term are intermediate terms used in the 
      															calculations of my primitive variables.*/      
      
      CCTK_REAL C2andP2 = crossB(time,z_dir)*crossB(time,z_dir) + plusA(time,z_dir)*plusA(time,z_dir);      
      CCTK_REAL sqr_gamma_det = pow(1-C2andP2*pow(cos(kk(time,z_dir)*t)*sin(kk(time,z_dir)*z_dir),2),.5);      
      CCTK_REAL temp =0;      
      
      for(int ii=1;ii<4;ii++){       							/*Now calculate the value of gamma_ij*u_i*u_j*/         
         for(int jj=1;jj<4;jj++){    							/*Store it in temp for use in finding uu0*/	    
           temp = temp + GG_contra_munu[ii][jj]*UU_mu[ii]*UU_mu[jj];         
         }      
      };               
      
      CCTK_REAL uu0 = (1+pow(temp,.5)/lapse();
      
      rho_star[vindex] = lapse()*sqr_gamma_det*rho*uu0;							//Duez 420 (34) from convservation of baryon numbers
      tau[vindex] = lapse()*lapse*sqr_gamma_det*TT_contra_munu[0][0]-rho_star;		//Duez 420 (37) energy variable
      
      S_x[vindex] = lapse()*sqr_gamma_det*TT_contra_munu[0][1];	//Duez 420 (35) momentum density variable
      S_y[vindex] = lapse()*sqr_gamma_det*TT_contra_munu[0][2];	//                ||
      S_z[vindex] = lapse()*sqr_gamma_det*TT_contra_munu[0][3];	//                ||
      
      bb_x[vindex] = sqr_gamma_det*BB_ofU[1];                  	//Duez 420 (27) no-monopole conatraint
      bb_y[vindex] = sqr_gamma_det*BB_ofU[2];                  	//                 ||
      bb_z[vindex] = sqr_gamma_det*BB_ofU[3];                  	//                 ||      
      
      vv_x[vindex] = 0;
      vv_y[vindex] = 0;
      vv_z[vindex] = 0;
      
      /*Initialize to zero for calculation*/ //Duez 420 (56) flow velocity of plasma normal to the observer        
      
      for(jj=1;jj<4;jj++){              							//Reset temp for use in vv
        vv_x[vindex] = vv_x[vindex]+GG_contra_munu[1][jj]*UU_mu[jj]/uu0 -shift();        
        vv_y[vindex] = vv_y[vindex]+GG_contra_munu[1][jj]*UU_mu[jj]/uu0 -shift();        
        vv_z[vindex] = vv_z[vindex]+GG_contra_munu[1][jj]*UU_mu[jj]/uu0 -shift();      
      }
      
      Ttt_mhd[vindex] = TT_contra_munu[0][0];
      Ttx_mhd[vindex] = TT_contra_munu[0][1];
      Tty_mhd[vindex] = TT_contra_munu[0][2];
      Ttz_mhd[vindex] = TT_contra_munu[0][3];
      Txx_mhd[vindex] = TT_contra_munu[1][1];
      Tyy_mhd[vindex] = TT_contra_munu[2][2];
      Tzz_mhd[vindex] = TT_contra_munu[3][3];
      Txy_mhd[vindex] = TT_contra_munu[1][2];
      Txz_mhd[vindex] = TT_contra_munu[1][3];
      Tyz_mhd[vindex] = TT_contra_munu[2][3];
      
      }//end of for (int i=istart; i<iend; i++)
    }//end of for (int j=jstart; j<jend; j++)
  }//end of for (int k=kstart; k<kend; k++)
}//end of extern "C" void Calc_Init(CCTK_ARGUMENTS)

void find_Tmunu(CCTK_REAL TT_contra_munu[4][4]); 
CCTK_REAL grav_wave(CCTK_REAL time, CCTK_REAL z_dir);
CCTK_REAL line_element(CCTK_REAL dx0, CCTK_REAL dx1, CCTK_REAL dx2, CCTK_REAL dx3, CCTK_REAL time); //dx_mu is determined elsewhere in the code
CCTK_REAL Y(CCTK_REAL time, CCTK_REAL exp_line_elem);                //The scale function
CCTK_REAL lapse();                                                   //coefficient of dx0^2
CCTK_REAL omicron(CCTK_REAL time, CCTK_REAL exp_line_elem);          //coefficient of dx1^2
CCTK_REAL iota(CCTK_REAL time, CCTK_REAL exp_line_elem);             //coefficient of dx2^2
CCTK_REAL upsilon(CCTK_REAL time, CCTK_REAL exp_line_elem);          //coefficient of dx3^2
CCTK_REAL shift();                                                   //the shift vector
CCTK_REAL plusA(CCTK_REAL time, CCTK_REAL z_dir, CCTK_REAL kk);      //Amplitude of the wave PLUS term
CCTK_REAL crossB(CCTK_REAL time, CCTK_REAL z_dir, CCTK_REAL kk);     //Amplitude of the wave CROSS term
CCTK_REAL kk(CCTK_REAL time, CCTK_REAL z_dir);                       //Wave number

                             
//*********************** code for GG_contr_munu ************************************************************
//GG_contr_munu = is the contravariant metric
CCTK_REAL GG_contra_munu[4][4]={ /* This initial contra metric is defined for t=0 and z=0*/
                            {-lapse(),                          0,                       0,             0},
                            { 0      , (omicron(0,qq)-plusA(0,0)),            -crossB(0,0),             0},
                            { 0      ,               -crossB(0,0), (iota(0,qq)+plusA(0,0)),             0},
                            { 0      ,                          0,                       0, upsilon(0,qq)} 
                           };

//***************************code for Tcontra_munu***********************************************************
void find_Tmunu(CCTK_REAL TT_contra_munu[4][4]){ //calculates full rank 2 Stress energy tensor
  for(int mu=0;mu<4;mu++){
    for(int nu=0;nu<4;nu++){
      TT_contra_munu[mu][nu]=( (rho*h+(1/(4*PI))*BB_ofU[mu]*BB_ofU[mu])*UU_mu[mu]*UU_mu[nu] 
                             + (PP+((1/(4*PI))*BB_ofU[mu]*BB_ofU[mu])/2)*GG_contra_munu[mu][nu]
                             - (1/(4*PI))*BB_ofU[mu]*BB_ofU[nu] );
	}//end of nu loop
  }//end of mu loop
}//end of void find_Tmunu
/*The plus and cross terms contains a dependence on the wave number but the wave number is a function so it is not explicitly passed in.*/
CCTK_REAL plusA(CCTK_REAL time, CCTK_REAL z_dir){                   //The plus term
 return 1.18e-4*sin(kk(time,z_dir)*z_dir)*cos(kk(time,z_dir)*time); //1.18e-4 amp from Duez 420 (79)
}

CCTK_REAL crossB(CCTK_REAL time, CCTK_REAL z_dir){                  //The cross term
 return 1.18e-4*sin(kk(time,z_dir)*z_dir)*cos(kk(time,z_dir)*time); //1.18e-4 amp from Duez 420 (79)}

CCTK_REAL kk(CCTK_REAL time, CCTK_REAL z_dir){                      //The wave number
 return 2*PI;                              //To ensure a computational domain covers two wavelengths 
}

CCTK_REAL lapse(){                         //The (coefficient of dx0)^2
 return 1;                                 //Minkowski Space
}

CCTK_REAL shift(){                          //The shift vector may be derived from the line element
 return 0;                                 //It is set equal to zero for the Minkowski metric
}

CCTK_REAL omicron(CCTK_REAL time, CCTK_REAL exp_lime_elem){    //The coefficient of dx1^2
 return 1*pow(Y(time,qq),2);                          //Minkowski Space
}

CCTK_REAL iota(CCTK_REAL time, CCTK_REAL exp_lime_elem){       //The coefficient of dx2^2
 return 1*pow(Y(time,qq),2);                          //Minkowski Space
}

CCTK_REAL upsilon(CCTK_REAL time, CCTK_REAL exp_lime_elem){    //The coefficient of dx3^2
 return 1*pow(Y(time,qq),2);                          //Minkowski Space
}

CCTK_REAL Y(CCTK_REAL time, CCTK_REAL exp_line_elem){  //This is the scaling function
 return pow(time,exp_line_elem);
}

CCTK_REAL line_element(CCTK_REAL time, CCTK_REAL dx0, CCTK_REAL dx1, CCTK_REAL dx2, CCTK_REAL dx3){
  CCTK_REAL ds_squared=0;
  ds_squared=-(lapse()*dx0)*(lapse()*dx0)             +omicron(time,qq)*(dx1)*(dx1)             +iota(time,qq)*(dx2)*(dx2)             +upsilon(time,qq)*(dx3)*(dx3);
 return ds_squared;
}
CCTK_REAL EEpsilon(){                 //specific internal energy
  return i_eps*Kilogram*Meter*Meter/(Second*Second);
}

/*
This is the gravity wave for perturbation in the plus and cross direction
*/
CCTK_REAL grav_wave(CCTK_REAL time,CCTK_REAL z_dir,CCTK_REAL kk){ 
 return plusA(time,z_dir,kk(time,z_dir))+crossB(time,z_dir,kk(time,z_dir));//linear addition of waves
}

/************************************************************************************************************
* int main is used for testing purposes only and should be commented out before incorporating this code in  *
* another program.JSH & CLB 22NOV2005                                                                       *
************************************************************************************************************/
/*
int main(void){
  CCTK_REAL TT_contra_munu[4][4];
  find_Tmunu(TT_contra_munu);
  cout<<TT_contra_munu[0][0]<<"\t\t"<<TT_contra_munu[0][1]<<"\t\t"<<TT_contra_munu[0][2]<<"\t\t"<<TT_contra_munu[0][3]<<'\n'
      <<TT_contra_munu[1][0]<<"\t\t"<<TT_contra_munu[1][1]<<"\t\t"<<TT_contra_munu[1][2]<<"\t\t"<<TT_contra_munu[1][3]<<'\n'
      <<TT_contra_munu[2][0]<<"\t\t"<<TT_contra_munu[2][1]<<"\t\t"<<TT_contra_munu[2][2]<<"\t\t"<<TT_contra_munu[2][3]<<'\n'
      <<TT_contra_munu[3][0]<<"\t\t"<<TT_contra_munu[3][1]<<"\t\t"<<TT_contra_munu[3][2]<<"\t\t"<<TT_contra_munu[3][3]<<'\n';
 return 1;
}*/
