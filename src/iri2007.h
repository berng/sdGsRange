#ifndef __IRI2007_H__
#define __IRI2007_H__
#define _FTRUE 1
#define _FFALSE 0

typedef int boolean;
typedef struct{				// true:1		 false:0			Recommended
boolean NeComputed;  	 	//Ne computed		 Ne not computed                     t
boolean TeTiComputed;    		//Te, Ti computed        Te, Ti not computed                 t	   
boolean NeNiComputed;	 	//Ne & Ni computed       Ni not computed                     t
boolean MagneticFieldByTable;	//B0 - Table option      B0 - Gulyaeva (1987)                t
boolean foF2_CCIR;		//foF2 - CCIR            foF2 - URSI                     false
boolean Ni_old;    		//Ni - DS-78 & DY-85     Ni - DS-95 & TTS-03             false
boolean Ne_tops;    		//Ne - Tops: f10.7<188   f10.7 unlimited                     t            
boolean foF2_model;    		//foF2 from model        foF2 or NmF2 - user input           t
boolean hmF2_model;    		//hmF2 from model        hmF2 or M3000F2 - user input        t
boolean TeStandard;    		//Te - Standard          Te - Using Te/Ne correlation        t
boolean NeStandard;    		//Ne - Standard Profile  Ne - Lay-function formalism         t
boolean MessagesToUnit6; 		//Messages to unit 6     no messages                         t
boolean foF1_model;    		//foF1 from model        foF1 or NmF1 - user input           t
boolean hmF1_model;    		//hmF1 from model        hmF1 - user input (only Lay version)t
boolean foE_model;    		//foE  from model        foE or NmE - user input             t
boolean hmE_model;    		//hmE  from model        hmE - user input                    t
boolean Rz12_from_file;  		//Rz12 from file         Rz12 - user input                   t
boolean IGRF_dip;  		//IGRF dip, magbr, modip old FIELDG using POGO68/10 for 1973 t
boolean F1_probability_model;	//F1 probability model   critical solar zenith angle (old)   t
boolean standard_F1;		//standard F1            standard F1 plus L condition        t
boolean IonDriftComputed;		//ion drift computed     ion drift not computed          false
boolean IonDensitiesInPercent;	//ion densities in %     ion densities in m-3                t
boolean Te_tops;			//Te_tops (Aeros,ISIS)   Te_topside (Intercosmos)        false
boolean Dregion_IRI95;		//D-region: IRI-95       Special: 3 D-region models          t
boolean F107D_from_AP_DAT;	//F107D from AP.DAT      F107D user input (oarr(41))         t
boolean foF2_storm_model;		//foF2 storm model       no storm updating                   t
boolean IG12_from_file;		//IG12 from file         IG12 - user input		     t
boolean spreadF_probability;	//spread-F probability 	 not computed                    false
boolean IRI01_topside;		//IRI01-topside          new options as def. by JF(30)   false
boolean IRI01_topside_corr;	//IRI01-topside corr.    NeQuick topside model   	 false 
//			     (IRI01_topside,IRI01_topside_corr) = (t,t) IRIold, 
//								  (f,t) IRIcor, 
//					    			  (f,f) NeQuick, 
//					    			  (t,f) TTS   
} iri_conf;

typedef struct{
       float NmF2_m3;        //user input for foF2/MHz or NmF2/m-3   (for iri_conf.foF2_model=false)
       float hmF2_km;	 //user input for hmF2/km or M(3000)F2   (for iri_conf.hmF2_model=false)
       float NmF1_m3;        //user input for foF1/MHz or NmF1/m-3   (.foF1_model=false)
       float hmF1_km;	 //user input for hmF1/km (.hmF1_model=false)
       float NmE_m3;         //user input for foE/MHz or NmE/m-3    (.foE_model=false)
       float hmE_km;		 //user input for hmE/km (.hmE_model=false)
       float NmD_m3;             
       float hmD_km;
       float h_half_km;
       float B0_km;
       float VALLEY_BASE_m3;
       float VALLEY_TOP_km;
       float TE_PEAK_K;          
       float TE_PEAK_height_km;
       float TE_MOD_300KM_K;     //user input for Ne(300km)/m-3. (for iri_conf.TeStandard=false) Use OARR()=-1 if one of these values is not available.
       float TE_MOD_400KM_K;   //user input for Ne(400km)/m-3. (for iri_conf.TeStandard=false) Use OARR()=-1 if one of these values is not available. If iri_conf.Te_tops=.false. then Ne(550km)/m-3.
       float TE_MOD_600KM_K;         
       float TE_MOD_1400KM_K;
       float TE_MOD_3000KM_K;   
       float TE_120KM_K;	//=TN=TI at 120km
       float TI_MOD_430KM;      
       float X_KM;		// WHERE TE=TI
       float SOL_ZENITH_ANG_DEG; 
       float SUN_DECLINATION_DEG;
       float DIP_deg;
       float DIP_LATITUDE_deg;
       float MODIFIED_DIP_LAT;  
       float DELA;
       float sunrise_dec_hours;
       float sunset_dec_hours;
       float ISEASON; 		//(1=spring) 
       float NSEASON; 		//(northern)
       float Rz12;             //user input for Rz12(for iri_conf.Rz12_from_file=false)
       float Covington_Index;
       float B1;
       float M_3000_F2;
       float TEC_m2;           
       float TEC_top_2_TEC100;	//TEC_top/TEC*100.
       float gind_IG12;       //gind(IG12), user input for IG12   (iri_conf.IG12_from_file=false)
       float F1_probability;	// (old)
       float F107_daily;       //user input for daily F10.7 index (
				    //      for
				    //     iri_conf.IonDriftComputed=true
				    //        or
				    //         iri_conf.Te_tops=false
				    //        or
				    //         iri_conf.Dregion_IRI95=false
				    //        or
				    //         iri_conf.F107D_from_AP_DAT=false
				    //	 
				    //	NOTE: if <0 then 12-month running mean is taken from internal file
				    //       )
       float c1; 		//(F1 shape)
       float daynr;
       float V_vert_eq;	//equatorial vertical ion drift in m/s
       float foF2_storm2foF2_quiet;// foF2_storm / foF2_quiet
       float F1_probability_without_L_condition;
       float F1_probability_with_L_condition_incl;
       float spreadF_occurrence_probability; //(Brazilian model)
//                # INPUT as well as OUTPUT parameter
//                $ special for IRIWeb (only place-holders)

    

}oarr_type;


typedef struct{
float Ne_m3;	// ELECTRON DENSITY/M-3
float Tn_K;	// NEUTRAL TEMPERATURE/K
float Ti_K;	// ION TEMPERATURE/K
float Te_K;	//ELECTRON TEMPERATURE/K
float O_comp;	// O+ ION DENSITY/% or /M-3 if jf(22)=f 
float H_comp;	//H+ ION DENSITY/% or /M-3 if jf(22)=f
float He_comp; //HE+ ION DENSITY/% or /M-3 if jf(22)=f
float O2_comp; //O2+ ION DENSITY/% or /M-3 if jf(22)=f
float NO_comp;	//NO+ ION DENSITY/% or /M-3 if jf(22)=f

//AND, IF iri_conf.Ni_old=.FALSE.:
float NiCluster; 	   //  CLUSTER IONS DEN/% or /M-3 if iri_conf.IonDensitiesInPercent=.false.
float Ni_Npos;	  	   //  N+ ION DENSITY/% or /M-3 if iri_conf.IonDensitiesInPercent=.false.
float NeDregion;  
float NeEregion; 	   //  D/E-region densities (Friedrich)
float StrongConditions;   //  1:7 Danilov for 60,65..90km; 8:14 for a
	                   //major Stratospheric Warming (SW=1) event; 15:21 
	                   //for strong Winter Anomaly (WA=1) conditions
float outf15; //free
float outf16; //free
float outf17; //free
float outf18; //free
float outf19; //free
float outf20; //free
}outf_type; //OUTF(1:20,1:100)
//typedef float outf2_type[20][100];
typedef float** outf2_type;

/*
C*****************************************************************
C*** THE ALTITUDE LIMITS ARE:  LOWER (DAY/NIGHT)  UPPER        ***
C***     ELECTRON DENSITY         60/80 KM       1000 KM       ***
C***     TEMPERATURES              120 KM        2500/3000 KM  ***
C***     ION DENSITIES             100 KM        1000 KM       ***
C*****************************************************************
C*****************************************************************
C*********            INTERNALLY                    **************
C*********       ALL ANGLES ARE IN DEGREE           **************
C*********       ALL DENSITIES ARE IN M-3           **************
C*********       ALL ALTITUDES ARE IN KM            **************
C*********     ALL TEMPERATURES ARE IN KELVIN       **************
C*********     ALL TIMES ARE IN DECIMAL HOURS       **************
C*****************************************************************
C*****************************************************************
C*****************************************************************
*/

/*
c !!!!! Subroutine INITIALIZE has to be called once before calling
c !!!!! iri_sub. This is already included in subroutine IRI_WEB which
c !!!!! calls iri_sub. 
C*****************************************************************
C********* INTERNATIONAL REFERENCE IONOSPHERE (IRI). *************
C*****************************************************************
C**************** ALL-IN-ONE SUBROUTINE  *************************
C*****************************************************************
*/

void initialize_();

void iri_sub_(iri_conf* JF,
	      int* JMAG,
	      float* ALATI,
	      float* ALONG,
	      int* IYYYY,
	      int* MMDD,
	      float* DHOUR, 
	      float* HEIBEG,
	      float* HEIEND,
	      float* HEISTP,
	      outf2_type OUTF,
	      oarr_type* OARR);

void iri_fast_(iri_conf* JF,
	      int* JMAG,
	      float* ALATI,
	      float* ALONG,
	      int* IYYYY,
	      int* MMDD,
	      float* DHOUR, 
	      float* HEIGHT,
	      void* OUTF,
	      oarr_type* OARR);



void iri_web_(int* jmag,
	      boolean* jf,
	      float* alati,
	      float* along,
	      int* iyyyy,
	      int* mmdd,
	      int* iut,
	      float* dhour,
	      float* height,
	      float* h_tec_max,
	      int* ivar,
	      float* vbeg,
	      float* vend,
	      float* vstp,
	      float* a,
	      float* b);
/*
c-----------------------------------------------------------------------        
c changes:
c       11/16/99 jf(30) instead of jf(17)
c
c-----------------------------------------------------------------------        
c input:   jmag,alati,along,iyyyy,mmdd,dhour  see IRI_SUB
c          height  height in km
c          h_tec_max  =0 no TEC otherwise upper boundary for integral
c          iut     =1 for UT       =0 for LT
c          ivar    =1      altitude
c                  =2,3    latitude,longitude
c                  =4,5,6  year,month,day
c                  =7      day of year
c                  =8      hour (UT or LT)
c          vbeg,vend,vstp  variable range (begin,end,step)
c output:  a       similar to outf in IRI_SUB
c          b       similar to oarr in IRI_SUB
c
c          numstp  number of steps; maximal 100
c-----------------------------------------------------------------------        
*/
#endif
