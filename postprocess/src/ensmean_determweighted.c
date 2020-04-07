# define _BSD_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <sys/types.h>
#include <stdint.h>
#include <ctype.h>
#include <hdf5.h>
#include <hdf5_hl.h>

#define DEFSETS 1000
#define UNIT_DBZ 0
#define UNIT_RATE 1

double ZR_A=223.0, ZR_B=1.53;


int H5get_variable_string(hid_t h5,char *datasetname, char *attrname, char *str);
uint8_t RfromdBZN(uint8_t dBZN);
uint8_t dBZNfromRi(uint16_t Ri, double gain, double offset);

int main(int argc, char *argv[])
{

   /* HDF5 variables */
  hid_t INH5,OUTH5,whatg,whereg,howg,datasetg,datag,datawhatg,attr;
  hid_t setwhatg,setwhereg,sethowg,dataset,plist,space,datatype,memtype,filetype;
  hsize_t dims[2],chunkdims[2]={100,100};
  size_t typesize,arrsize,databytes;
  
  int unit=UNIT_DBZ;
  int timesteps,max_leadtime;
  double determ_initweight=100.0; /* Initial weight of deterministic nowcast, percent of amount of members */
  /* E.g. if iw=100, the weight of deterministic is equal to all memebers together */
  double determ_weightspan=100.0; /* The weighting time span for deterministic nowcast, percent of max leadtime */
  /* E.g. if time is 200% the determ weight has dropped to half at the time of last timestep,
     if 50%, the weight drops to zero at the time of middle timestep. If negative, determ w is kept constant */  
  double determ_startw, determ_lapse;
  double determ_w; /* Actual weight of determ nowcast as a function of leadtime */
  int DETERM = 0, DETERM_W = 1; /* if deterministic "member" */
  double RdBZN[256],K=1.0;
  double xscale,yscale,*sumdata,fval,RATE,meanR;
  /* double LL_lon,LL_lat,UL_lon,UL_lat,UR_lon,UR_lat,LR_lon,LR_lat; */
  double gain=1.0,offset=0.0,fnodata=65535.0,fundetect=0.0,goffset;
  double norain_dBZ; /* from attribute /meta/configuration/NORAIN_VALUE, PPN internal for undetected echo */
  uint8_t dBZNRi[65536],undetect;
  double *countdata;
  long long int longtime;
  long xsize,ysize,N,Ri;
  uint16_t *membdata,*meandata,meanRi,nodata=65535,val;
  long dBZNi;
  uint8_t *dBZdata,dBZN, DBZ_DATA=0, RATE_DATA=0, MEANDBZ=0, MEANR=0;
  char varstr[200];
  char projdef[200];
  char timestamp[15]={0};
  char HDFoutfile[300];
  char unitstr[20];
  herr_t ret,status;
  int i,FIRST=1,len;
  H5T_class_t class_id;
  char datagname[100];
  char pgmname[100];
  char datasetname[100],validtimestr[50]={0};
  int membI=0,leadtI=0,leadtimes=DEFSETS,members=DEFSETS;
  FILE *PGMF;
  /*   char membstr[5]={0},leadtstr[5]={0}; */



  setbuf(stdout,NULL);
  H5Eset_auto(H5E_DEFAULT,NULL,NULL);

  /* get environment */
  {
    char *p;

    if(p=getenv("DETERM_INITWEIGHT")) determ_initweight = atof(p);
    if(p=getenv("DETERM_WEIGHTSPAN")) determ_weightspan = atof(p);
    if(p=getenv("GENERATE_MEANDBZ")) if(strcasecmp(p,"TRUE")==0) MEANDBZ=1;
    if(p=getenv("GENERATE_MEANR")) if(strcasecmp(p,"TRUE")==0) MEANR=1;
  }

  if(!(MEANR | MEANDBZ)) 
  {
    printf("GENERATE_MEANDBZ nor GENERATE_MEANR set TRUE, enable either!\n");
    return(1);
  }

  INH5=H5Fopen(argv[1], H5F_ACC_RDONLY, H5P_DEFAULT);
  datatype=H5T_STD_U16LE;

  if(H5get_variable_string(INH5,"/meta/configuration","ENSEMBLE_SIZE",varstr) >= 0)
     members = atoi(varstr);
  if(H5get_variable_string(INH5,"/meta/configuration","NUM_TIMESTEPS",varstr) >= 0)
     timesteps = atoi(varstr);
  if(H5get_variable_string(INH5,"/meta/configuration","MAX_LEADTIME",varstr) >= 0)
     max_leadtime = atoi(varstr);
  if(H5get_variable_string(INH5,"/meta/configuration","NORAIN_VALUE",varstr) >= 0)
     norain_dBZ = atoi(varstr);
  if(H5get_variable_string(INH5,"/meta/configuration","ZR_A",varstr) >= 0)
     ZR_A = atof(varstr);
  if(H5get_variable_string(INH5,"/meta/configuration","ZR_B",varstr) >= 0)
     ZR_B = atof(varstr);
  /* printf("%d\n",RitodBZN(1,1.0)); */

  determ_startw = 0.01*determ_initweight * (double)members;
  if(determ_weightspan <= 0.0) determ_lapse=0.0;
  else determ_lapse = 100.0*determ_startw/(determ_weightspan*(double)timesteps);
  /* determ_w = determ_startw - determ_lapse * (double)leadtI */ 

  /* R from dBZN LUT */
  {
    double R=0.0,B,rcW;
    int N;

    B = 0.05/ZR_B;
    rcW = -(3.2 + log10(ZR_A))/ZR_B;  

    RdBZN[0]=0;
    RdBZN[255]=-1e-6;
    for(N=1;N<255;N++) RdBZN[N]=pow(10.0,B*(double)N + rcW);
  }


  DETERM_W = 1;
  for(leadtI=0;leadtI<timesteps;leadtI++) 
  {
     determ_w = determ_startw - determ_lapse * (double)leadtI; /* weight of deterministic "member" */ 
     if(determ_w < 0.0) DETERM_W = 0; else printf("Determ weight = %.1f\n",determ_w);
     printf("Lead=%02d\nMembers: ",leadtI);
  
     DETERM=0;
     for(membI=0;membI<=members;membI++)
     {
         if(membI==members) DETERM = 1;

         if(DETERM) sprintf(datasetname,"/deterministic/leadtime-%02d",leadtI);
         else sprintf(datasetname,"/member-%02d/leadtime-%02d",membI,leadtI);

         if(FIRST) 
	 {
	   H5LTget_dataset_info(INH5, datasetname, dims, &class_id, &typesize );
	   /* printf("%ld %ld\n",xsize,ysize); */
	   xsize=(long)dims[1];
           ysize=(long)dims[0];
           arrsize=xsize*ysize;
           databytes=arrsize*typesize;
	   meandata=malloc(databytes);
           memset(meandata,255,databytes);
           dBZdata=malloc(arrsize);
           countdata=malloc(arrsize*sizeof(double));
           membdata=malloc(databytes);
           sumdata=malloc(arrsize*sizeof(double));
           H5LTget_attribute_double(INH5,datasetname,"gain",&gain);
           H5LTget_attribute_double(INH5,datasetname,"nodata",&fnodata);
           H5LTget_attribute_double(INH5,datasetname,"offset",&offset);
           goffset=gain/offset;
           undetect = (uint8_t)((norain_dBZ - offset)/gain);
	   /*           printf("Undetect %d\n",undetect); */
           len=H5get_variable_string(INH5,"meta","nowcast_units",unitstr);
	   /*           printf("varstr: %d %s\n",len,unitstr); */
           if(strstr(unitstr,"dBZ")) DBZ_DATA=1;
           if(strstr(unitstr,"rrate")) RATE_DATA=1;

           nodata=(uint16_t)(fnodata+1e-9);
	   for(Ri=0;Ri<=65535;Ri++) dBZNRi[Ri]=dBZNfromRi(Ri,gain,offset);
	   /*           printf("FIRST %d  gain %f, offset %f, DBZ_DATA %d  RATE_DATA %d\n",FIRST,gain,offset,DBZ_DATA,RATE_DATA); */
           FIRST=0;
 	 }

         ret=H5LTget_attribute_info(INH5, datasetname, "Valid for", dims, &class_id, &typesize);
	 if(class_id == H5T_INTEGER) 
	 {
	    ret=H5LTget_attribute_long_long(INH5, datasetname, "Valid for", &longtime);
	    /*
	    dateint=longtime/1000000;
	    timeint=longtime%1000000;
	    */
	    sprintf(timestamp,"%lld",longtime/100);
	 }

	 if(!membI)
	 { 
	    memset(sumdata,0,arrsize*sizeof(double));
	    memset(countdata,0,arrsize*sizeof(double));
	 }
	 H5LTread_dataset(INH5,datasetname,datatype,membdata);
	 memset(dBZdata,255,arrsize);
	 for(N=0;N<arrsize;N++)
	 {
	   val=membdata[N];
	   if(val==nodata) continue;
	   fval=(double)val;
	   if(RATE_DATA)
	   { 
	      RATE = fval*gain+offset;
	      dBZN=dBZNRi[val];
              if(val==undetect) dBZN=0;
	   }

	   if(DBZ_DATA)
	   { 
	      dBZNi=2*(fval*gain + offset)+64;
	      dBZN=(uint8_t)dBZNi;
	      if(dBZNi>=255) dBZN=255; 
	      if(dBZNi<0 || val==undetect) dBZN=0; 
	   }
	   dBZdata[N]=dBZN;

	   if(DBZ_DATA) RATE = RdBZN[dBZN];

	   if(DETERM)
	   { 
	      if(DETERM_W)
	      {
	         sumdata[N] += RATE * determ_w;
	         countdata[N] += determ_w;
	      }
	   }
	   else 
	   {
	      sumdata[N] += RATE;
	      countdata[N] += 1.0;
	   }
	 }
	 printf("%02d ",membI);

# if 0
	 if(!membI) /* Write member zero */
	 {
	    sprintf(pgmname,"membR_%s_%02d.pgm",timestamp,membI);
	    PGMF=fopen(pgmname,"w");
	    fprintf(PGMF,"P5\n%ld %ld\n65535\n",xsize,ysize);
	    swab(membdata,membdata,databytes);
	    fwrite(membdata,1,databytes,PGMF);
	    fclose(PGMF);

	    sprintf(pgmname,"memb0_dBZ_%s.pgm",timestamp);
	    PGMF=fopen(pgmname,"w");
	    fprintf(PGMF,"P5\n%ld %ld\n255\n",xsize,ysize);
	    fwrite(dBZdata,1,arrsize,PGMF);
	    fclose(PGMF);
	 }


      /* Write deterministic */
	 if(DETERM)
         {
            printf("Deterministic lead=%02d %s %s\n",leadtI,validtimestr,timestamp);
            sprintf(pgmname,"determdBZ_%s.pgm",timestamp);
            PGMF=fopen(pgmname,"w");
            fprintf(PGMF,"P5\n%ld %ld\n255\n",xsize,ysize);
            fwrite(dBZdata,1,arrsize,PGMF);
            fclose(PGMF);
         }
# endif

     }


     printf("\nMean lead=%02d %s %s\n",leadtI,validtimestr,timestamp);

     /* Get unperturbed */
     sprintf(datasetname,"/unperturbed/leadtime-%02d",leadtI);
     ret=H5LTget_attribute_info(INH5, datasetname, "Valid for", dims, &class_id, &typesize);
     if(ret >= 0) 
     {
        H5LTread_dataset(INH5,datasetname,datatype,membdata);
        memset(dBZdata,255,arrsize);
        for(N=0;N<arrsize;N++)
        {
           if(membdata[N]==nodata) continue;
	   if(RATE_DATA) dBZdata[N]=dBZNRi[membdata[N]];
           else dBZdata[N]=(uint8_t)(2*((double)membdata[N]*gain+offset)+64);
        }
        printf("Unperturbed lead=%02d %s %s\n",leadtI,validtimestr,timestamp);
        sprintf(pgmname,"unpertdBZ_%s.pgm",timestamp);
        PGMF=fopen(pgmname,"w");
        fprintf(PGMF,"P5\n%ld %ld\n255\n",xsize,ysize);
        fwrite(dBZdata,1,arrsize,PGMF);
        fclose(PGMF);
     }

     /* Calculate and write ensemble mean */
     for(N=0;N<arrsize;N++) 
     {
        meanR = sumdata[N]/countdata[N];
        if(meanR<655.36 && meanR>=0.0) meanRi = (uint16_t)(meanR*100);
        
        meandata[N]=meanRi;
        if(MEANDBZ)
	{
           dBZN = dBZNRi[meanRi];
           dBZdata[N]=dBZN;
	}
     }
     if(MEANR)
     {
        swab(meandata,meandata,databytes);
        sprintf(pgmname,"meanR_%s.pgm",timestamp);
        PGMF=fopen(pgmname,"w");
        fprintf(PGMF,"P5\n%ld %ld\n65535\n",xsize,ysize);
        fwrite(meandata,1,databytes,PGMF);
        fclose(PGMF);
     }

     if(MEANDBZ)
     {
        sprintf(pgmname,"meandBZ_%s.pgm",timestamp);
        PGMF=fopen(pgmname,"w");
        fprintf(PGMF,"P5\n%ld %ld\n255\n",xsize,ysize);
        fwrite(dBZdata,1,arrsize,PGMF);
        fclose(PGMF);
     }
  }

  
  free(membdata);
  free(meandata);
  free(sumdata);
  free(dBZdata);
  free(countdata);
  H5Fclose(INH5);

  return(0);
}

/*==================================================================================================*/

uint8_t dBZNfromRi(uint16_t Ri, double gain, double offset)
{
    double R,dBZ;
    int16_t dBZI;

    if(!Ri) dBZI=0;
    else
    {
       if(Ri==65535) dBZI=255;
       else
       {
	  R=gain*(double)Ri+offset;
	  R=0.01*(double)Ri;
          dBZ=10.0*log10(ZR_A*pow(R,ZR_B));
          if(dBZ > -32.0) dBZI=(int)(2.0*dBZ)+64; else dBZI=0;
          if(dBZI > 254) dBZI=254;
       }
    }

    /*    printf("%d=%d ",R,dBZI); */
    return((uint8_t)dBZI);
}


double  dBZNtoR(uint8_t dBZN)
{
  double R=0.0,B,rcW;

  if(dBZN)
  {
     B = 0.05/ZR_B;
     rcW = -(3.2 + log10(ZR_A))/ZR_B;  
     R=pow(10.0,B*dBZN + rcW);
  }
  return(R);
}

/*--------------------------------------------------------------------------------------------------*/

int H5get_variable_string(hid_t h5,char *datasetname, char *attrname, char *str)
{
   hid_t dataset,attr,memtype,space;
   hvl_t  rdata;             /* Pointer to vlen structures */
   char *ptr;
   int i,len;
   herr_t status;

   memset(str,0,strlen(str));
   dataset=H5Dopen(h5,datasetname,H5P_DEFAULT);
   if(dataset<0) dataset=H5Gopen(h5,datasetname,H5P_DEFAULT);
   if(dataset<0) return(-1);
   attr = H5Aopen(dataset, attrname, H5P_DEFAULT);

   space = H5Aget_space(attr);
   memtype = H5Tvlen_create(H5T_NATIVE_CHAR);
   status = H5Aread(attr, memtype, &rdata);
   if(status<0) return(-1);
   ptr = rdata.p;
   len = rdata.len;
   for (i=0; i<len; i++) str[i]=ptr[i];
   str[i]=0;
   status = H5Dvlen_reclaim (memtype, space, H5P_DEFAULT, &rdata);
   status = H5Aclose (attr);
   status = H5Dclose (dataset);
   status = H5Sclose (space);
   status = H5Tclose (memtype);
   return(len);
}



   
   


# if 0

     H5LTget_attribute_long(INH5,"where","xsize",&xsize); 
     H5LTget_attribute_long(INH5,"where","ysize",&ysize);
     H5LTget_attribute_double(INH5,"where","xscale",&xscale);          
     H5LTget_attribute_double(INH5,"where","yscale",&yscale);
     H5LTget_attribute_double(INH5,"where","LL_lon",&LL_lon);
     H5LTget_attribute_double(INH5,"where","LL_lat",&LL_lat);
     H5LTget_attribute_double(INH5,"where","UL_lon",&UL_lon);
     H5LTget_attribute_double(INH5,"where","UL_lat",&UL_lat);
     H5LTget_attribute_double(INH5,"where","UR_lon",&UR_lon);
     H5LTget_attribute_double(INH5,"where","UR_lat",&UR_lat);
     H5LTget_attribute_double(INH5,"where","LR_lon",&LR_lon);
     H5LTget_attribute_double(INH5,"where","LR_lat",&LR_lat);
     H5LTget_attribute_string(INH5,"where","projdef",projdef);
     H5Fclose(INH5);
  }

  arrsize=xsize*ysize;
  membsize=arrsize*sizeof(int); 
  dims[0]=(hsize_t)xsize;
  dims[1]=(hsize_t)ysize;
  datatype= H5T_STD_U8LE;







   sprintf(HDFoutfile,"%s/%s_%03d-%03d_PACC.h5",outdir,obstime,Accmins[p],Interval[p]);
    OUTH5 = H5Fcreate(HDFoutfile,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
    H5LTset_attribute_string(OUTH5,"/","Conventions","ODIM_H5/V2_1");

    dataspace=H5Screate_simple(2, dims, NULL);
    plist = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(plist, 2, chunkdims);
    H5Pset_deflate( plist, 6);

    whatg=H5Gcreate2(OUTH5,"what",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    whereg=H5Gcreate2(OUTH5,"where",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    howg=H5Gcreate2(OUTH5,"how",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

    H5LTset_attribute_string(OUTH5,"what","date",datadate);
    H5LTset_attribute_string(OUTH5,"what","time",datatime);
    H5LTset_attribute_string(OUTH5,"what","object","COMP");
    H5LTset_attribute_string(OUTH5,"what","source","ORG:EFKL");
    H5LTset_attribute_string(OUTH5,"what","version","H5rad 2.1");

    H5LTset_attribute_double(OUTH5,"where","LL_lon",&LL_lon,1);
    H5LTset_attribute_double(OUTH5,"where","LL_lat",&LL_lat,1);
    H5LTset_attribute_double(OUTH5,"where","UL_lon",&UL_lon,1);
    H5LTset_attribute_double(OUTH5,"where","UL_lat",&UL_lat,1);
    H5LTset_attribute_double(OUTH5,"where","UR_lon",&UR_lon,1);
    H5LTset_attribute_double(OUTH5,"where","UR_lat",&UR_lat,1);
    H5LTset_attribute_double(OUTH5,"where","LR_lon",&LR_lon,1);
    H5LTset_attribute_double(OUTH5,"where","LR_lat",&LR_lat,1);
    H5LTset_attribute_long(OUTH5,"where","xsize",&xsize,1);
    H5LTset_attribute_long(OUTH5,"where","ysize",&ysize,1);
    H5LTset_attribute_double(OUTH5,"where","xscale",&xscale,1);
    H5LTset_attribute_double(OUTH5,"where","yscale",&yscale,1);
    H5LTset_attribute_string(OUTH5,"where","projdef",projdef);

   for(f=0;f<Fields[p];f++)
    {
        /* Define dataset start and end times */
        for(i=0;i<2;i++)
        {
           Forinds[i]=AccInds[i][p][f];
           formins=Forinds[i]*5;
           secs=bsecs+formins*60;
           date_from_sec(fortime,secs);
           sprintf(fielddate[i],"%.8s",fortime);
           sprintf(fieldtime[i],"%.4s00",fortime+8);
        }

        /* HDF5 dataset attributes */
        sprintf(datasetname,"dataset%d",f+1);
        printf("dataset%d\n",f+1);
        datasetg=H5Gcreate2(OUTH5,datasetname,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        setwhatg=H5Gcreate2(datasetg,"what",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        H5LTset_attribute_string(datasetg,"what","product","COMP");
        H5LTset_attribute_string(datasetg,"what","quantity","PROB");
        H5LTset_attribute_double(datasetg,"what","gain",&gain,1);
        H5LTset_attribute_double(datasetg,"what","offset",&offset,1);
        H5LTset_attribute_double(datasetg,"what","nodata",&nodata,1);
        H5LTset_attribute_double(datasetg,"what","undetect",&undetect,1);
        H5LTset_attribute_string(datasetg,"what","startdate",fielddate[0]);
        H5LTset_attribute_string(datasetg,"what","starttime",fieldtime[0]);
        H5LTset_attribute_string(datasetg,"what","enddate",fielddate[1]);
        H5LTset_attribute_string(datasetg,"what","endtime",fieldtime[1]);

       /* Looping thru thresholds and writing data to dataset f for each threshold */
        for(thri=1;thri<=Thresholds[p];thri++)
        {
          sprintf(datagname,"data%ld",thri);        
          printf("\tdata%ld\n",thri);        
          datag=H5Gcreate2(datasetg,datagname,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
          datawhatg=H5Gcreate2(datag,"what",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
          H5LTset_attribute_long(datag,"what","threshold_id",&thri,1);
          H5LTset_attribute_double(datag,"what","threshold_value",&f_acclims[p][thri],1);
 
          /* create new empty dataset and attributes to destination file */ 
          dataset = H5Dcreate2(datag,"data", datatype, dataspace,
                    H5P_DEFAULT, plist, H5P_DEFAULT);
          H5LTset_attribute_string( datag, "data", "CLASS", "IMAGE");
          H5LTset_attribute_string( datag, "data", "IMAGE_VERSION", "1.2");
          H5Dwrite(dataset, datatype, H5S_ALL, dataspace, H5P_DEFAULT,&out_fracarr[thri][0]);
          H5Dclose(dataset);
          H5Gclose(datawhatg);
          H5Gclose(datag);
        }
        H5Gclose(setwhatg);
        H5Gclose(datasetg);
    }
    H5Gclose(whatg);
    H5Gclose(whereg);
    H5Gclose(howg);
  }
  H5Fclose(OUTH5);

# endif
