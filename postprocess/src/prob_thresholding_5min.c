#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdint.h>
#include <projects.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <gdal.h>
#include <cpl_string.h>
#include <ogr_srs_api.h>

#define ACCLIMS 32

/* static double f_acclims[TSTEPS][ACCLIMS],pdiv; */
static double *f_acclims[ACCLIMS],pdiv;
static int32_t Thresholds[100]={0};
static int Fields[100]={0},Periods=0,Interval[100],Accmins[100];
/* static int min5flags[TSTEPS]={0}, **AccInds[2]; */
static int *min5flags, **AccInds[2];
/* static int32_t *ravacc[TSTEPS][MEMBERS],acclims[TSTEPS][ACCLIMS],submems; */
static int32_t ***ravacc,*acclims[ACCLIMS],submems;
static PJ *projref;
static projUV SWgeo,SWcart,NEgeo,NEcart,NWgeo,NWcart,SEgeo,SEcart;
static char *Area;
static int members,timesteps;

static int compare(const void *i1, const void *i2) 
           {
            int32_t *ci1 = (int32_t *) i1;
            int32_t *ci2 = (int32_t *) i2;
            return (*ci1 < *ci2);
	   }

void get_fractiles_per_accs(int32_t *membarr, unsigned char *fracarr, int period);  
void date_from_sec(char *preddate, time_t secs);
time_t sec_from_date(char *stim);

int main(int argc, char *argv[])
{
  int32_t X,Y,A,A1,A0,R1,R2,i,m,*membarr_orig,arrsize,N,a,PASS,M;
  int32_t *membarr,count_zero,MaxThresholds=0,p,f;
  int32_t membsize,formins,Forinds[2];
  int GEOTIFF_OUT=0,PGM_OUT=0,NODETERM=0;
  int aI;
  float dR;
  long thri;
  unsigned char *fracarr,**out_fracarr,*monacc;
  char accname[200],outname[200]={0};
  char moname[300];
  char obstime[13],fortime[13],fmins[4],*outdir,*accdir;
  char datadate[9],datatime[7];
  char fielddate[2][13]={{0}},fieldtime[2][13]={{0}};
  char *accpref=NULL,def_accpref[10]="ACC";
  char timezone_str[100]="UTC";
  char source[200]={0},nodes[10000]={0};
  time_t bsecs,secs;
  FILE *OUTF;
  FILE *MONF;

  /* HDF5 variables */
  hid_t OUTH5,whatg,whereg,howg,datasetg,datag,datawhatg;
  hid_t setwhatg,setwhereg,sethowg,dataset,plist,dataspace,datatype;
  hsize_t dims[2],chunkdims[2]={100,100};
  double xscale,yscale;
  double LL_lon,LL_lat,UL_lon,UL_lat,UR_lon,UR_lat,LR_lon,LR_lat;
  double gain=1.0,offset=0.0,nodata=255.0,undetect=0.0;
  long xsize,ysize,k;
  char projdef[200];
  char HDFoutfile[300];
  herr_t ret;
  char datagname[100],datasetname[100];

  /* GDAL variables for GeoTIFF output */
  GDALDriverH hDriver;
  double adfGeoTransform[6];
  int GDT_datatype;
  char **papszOptions = NULL;
  OGRSpatialReferenceH hSRS;
  char *pszSRS_WKT = NULL;
  int Zlevel=6;
  char TIFFcompression[200]="DEFLATE",sZlevel[3]="6";
  int EPSG=0;


  setbuf(stdout,NULL);

  if(argc<7)
  {
    printf("Arguments are:\n1. Timestamp YYYYmmddHHMM\n");
    printf("2. Precipitation thresholds configuration file\n");
    printf("3. Directory of interpolated accumulation files (ravake-style *.dat)\n");
    printf("   e.g. prefix_201405021035-201405021515+280.dat\n");
    printf("4. Output directory\n");
    printf("5. Area domain\n");
    printf("6. Member count, including deterministic member if stored\n");
    printf("7. Timesteps of ppn output file\n");
    return(1);
  }

  sprintf(obstime,"%.12s",argv[1]);
  sprintf(datadate,"%.8s",obstime);
  sprintf(datatime,"%.4s",obstime+8); 
  bsecs=sec_from_date(obstime);

  accdir=argv[3]; /* RAVACC_201405021035-201405021515+280.dat */
  outdir=argv[4];
  Area=argv[5];
  members=atoi(argv[6]);
  timesteps=atoi(argv[7]);

  {
    char *ptr;

    if(ptr=getenv("PROB_EPSG")) EPSG=atoi(ptr);
    if(ptr=getenv("PROB_IGNORE_DETERM")) if(strcasecmp(ptr,"TRUE")==0)
       { NODETERM=1; printf("IGNORING DETERMINISTIC MEMBER\n"); }
    if(ptr=getenv("PROB_OUTPUT_GEOTIFF")) if(strcasecmp(ptr,"TRUE")==0) GEOTIFF_OUT=1;
    if(ptr=getenv("PROB_OUTPUT_PGM")) if(strcasecmp(ptr,"TRUE")==0) PGM_OUT=1;
    accpref=getenv("INTERP_NC_ACCPREF");
    if(!accpref)
    {
       printf("Set INTERP_NC_ACCPREF for prefix used in interpolation output.\n");
       return(2);
    }
  }

  /* Amount of members must include deterministic member if present in acc data, 
     but if it's ignored, then bypass the member #0 */
     members -= NODETERM;


  for(aI=0;aI<ACCLIMS;aI++)
  {
     f_acclims[aI]=calloc(timesteps,sizeof(double));
     acclims[aI]=calloc(timesteps,sizeof(int32_t));
  }
  min5flags=calloc(timesteps,sizeof(int));

  /* Geometry initialization */
  {
     char *hdrfile;
     hid_t hdrh5;

     hdrfile=getenv("ODIM_HDF5_TMPLFILE");
     if(!hdrfile) exit(1);

     H5Eset_auto(H5E_DEFAULT,NULL,NULL);
     hdrh5=H5Fopen(hdrfile, H5F_ACC_RDONLY, H5P_DEFAULT);

     H5LTget_attribute_long(hdrh5,"where","xsize",&xsize); 
     H5LTget_attribute_long(hdrh5,"where","ysize",&ysize);
     H5LTget_attribute_double(hdrh5,"where","xscale",&xscale);          
     H5LTget_attribute_double(hdrh5,"where","yscale",&yscale);
     H5LTget_attribute_double(hdrh5,"where","LL_lon",&LL_lon);
     H5LTget_attribute_double(hdrh5,"where","LL_lat",&LL_lat);
     H5LTget_attribute_double(hdrh5,"where","UL_lon",&UL_lon);
     H5LTget_attribute_double(hdrh5,"where","UL_lat",&UL_lat);
     H5LTget_attribute_double(hdrh5,"where","UR_lon",&UR_lon);
     H5LTget_attribute_double(hdrh5,"where","UR_lat",&UR_lat);
     H5LTget_attribute_double(hdrh5,"where","LR_lon",&LR_lon);
     H5LTget_attribute_double(hdrh5,"where","LR_lat",&LR_lat);
     H5LTget_attribute_string(hdrh5,"where","projdef",projdef);
     H5LTget_attribute_string(hdrh5,"how","nodes",nodes);
     H5LTget_attribute_string(hdrh5,"what","source",source);
     H5Fclose(hdrh5);
     projref=pj_init_plus(projdef);

     SWgeo.u=LL_lon*DEG_TO_RAD;
     SWgeo.v=LL_lat*DEG_TO_RAD;
     NEgeo.u=UR_lon*DEG_TO_RAD;
     NEgeo.v=UR_lat*DEG_TO_RAD;
     NWgeo.u=UL_lon*DEG_TO_RAD;
     NWgeo.v=UL_lat*DEG_TO_RAD;
     SEgeo.u=LR_lon*DEG_TO_RAD;
     SEgeo.v=LR_lat*DEG_TO_RAD;

     SWcart=pj_fwd(SWgeo,projref);
     NEcart=pj_fwd(NEgeo,projref);
     NWcart=pj_fwd(NWgeo,projref);
     SEcart=pj_fwd(SEgeo,projref);
  }

  arrsize=xsize*ysize;
  membsize=arrsize*sizeof(int32_t); 
  monacc=malloc(arrsize);
  dims[0]=(hsize_t)ysize;
  dims[1]=(hsize_t)xsize;
  datatype= H5T_STD_U8LE;


  /* GDAL initialization for GeoTIFF output */
  if(GEOTIFF_OUT)
  {
     adfGeoTransform[0] = NWcart.u;
     adfGeoTransform[3] = NWcart.v;
     adfGeoTransform[1] = xscale;
     adfGeoTransform[5] = -yscale;
     adfGeoTransform[2] = 0;
     adfGeoTransform[4] = 0;

     GDALAllRegister();

     GDT_datatype=1; /* uint8_t */
     hDriver = GDALGetDriverByName("GTiff");

     if(strcmp(TIFFcompression,"NONE")) 
     {
        papszOptions = CSLSetNameValue( papszOptions,"COMPRESS",TIFFcompression);
        if(Zlevel) papszOptions = CSLSetNameValue( papszOptions,"ZLEVEL",sZlevel);
     }
     hSRS = OSRNewSpatialReference( NULL );
     if(EPSG) OSRImportFromEPSG(hSRS,EPSG); 
     else OSRImportFromProj4(hSRS,projdef);
     OSRExportToWkt( hSRS, &pszSRS_WKT );
  }

  { /* parse accumulation thresholds and read the data */
    char thr[200],*thrline,*pc,accstr[10];
    char lims[3]={'+','/',':'},dums[3][10]={{0}};
    char fortime[13]={0};
    char accname[200];
    int thri,len,i,p,f,vals[3],period=0,formins,forind;
    int32_t rval;
    time_t secs,bsecs,forsecs,intsecs,accsecs;
    FILE *CFGF,*NC_ACCF;

    MaxThresholds=0;
    CFGF=fopen(argv[2],"r");
    while(1)
    {
      memset(thr,0,200);    
      fgets(thr,199,CFGF);
      if(feof(CFGF)) break;
      if(thr[0]=='#') continue;
      printf("Period %d config %s\n",period+1,thr);
      thrline=thr;
      for(i=0;i<3;i++)
      {
         pc=strchr(thrline,lims[i]);
         vals[i]=atoi(strncpy(dums[i],thrline,pc-thrline));
         /* printf("%d %d %s\n",i,vals[i],dums[i]); */
         thrline=pc+1;
      }

      thri=1;
      while(1)
      {
        pc=strchr(thrline,',');
        if(pc) len=pc-thrline; else len=strlen(thrline); 
        memset(accstr,0,10);
        strncpy(accstr,thrline,len);
        /* printf("%s\t",accstr); */
        f_acclims[thri][period]=atof(accstr);
        acclims[thri][period]=(int32_t)(100000.0000001*f_acclims[thri][period]);
        printf("%d %d %f %ld\n",thri,period,f_acclims[thri][period],acclims[thri][period]); 
        if(!pc) break;
        thrline=pc+1;
        thri++;       
      }
      Thresholds[period]=thri;
      if(thri>MaxThresholds) MaxThresholds=thri;
      Accmins[period]=vals[0];
      Interval[period]=vals[1];
      Fields[period]=vals[2];
      period++;
    }
    Periods=period;
  
    printf("Periods %d\n",Periods);     
    printf("Input accumulation files\n");     

    bsecs=sec_from_date(obstime);
    for(i=0;i<2;i++)
    {
       AccInds[i]=calloc(Periods,sizeof(int *));
       for(p=0;p<Periods;p++)
       { 
          intsecs=Interval[p]*60;
          accsecs=Accmins[p]*60;
       
          AccInds[i][p]=calloc(Fields[p],sizeof(int));
          for(f=0;f<Fields[p];f++)
          {
             if(!(i+f)) continue;
         
             forsecs=f*intsecs + i*accsecs;
             forind=forsecs/300; /* forecast index for 5-minute interval accumulations */
             min5flags[forind]=1; /* flag this forecast minute data to be read */
             AccInds[i][p][f]=forind; /* store the index for each period and field, begin and end */  
          }
       }
     }

     ravacc=calloc(timesteps,sizeof(int32_t **));
     for(i=0;i<timesteps;i++)
     {
       ravacc[i]=calloc(members,sizeof(int32_t *));
       if(min5flags[i]) /* the same data could be used in several datasets, so read only once */
       {
         formins=i*5;
         forsecs=formins*60;
         secs=bsecs+forsecs;
         date_from_sec(fortime,secs);
         sprintf(accname,"%s/%s_%s-%s+%03d_%s.dat",accdir,accpref,obstime,fortime,formins,Area);
         printf("%s\n",accname);
         NC_ACCF=fopen(accname,"r");
         /* reading accumulation data for each field */
         for(m=0;m<members;m++)
         {
           ravacc[i][m]=malloc(membsize); /* memory is allocated only for data to be read, others are NULL */
           fread(&ravacc[i][m][0],membsize,1,NC_ACCF);
           /* reading member #0 again if deterministic is ignored */
           if(!m && NODETERM) fread(&ravacc[i][m][0],membsize,1,NC_ACCF);
	   	   	   	   
	   /*
           sprintf(moname,"%s/MONITOR_%s-%s-%d.pgm",accdir,obstime,fortime,m);
           MONF=fopen(moname,"w");
           fprintf(MONF,"P5\n%ld %ld\n255\n",xsize,ysize);
	   for(k=0;k<arrsize;k++)
	   { 
	      rval=ravacc[i][m][k];
              if(rval>=0) monacc[k]=rval/1000;
              else monacc[k]=255;
           }
           fwrite(monacc,arrsize,1,MONF);
           fclose(MONF);
	   */
         }
         fclose(NC_ACCF);
       }
     }   
     
  }

  /* Memory allocations */
  membarr_orig=calloc(members,sizeof(int32_t));
  submems=members*10-9;
  pdiv=100.0/((double)submems-1);
  membarr=calloc(submems,sizeof(int32_t));
  out_fracarr=calloc(MaxThresholds+1,sizeof(unsigned char *));
  fracarr=malloc(MaxThresholds);
  for(a=0;a<=MaxThresholds;a++) out_fracarr[a]=malloc(arrsize);


  /* Looping thru different accumulation periods and writing one HDF5 file per time period */
  for(p=0;p<Periods;p++)
  {

    /* HDF5 file root attributes */
    sprintf(HDFoutfile,"%s/%s_%03d-%03d_AccProb_%s.h5",outdir,obstime,Accmins[p],Interval[p],Area);
    printf("Generating %s\n",HDFoutfile);
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
    H5LTset_attribute_string(OUTH5,"what","source",source);
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

    H5LTset_attribute_string(OUTH5,"how","nodes",nodes);

    /* Looping thru fields of accumulation period */
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
	/*        printf("dataset%d\n",f+1); */
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
        
        /* for(a=0;a<=MaxThresholds;a++) memset(&out_fracarr[a][0],255,arrsize); */
	/*
	for(Y=0;Y<ysize;Y++)
	{
	    for(X=0;X<xsize;X++)
	    {
		N=Y*xsize+X;
	*/
        for(N=0;N<arrsize;N++)
	  {
		PASS=0;
		count_zero=0;
                /* Looping thru members and calculating the accumulation for this field f */
       		for(m=0;m<members;m++)
		{
                  membarr_orig[m]=0;
		  if(!Forinds[0]) A=ravacc[Forinds[1]][m][N]; /* if first acc period */
		  else
		  { 
                    A1 = ravacc[Forinds[1]][m][N];
                    A0 = ravacc[Forinds[0]][m][N];
                    if(A1 <0 || A0 < 0) PASS = 1; else A=A1-A0;
		  }
		  if(A<0) PASS=1; 
                  if(PASS) break;
		  if(!A) count_zero++;
		  membarr_orig[m]=A;
		}

		if(PASS) /* If even one of members is outside the area, this pixel is also marked as outside one */
		{  
		    for(a=0;a<=Thresholds[p];a++) out_fracarr[a][N]=255;
		    continue; 
		}
                /* the 'no rain' probability field */          
		out_fracarr[0][N]=(100*count_zero)/members; 

                /* Sorting member values at each pixel */
		qsort(&membarr_orig[0],members,sizeof(int32_t),compare);

                /* interpolating 10 values between sorted members */
		memset(membarr,0,submems*sizeof(int32_t));
		for(m=1;m<members;m++)
		{
		    R1=membarr_orig[m-1];
                    if(!R1) break;
		    R2=membarr_orig[m];
		    dR=(float)(R1-R2)/10.0;
		    for(i=0;i<10;i++)
		    { 
			M=(m-1)*10+i;
			membarr[M]=R1-(int32_t)(dR*(float)i);
		    }
		}
		membarr[submems-1]=membarr_orig[members-1];
		/* 1008456 */
		/*
                if(f==1 && N==826516)  
		{
		   int thI;
                   
		   printf("%ld:",N);
		   for(m=0;m<members;m++) printf("%ld\t",membarr_orig[m]);
                   printf("\n\n");
		   for(M=0;M<submems;M++) printf("%ld ",membarr[M]);
                   printf("\n%d %d\n",p,Thresholds[p]);
                   for(thI=0;thI<=Thresholds[p];thI++) printf("%ld ",acclims[thI][p]);
                   printf("\n\n");
		}
		*/

                /* initialize probabilities to 100% */
		memset(fracarr,100,Thresholds[p]); 
		get_fractiles_per_accs(membarr,fracarr,p);
		for(a=1;a<=Thresholds[p];a++) out_fracarr[a][N]=fracarr[a-1];
	  }
	

        /* Looping thru thresholds and writing data to dataset f for each threshold */
        for(thri=1;thri<=Thresholds[p];thri++)
	{
	  sprintf(datagname,"data%ld",thri);        
	  /* printf("\tdata%ld\n",thri); */        
          datag=H5Gcreate2(datasetg,datagname,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
          datawhatg=H5Gcreate2(datag,"what",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
          H5LTset_attribute_long(datag,"what","threshold_id",&thri,1);
          H5LTset_attribute_double(datag,"what","threshold_value",&f_acclims[thri][p],1);
 
          /* create new empty dataset and attributes to destination file */ 
          dataset = H5Dcreate2(datag,"data", datatype, dataspace,
                    H5P_DEFAULT, plist, H5P_DEFAULT);
          H5LTset_attribute_string( datag, "data", "CLASS", "IMAGE");
          H5LTset_attribute_string( datag, "data", "IMAGE_VERSION", "1.2");
          H5Dwrite(dataset, datatype, H5S_ALL, dataspace, H5P_DEFAULT,&out_fracarr[thri][0]);
          H5Dclose(dataset);
          H5Gclose(datawhatg);
          H5Gclose(datag);

	  /* GeoTIFF and PGM outputs if needed */
          { 
             char fortime_start[15],fortime_end[15],origtime[15];
             double thres;
             uint8_t *data;


             thres=f_acclims[thri][p];
             sprintf(origtime,"%.8sT%.4sZ",obstime,obstime+8);
             sprintf(fortime_start,"%sT%.4sZ",fielddate[0],fieldtime[0]);
             sprintf(fortime_end,"%sT%.4sZ",fielddate[1],fieldtime[1]);
             data=&out_fracarr[thri][0];

             if(GEOTIFF_OUT)
 	     {    
                GDALDatasetH hDstDS;        
                GDALRasterBandH hBand;
                char **papszMetaData=NULL;
                char metatag[2500];
                char GTiffoutfile[200];

                /* Creating the GeoTIFF file */
                sprintf(GTiffoutfile,"%s/ExcProbAcc_%s_%s-%s_Thr=%.2f_Area=%s.tif",
                                     outdir,origtime,fortime_start,fortime_end,thres,Area);
                hDstDS = GDALCreate( hDriver, GTiffoutfile, xsize, ysize, 1, GDT_datatype, papszOptions );

                /* Insert metadata as XML-tags */
                sprintf(metatag,"<GDALMetadata>");
                sprintf(metatag,"%s<Item name=\"Observation time\" format=\"YYYYMMDDThhmmZ\">%s</Item>\n",metatag,obstime);
                sprintf(metatag,"%s<Item name=\"Forecast start time\" format=\"YYYYMMDDThhmmZ\">%s</Item>\n",metatag,fortime_start);
                sprintf(metatag,"%s<Item name=\"Forecast end time\" format=\"YYYYMMDDThhmmZ\">%s</Item>\n",metatag,fortime_end);
                sprintf(metatag,"%s<Item name=\"Time zone\">%s</Item>\n",metatag,timezone_str);

                sprintf(metatag,"%s<Item name=\"Quantity\" unit=\"%%\">Exceedance probability of precipitation accumulation</Item>\n",metatag);
                sprintf(metatag,"%s<Item name=\"Gain\">%f</Item>\n",metatag,gain);
                sprintf(metatag,"%s<Item name=\"Offset\">%f</Item>\n",metatag,offset);
                sprintf(metatag,"%s<Item name=\"Nodata\">%.0f</Item>\n",metatag,nodata);
                sprintf(metatag,"%s<Item name=\"Undetect\">%.0f</Item>\n",metatag,undetect);

                sprintf(metatag,"%s<Item name=\"Accumulation time\" unit=\"min\">%d</Item>\n",metatag,Accmins[p]);
                sprintf(metatag,"%s<Item name=\"Accumulation threshold\" unit=\"mm\">%.2f</Item>\n",metatag,thres);
                sprintf(metatag,"%s<Item name=\"Area name\">%s</Item>\n",metatag,Area);
                sprintf(metatag,"%s</GDALMetadata>",metatag);
                papszMetaData = CSLSetNameValue(papszMetaData,"GDAL_METADATA",metatag);
                GDALSetMetadata(hDstDS,papszMetaData,NULL);

                /* Insert geolocation information */
                GDALSetGeoTransform( hDstDS, adfGeoTransform );
                GDALSetProjection( hDstDS, pszSRS_WKT );

                /* Insert the data */
                hBand = GDALGetRasterBand( hDstDS, 1 );
                GDALRasterIO( hBand, GF_Write, 0, 0, xsize, ysize, 
                           data, xsize, ysize, GDT_datatype, 0, 0 );    

                GDALClose( hDstDS );
	     }
             
             if(PGM_OUT)  /* PGM output */
	     {
                char PGMoutfile[200];
                FILE *PGMF;
	  
                sprintf(PGMoutfile,"%s/ExcProbAcc_%.12s_%s%.4s-%s%.4s_Accmins=%d_Thr=%.2f_Area=%s.pgm",
                     outdir,obstime,fielddate[0],fieldtime[0],fielddate[1],fieldtime[1],Accmins[p],thres,Area);
                PGMF=fopen(PGMoutfile,"w");
                fprintf(PGMF,"P5\n# obstime %s%s\n",fielddate[0],fieldtime[0]);
                fprintf(PGMF,"# param Probability\n# projdef %s\n# AccumulationThreshold %.2f\n",projdef,thres);
                fprintf(PGMF,"# bottomleft %f %f\n# topright %f %f\n",LL_lon,LL_lat,UR_lon,UR_lat);
                fprintf(PGMF,"%ld %ld\n255\n",xsize,ysize);
                fwrite(data,xsize*ysize,1,PGMF);
                fclose(PGMF);
	     }
	  }
	}
        H5Gclose(setwhatg);
        H5Gclose(datasetg);
    }
    H5Gclose(whatg);
    H5Gclose(whereg);
    H5Gclose(howg);
  }
  H5Fclose(OUTH5);
  /*
  for(a=0;a<=Thresholds[0];a++)
  {
    sprintf(outname,"%s/%s-%s+%03d_Prob-of-acc=%04d.pgm",outdir,obstime,fortime,formins,(int)acclims[0][a]/100);
    OUTF=fopen(outname,"w");
    fprintf(OUTF,"P5\n%ld %ld\n255\n",xsize,ysize);
    fwrite(&out_fracarr[a][0],arrsize,1,OUTF);
    fclose(OUTF);
  }
  */
  
  if(GEOTIFF_OUT)
  {
     OSRDestroySpatialReference( hSRS );
     CPLFree( pszSRS_WKT );
     CSLDestroy( papszOptions );
  }

  for(a=0;a<=MaxThresholds;a++) free(out_fracarr[a]);
  free(out_fracarr);
  free(fracarr);
  for(i=0;i<timesteps;i++) 
  {
     for(m=0;m<members;m++)free(ravacc[i][m]);
     free(ravacc[i]);
  }
  free(ravacc);
  for(aI=0;aI<ACCLIMS;aI++) { free(f_acclims[aI]);  free(acclims[aI]); } 
  free(membarr);
  free(membarr_orig);
  printf("Thresholding ready\n");
  return(0);
}


void get_fractiles_per_accs(int32_t *membarr, unsigned char *fracarr, int period)  
{
  int32_t m,thI,R1,R2,Rmin;

  thI=Thresholds[period]-1; /* Threshold array index */
  for(m=1;m<submems;m++)
  {
    if(thI<0) break;
    Rmin=acclims[thI+1][period];
    R2=membarr[m];
    R1=membarr[m-1];
    /*    if(!R1) break; */
    if(Rmin>R1)
    {
       fracarr[thI]=0;
       thI--;
       continue;
    } 
    if(Rmin>=R2)
    {
      fracarr[thI]=m*pdiv;
      thI--;
      continue;
    }
  }
}



/*--------------------------------------------------------------------------------------------------*/

void date_from_sec(char *preddate,time_t secs)
{
   struct tm *Sdd;

   Sdd=gmtime(&secs);

   sprintf(preddate,"%d%02d%02d%02d%02d",Sdd->tm_year+1900,Sdd->tm_mon+1,
                                 Sdd->tm_mday,Sdd->tm_hour,
                                 Sdd->tm_min);
   return;
}

/*--------------------------------------------------------------------------------------------------*/

time_t sec_from_date(char *stim)
{
   struct tm Sdd;
   time_t secs;

   sscanf(stim,"%4d%2d%2d%2d%2d",&Sdd.tm_year,&Sdd.tm_mon,
                                    &Sdd.tm_mday,&Sdd.tm_hour,
                                    &Sdd.tm_min);
   Sdd.tm_year-=1900;
   Sdd.tm_mon--;
   Sdd.tm_sec=0;
   secs=mktime(&Sdd);
   return(secs);
}

