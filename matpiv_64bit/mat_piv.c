#include "mex.h"
#include "math.h"
#include "string.h"
#include "time.h"

#include "piv_util.h"
#include "matrix.h"
#include "fftw.h"


/* Global variables*/

short NROWS, NCOLS, NX, NY, NI, NJ;
int nri, nci;
IMG_DATA **I1, **I2, **I1_ssq, **I2_ssq;
double *xi, *yi, *uest, *vest, *u, *v, *c, *fout, *gout, *cout;

int i_weight_flag, N_sigma,shftmode,maxpass,subpix_maxpass;
float subpix_eps, min_corr_coef, Umin,Umax,Vmin,Vmax;
int i_GaussFilt_Flag, i_MedianFilt_Flag, OutC_flag;

fcomplex **h, **f, **g;
float **r;

/* #######################################################################*/
/* Parameters for MQD and PIV tracking ###################################*/
short NXMQ,NYMQ,NSRX,NSRY;  //Pattern size and search (tracking) radius
short NXmn,NYmn;            // Size of Qmn, Pmn, and Dmn
int Interrogation_Method;   // =0: correlation, =1: Tracking correlation, 
                            // =2: MQD
float **Qmn, **Pmn, **Dmn;  // Qmn is the cumsum of g2^2 at the tracking position
                            // Dmn is the quadratic difference
                            // Pmn is the cross-correlation
                            // Dmn = g1^2(a const) + Qmn - 2Pmn
/* #######################################################################*/

/* Parameters for Hart correlation*/
int Hart_flag;
float Hart_overlap_percentage;
int Hart_offset_x, Hart_offset_y;
float **HartR;

/* Parameters for sub-pixel estimation methods*/
int subpixel_method;
int subpixel_max_pass;
float subpixel_epsilon;

void PIV_Correlation();
void PIV_Correlation_Hart();
int Read_PIV_Parameters(mxArray * para);
float Sub_Sample_Img(float** I, float** I_ssq, int nrows, int ncols, fcomplex **f, int nx, int ny, float cx, float cy);
void PIV_MakePlan();
void PIV_End();

/* #######################################################################*/
/* functions for MQD ####################################################*/
void PIV_MQD();
void PIV_MQD_Hart();
float Sub_Img_MQD(float** I, float** I_ssq, int nrows, int ncols, fcomplex **f, int nx, int ny, int padx, int pady, float cx, float cy);
/* #######################################################################*/



/*This is the main function that act as the interface between MatLab and C routines*/

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int m, n, i, j;
	double * data;

	float utmp,vtmp,ctmp;
	float xii,yii,uee,vee;
	double time1,time2;

/* Check if the number of inputs is right*/
    if (nrhs!=7)
    {
        printf("MatPIV error: There must be be 7 input fields.\n");
        return;
    }
    
/* Check the formats of the input image pair, and read them in Matrics I1 and I2*/
   	if( !mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) )
	{
		printf("MatPIV Error: Input images are not DOUBLE type, please convert them to DOUBLE before calling MatPIV\n");
		return;
	}
    

	NROWS = mxGetM(prhs[0]);
	NCOLS = mxGetN(prhs[0]);
	if( NROWS!=mxGetM(prhs[1]) || NCOLS!=mxGetN(prhs[1]) )
	{
		printf("MatPIV Error: The dimensions of the image pair do not match\n");
		return;			
	}    
    I1=(IMG_DATA**) matrix(NROWS,NCOLS,sizeof(IMG_DATA));
    I2=(IMG_DATA**) matrix(NROWS,NCOLS,sizeof(IMG_DATA));
	data = mxGetPr(prhs[0]);    
    for(i=0;i<NROWS;i++)
        for(j=0;j<NCOLS;j++)
            I1[i][j]=data[j*NROWS+i];
	data = mxGetPr(prhs[1]);
    for(i=0;i<NROWS;i++)
        for(j=0;j<NCOLS;j++)
            I2[i][j]=data[j*NROWS+i];


    
  /* Get the grid matrice [xi, yi] and the estimated displacement field [ue, ve]*/
  /* Note: the dimension of xi,yi,ue and ve must be the same.*/
	nri=mxGetM(prhs[2]);
	nci=mxGetN(prhs[2]);    
    if(nri != mxGetM(prhs[3]) || nri != mxGetM(prhs[4]) || nri != mxGetM(prhs[5]) ||
       nci != mxGetN(prhs[3]) || nci != mxGetN(prhs[4]) || nci != mxGetN(prhs[5]) )
    {
        printf("MatPIV error: The input grid field and estimated velocity field are not consistent in size.\n");
        return;
    }
    if( !mxIsDouble(prhs[2]) || !mxIsDouble(prhs[3]) ||
        !mxIsDouble(prhs[4]) || !mxIsDouble(prhs[4]) )
    {
        printf("MatPIV error: The input grid field and estimated velocity field has to be double precision type.\n");
    }
	xi = mxGetPr(prhs[2]);
	yi = mxGetPr(prhs[3]);
	uest = mxGetPr(prhs[4]);
	vest = mxGetPr(prhs[5]);
     
/* Get PIV interrogation parameters (subwindow size, subpixel precision, filters, etc.)*/
/* All the parameters should be capsulated into a MatLab Structure */
    if( !mxIsStruct(prhs[6]) )
    {
        printf("MatPIV error: PIV parameters are not inputted as a Matlab Structure data.\n"); 
        return;
    }
    
    if( Read_PIV_Parameters(prhs[6]) == -1 )
    {
        return;
    }

/* Prepare space and map for output    */
	plhs[0]= mxCreateDoubleMatrix(nri,nci,mxREAL);
	plhs[1]= mxCreateDoubleMatrix(nri,nci,mxREAL);
	plhs[2]= mxCreateDoubleMatrix(nri,nci,mxREAL);
    

    OutC_flag=0;
    if(nlhs>3)
    {
        OutC_flag=1;
        if(Interrogation_Method==0)
            plhs[3]= mxCreateDoubleMatrix(NX*nri*nci,NY,mxREAL);
        else
            plhs[3]= mxCreateDoubleMatrix(NXmn*nri*nci,NYmn, mxREAL);
    }

	u = mxGetPr(plhs[0]);
	v = mxGetPr(plhs[1]);
    c = mxGetPr(plhs[2]);
    if (nlhs>3)
    {
        cout = mxGetPr(plhs[3]);
    }
    
    
/* Start PIV interrogation    */

	time1=(double)clock();

    PIV_MakePlan();

    if(Hart_flag==0)
    {
        if(Interrogation_Method==0)        PIV_Correlation();    
        else    PIV_MQD();
    }
    else
    {
        if(Interrogation_Method==0)        PIV_Correlation_Hart();    
        else    PIV_MQD_Hart();
    }
   
    PIV_End();
	time2=(double)clock();
    

	printf("Time cost = %f s (actually %f tickes)\n",(float)((time2-time1)/CLOCKS_PER_SEC), (float)(time2-time1));
  
	free_matrix(I1);
	free_matrix(I2);
}


/* This subroutine read out PIV parameters that were inputted as a MatLab Structure*/
/* If some parameters are not given by the input structure, default values will be given.*/
int Read_PIV_Parameters(mxArray * para)
{
    const char *fieldname;
    int nfields,i;
    mxArray *data;
    double *vel_thresh;

/*Default values for PIV interrogation parameters*/
    NX=32;  NY=32;  
    maxpass=1;
    min_corr_coef=0.0;
    Umin=-NX/2;  Umax=NX/2;     Vmin=-NY/2;  Vmax=NY/2;
    
    Interrogation_Method=0;
    NSRX=NX/2; NSRY=NY/2; 
    
    Hart_flag=0;
    Hart_overlap_percentage=0.5;
    
    subpixel_method=1;
    subpixel_max_pass=3;
    subpixel_epsilon=0.001;
    
    nfields=mxGetNumberOfFields(para);
    
    for (i=0;i<nfields;i++)
    {        
        fieldname=mxGetFieldNameByNumber(para,i);
        data=mxGetFieldByNumber(para,0,i);
    
    /* Interrogation subwindow size, NX and NY  */ 
        if( strcmpi(fieldname,"nx")==0 )
        {
            NX=(int)(mxGetScalar(data));
            NSRX=NX/2;
        }
        else if ( strcmpi(fieldname,"ny")==0 )
        {
            NY=(int)(mxGetScalar(data));  
            NSRY=NY/2; 
        }
    /* Maximal number of interations (for dynamic subwindow shifting)    */ 
        else if( strcmpi(fieldname,"max_pass")==0 )
        {
            maxpass=(int)(mxGetScalar(data));
        }
    /* Threshold value for the cross-correlation coeffecient    */
        else if( strcmpi(fieldname,"Min_Corr_Coef")==0 )
        {
            min_corr_coef=(float)(mxGetScalar(data));
        }       
    /* Valid displacement measurement range (w.r.t. uest and vest)    */
        else if( strcmpi(fieldname,"Vel_Range")==0 )
        {
            if(mxGetM(data)*mxGetN(data)!=4)
            {
                printf("MatPIV Error: Threshold for velocity range (Vel_Range) must be a vecotr of 4 elements.\n");
                return -1;
            }
            vel_thresh=mxGetPr(data);
            Umin=vel_thresh[0];     Umax=vel_thresh[1];
            Vmin=vel_thresh[2];     Vmax=vel_thresh[3];
            if(Umin>Umax || Vmin>Vmax)
            {
                printf("MatPIV Error: For the given velocity range, Umin>Umax or Vmin>Vmax.\n");
                return -1;
            }
        }
    /* Switch of MQD method*/
        else if( strcmpi(fieldname,"method")==0 )
        {
            Interrogation_Method = (int) (mxGetScalar(data));
            if(Interrogation_Method<0)  Interrogation_Method=0;
            if(Interrogation_Method>2)  Interrogation_Method=2;
        }
    /* Parameters for MQD method*/
        else if( strcmpi(fieldname,"PIV_Tracking_radius_x")==0 )
        {
            NSRX = (int)(mxGetScalar(data));
        }
        else if( strcmpi(fieldname,"PIV_Tracking_radius_y")==0 )
        {
            NSRY = (int)(mxGetScalar(data));
        }
        else if( strcmpi(fieldname,"Hart_correlation")==0 )
        {
            Hart_flag=(int)(mxGetScalar(data));
            Hart_flag >0 ? 1 : 0;
        }
        else if( strcmpi(fieldname,"Hart_overlap")==0 )
        {
            Hart_overlap_percentage=mxGetScalar(data);
            if(Hart_overlap_percentage>1.0 || Hart_overlap_percentage <0.0)
            {
                printf("MatPIV Error: Subwindow overlap for Hart correlation must be 0.0~1.0");
                return -1;
            }
        }
        else if( strcmpi(fieldname,"Sub_Pixel_method")==0 )
        {
            subpixel_method = (int) (mxGetScalar(data));
            if(subpixel_method<0)   subpixel_method=0;
            if(subpixel_method>2)   subpixel_method=2;
        }
        else
        {
            printf("MatPIV Error: PIV Parameter %s is not a valid field name.\n",fieldname);
            return -1;
        }        
    }
    
    Hart_offset_x=(int)((1.0-Hart_overlap_percentage)*NX*0.5);
    Hart_offset_y=(int)((1.0-Hart_overlap_percentage)*NY*0.5);
 
    NXMQ=NX+NSRX*2;
    NYMQ=NY+NSRY*2;
    
    NXmn=NSRX*2+1;
    NYmn=NSRY*2+1;
    
    return 0;
}

void PIV_MakePlan()
{
    int i,j;
    int nxh,nyh;
    
    if(Interrogation_Method==0)
    {
        h =(fcomplex**) matrix(NY, NX, sizeof(fcomplex));
        r =(float**) matrix(NY, NX, sizeof(float));
        f =(fcomplex**) matrix(NY, NX, sizeof(fcomplex));
        g =(fcomplex**) matrix(NY, NX, sizeof(fcomplex));
        HartR= (float**) matrix(NY, NX, sizeof(float));
    }
    else
    {
        h =(fcomplex**) matrix(NYMQ, NXMQ, sizeof(fcomplex));
        r =(float**) matrix(NYMQ, NXMQ, sizeof(float));
        f =(fcomplex**) matrix(NYMQ, NXMQ, sizeof(fcomplex));
        g =(fcomplex**) matrix(NYMQ, NXMQ, sizeof(fcomplex));
        Qmn = (float**) matrix(NYmn, NXmn, sizeof(float));
        Pmn = (float**) matrix(NYmn, NXmn, sizeof(float));
        Dmn = (float**) matrix(NYmn, NXmn, sizeof(float));
        HartR= (float**) matrix(NYMQ, NXMQ, sizeof(float));
    }
    
    I1_ssq = (IMG_DATA**) matrix(NROWS+1,NCOLS+1,sizeof(IMG_DATA));
    I2_ssq = (IMG_DATA**) matrix(NROWS+1,NCOLS+1,sizeof(IMG_DATA));   
    I1_ssq[0][0]=0.0;   I2_ssq[0][0]=0.0;
    for(i=0;i<NROWS;i++)
    {
        I1_ssq[i+1][0]=0.0;
        I2_ssq[i+1][0]=0.0;
        for(j=0;j<NCOLS;j++)
        {
            I1_ssq[i+1][j+1]=I1_ssq[i+1][j]+I1[i][j]*I1[i][j];
            I2_ssq[i+1][j+1]=I2_ssq[i+1][j]+I2[i][j]*I2[i][j];
        }
    }
    for (j=0;j<NCOLS;j++)
    {
        I1_ssq[0][j+1]=0.0;
        I2_ssq[0][j+1]=0.0;
        for(i=0;i<NROWS;i++)
        {
            I1_ssq[i+1][j+1] += I1_ssq[i][j+1];
            I2_ssq[i+1][j+1] += I2_ssq[i][j+1];
        }
    }
}

float Sub_Sample_Img(float** I, float ** I_ssq, int nrows, int ncols, fcomplex **f, int nx, int ny, float cx, float cy)
{
	int i,j,ii,jj,si,ei,sj,ej;
	int sx,sy;
    float ssq,factor;
    float Left,Upper,LU;

    sx=(int)floor(cx+0.499999)-nx/2;
	sy=(int)floor(cy+0.499999)-ny/2;
    
    memset(*f,0,sizeof(fcomplex)*nx*ny);
    si=(sy>=0 ? 0 : -sy);
    sj=(sx>=0 ? 0 : -sx);
    ei=((ny+sy)<=nrows ? ny: nrows-sy);
    ej=((nx+sx)<=ncols ? nx: ncols-sx);
    
    Left = I_ssq[ei+sy][sj+sx];
    Upper = I_ssq[si+sy][ej+sx];
    LU = I_ssq[si+sy][sj+sx];    
    ssq = I_ssq[ei+sy][ej+sx] - Left - Upper + LU;
            
	for(i=si;i<ei;i++)
	{
		ii=i+sy;
		for(j=sj;j<ej;j++)
		{
			jj=j+sx;			
            f[i][j].r=I[ii][jj];
		}
	}
    
    return ssq;
}

void PIV_End()
{
    free_matrix(h);
    free_matrix(r);
    free_matrix(f);
    free_matrix(g);
    free_matrix(HartR);
    
    if(Interrogation_Method>0)
    {
        free_matrix(Qmn);
        free_matrix(Pmn);
        free_matrix(Dmn);
    }
    
    free_matrix(I1_ssq);
    free_matrix(I2_ssq);
}

void PIV_Correlation()
{
    int i, j, ii, jj, npass;
    float ex,ey,ssx,ssy;
    float maxcc=0.5, cc_first=0.0;
    int IP=0, JP=0, success_flag;
    float ssq1,ssq2,factor;
    
    for (i=0;i<nri*nci;i++)
    {
        ex=floor(uest[i]);
        ey=-floor(vest[i]);
        
        success_flag=0;
        npass=0;
        do{
            ssq1=Sub_Sample_Img(I1,I1_ssq,NROWS,NCOLS,f,NX,NY,xi[i]-ex/2,yi[i]-ey/2);
            ssq2=Sub_Sample_Img(I2,I2_ssq,NROWS,NCOLS,g,NX,NY,xi[i]+ex/2,yi[i]+ey/2);
            factor=(ssq1+ssq2)/2.0; 
            factor= (factor>0? 1.0/factor : 1.0);
            Cross_CR(f,g,factor,NY,NX,r,h,&maxcc,&IP,&JP);
            if(npass==0) cc_first=maxcc;   
            ex+=JP;
            ey+=IP;
            npass++;
        }while(npass<maxpass & (IP!=0||JP!=0));
        
        if(npass<maxpass||maxcc>=cc_first)   success_flag=1;
        else    success_flag=0;
                
        if(success_flag==1 && abs(IP)<NY/3.0 && abs(JP)<NX/3 
           && maxcc>min_corr_coef && ex>=Umin && ex<=Umax && -ey>=Vmin && -ey<=Vmax)
        {
            if(subpixel_method==1)
                Iter_Peak_Fit(r,NY,NX,(NY>>1)+IP,(NX>>1)+JP,&ssx,&ssy);
            else
                Peak_Fit(r,(NY>>1)+IP,(NX>>1)+JP,&ssx,&ssy);
            if(ssx!=2000. && ssy!=2000.)
            {
                u[i]=ex+ssx;    v[i]=-ey-ssy;   c[i]=maxcc;
            }
            else
            {
                u[i]=0.0;   v[i]=0.0;   c[i]=0.0;
            }
        }        
        else
        {
            u[i]=0.0;   v[i]=0.0;   c[i]=0.0;
        }
        
        if(OutC_flag>0)
        {
            for(ii=0;ii<NY;ii++)
            {
                for(jj=0;jj<NX;jj++)
                {
                    cout[ii*NX*nci*nri + i*NX+jj]=r[ii][jj];
                }
            }
        }
    } 

}

void PIV_Correlation_Hart()
{
    int i, j, ii, jj, npass;
    float ex,ey,ssx,ssy;
    float maxcc=0.5, cc_first=0.0;
    int IP=0, JP=0, success_flag;
    float ssq1,ssq2,factor;
    
    for (i=0;i<nri*nci;i++)
    {
        ex=floor(uest[i]);
        ey=-floor(vest[i]);
        
        success_flag=0;
        npass=0;
        do{
            ssq1=Sub_Sample_Img(I1,I1_ssq,NROWS,NCOLS,f,NX,NY,xi[i]-ex/2,yi[i]-ey/2);
            ssq2=Sub_Sample_Img(I2,I2_ssq,NROWS,NCOLS,g,NX,NY,xi[i]+ex/2,yi[i]+ey/2);
            factor=(ssq1+ssq2)/2.0; 
            factor= (factor>0? 1.0/factor : 1.0);
            Cross_CR(f,g,factor,NY,NX,r,h,&maxcc,&IP,&JP);
            if(npass==0) cc_first=maxcc;   
            ex+=JP;
            ey+=IP;
            npass++;
        }while(npass<maxpass & (IP!=0||JP!=0));
        
        if(npass<maxpass||maxcc>=cc_first)   success_flag=1;
        else    success_flag=0;
                
        if(success_flag==1 && abs(IP)<NY/3.0 && abs(JP)<NX/3 
           && maxcc>min_corr_coef && ex>=Umin && ex<=Umax && -ey>=Vmin && -ey<=Vmax)
        {
            Iter_Peak_Fit(r,NY,NX,(NY>>1)+IP,(NX>>1)+JP,&ssx,&ssy);
            if(ssx!=2000. && ssy!=2000.)
            {
                u[i]=ex+ssx;    v[i]=-ey-ssy;   c[i]=maxcc;
            }
            else
            {
                u[i]=0.0;   v[i]=0.0;   c[i]=0.0;
            }
        }        
        else
        {
            u[i]=0.0;   v[i]=0.0;   c[i]=0.0;
        }
        
        if(OutC_flag>0)
        {
            for(ii=0;ii<NY;ii++)
            {
                for(jj=0;jj<NX;jj++)
                {
                    cout[ii*NX*nci*nri + i*NX+jj]=r[ii][jj];
                }
            }
        }
    } 

}


float Sub_Img_MQD(float** I, float** I_ssq, int nrows, int ncols, fcomplex **f, int nx, int ny, int nnx, int nny, float cx, float cy)
{
	int i,j,ii,jj,si,ei,sj,ej;
	int sx,sy;
    int lux,luy;
    
    float Left,Upper,LU,ssq,factor;
    
    memset(*f,0,sizeof(fcomplex)*nx*ny);
    
    lux=(nx-nnx)/2;
    luy=(ny-nny)/2;
    
    sx=(int)floor(cx+0.499999)-nnx/2;
	sy=(int)floor(cy+0.499999)-nny/2;
    
    si=(sy>=0 ? 0 : -sy);
    sj=(sx>=0 ? 0 : -sx);
    ei=((nny+sy)<=nrows ? nny: nrows-sy);
    ej=((nnx+sx)<=ncols ? nnx: ncols-sx); 

    Left = I_ssq[ei+sy][sj+sx];
    Upper = I_ssq[si+sy][ej+sx];
    LU = I_ssq[si+sy][sj+sx];    
    ssq = I_ssq[ei+sy][ej+sx] - Left - Upper + LU;
    
        
	for(i=si;i<ei;i++)
	{
		ii=i+sy;
		for(j=sj;j<ej;j++)
		{
			jj=j+sx;			
            f[i+luy][j+lux].r=I[ii][jj];
		}
	}
    
    return ssq;
}

void Sub_Q(float **Is2, int nrows, int ncols, float **Qmn, int nx, int ny, int nnx, int nny, float cx, float cy)
{
    int i,j,ii,jj,si,ei,sj,ej,mnx,mny;
	int sx,sy;
        
    mnx=nx-nnx+1;
    mny=ny-nny+1;
    
    sx=(int)floor(cx+0.499999)-nx/2;
	sy=(int)floor(cy+0.499999)-ny/2;    

    if(sx>=0 && sy>=0 && sx+nx<=ncols && sy+ny<=nrows)
    {
        
        for(i=0;i<mny;i++)
        {
            for(j=0;j<mnx;j++)
            {
                Qmn[i][j]=Is2[sy+i+nny][sx+j+nnx] - Is2[sy+i][sx+j+nnx] 
                          -Is2[sy+i+nny][sx+j] + Is2[sy+i][sx+j];
            }
        }
    }    
    else
    {
        for(i=0;i<mny;i++)
        {
            si= sy+i;   ei=sy+i+nny;
            si=(si+abs(si))/2;
            si=nrows-((nrows-si)+abs(nrows-si))/2;
            ei=(ei+abs(ei))/2;
            ei=nrows-((nrows-ei)+abs(nrows-ei))/2;
            for(j=0;j<mnx;j++)
            {
                sj=sx+j;    ej=sx+j+nnx;                
                sj=(sj+abs(sj))/2;
                sj=ncols-((ncols-sj)+abs(ncols-sj))/2;
                ej=(ej+abs(ej))/2;
                ej=ncols-((ncols-ej)+abs(ncols-ej))/2;
                Qmn[i][j]=Is2[ei][ej]-Is2[si][ej]-Is2[ei][sj]+Is2[si][sj];
            }
        }
    }
    
}

void PIV_MQD()
{
    int i, j, ii, jj, npass;
    float ex,ey,ssx,ssy;
    float maxcc=0.5, cc_first=0.0;
    int IP=0, JP=0, success_flag;
    float ssq1, ssq2, factor;
    
    for (i=0;i<nri*nci;i++)
    {
        ex=floor(uest[i]);
        ey=-floor(vest[i]);
        
        success_flag=0;
        npass=0;
        do{
            ssq1=Sub_Img_MQD(I1,I1_ssq,NROWS,NCOLS,f,NXMQ,NYMQ,NX,NY,xi[i]-ex/2.,yi[i]-ey/2.);
            ssq2=Sub_Img_MQD(I2,I2_ssq,NROWS,NCOLS,g,NXMQ,NYMQ,NXMQ,NYMQ,xi[i]+ex/2.,yi[i]+ey/2.);
            factor=(ssq1+ssq2)/2.0; 
            factor= (factor>0? 1.0/factor : 1.0);
            if(Interrogation_Method==2)
                Sub_Q(I2_ssq,NROWS,NCOLS,Qmn,NXMQ,NYMQ,NX,NY,xi[i]+ex/2.,yi[i]+ey/2.);        
            Cross_CR(f,g,factor,NYMQ,NXMQ,r,h,&maxcc,&IP,&JP);
            maxcc=0.0;
            for(ii=0;ii<NYmn;ii++)
            {
                for(jj=0;jj<NXmn;jj++)
                {
                    Pmn[ii][jj]=r[ii+NYMQ/2-NSRY][jj+NXMQ/2-NSRX];
                    Dmn[ii][jj]=Pmn[ii][jj];
                    if(Interrogation_Method==2)
                        Dmn[ii][jj]-=(Qmn[ii][jj]*factor*0.5+ssq1*factor*0.5-1.0);
        
                    if(Dmn[ii][jj]>maxcc)
                    {
                        maxcc=Dmn[ii][jj];
                        IP=ii-NYmn/2;
                        JP=jj-NXmn/2;
                    }
                }
            }
            if(npass==0)    cc_first=maxcc;
            ex+=JP;
            ey+=IP;
            npass++;
        }while(npass<maxpass & (IP!=0||JP!=0));
        
        if(npass<maxpass||maxcc>=cc_first*0.9)   success_flag=1;
        else    success_flag=0;      
        
        if(success_flag==1 && maxcc>min_corr_coef 
            &&ex>=Umin && ex<=Umax && -ey>=Vmin && -ey<=Vmax)
        {
            if(subpixel_method==1)
                Iter_Peak_Fit(Dmn, NYmn, NXmn, (NYmn>>1)+IP,(NXmn>>1)+JP,&ssx,&ssy);
            else
                Peak_Fit(Dmn, (NYmn>>1)+IP,(NXmn>>1)+JP,&ssx,&ssy);
                
            if(ssx!=2000. && ssy!=2000.)
            {
                u[i]=ex+ssx;    v[i]=-ey-ssy;   c[i]=maxcc;
            }
            else
            {
                u[i]=0.0;   v[i]=0.0;   c[i]=0.0;
            }
        }        
        else
        {
            u[i]=0.0;   v[i]=0.0;   c[i]=0.0;
        }
        
        if(OutC_flag>0)
        {
            for(ii=0;ii<NYmn;ii++)
            {
                for(jj=0;jj<NXmn;jj++)
                {
                    cout[ii*NXmn*nci*nri + i*NXmn+jj]=Dmn[ii][jj];
                }
            }
        }
    }
    
}


/*Currently does not support MQD, only PIV_tracking*/
void PIV_MQD_Hart()
{
    int i, j, ii, jj, npass;
    float ex,ey,ssx,ssy;
    float maxcc=0.5, cc_first=0.0, mincc;
    int IP=0, JP=0, success_flag;
    float ssq1,ssq2,factor;
    
    mincc=min_corr_coef*min_corr_coef*min_corr_coef*min_corr_coef;
    
    for (i=0;i<nri*nci;i++)
    {
        ex=floor(uest[i]);
        ey=-floor(vest[i]);
        
        success_flag=0;
        npass=0;
        do{
            ssq1=Sub_Img_MQD(I1,I1_ssq,NROWS,NCOLS,f,NXMQ,NYMQ,NX,NY,xi[i]-ex/2.-Hart_offset_x,yi[i]-ey/2.-Hart_offset_y);
            ssq2=Sub_Img_MQD(I2,I2_ssq,NROWS,NCOLS,g,NXMQ,NYMQ,NXMQ,NYMQ,xi[i]+ex/2.-Hart_offset_x,yi[i]+ey/2.-Hart_offset_y);
            factor=(ssq1+ssq2)/2.0; 
            factor= (factor>0? 1.0/factor : 1.0);
            Cross_CR(f,g,factor,NYMQ,NXMQ,r,h,&maxcc,&IP,&JP);
            for(ii=0;ii<NYMQ;ii++)
                for(jj=0;jj<NXMQ;jj++)
                    HartR[ii][jj]=r[ii][jj];
            
            ssq1=Sub_Img_MQD(I1,I1_ssq,NROWS,NCOLS,f,NXMQ,NYMQ,NX,NY,xi[i]-ex/2.+Hart_offset_x,yi[i]-ey/2.-Hart_offset_y);
            ssq2=Sub_Img_MQD(I2,I2_ssq,NROWS,NCOLS,g,NXMQ,NYMQ,NXMQ,NYMQ,xi[i]+ex/2.+Hart_offset_x,yi[i]+ey/2.-Hart_offset_y);
            factor=(ssq1+ssq2)/2.0; 
            factor= (factor>0? 1.0/factor : 1.0);
            Cross_CR(f,g,factor,NYMQ,NXMQ,r,h,&maxcc,&IP,&JP);
            for(ii=0;ii<NYMQ;ii++)
                for(jj=0;jj<NXMQ;jj++)
                    HartR[ii][jj]*=r[ii][jj];
            
            ssq1=Sub_Img_MQD(I1,I1_ssq,NROWS,NCOLS,f,NXMQ,NYMQ,NX,NY,xi[i]-ex/2.-Hart_offset_x,yi[i]-ey/2.+Hart_offset_y);
            ssq2=Sub_Img_MQD(I2,I2_ssq,NROWS,NCOLS,g,NXMQ,NYMQ,NXMQ,NYMQ,xi[i]+ex/2.-Hart_offset_x,yi[i]+ey/2.+Hart_offset_y);
            factor=(ssq1+ssq2)/2.0; 
            factor= (factor>0? 1.0/factor : 1.0);
            Cross_CR(f,g,factor,NYMQ,NXMQ,r,h,&maxcc,&IP,&JP);
            for(ii=0;ii<NYMQ;ii++)
                for(jj=0;jj<NXMQ;jj++)
                    HartR[ii][jj]*=r[ii][jj];
            
            ssq1=Sub_Img_MQD(I1,I1_ssq,NROWS,NCOLS,f,NXMQ,NYMQ,NX,NY,xi[i]-ex/2.+Hart_offset_x,yi[i]-ey/2.+Hart_offset_y);
            ssq2=Sub_Img_MQD(I2,I2_ssq,NROWS,NCOLS,g,NXMQ,NYMQ,NXMQ,NYMQ,xi[i]+ex/2.+Hart_offset_x,yi[i]+ey/2.+Hart_offset_y);
            factor=(ssq1+ssq2)/2.0; 
            factor= (factor>0? 1.0/factor : 1.0);
            Cross_CR(f,g,factor,NYMQ,NXMQ,r,h,&maxcc,&IP,&JP);
            for(ii=0;ii<NYMQ;ii++)
                for(jj=0;jj<NXMQ;jj++)
                    HartR[ii][jj]*=r[ii][jj];
            
            maxcc=0.0;
            for(ii=0;ii<NYmn;ii++)
            {
                for(jj=0;jj<NXmn;jj++)
                {
                    Pmn[ii][jj]=HartR[ii+NYMQ/2-NSRY][jj+NXMQ/2-NSRX];
                    Dmn[ii][jj]=Pmn[ii][jj];
        
                    if(Dmn[ii][jj]>maxcc)
                    {
                        maxcc=Dmn[ii][jj];
                        IP=ii-NYmn/2;
                        JP=jj-NXmn/2;
                    }
                }
            }
            if(npass==0)    cc_first=maxcc;
            ex+=JP;
            ey+=IP;
            npass++;
        }while(npass<maxpass & (IP!=0||JP!=0));
        
        if(npass<maxpass||maxcc>=cc_first*0.9)   success_flag=1;
        else    success_flag=0;      
        
        if(success_flag==1 && maxcc>mincc 
            &&ex>=Umin && ex<=Umax && -ey>=Vmin && -ey<=Vmax)
        {
            if(subpixel_method==1)
                Iter_Peak_Fit(Dmn, NYmn, NXmn, (NYmn>>1)+IP,(NXmn>>1)+JP,&ssx,&ssy);
            else
                Peak_Fit(Dmn, (NYmn>>1)+IP,(NXmn>>1)+JP,&ssx,&ssy);
                
            if(ssx!=2000. && ssy!=2000.)
            {
                u[i]=ex+ssx;    v[i]=-ey-ssy;   c[i]=maxcc;
            }
            else
            {
                u[i]=0.0;   v[i]=0.0;   c[i]=0.0;
            }
        }        
        else
        {
            u[i]=0.0;   v[i]=0.0;   c[i]=0.0;
        }
        
        if(OutC_flag>0)
        {
            for(ii=0;ii<NYmn;ii++)
            {
                for(jj=0;jj<NXmn;jj++)
                {
                    cout[ii*NXmn*nci*nri + i*NXmn+jj]=Dmn[ii][jj];
                }
            }
        }
    }
    
}

