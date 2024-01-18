#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "piv_util.h"
#include "fftw.h"

/*1111111111111111111111111111111111111111111111111111111111111111111111111111111111*/

/* Calculate the cross correlation of two sub-images f and g, both with dimension (ni, nj);          */
/* (1). h is their co-spectrum, and r is the resultant correlation plane                             */
/* (2). The returned correlation r plane is normalized to 0~1.0, i.e. r=1 means perfect correlation  */
/* (3). The returned correlation plane has been arranged such that the correlation for zero          */
/*      displacement is located at the center of the matrix. i.e. center index =floor(N/2), where    */
/*      indices range from 0 to N-1 (C convention). The resultant corrletion plane r is identical to */
/*      r=imfilter(f,g,'circular')/sqrt( sum(f(:).^2) * sum(g(:).^2) ) im MatLab                     */
/* (4). The maximal correlation coefficent and the relative location (to the center of the matrix)   */
/*      is also returned such that the displacement of partern from f to g is IP and JP (in y and x) */


void Cross_CR(fcomplex **f, fcomplex **g, float maxcorr, int ni, int nj, float **r, fcomplex **h, float*maxcc, int*IP, int *JP)
{
	int i,j,nci,ncj,ii,jj;
	float factor;
    fcomplex **ff,**gg;
    fftwnd_plan p;
    
    factor=1.0/(float)ni/(float)nj*maxcorr;
    nci = (ni>>1);
    ncj = (nj>>1);
    
    ff=(fcomplex**)matrix(ni,nj,sizeof(fcomplex));
    gg=(fcomplex**)matrix(ni,nj,sizeof(fcomplex));    
    
    p=fftw2d_create_plan(ni,nj,FFTW_FORWARD,FFTW_ESTIMATE);
	fftwnd_one(p,*f,*ff);
    fftwnd_one(p,*g,*gg);
    fftwnd_destroy_plan(p);   
	
	for(i=0;i<ni;i++)
	{
		for(j=0;j<nj;j++)
		{
            h[i][j]=Cmul( Conjg(ff[i][j]), gg[i][j] );
		}
	}	
    
    p=fftw2d_create_plan(ni,nj,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftwnd_one(p,*h,*ff);
    fftwnd_destroy_plan(p);

    *maxcc=0.0;
    *IP=0;  *JP=0;
	for(i=0;i<ni;i++)
	{
        ii=(i+nci)%ni;
		for(j=0;j<nj;j++)
		{
            jj=(j+ncj)%nj;
            r[ii][jj]=ff[i][j].r* factor;
            if(*maxcc<r[ii][jj])
            {
                *maxcc=r[ii][jj];
                *IP=ii;
                *JP=jj;
            }
		}
	}
    
    *IP=*IP-(ni>>1);
    *JP=*JP-(nj>>1);

    free_matrix(ff);
    free_matrix(gg);	
}
/*1111111111111111111111111111111111111111111111111111111111111111111111111111111111*/

/*2222222222222222222222222222222222222222222222222222222222222222222222222222222222*/

/*  Allocate a 2-D matrix with size [nrow][ncol], and            */
/*  the sizeof(DATATYPE) = nbytes                                 */
void ** matrix(long nrow, long ncol, size_t nbytes)
{
	BYTE ** m;
	long i;
	
	m=malloc( sizeof(void*)*nrow );
	if ( !m ) return NULL;
	
	m[0]=malloc( (size_t)(nrow*ncol)*nbytes );
	if ( !m[0] ) 
	{
		free(m);
		return NULL;
	}
	
    for ( i=1; i<nrow; i++ )
	    m[i] = m[i-1] + ncol*nbytes;
	
    
    return (void**) m;

}

void free_matrix(void** m)
{
	free(m[0]);
	free(m);
}

void ** vec_mtx(void* data, long nrow, long ncol, size_t nbytes)
{
	BYTE ** m;

	long i;
	
	m=malloc( sizeof(void*)*nrow );
	if ( !m ) return NULL;
	
	m[0]=data;
	
    for ( i=1; i<nrow; i++ )
	    m[i] = m[i-1] + ncol*nbytes;
	
    
    return (void**) m;
}


fcomplex Cadd(fcomplex a, fcomplex b)
{
	fcomplex c;
	c.r=a.r+b.r;
	c.i=a.i+b.i;
	return c;
}

fcomplex Csub(fcomplex a, fcomplex b)
{
	fcomplex c;
	c.r=a.r-b.r;
	c.i=a.i-b.i;
	return c;
}


fcomplex Cmul(fcomplex a, fcomplex b)
{
	fcomplex c;
	c.r=a.r*b.r-a.i*b.i;
	c.i=a.i*b.r+a.r*b.i;
	return c;
}

fcomplex Complex(float re, float im)
{
	fcomplex c;
	c.r=re;
	c.i=im;
	return c;
}

fcomplex Conjg(fcomplex z)
{
	fcomplex c;
	c.r=z.r;
	c.i = -z.i;
	return c;
}

fcomplex Cdiv(fcomplex a, fcomplex b)
{
	fcomplex c;
	float r,den;
	if (fabs(b.r) >= fabs(b.i)) {
		r=b.i/b.r;
		den=b.r+r*b.i;
		c.r=(a.r+r*a.i)/den;
		c.i=(a.i-r*a.r)/den;
	} else {
		r=b.r/b.i;
		den=b.i+r*b.r;
		c.r=(a.r*r+a.i)/den;
		c.i=(a.i*r-a.r)/den;
	}
	return c;
}

float Cabs(fcomplex z)
{
	float x,y,ans,temp;
	x=fabs(z.r);
	y=fabs(z.i);
	if (x == 0.0)
		ans=y;
	else if (y == 0.0)
		ans=x;
	else if (x > y) {
		temp=y/x;
		ans=x*sqrt(1.0+temp*temp);
	} else {
		temp=x/y;
		ans=y*sqrt(1.0+temp*temp);
	}
	return ans;
}

float Cabs2(fcomplex z)
{
	return z.r*z.r+z.i*z.i;
}

fcomplex Csqrt(fcomplex z)
{
	fcomplex c;
	float x,y,w,r;
	if ((z.r == 0.0) && (z.i == 0.0)) {
		c.r=0.0;
		c.i=0.0;
		return c;
	} else {
		x=fabs(z.r);
		y=fabs(z.i);
		if (x >= y) {
			r=y/x;
			w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
		} else {
			r=x/y;
			w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
		}
		if (z.r >= 0.0) {
			c.r=w;
			c.i=z.i/(2.0*w);
		} else {
			c.i=(z.i >= 0) ? w : -w;
			c.r=z.i/(2.0*c.i);
		}
		return c;
	}
}

fcomplex RCmul(float x, fcomplex a)
{
	fcomplex c;
	c.r=x*a.r;
	c.i=x*a.i;
	return c;
}

fcomplex EXPComplex(double a)
{
	fcomplex c;
	c.r=(float)cos(a);
	c.i=(float)sin(a);
	return c;
}
/*2222222222222222222222222222222222222222222222222222222222222222222222222222222222*/


/*3333333333333333333333333333333333333333333333333333333333333333333333333333333333*/
/* PIV interrogation tools*/

void Peak_Fit2(float** cor, int IP, int JP, float *xp, float *yp)
{
	float a11,a22,a12,b1,b2;
    
   	if( !(cor[IP][JP-1]<cor[IP][JP] || cor[IP][JP]>cor[IP][JP+1]) )
    {
        *xp=2000;
        return;
    }
   	if( !(cor[IP-1][JP]<cor[IP][JP] || cor[IP][JP]>cor[IP+1][JP]) )
    {
        *yp=2000;
        return;
    }
 
    
    a11=cor[IP][JP+1]-2*cor[IP][JP]+cor[IP][JP-1];
    a22=cor[IP+1][JP]-2*cor[IP][JP]+cor[IP-1][JP];
    a12=(cor[IP+1][JP+1]-cor[IP-1][JP+1]+cor[IP+1][JP-1]-cor[IP-1][JP-1])*0.25;
    b1=0.5*(cor[IP][JP-1]-cor[IP][JP+1]);
    b2=0.5*(cor[IP-1][JP]-cor[IP+1][JP]);
    
    *xp=(a22*b1-a12*b2)/(a11*a22-a12*a12);
    *yp=(a11*b2-a12*b1)/(a11*a22-a12*a12);
    
}


void Peak_Fit(float** cor, int IP, int JP, float *xp, float *yp)
{
	float y1,y2,y3;
	float A,C;
    
	y1=cor[ IP ][ JP-1 ];
	y2=cor[ IP ][ JP ] ;
	y3=cor[ IP ][ JP+1 ];
	
	if( y1<y2 || y2>y3 )
	{
		A=log(y1/y2);
		C=log(y2/y3);
		*xp=(C+A)/(2*(A-C));
	}
	else
	{
		*xp=2000.0;
	}

	y1=cor[ IP-1 ][ JP ];
	y2=cor[ IP ][ JP ] ;
	y3=cor[ IP+1 ][ JP ];
	
	if( y1<y2 || y2>y3 )
	{
		A=log(y1/y2);
        C=log(y2/y3);
		*yp=(C+A)/(2*(A-C));
	}
	else
	{
		*yp=2000.0;
	}

}

void Iter_Peak_Fit(float** cor, int ni, int nj, int IP, int JP, float *xp, float *yp)
{
	float y1,y2,y3,A,C,ex,ey,ex1,ey1;
    int i,j,i1,i2,i3,j1,j2,j3;
    int mpass=7, mccp, eps=0.01, np;
    fcomplex *hhx, *hhy, *rrx, *rry;
    fftw_plan px,py;
    int *idi, *idj;
        
    rrx=(fcomplex*)malloc(nj*sizeof(fcomplex));
    rry=(fcomplex*)malloc(ni*sizeof(fcomplex));
    hhx=(fcomplex*)malloc(nj*sizeof(fcomplex));
    hhy=(fcomplex*)malloc(ni*sizeof(fcomplex));
    idi=(int*)malloc(ni*sizeof(int));
    idj=(int*)malloc(nj*sizeof(int));
    
    for (i=0;i<ni;i++)  
    {
        rry[i]=Complex(cor[i][JP],0.0);
        hhy[i]=Complex(cor[i][JP],0.0);
        idi[i] = (i<ni/2 ? i : i-ni);
    }
    py = fftw_create_plan(ni, FFTW_FORWARD, FFTW_ESTIMATE|FFTW_IN_PLACE);
    fftw_one(py, hhy, NULL);
    fftw_destroy_plan(py);    
    for (j=0;j<nj;j++)  
    {
        rrx[j]=Complex(cor[IP][j],0.0);
        hhx[j]=Complex(cor[IP][j],0.0);
        idj[j] = (j<nj/2 ? j : j-nj);
    }
    px = fftw_create_plan(nj, FFTW_FORWARD, FFTW_ESTIMATE|FFTW_IN_PLACE);
    fftw_one(px, hhx, NULL);
    fftw_destroy_plan(px);    
    
  /* Sub-pixel fitting in x- direction*/
    j1= ( JP==0 ? nj-1: JP-1);    j2= JP;    j3= ( JP==nj-1 ? 0: JP+1);    
    px = fftw_create_plan(nj, FFTW_BACKWARD, FFTW_ESTIMATE|FFTW_IN_PLACE);
    np=0;   *xp=0.0;    ex1=0.0;
    do {
        y1=rrx[j1].r; y2=rrx[j2].r; y3=rrx[j3].r;        
        if(y1>y2 || y3>y2 ) 
        {
            ex1=2000.;
            break;
        }
        A=log(y1/y2);   C=log(y2/y3);    ex=(C+A)/2./(A-C);
        if(np==0)   ex1=ex;
        *xp+=ex;
        for (j=0;j<nj;j++)
            rrx[j]=Cmul(hhx[j], EXPComplex(6.28318530717959*(double)(idj[j]*(*xp)/nj)) );
        fftw_one(px, rrx, NULL);
        np++;
    }while(np<mpass && fabs(ex)>eps);
    fftw_destroy_plan(px);
    if(fabs(*xp) >1 )    *xp=ex1;
    
  /* Sub-pixel fitting in y- direction*/
    i1= ( IP==0 ? ni-1: IP-1);    i2= IP;    i3= ( IP==ni-1 ? 0: IP+1);
    py = fftw_create_plan(ni, FFTW_BACKWARD, FFTW_ESTIMATE|FFTW_IN_PLACE);
    np=0;   *yp=0.0;    ey1=0.0;
    do {
        y1=rry[i1].r; y2=rry[i2].r; y3=rry[i3].r;
        if(y1>y2 || y3>y2 ) 
        {
            ey1=2000.;
            break;
        }
        A=log(y1/y2);   C=log(y2/y3);        ey=(C+A)/2./(A-C);       
        if(np==0)   ey1=ey;
        *yp+=ey;
        for (i=0;i<ni;i++)
            rry[i]=Cmul(hhy[i], EXPComplex(6.28318530717959*(double)(idi[i]*(*yp)/ni)) );
        fftw_one(py, rry, NULL);
        np++;
    }while(np<mpass && fabs(ey)>eps);
    fftw_destroy_plan(py);
    if(fabs(*yp) >1 )    *yp=ey1;
    
    free(rrx);
    free(rry);
    free(hhx);
    free(hhy);
    free(idi);
    free(idj);
}

int FFT_Peak_Fit(fcomplex **hh, int ni, int nj, float ex, float ey, float *xp, float *yp, int mpass, float eps)
{
	int i,j,npass=0;
	int ii,jj,IP,JP,imax,jmax;
	fcomplex temp,**ff;
	float maxcc, xxp, yyp, **rr;
    fftwnd_plan p;
    float exx,eyy;
    int nci,ncj;
    
    exx=ex;
    eyy=ey;
    
    nci = (ni>>1);
    ncj = (nj>>1);
	
    ff=(fcomplex**)matrix(ni,nj,sizeof(fcomplex));
    rr=(float**)matrix(ni,nj,sizeof(float));
    
	do{
		for(i=0;i<ni;i++)
		{
			for(j=0;j<nj;j++)
			{
				ii = (i<ni/2 ? i : i-ni);
				jj = (j<nj/2 ? j : j-nj);
                temp = EXPComplex(6.28318530717959*(double)(ii*ey/ni+jj*ex/nj));                
				ff[i][j] = Cmul(hh[i][j], temp);
			}
		}
		
        p=fftw2d_create_plan(ni,nj,FFTW_BACKWARD,FFTW_ESTIMATE|FFTW_IN_PLACE);
        fftwnd_one(p,*ff,NULL);
        fftwnd_destroy_plan(p);
	
        maxcc=0.0;
        IP=0;  JP=0;
        for(i=0;i<ni;i++)
        {
            ii=(i+nci)%ni;
            for(j=0;j<nj;j++)
            {
                jj=(j+ncj)%nj;
                rr[ii][jj]=ff[i][j].r;
                if(maxcc<rr[ii][jj])
                {
                    maxcc=rr[ii][jj];
                    IP=ii;
                    JP=jj;
                }
            }
        }
        if(IP==0||JP==0||IP==ni-1||JP==nj-1)
        {
            npass=mpass;
            *xp=exx;    *yp=eyy;
            free_matrix(rr);
            free_matrix(ff);
            return;
        }
		Peak_Fit(rr, IP, JP, &xxp, &yyp);
		
		ex +=xxp;
		ey +=yyp;
		npass ++;
		
	}while(npass<mpass && (fabs(xxp)>eps || fabs(yyp) >eps));
	
	if(npass < mpass)
	{
		*xp =ex;
		*yp =ey;
	}	
	else
	{
		*xp =exx;
		*yp =eyy;
	}
    
    free_matrix(rr);
    free_matrix(ff);
}

/*3333333333333333333333333333333333333333333333333333333333333333333333333333333333*/

