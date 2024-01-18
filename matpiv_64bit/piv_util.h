/* The data type of the image pixel. Usually, camera is 16 bit (unsigned short). */
/* But this can be changed according to the image type (sometimes we might       */
/* use 8 bit camera, then it should be unsigned char).                           */

#define IMG_DATA float
#define BYTE unsigned char
typedef struct FCOMPLEX {float r,i;} fcomplex;

/* Calculate the cross correlation of two sub-images f and g, both with dimension (ni, nj);          */
/* (1). h is their co-spectrum, and r is the resultant correlation plane                             */
/* (2). The returned correlation r plane is normalized to 0~1.0, i.e. r=1 means perfect correlation  */
/* (3). The returned correlation plane has been arranged such that the correlation for zero          */
/*      displacement is located at the center of the matrix. i.e. center index =floor(N/2), where    */
/*      indices range from 0 to N-1 (C convention). The resultant corrletion plane r is identical to */
/*      r=imfilter(f,g,'circular')/sqrt( sum(f(:).^2) * sum(g(:).^2) ) im MatLab                     */
/* (4). The maximal correlation coefficent and the relative location (to the center of the matrix)   */
/*      is also returned such that the displacement of partern from f to g is IP and JP (in y and x) */
void Cross_CR(fcomplex **f, fcomplex **g, float maxcorr, int ni, int nj, float **r, fcomplex **h, float*maxcc, int*IP, int *JP);


/* 2D matrix allocation and deallocation in C */
void ** matrix(long nrow, long ncol, size_t nbytes);
void free_matrix(void** m);
void ** vec_mtx(void* data, long nrow, long ncol, size_t nbytes);


/**************************************************************************/
/*The followings are definitions and functions related to complex calculus*/


fcomplex Cadd(fcomplex a, fcomplex b);
/*Returns the complex sum of two complex numbers.*/

fcomplex Csub(fcomplex a, fcomplex b);
/*Returns the complex difference of two complex numbers.*/

fcomplex Cmul(fcomplex a, fcomplex b);
/*Returns the complex product of two complex numbers.*/

fcomplex Cdiv(fcomplex a, fcomplex b);
/*Returns the complex quotient of two complex numbers.*/

fcomplex Csqrt(fcomplex z);
/*Returns the complex square root of a complex number.*/

fcomplex Conjg(fcomplex z);
/*Returns the complex conjugate of a complex number.*/

float Cabs(fcomplex z);
float Cabs2(fcomplex z);
/*Returns the absolute value (modulus)of a complex number.*/

fcomplex Complex(float re, float im);
/*Returns a complex umber with speci .ed real a d imaginary parts.*/

fcomplex RCmul(float x, fcomplex a);
/*Returns the complex product of a real number and a complex number.*/

fcomplex EXPComplex(double a);
/*Returns exp(i*a), where i=sqrt(-1)*/

/**************************************************************************/


/* PIV interrogation tools*/
void Peak_Fit(float** cor, int IP, int JP, float *xp, float *yp);
void Peak_Fit2(float** cor, int IP, int JP, float *xp, float *yp);
void Iter_Peak_Fit(float** cor, int ni, int nj, int IP, int JP, float *xp, float *yp);

int FFT_Peak_Fit(fcomplex **hh, int ni, int nj, float ex, float ey, float *xp, float *yp, int mpass, float eps);


