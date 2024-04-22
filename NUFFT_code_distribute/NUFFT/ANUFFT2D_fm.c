#include <math.h>
#include <stdio.h>


void agridlut2d(sk_real,sk_imag,K,N,J,kx1,ky1,
	sg_real,sg_imag,fnx_real,fnx_imag,fny_real,fny_imag,Ofactor,len,dcf)

double *sg_real;		
double *sg_imag; 	
int K,N,J;		
double *kx1;		
double *ky1;
double *sk_real;	
double *sk_imag;	  	 

double *fnx_real;	
double *fnx_imag;
double *fny_real;	
double *fny_imag;	
int Ofactor;
int len;
double *dcf;
		
{
int i;		  
int minx,maxx,miny,maxy;	
double kxvalue,kyvalue;
double x,y;
double xvalue,yvalue;
int floorx,floory;
double interpx_real,interpy_real,interp_real;
double interpx_imag,interpy_imag,interp_imag;
int xtemp,ytemp;
double resx,resy;
double e;
int middle;
middle=Ofactor*J/2+5;


for (i=0;i<K*K;i++)
       {
*(sg_real+i) = 0;
*(sg_imag+i) = 0;
       }

    for(i=1;i<=len;i++)
   
   {
       kxvalue = *(kx1+i-1)+1; 
       
       kyvalue=*(ky1+i-1)+1;
       
       minx = floor(kxvalue-(double)J/2);
       maxx =ceil(kxvalue+(double)J/2);
       miny = floor(kyvalue-(double)J/2);
       maxy =ceil(kyvalue+(double)J/2);
       for (x=minx;x<=maxx;x++)
    
          {  
              
              
              xvalue = fabs(kxvalue - x);
              e=Ofactor;
             if(xvalue>(double)J/2+1/e) 
              {interpx_real = 0;
             interpx_imag = 0;}
           else
            { 
                floorx = floor(xvalue*Ofactor);
                e=floorx;
              resx = (xvalue-e/(double)Ofactor)*(double)Ofactor;
              if (x>kxvalue||x==kxvalue)
              
              {interpx_real = *(fnx_real+middle+floorx)*(1-resx)+*(fnx_real+middle+floorx+1)*resx;
               interpx_imag = *(fnx_imag+middle+floorx)*(1-resx)+*(fnx_imag+middle+floorx+1)*resx;
               }
              else
                  
                   {interpx_real = *(fnx_real+middle-floorx)*(1-resx)+*(fnx_real+middle-floorx-1)*resx;
               interpx_imag = *(fnx_imag+middle-floorx)*(1-resx)+*(fnx_imag+middle-floorx-1)*resx;
            
            }
            }
            xtemp = x;
            if(xtemp<=0) xtemp = xtemp + K;
            else
                if(xtemp>K) xtemp = xtemp - K;
          
         for (y=miny;y<=maxy;y++)
          {  yvalue = fabs(kyvalue - y); 
          
           if(yvalue>(double)J/2+1/(double)Ofactor) 
           { interpy_real = 0;interpy_imag = 0;}
           else
          {
               floory = floor(yvalue*Ofactor);
              e=floory;
              resy = (yvalue-(double)floory/(double)Ofactor)*(double)Ofactor;
              if (y>kyvalue||y==kyvalue)
              
              {interpy_real = *(fny_real+middle+floory)*(1-resy)+*(fny_real+middle+floory+1)*resy;
               interpy_imag = *(fny_imag+middle+floory)*(1-resy)+*(fny_imag+middle+floory+1)*resy;
               }
              else
                   {interpy_real = *(fny_real+middle-floory)*(1-resy)+*(fny_real+middle-floory-1)*resy;
               interpy_imag = *(fny_imag+middle-floory)*(1-resy)+*(fny_imag+middle-floory-1)*resy;
            
            }
          }
            ytemp = y;
            
            //if(ytemp>0.00)printf("y=%d\n",ytemp);
            if(ytemp<=0) {ytemp = ytemp + K;} 
            
            else
                if(ytemp>K) ytemp = ytemp-K; 
             
            interp_real=interpx_real*interpy_real-interpx_imag*interpy_imag;
            interp_imag=interpx_imag*interpy_real+interpy_imag*interpx_real;
            
            
               
          *(sg_real+(xtemp-1)+(ytemp-1)*K)+=*(sk_real+i-1)*interp_real**(dcf+i-1)-*(sk_imag+i-1)*interp_imag**(dcf+i-1);
            
          *(sg_imag+(xtemp-1)+(ytemp-1)*K)+=*(sk_imag+i-1)*interp_real**(dcf+i-1)+*(sk_real+i-1)*interp_imag**(dcf+i-1);
         
            //printf("g=%f",*(sg_real+(xtemp-1)+(ytemp-1)*K));
            
           
        
            
         }
       }     
         
    }  


}




