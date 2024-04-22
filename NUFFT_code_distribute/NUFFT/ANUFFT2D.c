#include <math.h>
#include <stdio.h>


void gridlut(sk_real,sk_imag,K,N,J,kx1,ky1,
	sg_real,sg_imag,fn,Ofactor,len,dcf)

double *sg_real;		
double *sg_imag; 	
int K,N,J;		
double *kx1;		
double *ky1;
double *sk_real;	
double *sk_imag;	  	 

double *fn;	
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
double interpx,interpy,interp;
int xtemp,ytemp;
double resx,resy;
double e;



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
              interpx = 0;
           else
            { 
                floorx = floor(xvalue*Ofactor);
                e=floorx;
              resx = (xvalue-e/(double)Ofactor)*(double)Ofactor;
               interpx = *(fn+floorx)*(1-resx)+*(fn+floorx+1)*resx;
            }
            xtemp = x;
            if(xtemp<=0) xtemp = xtemp + K;
            else
                if(xtemp>K) xtemp = xtemp - K;
          
         for (y=miny;y<=maxy;y++)
          {  yvalue = fabs(kyvalue - y); 
          
           if(yvalue>(double)J/2+1/(double)Ofactor) 
             interpy = 0;
           else
          {
               floory = floor(yvalue*Ofactor);
              e=floory;
              resy = (yvalue-(double)floory/(double)Ofactor)*(double)Ofactor;
              interpy = *(fn+floory)*(1-resy)+*(fn+floory+1)*resy;
          }
            ytemp = y;
            
            //if(ytemp>0.00)printf("y=%d\n",ytemp);
            if(ytemp<=0) {ytemp = ytemp + K;} 
            
            else
                if(ytemp>K) ytemp = ytemp-K; 
             
            interp=interpx*interpy;
            
            
            *(sg_real+(xtemp-1)+(ytemp-1)*K)+=*(sk_real+i-1)*interp**(dcf+i-1);
            
          *(sg_imag+(xtemp-1)+(ytemp-1)*K)+=*(sk_imag+i-1)*interp**(dcf+i-1);
        
            
         }
       }     
         
    }  


}




