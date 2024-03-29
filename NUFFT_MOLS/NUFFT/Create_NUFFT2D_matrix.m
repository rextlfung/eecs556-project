function S=Create_NUFFT2D_matrix(kloc,fn,J,K,N,Ofactor)

kloc=(kloc-1-1i)*K/N;

len=length(kloc);
S = sparse(len,K^2);

for i=1:1:len
      
       kxvalue =real(kloc(i))+1;
       kyvalue=imag(kloc(i))+1;
      
       minx = floor(kxvalue-J/2);
       maxx =ceil(kxvalue+J/2);
       miny = floor(kyvalue-J/2);
       maxy =ceil(kyvalue+J/2);
       for x=minx:1:maxx
    
              xvalue = abs(kxvalue - x);
         
             if xvalue>J/2-1/Ofactor
              interpx = 0;
             else
            
                floorx = floor(xvalue*Ofactor);
                
              resx = (xvalue-floorx/Ofactor)*Ofactor;
              
               interpx = fn(floorx+1)*(1-resx)+fn(floorx+2)*resx;
             end
             
            xtemp = x;
            if xtemp<=0
                xtemp = xtemp + K;
            else
                if(xtemp>K)
                    xtemp = xtemp - K;
                end
            end
          
         for y=miny:1:maxy
           yvalue = abs(kyvalue - y); 
          
           if yvalue>J/2-1/Ofactor 
             interpy = 0;
           else
          
               floory = floor(yvalue*Ofactor);
             
              resy = (yvalue-floory/Ofactor)*Ofactor;
              interpy = fn(floory+1)*(1-resy)+fn(floory+2)*resy;
           end
            ytemp = y;
       
            if ytemp<=0
                ytemp = ytemp + K;
            else
                if ytemp>K 
                    ytemp = ytemp-K; 
                end
            end
             
            interp=interpx*interpy;
            S(i, xtemp +(ytemp-1)*K)=interp;
            
            
         end
       end
end
end