classdef NUFFT2D
    
    properties
        S
        J
        K
        N
        Ofactor
        prefilter_2D
        St
    end
    
    methods
        function NFT = NUFFT2D(kloc, fn, prefilter_2D, J,K,N,Ofactor);
           NFT.S = Create_NUFFT2D_matrix_fast(kloc, transpose(fn), J,K,N,Ofactor);
           %NFT.S = Create_NUFFT2D_matrix_int(kloc, fn, J,K,N,Ofactor);
           
           NFT.J = J;
           NFT.K = K;
           NFT.N = N;
           NFT.Ofactor = Ofactor;
           NFT.prefilter_2D = prefilter_2D;   
           NFT.St = transpose(NFT.S);
        end
        
        function data = A(NFT,x)
            paddedimage = zeros(NFT.K,NFT.K);
            diff=NFT.K-NFT.N;
            paddedimage(diff/2+1:end-diff/2,diff/2+1:end-diff/2) = x;
            paddedimage = paddedimage.*NFT.prefilter_2D;    
            KFFTimag=fftshift(fft2(fftshift(paddedimage)));
            data = NFT.S*KFFTimag(:);
        end
        
        function data = At(NFT,x)
            gdata = reshape(NFT.St*x,NFT.K,NFT.K);
            gdata = fftshift(ifft2(fftshift(gdata)))*NFT.K*NFT.K;
            gdata=gdata.*prefilter_2D;
            diff=NFT.K-NFT.N;
            Ax =gdata(diff/2+1:end-diff/2,diff/2+1:end-diff/2); 
        end
        
            
            
    end
end


function S=Create_NUFFT2D_matrix_int(kloc,fn,J,K,N,Ofactor)
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


function S=Create_NUFFT2D_matrix_fast(kloc,fn,J,K,N,Ofactor)
kloc=(kloc+N/2+1i*N/2)*K/N;

len=length(kloc);

S = sparse(len,K^2);

for i=1:1:length(kloc)
      
       imgx = zeros(1,K);
       imgy = zeros(K,1);
        
       kxvalue =real(kloc(i));
       kyvalue=imag(kloc(i));
      
       xvalues = [floor(kxvalue-J/2):ceil(kxvalue+J/2)];
       yvalues = [floor(kyvalue-J/2):ceil(kyvalue+J/2)];
       
       ixvalue = abs(kxvalue - xvalues);
       iyvalue = abs(kyvalue - yvalues);

       index = find(ixvalue < J/2-1/Ofactor);
       interpx = zeros(size(ixvalue));
       floorx = floor(ixvalue*Ofactor);
       resx = (ixvalue-floorx./Ofactor)*Ofactor;
       interpx(index) = fn(floorx(index)+1).*(1-resx(index)) + fn(floorx(index)+2).*resx(index);
       
       index = find(iyvalue < J/2-1/Ofactor);
       interpy = zeros(size(iyvalue));
       floory = floor(iyvalue*Ofactor);
       resy = (iyvalue-floory./Ofactor)*Ofactor;
       interpy(index) = fn(floory(index)+1).*(1-resy(index)) + fn(floory(index)+2).*resy(index);
  
       xtemp = xvalues; xtemp = xtemp + K*(xtemp<0);xtemp = xtemp - K*(xtemp>=K);
       ytemp = yvalues; ytemp = ytemp + K*(ytemp<0);ytemp = ytemp - K*(ytemp>=K);

       imgx(xtemp+1) = interpx;
       imgy(ytemp+1) = interpy;
       img = imgy*imgx;
       S(i,:) = sparse(img(:));
    end
end