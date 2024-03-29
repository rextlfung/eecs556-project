function kloc=getpolar(L,N)

thc = linspace(0, pi-pi/L, L);

kloc=[];
n=1;
% full mask
for ll = 1:L

	if ((thc(ll) <= pi/4) | (thc(ll) > 3*pi/4))
		yr = (tan(thc(ll))*(-N/2+1:N/2-1))+N/2+1;
    	 nn = 1:N-1;
      	kloc=[kloc,nn+1+yr*i];
        %n=n+1;
      
  else 
		xc = (cot(thc(ll))*(-N/2+1:N/2-1))+N/2+1;
		 nn = 1:N-1;
			kloc=[kloc,xc+nn*i+i];
            n=n+1;
		
	end

end


end