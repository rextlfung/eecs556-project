function x=INFT_n(X,A2,B2,N)
x=[];

for i=1:N
    AX(:,i) = X(:).*A2(:,i);
end
   x=AX.'*B2;

          
         
end 