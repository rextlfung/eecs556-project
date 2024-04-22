function X=NFT_n(x,N,A1,B1,K)

X=[];
for u=1:1:K
    X(u)=sum(sum(x.*(reshape(A1(u,:),N,1)*B1(u,:))));
    %X(u)=A(u,:)*(x*(B(u,:))');
end


end