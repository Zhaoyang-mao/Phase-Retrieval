function Xd=convert2diag(X)

% This function reorder the matrix X by diagonal;
% Written by Tamir Bendory; last update 28/08/2016.

%% -------------------Initialization

N=size(X,1);
Xd=zeros(size(X)); t1=Xd; t2=Xd;


%% --------------------------------------

for ii=1:N
   
    ind1=1:N-ii+1;    ind2=ii-1;
    
    t1=diag(X(ind1,ii),ii-1);
    
    if ind2~=0
    t2=diag(X(N-ind2+1:N,ii),-N+ii-1);
    end
    
   Xd=Xd+t1+t2;
   
end