function [d]=DisplacedWeightedTotalVariation(d1,d2,dirc,dx,dy,l) #dirc=1: horizontal shifting, dirc=2: vertical shifting
    dist1=squeeze(d1{1});
    dist2=squeeze(d2{1});
    [m,n]=size(dist1);
    global curr;
    global W;
    if size(W,1)!=m || curr(1,1)!=d1{2} || curr(2,1)!=d1{3}  
        if size(W,1)!=m
            W=zeros(m,n);
        end
        curr(1,1)=d1{2};
        curr(2,1)=d1{3};
        for i=1:m
            for j=1:n
                W(i,j)=sqrt(((i-curr(1,1))*dx)^2+((j-curr(2,1))*dy)^2);
            end
        end
    end
    assert(size(dist1)==size(dist2),"nonconformant arguments in DisplacedTotalVariation");
    if nargin == 5
        l=1;
    end
    if dirc==1
        d=sum(sum(abs(dist1(:,1:n-1)-dist2(:,2:n)).^l.*W(:,1:n-1)));
    elseif
        d=sum(sum(abs(dist1(1:m-1,:)-dist2(2:m,:)).^l.*W(1:m-1,:)));
    end
end
global curr;
curr=zeros(2,1);
global W;
W=zeros(1,1);