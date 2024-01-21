function M2_spar=right_mat(miu,m,n)

count=1;
temp=5*(m-2)*(n-2);  % NO. of non-zero entries in the sparse matrix
temp=zeros(temp,3);

for i=1:m-2
    for j=1:n-2
        index=(i-1)*(n-2)+j; % queue of mesh points       
            
            temp(count,1)=index;
            temp(count,2)=index+2*i-1;
            temp(count,3)=0.5*miu;
            count=count+1;
            
            temp(count,1)=index;
            temp(count,2)=n*i+j;
            temp(count,3)=0.5*miu;
            count=count+1;
            
            temp(count,1)=index;
            temp(count,2)=n*i+j+1;
            temp(count,3)=1-2*miu;
            count=count+1;
            
            temp(count,1)=index;
            temp(count,2)=n*i+j+2;
            temp(count,3)=0.5*miu;
            count=count+1;
            
            temp(count,1)=index;
            temp(count,2)=n*(i+1)+j+1;
            temp(count,3)=0.5*miu;
            count=count+1;         
      
    end
end

temp=temp(1:count-1,:);
temp(count,:)=[(m-2)*(n-2),m*n,0];
M2_spar=spconvert(temp);