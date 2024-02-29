function [d]=MatTrajMetric(t1,t2,K,DistriMetric,Decoder)
    s=0;
    for i=1:K
        if nargin==5
            s+=DistriMetric(Decoder(t1,i),Decoder(t2,i))^2;
        elseif nargin==4
            s+=DistriMetric(t1(i,:,:),t2(i,:,:))^2;
        else
            s+=sum(sum((t1(i,:,:)-t2(i,:,:)).^2));
        end
    end
    d=sqrt(s/K);