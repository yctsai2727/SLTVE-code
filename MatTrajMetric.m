function [d]=MatTrajMetric(t1,t2,DistriMetric,Decoder)
    t1=squeeze(t1);
    t2=squeeze(t2);
    [K,m,n]=size(t1);
    assert(size(t1)==size(t2),"nonconformant arguments in MatTrajMetric");
    s=0;
    for i=1:K
        if nargin==4
            s+=DistriMetric(Decoder(t1(i,:,:)),Decoder(t2(i,:,:)))^2;
        elseif nargin==3
            s+=DistriMetric(t1(i,:,:),t2(i,:,:))^2;
        else
            s+=norm(t1(i,:,:),t2(i,:,:),2)^2;
        end
    end
    d=sqrt(s/K);