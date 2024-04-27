function [d]=femd(d1,d2,level)
    d1=squeeze(d1);
    d2=squeeze(d2);
    assert(size(d1)==size(d2),"nonconformant arguments in femd");
    [m,n]=size(d1);
    if nargin==2
        level=4;
    end
    c=fwt2(d1-d2,'sym5',level);
    c=abs(c);
    [m_c,n_c] = size(c);
    m_cut = floor(m_c/2);
    n_cut = floor(n_c/2);
    d = (2^(-2*level))*(sum(sum(c(1:m_cut,n_cut+1:n_c)))+sum(sum(c(m_cut+1:m_c,:))));
end