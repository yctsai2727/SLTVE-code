function [A]=LowRankDecoder(M,r,m,n)
    M=squeeze(M);
    [a,b]=size(M);
    assert(a==r&&b==m+n+1,"invalid format of input matrix in LowRankDecoder");
    S=diag(M(:,1));
    U=M(:,2:m+1)';
    V=M(:,m+2:m+n+1);
    A=U*S*V;
