function [e]=CovarianceEigen(d1,d2)
    d1=squeeze(d1);
    d2=squeeze(d2);
    assert(size(d1)==size(d2),"nonconformant arguments in CovarianceEigen");
    