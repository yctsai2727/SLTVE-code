% Written by Shing-yu Leung
% Purpose: 5th order WENO reconstruction to determine u+,u-,v+,v- for HJ equation
%	phi_t+H(phi,D phi)=0
% Input:
%	phi: the function
%	xm1: phi_{i-1,j}
%	xp1: phi_{i+1,j}
%	ym1: phi_{i,j-1}
%	yp1: phi_{i,j+1}
%	all above matrices are with same size
% Output:
%	u+ (up), u- (um), v+ (vp), v- (vm)

function [up,um,vp,vm]=WENO2(phi,xm1,xp1,ym1,yp1,dx,dy)

[m,n]=size(phi);

xm2=[xm1(1,:);xm1(1:m-1,:)];
xm3=[xm2(1,:);xm2(1:m-1,:)];
xp2=[xp1(2:m,:);xp1(m,:)];
xp3=[xp2(2:m,:);xp2(m,:)];
ym2=[ym1(:,1),ym1(:,1:n-1)];
ym3=[ym2(:,1),ym2(:,1:n-1)];
yp2=[yp1(:,2:n),yp1(:,n)]; 
yp3=[yp2(:,2:n),yp2(:,n)]; 

order=5;

if (order==1)
    up=(xp1-phi)/dx;
    um=(phi-xm1)/dx;
    vp=(yp1-phi)/dy;
    vm=(phi-ym1)/dy;
    
    return;
end

% x-direction

v1=(xm2-xm3)/dx;
v2=(xm1-xm2)/dx;
v3=(phi-xm1)/dx;
v4=(xp1-phi)/dx;
v5=(xp2-xp1)/dx;

phix1=v1/3-7/6*v2+11/6*v3;
phix2=-v2/6+5/6*v3+v4/3;
phix3=v3/3+5/6*v4-v5/6;

s1=13/12*(v1-2*v2+v3).^2+1/4*(v1-4*v2+3*v3).^2;
s2=13/12*(v2-2*v3+v4).^2+1/4*(v2-v4).^2;
s3=13/12*(v3-2*v4+v5).^2+1/4*(3*v3-4*v4+v5).^2;

tol=1e-6*max(max(max(max(v1.^2,v2.^2),v3.^2),v4.^2),v5.^2)+eps;

alpha1=0.1./(s1+tol).^2;
alpha2=0.6./(s2+tol).^2;
alpha3=0.3./(s3+tol).^2;

alpha1=[zeros(2,n);alpha1(3:m,:)];
alpha2=[zeros(2,n);alpha2(3:m-1,:);zeros(1,n)];
alpha3=[alpha3(1:m-1,:);zeros(1,n)];

um=alpha1.*phix1+alpha2.*phix2+alpha3.*phix3;
um=um./(alpha1+alpha2+alpha3);

w1=(xp3-xp2)/dx;
w2=v5;
w3=v4;
w4=v3;
w5=v2;

v1=w1;
v2=w2;
v3=w3;
v4=w4;
v5=w5;

phix1=v1/3-7/6*v2+11/6*v3;
phix2=-v2/6+5/6*v3+v4/3;
phix3=v3/3+5/6*v4-v5/6;

s1=13/12*(v1-2*v2+v3).^2+1/4*(v1-4*v2+3*v3).^2;
s2=13/12*(v2-2*v3+v4).^2+1/4*(v2-v4).^2;
s3=13/12*(v3-2*v4+v5).^2+1/4*(3*v3-4*v4+v5).^2;

tol=1e-6*max(max(max(max(v1.^2,v2.^2),v3.^2),v4.^2),v5.^2)+eps;

alpha1=0.1./(s1+tol).^2;
alpha2=0.6./(s2+tol).^2;
alpha3=0.3./(s3+tol).^2;

alpha1=[alpha1(1:m-2,:);zeros(2,n)];
alpha2=[zeros(1,n);alpha2(2:m-2,:);zeros(2,n)];
alpha3=[zeros(1,n);alpha3(2:m,:)];

up=alpha1.*phix1+alpha2.*phix2+alpha3.*phix3;
up=up./(alpha1+alpha2+alpha3);

% y-direction

v1=(ym2-ym3)/dy;
v2=(ym1-ym2)/dy;
v3=(phi-ym1)/dy;
v4=(yp1-phi)/dy;
v5=(yp2-yp1)/dy;

phix1=v1/3-7/6*v2+11/6*v3;
phix2=-v2/6+5/6*v3+v4/3;
phix3=v3/3+5/6*v4-v5/6;

s1=13/12*(v1-2*v2+v3).^2+1/4*(v1-4*v2+3*v3).^2;
s2=13/12*(v2-2*v3+v4).^2+1/4*(v2-v4).^2;
s3=13/12*(v3-2*v4+v5).^2+1/4*(3*v3-4*v4+v5).^2;

tol=1e-6*max(max(max(max(v1.^2,v2.^2),v3.^2),v4.^2),v5.^2)+eps;

alpha1=0.1./(s1+tol).^2;
alpha2=0.6./(s2+tol).^2;
alpha3=0.3./(s3+tol).^2;

alpha1=[zeros(m,2),alpha1(:,3:n)];
alpha2=[zeros(m,2),alpha2(:,3:n-1),zeros(m,1)];
alpha3=[alpha3(:,1:n-1),zeros(m,1)];

vm=alpha1.*phix1+alpha2.*phix2+alpha3.*phix3;
vm=vm./(alpha1+alpha2+alpha3);

w1=(yp3-yp2)/dy;
w2=v5;
w3=v4;
w4=v3;
w5=v2;

v1=w1;
v2=w2;
v3=w3;
v4=w4;
v5=w5;

phix1=v1/3-7/6*v2+11/6*v3;
phix2=-v2/6+5/6*v3+v4/3;
phix3=v3/3+5/6*v4-v5/6;

s1=13/12*(v1-2*v2+v3).^2+1/4*(v1-4*v2+3*v3).^2;
s2=13/12*(v2-2*v3+v4).^2+1/4*(v2-v4).^2;
s3=13/12*(v3-2*v4+v5).^2+1/4*(3*v3-4*v4+v5).^2;

tol=1e-6*max(max(max(max(v1.^2,v2.^2),v3.^2),v4.^2),v5.^2)+eps;

alpha1=0.1./(s1+tol).^2;
alpha2=0.6./(s2+tol).^2;
alpha3=0.3./(s3+tol).^2;

alpha1=[alpha1(:,1:n-2),zeros(m,2)];
alpha2=[zeros(m,1),alpha2(:,2:n-2),zeros(m,2)];
alpha3=[zeros(m,1),alpha3(:,2:n)];

vp=alpha1.*phix1+alpha2.*phix2+alpha3.*phix3;
vp=vp./(alpha1+alpha2+alpha3);
