function [u_field,v_field]=DataPreprocessing(u_2014,v_2014);
% Both u_2014 and v_2014 are 4D datasets. The first and second dimension correspond to
% the longitude and latitude, respectively, both with an equal interval of 1/3 degree. The third dimension can be simply set to be constant 1.
% The last dimension corresponds to the date in 2014, with an equal
% interval of 5 days.

%%%%%% focus on E180-E230 & S17-N8 %%%%%%%%%%%%%
for k=1:72
   u(:,:,k)=u_2014(481:631,265:-1:190,1,k);
   v(:,:,k)=v_2014(481:631,265:-1:190,1,k);
   u_temp(:,:,5*k-4)=u(:,:,k);
   v_temp(:,:,5*k-4)=v(:,:,k);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%Interpolate the velocity data to obtain data for each day %%%
for k=1:71
   u_temp(:,:,5*k-3)=0.8*u_temp(:,:,5*k-4)+0.2*u_temp(:,:,5*k+1); %linear interpolation 
   u_temp(:,:,5*k-2)=0.6*u_temp(:,:,5*k-4)+0.4*u_temp(:,:,5*k+1);
   u_temp(:,:,5*k-1)=0.4*u_temp(:,:,5*k-4)+0.6*u_temp(:,:,5*k+1);
   u_temp(:,:,5*k)  =0.2*u_temp(:,:,5*k-4)+0.8*u_temp(:,:,5*k+1);
   
   v_temp(:,:,5*k-3)=0.8*v_temp(:,:,5*k-4)+0.2*v_temp(:,:,5*k+1); %linear interpolation 
   v_temp(:,:,5*k-2)=0.6*v_temp(:,:,5*k-4)+0.4*v_temp(:,:,5*k+1);
   v_temp(:,:,5*k-1)=0.4*v_temp(:,:,5*k-4)+0.6*v_temp(:,:,5*k+1);
   v_temp(:,:,5*k)  =0.2*v_temp(:,:,5*k-4)+0.8*v_temp(:,:,5*k+1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%Interpolate the velocity data to obtain data for every half day%%%%%%%%%%
for k=1:356
   u_temp1(:,:,2*k-1)=u_temp(:,:,k);
   v_temp1(:,:,2*k-1)=v_temp(:,:,k);
end
for k=1:355
   u_temp1(:,:,2*k)=(u_temp1(:,:,2*k+1)+u_temp1(:,:,2*k-1))/2;
   v_temp1(:,:,2*k)=(v_temp1(:,:,2*k+1)+v_temp1(:,:,2*k-1))/2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%Extend the velocity field to four times a day%%%%%%%%%%
for k=1:711
   u_temp2(:,:,2*k-1)=u_temp1(:,:,k);
   v_temp2(:,:,2*k-1)=v_temp1(:,:,k);
end
for k=1:710
   u_temp2(:,:,2*k)=(u_temp2(:,:,2*k+1)+u_temp2(:,:,2*k-1))/2;
   v_temp2(:,:,2*k)=(v_temp2(:,:,2*k+1)+v_temp2(:,:,2*k-1))/2;
end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%Interpolate to obtain data with equal longitude and latitude interval of 1/6 degree%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n,k]=size(u_temp2);
u_temp3(1:2:2*m-1,1:n,1:k)=u_temp2;
v_temp3(1:2:2*m-1,1:n,1:k)=v_temp2;
u_temp3(2:2:2*m-2,:,:)=(u_temp3(1:2:2*m-3,:,:)+u_temp3(3:2:2*m-1,:,:))/2;
v_temp3(2:2:2*m-2,:,:)=(v_temp3(1:2:2*m-3,:,:)+v_temp3(3:2:2*m-1,:,:))/2;

u_temp4(1:2*m-1,1:2:2*n-1,1:k)=u_temp3;
v_temp4(1:2*m-1,1:2:2*n-1,1:k)=v_temp3;
u_temp4(:,2:2:2*n-2,:)=(u_temp4(:,1:2:2*n-3,:)+u_temp4(:,3:2:2*n-1,:))/2;
v_temp4(:,2:2:2*n-2,:)=(v_temp4(:,1:2:2*n-3,:)+v_temp4(:,3:2:2*n-1,:))/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%Interpolate to obtain data with equal longitude and latitude interval of 1/12 degre%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [m,n,k]=size(u_temp4);
% u_temp5(1:2:2*m-1,1:n,1:k)=u_temp4;
% v_temp5(1:2:2*m-1,1:n,1:k)=v_temp4;
% u_temp5(2:2:2*m-2,:,:)=(u_temp5(1:2:2*m-3,:,:)+u_temp5(3:2:2*m-1,:,:))/2;
% v_temp5(2:2:2*m-2,:,:)=(v_temp5(1:2:2*m-3,:,:)+v_temp5(3:2:2*m-1,:,:))/2;

% u_temp6(1:2*m-1,1:2:2*n-1,1:k)=u_temp5;
% v_temp6(1:2*m-1,1:2:2*n-1,1:k)=v_temp5;
% u_temp6(:,2:2:2*n-2,:)=(u_temp6(:,1:2:2*n-3,:)+u_temp6(:,3:2:2*n-1,:))/2;
% v_temp6(:,2:2:2*n-2,:)=(v_temp6(:,1:2:2*n-3,:)+v_temp6(:,3:2:2*n-1,:))/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%Convert the unit of data as 5 degrees per 10 days%%%%%%%%%%%%%%%%%%%%%%%%

Radius=6.371*10^6;
u_field=u_temp4*24*3600*360/pi/Radius/2;
v_field=v_temp4*24*3600*360/pi/Radius/2;

u_field=u_field(:,:,1:250);  % Consider data of the first 50 days
v_field=v_field(:,:,1:250);




