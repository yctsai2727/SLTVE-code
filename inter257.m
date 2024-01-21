  

   [xx,yy]=meshgrid(xmin:dx/2:xmax,ymin:dy/2:ymax);
     xx=xx';
     yy=yy';
    FTLE_257=interp2(x',y',FTLEd',xx',yy','cubic')'; 
    

%     filename=strcat('FlowMap_',num2str(m),'_',num2str(n),'_',num2str(D0*10000));
%     print('-dpsc2',filename)
%     print('-djpeg',filename)
%     
%     filename=strcat('FTLEDiffusion_',num2str(m),'_',num2str(n),'.mat');
%     save(filename)
