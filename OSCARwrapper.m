function [velo] = OSCARwrapper(dt)
    global Trigger = false;
    global U = [];
    global V = [];

    if !Trigger
        ['Loading and Processing Data...']
        addpath('./OSCAR files/');

        load u_2014.mat;
        load v_2014.mat;

        [U,V]=DataPreprocessing(u_2014,v_2014);
        Trigger = true;
    end
    velo = @(x,y,t,tf) Field(round(2*t/dt)+1);
end

function [u,v] = Field(k)
    global U;
    global V;
    u = U(:,:,k);
    v = V(:,:,k);
end
