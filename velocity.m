function [func]=velocity(c)
  if c==1
    func=@velo1;
  elseif c==2
    func=@velo2;
  else
    func=@iden;
  end

function [u,v]=velo1(x,y,t,tf)
  eps=0.1; 
  omega=2*pi/10;
  A=0.1;
  a=eps*sin(omega*t);
  b=1-2*a;
  f=a*x.^2+b*x;
  u=-pi*A*sin(pi*f).*cos(pi*y);
  v=pi*A*cos(pi*f).*sin(pi*y).*(2*a*x+b);

function [u,v]=velo2(x,y,t,tf)
  eps=0.01;
  u=y.*(abs(sqrt(x.*x+y.*y)-0.5)<eps);
  v=-x.*(abs(sqrt(x.*x+y.*y)-0.5)<eps);

function [u,v]=iden(x,y,t,tf)
  u=0;
  v=0;