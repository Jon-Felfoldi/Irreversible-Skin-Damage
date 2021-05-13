function tc = TwoHeat_Tc_V5_new(Nu,Tc,g,z)
%Numerical solution of the Heat Equation with Critical Time

IC = zeros(100,1); %Initial Condition
D1=1;
NU=Nu;
G=g;
Z=z;
N=length(IC);

tend=100000000;

X=linspace(0,1,N);
dx1=X(2)-X(1);

options =odeset('AbsTol',1e-6,'RelTol',1e-6,'Events',@(t,y)TcHit(t,y,Tc,N,Z,dx1));
sol = ode23s(@(t,y) heat2(t,y,dx1,D1,N,NU,G,Z),[0 tend], IC,options);

tc=sol.xe;
if isempty(tc)
    tc=tend;
end
% figure(10)
% surf(X,sol.x,sol.y','Edgecolor','none')
% view(0,90)
end

function dy = heat2(~,y,dx1,D1,N,NU,G,Z)
dy=zeros(N,1);
y(1)=(-y(3)+4*y(2)+2*dx1*NU)/(3+2*dx1*NU);
y(N) =(4*y(N-1)-y(N-2))/(3+2*dx1*Z);

for i=2:N-1
    dy(i) = D1* 1/dx1^2 * (y(i+1)-2*y(i)+y(i-1))-G*(y(i));
end

end

function [position,isterminal,direction] = TcHit(~,y,Tc,N,Z,dx1)
y(N) =(4*y(N-1)-y(N-2))/(3+2*dx1*Z);
position = y(N)-Tc;% The value that we want to be zero
isterminal = 1;  % Halt integration 
direction = 0;   % The zero can be approached from either direction
end