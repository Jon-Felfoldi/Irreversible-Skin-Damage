function tc = TwoHeat_Tc_V3(D2s,a,Nu,Tc)
%Numerical solution of the Heat Equation with Critical Time

IC = zeros(100,1); %Initial Condition
l1=a;
D1=1;
D2=D2s;
NU=Nu;
N=length(IC);

tend=10;

X=linspace(0,1,N);
dx1=X(2)-X(1);

[~,M]=min(abs(X-l1));

options =odeset('AbsTol',1e-4,'RelTol',1e-4,'Events',@(t,y)TcHit(t,y,Tc,M));
sol = ode23s(@(t,y) heat2(t,y,dx1,D1,D2,M,N,NU),[0 tend], IC,options);

tc=sol.xe;
if isempty(tc)
    tc=tend;
end

figure(10)
surf(X,sol.x,sol.y','Edgecolor','none')
view(0,90)

end


function dy = heat2(~,y,dx1,D1,D2,M,N,NU)
dy=zeros(N,1);
y(1)=(-y(3)+4*y(2)+2*dx1*NU)/(3+2*dx1*NU);
 y(N)=4/3*y(N-1)-1/3*y(N-2);

for i=2:M-1
    dy(i) = D1* 1/dx1^2 * (y(i+1)-2*y(i)+y(i-1));
end

dy(M) = 1/dx1^2*(D2*(-1.5*y(M)+2*y(M+1)-0.5*y(M+2))-D1*(0.5*y(M-2)-2*y(M-1)+1.5*y(M)));
 
 for i=(M+1):(N-1)
     dy(i) = D2/dx1^2 * (y(i+1)-2*y(i)+y(i-1)) ;
 end

end

function [position,isterminal,direction] = TcHit(~,y,Tc,M)
position = y(M)-Tc;% The value that we want to be zero
isterminal = 1;  % Halt integration 
direction = 0;   % The zero can be approached from either direction
end