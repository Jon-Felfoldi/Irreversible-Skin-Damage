function [e,R,tc_dim] = ScaldMin2(d1, Nu, Tc, G, Z)
X=[d1,Nu,Tc,G,Z];

R = load("ScaldData.mat");

R=R.X;

if min(X)<0 || Tc>min(R(:,2)) || Tc<37
    e=Inf;
else
    [n,~] = size(R);
    tc_dim = zeros(n,1);
    
    l1 = 1.6;
    
    for i = 1:n
        
        Tc_nd = (Tc-37)/(R(i,2)-37);
        
        %account for G = gamma*L^2/D in non-dimensional model
        
        tc_nd = TwoHeat_Tc_V5_new(Nu,Tc_nd,G*(l1)^2/d1,Z);
        
        tc_nd = max(tc_nd);
        tc_dim(i,1) = (tc_nd*(l1)^2)/d1;
        
    end
    
%     R1 = R;
%     R2 = R;
%     R1(:,1) = R(:,1) - 1800;, R1(:,2),R1(:,1),'g+',R2(:,2),R2(:,1),'b+'
%     R2(:,1) = R(:,1) + 1800;
    
    figure(90)
    semilogy(R(:,2),R(:,1),'r+',R(:,2),tc_dim,'b-o')

    e = sum( ((R(:,1)-tc_dim(:,1))./R(:,1)).^2  ); %relative square error
    
end
end