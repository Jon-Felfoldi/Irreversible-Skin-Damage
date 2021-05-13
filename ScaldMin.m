function e = ScaldMin(Dh, Nu, Tc, D1,G)
X=[Dh,Nu,Tc,D1,G];

R = load("ScaldData.mat");

R=R.X;

if min(X)<0 || Tc>min(R(:,2)) || Tc<37
    e=Inf;
else
    [n,~] = size(R);
    tc_dim = zeros(n,1);
    
    l1 = 1.6;
    l2 = 2.4;
    a = l1/(l1+l2);
    
    for i = 1:n
        
        Tc_nd = (Tc-37)/(R(i,2)-37);
        
        %account for G = gamma*L^2/D in non-dimensional model
        
        tc_nd = TwoHeat_Tc_V4_new(Dh,a,Nu,Tc_nd,G*(l1+l2)^2/D1);
        
        tc_nd = max(tc_nd);
        tc_dim(i,1) = (tc_nd*(l1+l2)^2)/D1;
        
    end
    
    figure(100)
    semilogy(R(:,2),R(:,1),'r+',R(:,2),tc_dim,'b-o')
    %   
    e = sum(((R(:,1)-tc_dim(:,1))./R(:,1)).^2  ); %relative square error
    
end
end