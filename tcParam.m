function tcParam

%Code used for sensitivity analysis

D1s=0.09; %Skin Layer diffusivity
D2s=.1*D1s; %Fat Layer diffusivity

Dh=D2s/D1s; %Non-Dimensional diffusivity for fat layer

l1=1.6; %Length of skin layer
l2=4-l1; %Length of fat layer
a=l1/(l1+l2); %Ratio of skin layer to total length/Interface Position

NU=3; %Nusselt number [Convective heat transfer coeff=3]

Tc=0.538; %By Moritz' Paper [Irreversible Damage at 44C]

X=linspace(0,1,100)*16/0.09;

for i=1:length(X)
    
%tc(i) = TwoHeat_Tc_V3(0.1,X(i),1,0.2121212121212);
tc(i) = TwoHeat_Tc_V4_new(0.1,0.4,1,0.2121212121212,X(i));
%tc(i) = TwoHeat_Tc_V5_new(1,0.2121212121212,1*10^(-7)*16/0.1,X(i));
%tcdim(i) = (16/D1s)*tc(i);
%Xdim(i) = 13*X(i)+37;

end

plot(X,tc);
xlabel('Nu');
ylabel('Critical Time','Interpreter','latex');
title('Relationship between the Critical Time and the Nusselt Number','Interpreter','latex');
set(gca, 'FontSize', 20)
