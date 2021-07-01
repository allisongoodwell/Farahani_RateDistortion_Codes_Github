% 04292021 
% Differences between the entropy of the quantized forcing data 
% of Ta, U and VPD into N bins using fixed binning and Lloyd algorithm.
% (diffentropy_function.m)
% compare them with maximum possible entropy H_{max} for each level of quantization
% Figure 4

clear
close all

load Forcing
Forcing = table2array(Forcing);

%% input variable 
Ta = Forcing(:,1);
VPD = Forcing(:,5);
U = Forcing(:,6);

[Hx_Ta,  Hx_Ta_Q] = diffentropy_function(Ta) ;
[Hx_VPD,  Hx_VPD_Q] = diffentropy_function(VPD) ;
[Hx_U,  Hx_U_Q] = diffentropy_function(U) ;

for n=1:4
    H_max(:,n) = log2(n+1);  %maximum possible entropy H_{max} for each level of quantization
end

%%
Ta_U_VPD = [Ta U VPD];
for j=1:4
[pdf, Coords]= compute_pdf(Ta_U_VPD,j+1);
info_Model = compute_info_measures(pdf);
H3D_Model(:,j) = info_Model.H3D;
end
%%

X_L = categorical({'N=2' 'N=3' 'N=4' 'N=5' });
X_L = reordercats(X_L,{'N=2' 'N=3' 'N=4' 'N=5'});

X=[2:1:5];
figure (1)
scatter(X-0.2,Hx_Ta_Q, 'r','filled');
grid on
hold on
scatter(X,Hx_U_Q, 'b','filled')
scatter(X+0.2,Hx_VPD_Q,'g','filled')
scatter(X, H_max, 'k','filled');
scatter(X-0.2,Hx_Ta,'r','filled');
scatter(X+0.2,Hx_VPD, 'g','filled')
scatter(X,Hx_U, 'b','filled')
set(gca,'XTick',[2 3 4 5 ],'FontSize',10);
ylim([Hx_U(1,1) 2.4])
h= set(gca,'FontSize',10);
legendstr={'Ta','U','VPD','H_{max}','',''};
legend(legendstr(1:4),'Location','northwest')
ylabel('Entropy (bits)', 'FontSize',13)
xlabel('Level of quantization', 'FontSize',13)
hold off
