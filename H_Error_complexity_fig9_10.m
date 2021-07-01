%07192020   
% calculate Entropy and Forcing complexity  level for each case 
% X = Entropy(Ta,U,VPD)./ Entropy(Ta,U,VPD)Model
% Comparison between the model and all quantized cases for (N = [2, 3, 4, 5]) based on their forcing data complexity, 
% and 3-hour averaged modeled diurnal error, RMSE, at 18:00-21:00 for{LE}

%Comparison between full model and quantized cases (N_{T} =2, N_{U}=2 and N_{VPD}=2 levels) based on
%their forcing data complexity and 3-hour averaged model diurnal RMSE for {LE},{SH}, and {Fc} for three times the day
% at 3:00-6:00 am, at 9:00-12:00 and at 18:00-21:00.

% use Error_diurnal function

clear
close all

%% input variable quantization for n = 2:5 level of quantization
% load Forcing
% Forcing = table2array(Forcing);

% Ta = Forcing(:,1);
% Ta_Q =  Quantization_function(Ta);
 
% VPD = Forcing(:,5);
% VPD_Q =  Quantization_function(VPD);
  
% U = Forcing(:,6);
% U_Q =  Quantization_function(U);

% AEG: this takes a while, going to save these variables for future runs...
% save('quantdata.mat','Ta','Ta_Q','VPD','VPD_Q','U','U_Q')
load quantdata

%% calculate Entropy for the whole model

Ta_U_VPD = [Ta U VPD];
N_Model= 15;

[pdf, Coords]= compute_pdf(Ta_U_VPD,N_Model);
info_Model = compute_info_measures(pdf);
H3D_Model = info_Model.H3D;

%% case 1: Ta only
%case 1: H(Ta_2, U_mod, VPD mod), H(Ta_3, U_mod, VPD_mod)... (since only Ta is changing)
m =1;
for i =1:size(Ta_Q,2)
    case_1(:,:,m) = [Ta_Q(:,i) U VPD];
    m = m+1;
end

for i =1:size(case_1,3)
    [pdf]= compute_pdf(case_1(:,:,i),N_Model);
    info_case1 = compute_info_measures(pdf);
    H3D_case1(:,i) = info_case1.H3D;
    X_1(:,i) = H3D_case1(:,i)./H3D_Model; % Forcing complexity  level:X
end

%% case 2: U only
m =1;
for i =1:size(U_Q,2)
    case_2(:,:,m) = [Ta U_Q(:,i) VPD];
    m = m+1;
end

for i =1:size(case_2,3)
    [pdf]= compute_pdf(case_2(:,:,i),N_Model);
    info_case2 = compute_info_measures(pdf);
    H3D_case2(:,i) = info_case2.H3D;
    X_2(:,i) = H3D_case2(:,i)./H3D_Model;
end

%% case 3: VPD only
m =1;
for i =1:size(VPD_Q,2)
    case_3(:,:,m) = [Ta U VPD_Q(:,i)];
    m = m+1;
end

for i =1:size(case_3,3)
    [pdf]= compute_pdf(case_3(:,:,i),N_Model);
    info_case3 = compute_info_measures(pdf);
    H3D_case3(:,i) = info_case3.H3D;
    X_3(:,i) = H3D_case3(:,i)./H3D_Model;
end

%% case 4: Ta + U
m =1;
for i =1:size(Ta_Q,2)
    for j =1:size(U_Q,2)
        case_4(:,:,m) = [Ta_Q(:,i) U_Q(:,j) VPD];
        m = m+1;        
    end
end

for i =1:size(case_4,3)
    [pdf]= compute_pdf(case_4(:,:,i),N_Model);
    info_case4 = compute_info_measures(pdf);
    H3D_case4(:,i) = info_case4.H3D;
    X_4(:,i) = H3D_case4(:,i)./H3D_Model;
end

%% case 5: Ta + VPD 
m =1;
for i =1:size(Ta_Q,2)
    for j =1:size(VPD_Q,2)
        case_5(:,:,m) = [Ta_Q(:,i) U VPD_Q(:,j)];
        m = m+1;        
    end
end

for i =1:size(case_5,3)
    [pdf]= compute_pdf(case_5(:,:,i),N_Model);
    info_case5 = compute_info_measures(pdf);
    H3D_case5(:,i) = info_case5.H3D;
    X_5(:,i) = H3D_case5(:,i)./H3D_Model;
end

%% case 6: U + VPD 
m =1;
for i =1:size(U_Q,2)
    for j =1:size(VPD_Q,2)
        case_6(:,:,m) = [Ta U_Q(:,i) VPD_Q(:,j)];
        m = m+1;
    end
end
for i =1:size(case_6,3)
    [pdf]= compute_pdf(case_6(:,:,i),N_Model);
    info_case6 = compute_info_measures(pdf);
    H3D_case6(:,i) = info_case6.H3D;
    X_6(:,i) = H3D_case6(:,i)./H3D_Model;
end

%% case 7: Ta + U + VPD
m =1;
for i =1:size(Ta_Q,2)
    for j =1:size(U_Q,2)
        for k =1:size(VPD_Q,2)
            case_7(:,:,m) = [Ta_Q(:,i) U_Q(:,j) VPD_Q(:,k)];
            m = m+1;
        end
    end
end

for i =1:size(case_7,3)
    [pdf]= compute_pdf(case_7(:,:,i),N_Model);
    info_case7 = compute_info_measures(pdf);
    H3D_case7(:,i) = info_case7.H3D;
    X_7(:,i) = H3D_case7(:,i)./H3D_Model;
end

%% Model and  Quantized diurnal error
[MQ_error_diurnal_1,Model_error_diurnal] = Error_diurnal('Results_addedfraction/Ta only');
MQ_error_diurnal_2  = Error_diurnal('Results_addedfraction/U only');
MQ_error_diurnal_3  = Error_diurnal('Results_addedfraction/VPD only');
MQ_error_diurnal_4  = Error_diurnal('Results_addedfraction/Ta_U');
MQ_error_diurnal_5  = Error_diurnal('Results_addedfraction/Ta_VPD');
MQ_error_diurnal_6 = Error_diurnal('Results_addedfraction/VPD_U');
MQ_error_diurnal_7  = Error_diurnal('Results_addedfraction/Ta_U_VPD');


LE_MQ_error_diurnal_1(:,:) = MQ_error_diurnal_1(:,1,:);
LE_MQ_error_diurnal_2(:,:) = MQ_error_diurnal_2(:,1,:);
LE_MQ_error_diurnal_3(:,:) = MQ_error_diurnal_3(:,1,:);
LE_MQ_error_diurnal_4(:,:) = MQ_error_diurnal_4(:,1,:);
LE_MQ_error_diurnal_5(:,:) = MQ_error_diurnal_5(:,1,:);
LE_MQ_error_diurnal_6(:,:) = MQ_error_diurnal_6(:,1,:);
LE_MQ_error_diurnal_7(:,:) = MQ_error_diurnal_7(:,1,:);

H_MQ_error_diurnal_1(:,:) = MQ_error_diurnal_1(:,2,:);
H_MQ_error_diurnal_2(:,:) = MQ_error_diurnal_2(:,2,:);
H_MQ_error_diurnal_3(:,:) = MQ_error_diurnal_3(:,2,:);
H_MQ_error_diurnal_4(:,:) = MQ_error_diurnal_4(:,2,:);
H_MQ_error_diurnal_5(:,:) = MQ_error_diurnal_5(:,2,:);
H_MQ_error_diurnal_6(:,:) = MQ_error_diurnal_6(:,2,:);
H_MQ_error_diurnal_7(:,:) = MQ_error_diurnal_7(:,2,:);

Fc_MQ_error_diurnal_1(:,:) = MQ_error_diurnal_1(:,3,:);
Fc_MQ_error_diurnal_2(:,:) = MQ_error_diurnal_2(:,3,:);
Fc_MQ_error_diurnal_3(:,:) = MQ_error_diurnal_3(:,3,:);
Fc_MQ_error_diurnal_4(:,:) = MQ_error_diurnal_4(:,3,:);
Fc_MQ_error_diurnal_5(:,:) = MQ_error_diurnal_5(:,3,:);
Fc_MQ_error_diurnal_6(:,:) = MQ_error_diurnal_6(:,3,:);
Fc_MQ_error_diurnal_7(:,:) = MQ_error_diurnal_7(:,3,:);
%% 3hour averaged
LE_Model_error_3hr= mean(reshape(Model_error_diurnal(:,1), 12, []));
H_Model_error_3hr= mean(reshape(Model_error_diurnal(:,2), 12, []));
Fc_Model_error_3hr= mean(reshape(Model_error_diurnal(:,3), 12, []));

for i=1:4
LE_MQ_error_3hr_1(:,i)= mean(reshape(LE_MQ_error_diurnal_1(:,i), 12, []));
LE_MQ_error_3hr_2(:,i)= mean(reshape(LE_MQ_error_diurnal_2(:,i), 12, []));
LE_MQ_error_3hr_3(:,i)= mean(reshape(LE_MQ_error_diurnal_3(:,i), 12, []));

H_MQ_error_3hr_1(:,i)= mean(reshape(H_MQ_error_diurnal_1(:,i), 12, []));
H_MQ_error_3hr_2(:,i)= mean(reshape(H_MQ_error_diurnal_2(:,i), 12, []));
H_MQ_error_3hr_3(:,i)= mean(reshape(H_MQ_error_diurnal_3(:,i), 12, []));

Fc_MQ_error_3hr_1(:,i)= mean(reshape(Fc_MQ_error_diurnal_1(:,i), 12, []));
Fc_MQ_error_3hr_2(:,i)= mean(reshape(Fc_MQ_error_diurnal_2(:,i), 12, []));
Fc_MQ_error_3hr_3(:,i)= mean(reshape(Fc_MQ_error_diurnal_3(:,i), 12, []));
end
for i=1:16
LE_MQ_error_3hr_4(:,i)= mean(reshape(LE_MQ_error_diurnal_4(:,i), 12, []));
LE_MQ_error_3hr_5(:,i)= mean(reshape(LE_MQ_error_diurnal_5(:,i), 12, []));
LE_MQ_error_3hr_6(:,i)= mean(reshape(LE_MQ_error_diurnal_6(:,i), 12, []));

H_MQ_error_3hr_4(:,i)= mean(reshape(H_MQ_error_diurnal_4(:,i), 12, []));
H_MQ_error_3hr_5(:,i)= mean(reshape(H_MQ_error_diurnal_5(:,i), 12, []));
H_MQ_error_3hr_6(:,i)= mean(reshape(H_MQ_error_diurnal_6(:,i), 12, []));

Fc_MQ_error_3hr_4(:,i)= mean(reshape(Fc_MQ_error_diurnal_4(:,i), 12, []));
Fc_MQ_error_3hr_5(:,i)= mean(reshape(Fc_MQ_error_diurnal_5(:,i), 12, []));
Fc_MQ_error_3hr_6(:,i)= mean(reshape(Fc_MQ_error_diurnal_6(:,i), 12, []));
end
for i=1:64
LE_MQ_error_3hr_7(:,i)= mean(reshape(LE_MQ_error_diurnal_7(:,i), 12, []));
H_MQ_error_3hr_7(:,i)= mean(reshape(H_MQ_error_diurnal_7(:,i), 12, []));
Fc_MQ_error_3hr_7(:,i)= mean(reshape(Fc_MQ_error_diurnal_7(:,i), 12, []));
end

%% some plotting
%0-3, 3-6(2),6-9,9-12(4),12-3,3-6,6-9(7),9-12
figure (2)
subplot(3,3,1)
plot(X_1(1,1),LE_MQ_error_3hr_1(2,1),'b*','LineWidth',2,'MarkerSize',8)
hold on
plot(X_2(1,1),LE_MQ_error_3hr_2(2,1),'r*','LineWidth',2,'MarkerSize',8)
plot(X_3(1,1),LE_MQ_error_3hr_3(2,1),'k*','LineWidth',2,'MarkerSize',8)
plot(X_4(1,1),LE_MQ_error_3hr_4(2,1),'bo','LineWidth',2,'MarkerSize',8)
plot(X_5(1,1),LE_MQ_error_3hr_5(2,1),'ro','LineWidth',2,'MarkerSize',8)
plot(X_6(1,1),LE_MQ_error_3hr_6(2,1),'ko','LineWidth',2,'MarkerSize',8)
plot(X_7(1,1),LE_MQ_error_3hr_7(2,1),'gs','LineWidth',2,'MarkerSize',8)
plot(1,LE_Model_error_3hr(:,2),'p','MarkerSize',15,'MarkerFaceColor','y')
ylabel('LE (W/m^2)', 'FontSize',12)
% ylim([0 2])
title('3-hr averaged RMSE at 3:00-6:00', 'FontSize',10)

subplot(3,3,4)
plot(X_1(1,1),H_MQ_error_3hr_1(2,1),'b*','LineWidth',2,'MarkerSize',8)
hold on
plot(X_2(1,1),H_MQ_error_3hr_2(2,1),'r*','LineWidth',2,'MarkerSize',8)
plot(X_3(1,1),H_MQ_error_3hr_3(2,1),'k*','LineWidth',2,'MarkerSize',8)
plot(X_4(1,1),H_MQ_error_3hr_4(2,1),'bo','LineWidth',2,'MarkerSize',8)
plot(X_5(1,1),H_MQ_error_3hr_5(2,1),'ro','LineWidth',2,'MarkerSize',8)
plot(X_6(1,1),H_MQ_error_3hr_6(2,1),'ko','LineWidth',2,'MarkerSize',8)
plot(X_7(1,1),H_MQ_error_3hr_7(2,1),'gs','LineWidth',2,'MarkerSize',8)
plot(1,H_Model_error_3hr(:,2),'p','MarkerSize',15,'MarkerFaceColor','y')
ylabel('SH (W/m^2)', 'FontSize',12)
% ylim([63 158])

subplot(3,3,7)
plot(X_1(1,1),Fc_MQ_error_3hr_1(2,1),'b*','LineWidth',2,'MarkerSize',8)
hold on
plot(X_2(1,1),Fc_MQ_error_3hr_2(2,1),'r*','LineWidth',2,'MarkerSize',8)
plot(X_3(1,1),Fc_MQ_error_3hr_3(2,1),'k*','LineWidth',2,'MarkerSize',8)
plot(X_4(1,1),Fc_MQ_error_3hr_4(2,1),'bo','LineWidth',2,'MarkerSize',8)
plot(X_5(1,1),Fc_MQ_error_3hr_5(2,1),'ro','LineWidth',2,'MarkerSize',8)
plot(X_6(1,1),Fc_MQ_error_3hr_6(2,1),'ko','LineWidth',2,'MarkerSize',8)
plot(X_7(1,1),Fc_MQ_error_3hr_7(2,1),'gs','LineWidth',2,'MarkerSize',8)
plot(1,Fc_Model_error_3hr(:,2),'p','MarkerSize',15,'MarkerFaceColor','y')
xlabel('Forcing data complexity (C_m)', 'FontSize',10) 
ylabel('Fc (\mumol/m^2s)', 'FontSize',12)
% ylim([9 11])

subplot(3,3,2)
plot(X_1(1,1),LE_MQ_error_3hr_1(4,1),'b*','LineWidth',2,'MarkerSize',8)
hold on
plot(X_2(1,1),LE_MQ_error_3hr_2(4,1),'r*','LineWidth',2,'MarkerSize',8)
plot(X_3(1,1),LE_MQ_error_3hr_3(4,1),'k*','LineWidth',2,'MarkerSize',8)
plot(X_4(1,1),LE_MQ_error_3hr_4(4,1),'bo','LineWidth',2,'MarkerSize',8)
plot(X_5(1,1),LE_MQ_error_3hr_5(4,1),'ro','LineWidth',2,'MarkerSize',8)
plot(X_6(1,1),LE_MQ_error_3hr_6(4,1),'ko','LineWidth',2,'MarkerSize',8)
plot(X_7(1,1),LE_MQ_error_3hr_7(4,1),'gs','LineWidth',2,'MarkerSize',8)
plot(1,LE_Model_error_3hr(:,4),'p','MarkerSize',15,'MarkerFaceColor','y')
% ylim([122 146])
title('3-hr averaged RMSE at 9:00-12:00', 'FontSize',10)

subplot(3,3,5)
plot(X_1(1,1),H_MQ_error_3hr_1(4,1),'b*','LineWidth',2,'MarkerSize',8)
hold on
plot(X_2(1,1),H_MQ_error_3hr_2(4,1),'r*','LineWidth',2,'MarkerSize',8)
plot(X_3(1,1),H_MQ_error_3hr_3(4,1),'k*','LineWidth',2,'MarkerSize',8)
plot(X_4(1,1),H_MQ_error_3hr_4(4,1),'bo','LineWidth',2,'MarkerSize',8)
plot(X_5(1,1),H_MQ_error_3hr_5(4,1),'ro','LineWidth',2,'MarkerSize',8)
plot(X_6(1,1),H_MQ_error_3hr_6(4,1),'ko','LineWidth',2,'MarkerSize',8)
plot(X_7(1,1),H_MQ_error_3hr_7(4,1),'gs','LineWidth',2,'MarkerSize',8)
plot(1,H_Model_error_3hr(:,4),'p','MarkerSize',15,'MarkerFaceColor','y')
% ylim([63 158])

subplot(3,3,8)
plot(X_1(1,1),Fc_MQ_error_3hr_1(4,1),'b*','LineWidth',2,'MarkerSize',8)
hold on
plot(X_2(1,1),Fc_MQ_error_3hr_2(4,1),'r*','LineWidth',2,'MarkerSize',8)
plot(X_3(1,1),Fc_MQ_error_3hr_3(4,1),'k*','LineWidth',2,'MarkerSize',8)
plot(X_4(1,1),Fc_MQ_error_3hr_4(4,1),'bo','LineWidth',2,'MarkerSize',8)
plot(X_5(1,1),Fc_MQ_error_3hr_5(4,1),'ro','LineWidth',2,'MarkerSize',8)
plot(X_6(1,1),Fc_MQ_error_3hr_6(4,1),'ko','LineWidth',2,'MarkerSize',8)
plot(X_7(1,1),Fc_MQ_error_3hr_7(4,1),'gs','LineWidth',2,'MarkerSize',8)
plot(1,Fc_Model_error_3hr(:,4),'p','MarkerSize',15,'MarkerFaceColor','y')
xlabel('Forcing data complexity (C_m)', 'FontSize',10) 
% ylim([9 11])

subplot(3,3,3)
plot(X_1(1,1),LE_MQ_error_3hr_1(7,1),'b*','LineWidth',2,'MarkerSize',8)
hold on
plot(X_2(1,1),LE_MQ_error_3hr_2(7,1),'r*','LineWidth',2,'MarkerSize',8)
plot(X_3(1,1),LE_MQ_error_3hr_3(7,1),'k*','LineWidth',2,'MarkerSize',8)
plot(X_4(1,1),LE_MQ_error_3hr_4(7,1),'bo','LineWidth',2,'MarkerSize',8)
plot(X_5(1,1),LE_MQ_error_3hr_5(7,1),'ro','LineWidth',2,'MarkerSize',8)
plot(X_6(1,1),LE_MQ_error_3hr_6(7,1),'ko','LineWidth',2,'MarkerSize',8)
plot(X_7(1,1),LE_MQ_error_3hr_7(7,1),'gs','LineWidth',2,'MarkerSize',8)
plot(1,LE_Model_error_3hr(:,7),'p','MarkerSize',15,'MarkerFaceColor','y')
legend('[$\hat{Ta}_2$]','$[\hat{U}_2]$','$[\hat{VPD}_2]$','$[\hat{Ta}_2,\hat{U}_2]$',...
    '$[\hat{Ta}_2,\hat{VPD}_2]$','$[\hat{U}_2,\hat{VPD}_2]$','$[\hat{Ta}_2,\hat{U}_2,\hat{VPD}_2]$',...
    'full','fontsize',10,'Interpreter','latex','orientation','horizontal')
% ylim([122 146])
title('3-hr averaged RMSE at 18:00-21:00', 'FontSize',10)

subplot(3,3,6)
plot(X_1(1,1),H_MQ_error_3hr_1(7,1),'b*','LineWidth',2,'MarkerSize',8)
hold on
plot(X_2(1,1),H_MQ_error_3hr_2(7,1),'r*','LineWidth',2,'MarkerSize',8)
plot(X_3(1,1),H_MQ_error_3hr_3(7,1),'k*','LineWidth',2,'MarkerSize',8)
plot(X_4(1,1),H_MQ_error_3hr_4(7,1),'bo','LineWidth',2,'MarkerSize',8)
plot(X_5(1,1),H_MQ_error_3hr_5(7,1),'ro','LineWidth',2,'MarkerSize',8)
plot(X_6(1,1),H_MQ_error_3hr_6(7,1),'ko','LineWidth',2,'MarkerSize',8)
plot(X_7(1,1),H_MQ_error_3hr_7(7,1),'gs','LineWidth',2,'MarkerSize',8)
plot(1,H_Model_error_3hr(:,7),'p','MarkerSize',15,'MarkerFaceColor','y')
% ylim([63 158])

subplot(3,3,9)
plot(X_1(1,1),Fc_MQ_error_3hr_1(7,1),'b*','LineWidth',2,'MarkerSize',8)
hold on
plot(X_2(1,1),Fc_MQ_error_3hr_2(7,1),'r*','LineWidth',2,'MarkerSize',8)
plot(X_3(1,1),Fc_MQ_error_3hr_3(7,1),'k*','LineWidth',2,'MarkerSize',8)
plot(X_4(1,1),Fc_MQ_error_3hr_4(7,1),'bo','LineWidth',2,'MarkerSize',8)
plot(X_5(1,1),Fc_MQ_error_3hr_5(7,1),'ro','LineWidth',2,'MarkerSize',8)
plot(X_6(1,1),Fc_MQ_error_3hr_6(7,1),'ko','LineWidth',2,'MarkerSize',8)
plot(X_7(1,1),Fc_MQ_error_3hr_7(7,1),'gs','LineWidth',2,'MarkerSize',8)
plot(1,Fc_Model_error_3hr(:,7),'p','MarkerSize',15,'MarkerFaceColor','y')
xlabel('Forcing data complexity (C_m)', 'FontSize',10) 
%%
figure (1)
plot(X_1,LE_MQ_error_3hr_1(7,:),'b*','MarkerSize',8,'LineWidth',1.5)
hold on
set(gca,'FontSize',15)        
plot(X_2,LE_MQ_error_3hr_2(7,:),'r*','MarkerSize',8,'LineWidth',1.5)
plot(X_3,LE_MQ_error_3hr_3(7,:),'k*','MarkerSize',8,'LineWidth',1.5)
plot(X_4,LE_MQ_error_3hr_4(7,:),'bo','MarkerSize',8,'LineWidth',1.5)
plot(X_5,LE_MQ_error_3hr_5(7,:),'ro','MarkerSize',8,'LineWidth',1.5)
plot(X_6,LE_MQ_error_3hr_6(7,:),'ko','MarkerSize',8,'LineWidth',1.5)
plot(X_7,LE_MQ_error_3hr_7(7,:),'gs','MarkerSize',8,'LineWidth',1.5)
plot(1,LE_Model_error_3hr(:,7),'p','MarkerSize',25,'MarkerFaceColor','y','LineWidth',2)
legend('[$\hat{Ta}$]','$[\hat{U}]$','$[\hat{VPD}]$','$[\hat{Ta},\hat{U}]$',...
    '$[\hat{Ta},\hat{VPD}]$','$[\hat{U},\hat{VPD}]$','$[\hat{Ta},\hat{U},\hat{VPD}]$',...
    'full','fontsize',12,'Interpreter','latex','orientation','horizontal')
xlabel('Forcing data complexity (C_m)', 'FontSize',15) 
ylabel('3 hr averaged RMSE for LE (W/m^2) at 18:00-21:00', 'FontSize',15)



