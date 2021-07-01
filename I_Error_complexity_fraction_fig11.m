%06102021   Error-Complexity
% Calculate Entropy and Forcing complexity level for each case 
% X = C_m = H(Ta,U,VPD_Quantized./ H(Ta,U,VPD)_Model
% omparison between all quantized cases for (N =[2, 3, 4, 5]) based on their forcing data complexity, 
% C_m, and the fraction of 96 time steps of day where quantized model performs better for 
% {LE}, {SH}, and {Fc}.

% use Error_diurnal function


clear
close all
%% calculate Entropy for the whole model
% load quantized input variable for n = 2:5 levels:
load quantdata
% original forcing data
Ta_U_VPD = [Ta U VPD];
N_Model= 15; %set it greater than 5 levels.

[pdf, Coords]= compute_pdf(Ta_U_VPD,N_Model);
info_Model = compute_info_measures(pdf);
H3D_Model = info_Model.H3D; %compute 3D entropy for full model

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
LE_Model_error_diurnal= Model_error_diurnal(:,1);

H_MQ_error_diurnal_1(:,:) = MQ_error_diurnal_1(:,2,:);
H_MQ_error_diurnal_2(:,:) = MQ_error_diurnal_2(:,2,:);
H_MQ_error_diurnal_3(:,:) = MQ_error_diurnal_3(:,2,:);
H_MQ_error_diurnal_4(:,:) = MQ_error_diurnal_4(:,2,:);
H_MQ_error_diurnal_5(:,:) = MQ_error_diurnal_5(:,2,:);
H_MQ_error_diurnal_6(:,:) = MQ_error_diurnal_6(:,2,:);
H_MQ_error_diurnal_7(:,:) = MQ_error_diurnal_7(:,2,:);
H_Model_error_diurnal= Model_error_diurnal(:,2);

Fc_MQ_error_diurnal_1(:,:) = MQ_error_diurnal_1(:,3,:);
Fc_MQ_error_diurnal_2(:,:) = MQ_error_diurnal_2(:,3,:);
Fc_MQ_error_diurnal_3(:,:) = MQ_error_diurnal_3(:,3,:);
Fc_MQ_error_diurnal_4(:,:) = MQ_error_diurnal_4(:,3,:);
Fc_MQ_error_diurnal_5(:,:) = MQ_error_diurnal_5(:,3,:);
Fc_MQ_error_diurnal_6(:,:) = MQ_error_diurnal_6(:,3,:);
Fc_MQ_error_diurnal_7(:,:) = MQ_error_diurnal_7(:,3,:);
Fc_Model_error_diurnal= Model_error_diurnal(:,3);
%% Label each case for each time step
%"worse" == 0 and "better" ==1;
LE_Label_1 = zeros(96,size(X_1,2)) ;
H_Label_1 = zeros(96,size(X_1,2)) ;
Fc_Label_1 = zeros(96,size(X_1,2)) ;
for i = 1:size(X_1,2)
    for j= 1:96
        if LE_MQ_error_diurnal_1(j,i)>LE_Model_error_diurnal(j,1)
            LE_Label_1(j,i)= 0;
        else
            LE_Label_1(j,i) = 1;
        end
        
        if H_MQ_error_diurnal_1(j,i)>H_Model_error_diurnal(j,1)
            H_Label_1(j,i)= 0;
        else
            H_Label_1(j,i) = 1;
        end
        
        if Fc_MQ_error_diurnal_1(j,i)>Fc_Model_error_diurnal(j,1)
            Fc_Label_1(j,i)= 0;
        else
            Fc_Label_1(j,i) = 1;
        end
    end
    LE_frac_better_1(:,i) = sum(LE_Label_1(:,i) == 1)./96;
    H_frac_better_1(:,i) = sum(H_Label_1(:,i) == 1)./96;
    Fc_frac_better_1(:,i) = sum(Fc_Label_1(:,i) == 1)./96;
end

LE_Label_2 = zeros(96,size(X_2,2)) ;
H_Label_2 = zeros(96,size(X_2,2)) ;
Fc_Label_2 = zeros(96,size(X_2,2)) ;
for i = 1:size(X_2,2)
    for j= 1:96
        if LE_MQ_error_diurnal_2(j,i)>LE_Model_error_diurnal(j,1)
            LE_Label_2(j,i)= 0;
        else
            LE_Label_2(j,i) = 1;
        end
        
        if H_MQ_error_diurnal_2(j,i)>H_Model_error_diurnal(j,1)
            H_Label_2(j,i)= 0;
        else
            H_Label_2(j,i) = 1;
        end
        
        if Fc_MQ_error_diurnal_2(j,i)>Fc_Model_error_diurnal(j,1)
            Fc_Label_2(j,i)= 0;
        else
            Fc_Label_2(j,i) = 1;
        end
    end
    LE_frac_better_2(:,i) = sum(LE_Label_2(:,i) == 1)./96;
    H_frac_better_2(:,i) = sum(H_Label_2(:,i) == 1)./96;
    Fc_frac_better_2(:,i) = sum(Fc_Label_2(:,i) == 1)./96;
end

LE_Label_3 = zeros(96,size(X_3,2)) ;
H_Label_3 = zeros(96,size(X_3,2)) ;
Fc_Label_3 = zeros(96,size(X_3,2)) ;
for i = 1:size(X_3,2)
    for j= 1:96
        if LE_MQ_error_diurnal_3(j,i)>LE_Model_error_diurnal(j,1)
            LE_Label_3(j,i)= 0;
        else
            LE_Label_3(j,i) = 1;
        end
        
        if H_MQ_error_diurnal_3(j,i)>H_Model_error_diurnal(j,1)
            H_Label_3(j,i)= 0;
        else
            H_Label_3(j,i) = 1;
        end
        
        if Fc_MQ_error_diurnal_3(j,i)>Fc_Model_error_diurnal(j,1)
            Fc_Label_3(j,i)= 0;
        else
            Fc_Label_3(j,i) = 1;
        end
    end
    LE_frac_better_3(:,i) = sum(LE_Label_3(:,i) == 1)./96;
    H_frac_better_3(:,i) = sum(H_Label_3(:,i) == 1)./96;
    Fc_frac_better_3(:,i) = sum(Fc_Label_3(:,i) == 1)./96;
end

LE_Label_4 = zeros(96,size(X_4,2)) ;
H_Label_4 = zeros(96,size(X_4,2)) ;
Fc_Label_4 = zeros(96,size(X_4,2)) ;
for i = 1:size(X_4,2)
    for j= 1:96
        if LE_MQ_error_diurnal_4(j,i)>LE_Model_error_diurnal(j,1)
            LE_Label_4(j,i)= 0;
        else
            LE_Label_4(j,i) = 1;
        end
        
        if H_MQ_error_diurnal_4(j,i)>H_Model_error_diurnal(j,1)
            H_Label_4(j,i)= 0;
        else
            H_Label_4(j,i) = 1;
        end
        
        if Fc_MQ_error_diurnal_4(j,i)>Fc_Model_error_diurnal(j,1)
            Fc_Label_4(j,i)= 0;
        else
            Fc_Label_4(j,i) = 1;
        end
    end
    LE_frac_better_4(:,i) = sum(LE_Label_4(:,i) == 1)./96;
    H_frac_better_4(:,i) = sum(H_Label_4(:,i) == 1)./96;
    Fc_frac_better_4(:,i) = sum(Fc_Label_4(:,i) == 1)./96;
end

LE_Label_5 = zeros(96,size(X_5,2)) ;
H_Label_5 = zeros(96,size(X_5,2)) ;
Fc_Label_5 = zeros(96,size(X_5,2)) ;
for i = 1:size(X_5,2)
    for j= 1:96
        if LE_MQ_error_diurnal_5(j,i)>LE_Model_error_diurnal(j,1)
            LE_Label_5(j,i)= 0;
        else
            LE_Label_5(j,i) = 1;
        end
        
        if H_MQ_error_diurnal_5(j,i)>H_Model_error_diurnal(j,1)
            H_Label_5(j,i)= 0;
        else
            H_Label_5(j,i) = 1;
        end
        
        if Fc_MQ_error_diurnal_5(j,i)>Fc_Model_error_diurnal(j,1)
            Fc_Label_5(j,i)= 0;
        else
            Fc_Label_5(j,i) = 1;
        end
    end
    LE_frac_better_5(:,i) = sum(LE_Label_5(:,i) == 1)./96;
    H_frac_better_5(:,i) = sum(H_Label_5(:,i) == 1)./96;
    Fc_frac_better_5(:,i) = sum(Fc_Label_5(:,i) == 1)./96;
end

LE_Label_6 = zeros(96,size(X_6,2)) ;
H_Label_6 = zeros(96,size(X_6,2)) ;
Fc_Label_6 = zeros(96,size(X_6,2)) ;
for i = 1:size(X_6,2)
    for j= 1:96
        if LE_MQ_error_diurnal_6(j,i)>LE_Model_error_diurnal(j,1)
            LE_Label_6(j,i)= 0;
        else
            LE_Label_6(j,i) = 1;
        end
        
        if H_MQ_error_diurnal_6(j,i)>H_Model_error_diurnal(j,1)
            H_Label_6(j,i)= 0;
        else
            H_Label_6(j,i) = 1;
        end
        
        if Fc_MQ_error_diurnal_6(j,i)>Fc_Model_error_diurnal(j,1)
            Fc_Label_6(j,i)= 0;
        else
            Fc_Label_6(j,i) = 1;
        end
    end
    LE_frac_better_6(:,i) = sum(LE_Label_6(:,i) == 1)./96;
    H_frac_better_6(:,i) = sum(H_Label_6(:,i) == 1)./96;
    Fc_frac_better_6(:,i) = sum(Fc_Label_6(:,i) == 1)./96;
end

LE_Label_7 = zeros(96,size(X_7,2)) ;
H_Label_7 = zeros(96,size(X_7,2)) ;
Fc_Label_7 = zeros(96,size(X_7,2)) ;
for i = 1:size(X_7,2)
    for j= 1:96
        if LE_MQ_error_diurnal_7(j,i)>LE_Model_error_diurnal(j,1)
            LE_Label_7(j,i)= 0;
        else
            LE_Label_7(j,i) = 1;
        end
        
        if H_MQ_error_diurnal_7(j,i)>H_Model_error_diurnal(j,1)
            H_Label_7(j,i)= 0;
        else
            H_Label_7(j,i) = 1;
        end
        
        if Fc_MQ_error_diurnal_7(j,i)>Fc_Model_error_diurnal(j,1)
            Fc_Label_7(j,i)= 0;
        else
            Fc_Label_7(j,i) = 1;
        end
    end
    LE_frac_better_7(:,i) = sum(LE_Label_7(:,i) == 1)./96;
    H_frac_better_7(:,i) = sum(H_Label_7(:,i) == 1)./96;
    Fc_frac_better_7(:,i) = sum(Fc_Label_7(:,i) == 1)./96;
end
figure (1)
subplot(1,3,1)
plot(X_1,LE_frac_better_1,'b*','LineWidth',1.5,'MarkerSize',6)
hold on 
plot(X_2,LE_frac_better_2,'r*','LineWidth',1.5,'MarkerSize',6)
plot(X_3,LE_frac_better_3,'k*','LineWidth',1.5,'MarkerSize',6)
plot(X_4,LE_frac_better_4,'bo','LineWidth',1.5,'MarkerSize',6)
plot(X_5,LE_frac_better_5,'ro','LineWidth',1.5,'MarkerSize',6)
plot(X_6,LE_frac_better_6,'ko','LineWidth',1.5,'MarkerSize',6)
plot(X_7,LE_frac_better_7,'gs','LineWidth',1.5,'MarkerSize',6)
title('LE', 'FontSize',12)
ylabel('Fraction of time steps of day where quantized model performs better', 'FontSize',12)
ylim([0 1]);

subplot(1,3,2)
plot(X_1,H_frac_better_1,'b*','LineWidth',1.5,'MarkerSize',6)
hold on 
plot(X_2,H_frac_better_2,'r*','LineWidth',1.5,'MarkerSize',6)
plot(X_3,H_frac_better_3,'k*','LineWidth',1.5,'MarkerSize',6)
plot(X_4,H_frac_better_4,'bo','LineWidth',1.5,'MarkerSize',6)
plot(X_5,H_frac_better_5,'ro','LineWidth',1.5,'MarkerSize',6)
plot(X_6,H_frac_better_6,'ko','LineWidth',1.5,'MarkerSize',6)
plot(X_7,H_frac_better_7,'gs','LineWidth',1.5,'MarkerSize',6)
xlabel('Forcing data complexity (C_m)', 'FontSize',12) 
title('SH', 'FontSize',12)
ylim([0 1]);

subplot(1,3,3)
plot(X_1,Fc_frac_better_1,'b*','LineWidth',1.5,'MarkerSize',6)
hold on 
plot(X_2,Fc_frac_better_2,'r*','LineWidth',1.5,'MarkerSize',6)
plot(X_3,Fc_frac_better_3,'k*','LineWidth',1.5,'MarkerSize',6)
plot(X_4,Fc_frac_better_4,'bo','LineWidth',1.5,'MarkerSize',6)
plot(X_5,Fc_frac_better_5,'ro','LineWidth',1.5,'MarkerSize',6)
plot(X_6,Fc_frac_better_6,'ko','LineWidth',1.5,'MarkerSize',6)
plot(X_7,Fc_frac_better_7,'gs','LineWidth',1.5,'MarkerSize',6)
legend('[$\hat{Ta}$]','$[\hat{U}]$','$[\hat{VPD}]$','$[\hat{Ta},\hat{U}]$',...
    '$[\hat{Ta},\hat{VPD}]$','$[\hat{U},\hat{VPD}]$','$[\hat{Ta},\hat{U},\hat{VPD}]$',...
    'full','fontsize',8,'Interpreter','latex')
title('Fc', 'FontSize',12)
ylim([0 1]);
