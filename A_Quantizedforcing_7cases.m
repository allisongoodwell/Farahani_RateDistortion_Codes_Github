%06012020 n_level Qunatization, Tempreture and wind speed, VPD
% FILE 1
% We define seven model cases based on individual or joint quantization 
% of the input variables (Ta, U, and VPD)
% Save Simplified forcing data for each level of quantization inorder to use in running the MLCan model

% need Quantization_function to quantize forcing variable
clear
close all

load Forcing
Forcing = table2array(Forcing);

%% input variable quantization for n = 2:5 levels of quantization
Ta = Forcing(:,1);
Ta_Q =  Quantization_function(Ta);

VPD = Forcing(:,5);
VPD_Q =  Quantization_function(VPD);

U = Forcing(:,6);
U_Q =  Quantization_function(U);

%% Save Simplified forcing data for each level of quantization
%case 1: Ta only
for n_threshold_Ta = 1:4
    load Maize_Forcing_Goosecreek2018
    fprintf('save threshold_Ta case: %d \n', n_threshold_Ta)
    Ta_crop(14305:24000,:)= Ta_Q(:,n_threshold_Ta);
    save (sprintf('Quantized_Forcing_Maize/Ta_%d_Q_Maize_Goosecreek_Forcing2018.mat',n_threshold_Ta));
end
%case 2: U only
for n_threshold_U = 1:4
    load Maize_Forcing_Goosecreek2018
    fprintf('save threshold_U case: %d \n', n_threshold_U)
    U_crop(14305:24000,:)= U_Q(:,n_threshold_U);
    save (sprintf('Quantized_Forcing_Maize/U_%d_Q_Maize_Goosecreek_Forcing2018.mat',n_threshold_U));
end
%case 3: VPD only
for n_threshold_VPD = 1:4
    load Maize_Forcing_Goosecreek2018
    fprintf('save threshold_VPD case: %d \n', n_threshold_VPD)
    VPD_crop(14305:24000,:)= VPD_Q(:,n_threshold_VPD);
    save (sprintf('Quantized_Forcing_Maize/VPD_%d_Q_Maize_Goosecreek_Forcing2018.mat',n_threshold_VPD));
end
%case 4: Ta + U 
for n_threshold_Ta = 1:4
    load Maize_Forcing_Goosecreek2018
    fprintf('save threshold_Ta_U case: %d \n', n_threshold_Ta)
    Ta_crop(14305:24000,:)= Ta_Q(:,n_threshold_Ta);
    for n_threshold_U = 1:4
        U_crop(14305:24000,:)= U_Q(:,n_threshold_U);
        save (sprintf('Quantized_Forcing_Maize/Ta_%d_U_%d_Q_Maize_Goosecreek_Forcing2018.mat',n_threshold_Ta,n_threshold_U));
    end
end
%case 5: Ta + VPD 
for n_threshold_Ta = 1:4
    load Maize_Forcing_Goosecreek2018
    fprintf('save threshold_Ta_VPD case: %d \n', n_threshold_Ta)
    Ta_crop(14305:24000,:)= Ta_Q(:,n_threshold_Ta);
    for n_threshold_VPD = 1:4
        VPD_crop(14305:24000,:)= VPD_Q(:,n_threshold_VPD);
        save (sprintf('Quantized_Forcing_Maize/Ta_%d_VPD_%d_Q_Maize_Goosecreek_Forcing2018.mat',n_threshold_Ta,n_threshold_VPD));
    end
end
%case 6: U + VPD 
for n_threshold_U = 1:4
    load Maize_Forcing_Goosecreek2018
    fprintf('save threshold_U_VPD case: %d \n', n_threshold_U)
    U_crop(14305:24000,:)= U_Q(:,n_threshold_U);
    for n_threshold_VPD = 1:4
        VPD_crop(14305:24000,:)= VPD_Q(:,n_threshold_VPD);
        save (sprintf('Quantized_Forcing_Maize/U_%d_VPD_%d_Q_Maize_Goosecreek_Forcing2018.mat',n_threshold_U,n_threshold_VPD));
    end
end
%case 7: Ta + VPD + U
for n_threshold_Ta = 1:4
    load Maize_Forcing_Goosecreek2018
    fprintf('save threshold_Ta_VPD_U case: %d \n', n_threshold_Ta)
    Ta_crop(14305:24000,:)= Ta_Q(:,n_threshold_Ta);
    for n_threshold_VPD = 1:4
        VPD_crop(14305:24000,:)= VPD_Q(:,n_threshold_VPD);
        for n_threshold_U = 1:4
            U_crop(14305:24000,:)= U_Q(:,n_threshold_U);
            save (sprintf('Quantized_Forcing_Maize/Ta_%d_VPD_%d_U_%d_Q_Maize_Goosecreek_Forcing2018.mat',n_threshold_Ta,n_threshold_VPD,n_threshold_U));
        end
    end
end