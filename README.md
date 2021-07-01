# Farahani_RateDistortion_Codes_Github
 repository for Farahani et al, submitted July 2021

Matlab codes for “Evaluating ecohydrological model sensitivity to input variability with an information theory-based approach” by: M. A. Farahani, A. Vahid, A.E. Goodwell, submitted July 2021 to Environmental Modeling and Software

List of codes:
A_Quantizedforcing_7cases: We define seven model cases based on individual or joint quantization of the input variables (Ta, U, and VPD). This code simplified forcing data for each level of quantization in order to use in running the MLCan model. 
B_NLQ_fig2: quantize air temperature to N=2, 3, 4, and 5 levels in order to find the representation points and final decision threshold correspond to each level of quantization. This also computes representation points calculated with the fixed binning method (Figure 2).
C_simplified_Ta_fig3: This code illustrates original Ta forcing data and simplified Ta to N=2, 3, 4, and 5 levels of quantization for growing season (Figure 3).
D_diffentropy_fig4: Compute the differences between the entropy of the quantized forcing data of Ta, U and VPD into N bins using fixed binning and Lloyd algorithm and compare them with maximum possible entropy Hmax for each level of quantization (Figure 4).
E_Error_diurnal_plot_fig6: This code computes diurnal quantized model, diurnal full model and diurnal observation results for LE, SH, and Fc in order to compare them. It also computes the diurnal cycle error (RMSE) of full model and quantized model and the difference between quantized model (Case 1) and full model RMSE (Figure6).
F_Error_diurnal_plot_fig7:  This code computes the quantized versus full model performance for all cases at different time of day, defined as the quantized diurnal error divided by full model diurnal error (Figure 7).
G_Error_20w_plot_fig8: This code compute the difference between quantized and full model error (DeltaRMSE,t) of latent heat within five 20-day time windows through the growing season (Figure 8).
H_Error_complexity_fig9_10: First, this code calculates Forcing data complexity (Cm) for each case in order to do comparison between the model and all quantized cases for (N = 2, 3, 4, 5) based on their forcing data complexity, and 3-hour averaged modeled diurnal error, RMSE, at 18:00-21:00 for LE (Figure 9). Second, compares between full model and quantized cases for N=2 levels of quantization based on their forcing data complexity and 3-hour averaged model diurnal RMSE for LE, SH, and Fc for three times the day at 3:00-6:00 am, at 9:00-12:00 and at 18:00-21:00 (Figure 10).
I_Error_complexity_fraction_fig11: First, this code calculates Forcing data complexity (Cm) for each case in order to do comparison between the model and all quantized cases for (N = 2, 3, 4, 5) based on their forcing data complexity and the fraction of 96-time steps of day where quantized model performs better for LE, SH and Fc.


Functions:
compute_info_measures: used in “H_Error_complexity_fig9_10”, “I_Error_complexity_fraction_fig11”
compute_pdf: used in “B_NLQ_fig2”, “H_Error_complexity_fig9_10”, “I_Error_complexity_fraction_fig11”
cornsoy_fractions
diffentropy_function: used in “D_diffentropy_fig4” in order to compute the differences between the entropy of the quantized forcing data of Ta, U and VPD into N bins using fixed binning and Lloyd algorithm.
Error_diurnal: used in “H_Error_complexity_fig9_10”, “I_Error_complexity_fraction_fig11”
Quantization_function: used in “A_Quantizedforcing_7cases” to quantize forcing variable.
Model result folder:
“Results_addedfraction” folder: Contain all corn and soy MLCan  models results combined after 
quantization (quantized mode results) for all 124 subcases.

Others (data files):
Directmeasurment.mat: flux tower observation measurements
cornsoy_fractions.mat: fractions of corn and soybean fields according to the flux footprint model
Forcing.mat: input forcing variable for growing season
Maize_Forcing_Goosecreek2018.mat: input forcing variable to run the MLCan
Maize_result_Model.mat: the full MLCan result for corn which forced original forcing data
soybean_3.9.2020.mat: the full MLCan result for soybean which forced original forcing data
Quantdata.mat: quantizaed input forcing variable for growing season
Ta_crop: air temperature measurement for the whole year
Ta.mat: air temperature measurement for growing season
VPD.mat: vapor pressure dificit measurement for growing season
U.mat: wind speed measurement for growing season



