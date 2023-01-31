Mechanotransduction model

The process of adhesion formation, maturation and disassembly are modelled with three major proteins (integrins, talin and vinculin) using differential equations. 

Description of the files:

ModelPipeline3.m - the main model file. Dependencies: nicePlot.m and DE_definitions11.m

SensPipeline3.m - Code to run sensitivity analysis on the model. Dependencies:  SensitivityAnalysisV3.m and testV3.m. 

PlottingPipeline.m - Code to plot the sensitivity analysis results. Dependencies: nicePlot.m

DE_definitions11.m - File with all the differential equations. 

calib_catchbond_opt.m - Code to optimize parameters for force dependent catch and slip bond rupture rates from calibrationData.

sseval_catchbond.m - function evaluating sum of square errors 

sseval_slip.m - function evaluating sum of square errors 

SensitivityAnalysisV3 - initializes variables that change on each loop of the sensitivity analysis

testV3 - runs the euler integration loop 

nicePlot.m - makes nice plots from not-so-nice defualt MATLAB plots. 

--

For running simulations with the baseline parameters, ModelPipeline3 can be used. Change the path to your choice of save directory based on where the figures need to be saved. 

Baseline figures - Figs 3B, 4, S2, S3, S4, S6 A-D, S7 A, S11. 

Fig 3A - Uncomment lines 135-172 of DE_definitions11.m - essentially setting those reaction rates to 0

Fig S5 - Comment out line 450 of ModelPipeline3.m 

Fig S6 E-H - Set TDRM_bool = 0 in line 268 of ModelPipeline3.m

Fig S7 B - Set TDRM_bool = 0 in line 268 of ModelPipeline3.m

Fig S8 - set talin_refold_factor = 0.2 on line 278 of ModelPipeline3.m

Fig S12 - set order1 = 1, and order2 = 1 on lines 179 and 180 of DE_definitions11.m

Fig S13 - set cat=1 in line 271 of ModelPipeline3.m 

--
To run the catch/slip bond parameter calibration, calib_catchbond_opt.m should be used. It requires the path to the folder with calibration data - change this path in the code. 
