# Master thesis - Yosef Emshani

## General

### "2_04_P_real_cropped_single"

"/DICOM" has Gold_Atlas 2_04_P real CT data.

"/Data" has results of MC simulations and developed proton transfer algorithm.

"/Python":
* "CT_Faker.py": Fakes CT slices and HU values.
* "everything.py": Transforms DICOM to SPR values using HLUT one can define using function.
* alternative: "SPR.py" imports CT and outputs SPR values using HLUT.
* "SliceFaker.py" manipulates Slice Thickness, Pixel data etc. of CT slices.
* "SlicesToBone.py" and "SlicesToWater.py" manipulates HU values of CT slice. You can import any HU value you desire.
* "calculator.py" calculates r = (m_p * v)/(q_e * B) for vacuum analysis.
* "combinedBeam.py" imports CT slice, trajectory calculated with proton transfer algorithm and MC dose distribution and outputs comparison plot.
* "topasBeam.py" only plots out MC trajectory through dose distribution heatmap.
* "voiTarget.py" analyzes matRad CTV, matRad stopping positions and imported algorithm stopping positions.
* "yshift.py" calculates Gaussian fit, R80 of MC simulation and trajectory lengths for both B = 0 T and B > 0 T with analyzed relative dose needed for both of them to be equal.

"/TOPAS":
* Use slice.txt for used beam settings and MC simulations of various cases.

### "matRad-master"
* "cihangir.m" analyzes dose influence matrix after TP creation and outputs R80 x-y positions in mm.
* "ProtonSimulation_half.m" is the code of calculating the proton transfer with magnetic fields.
* "run_ProtonSimulation_half.m" executes prior mentioned code with input of initial position, energy, SPR map, magnetic field, etc.
* "VacuumSimulation.m" calculates vacuum case proton transfer.
* "run_VacuumSimulation.m" executes prior mentioned code for analyzing vacuum case.
* "matRad_generateStf.m" is matRad initial given code. Analyze "voiTarget" here for CTV.
* "matRad_CT.m" is source code of generating a matRad TP, here using "2_04_P_CT_CST.mat".
* "run_LoopFinal.m" compares initial proton stopping positions calculated using the developed proton transfer algorithm and with matRad TP.
* "run_Gradient.m" uses gradient descent and outputs final positions.

## Tutorial

* 1. Prepare CT dataset and analyze its' HU values with conversion into an SPR map using SPR.py or everything.py.
* 2. Import SPR.csv into "run_ProtonSimulation.m" and give initial energy, etc. You receive "trajectory.csv"
* 3. You could transform those results into mm, since initial result should be printed out in cm.
* 4. Now you have your analytical result of proton stopping positions.
* 5. For MC results run "/TOPAS/slice.txt" and then "yshift.py" for example for R80 analysis using TOPAS_Beam.csv (output of slice.txt) dose distribution.
* 6. You can compare both using "combinedBeam.py".
* 7. Use "yshift.py" for B = 0 T and B > 0 T analysis of stopping positions of MC. Thus, you can now compare algorithm & MC.
* 8. Now for matRad TP use "matRad_CT.m".
* 9. Analyze comparison using "run_LoopFinal.m".
* 10. Run gradient descent algorithm using "run_Gradient.m".
* 11. Visualize results of matRad & algorithm comparison using "voiTarget.py".