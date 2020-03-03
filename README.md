# Human_all_active_models_EAP

1. Download morphology from allsdk
   To run the code, we need to install allensdk(https://alleninstitute.github.io/AllenSDK/install.html).
   In the "codes" folder: Yina_download_morph_step1.ipynb

2. Calculate rotation angels for simulations
   In the "codes" folder: run Yina_cal_rotation_angle_step2.ipynb to calculate the rotation angle for the cell to make sure that the apical dendrites ascend toward the pia in the simulation. 

3. Run simulations
   To run the simulations, we have two steps to do:
   1) Install bmtk and create an enviroment(for example, we named it "bmtk_ateam"). 
      Please follow the instructions in this link (https://github.com/AllenInstitute/bmtk).
   2) Compile the modfiles
      After installing bmtk, run the command: source activate bmtk_ateam
      In the "examples"->"biophys_components"->"mechanisms", run the command:
	  nrnivmodl modfiles/
	  It will create a new folder named "x86_64". 
      If you have trouble in this step, make sure you have deleted "x86_84" before compling the modfiles.

   In the folder "571654895_example", run the simulations using the following commands:
 	  source activate bmtk_ateam
	  python build_network.py
	  python run_bionet.py
	
4. Run EAP analysis:
   In the "codes" folder: run Yina_EAP_analysis_step4.ipynb to calculate and save extracellular action potential (EAP) of each electrode recordings for further analysis. 


In the "assets" folder, we have included the .json files and morphologies for all the models. We also included the templates for simulations, use build_network_in.py for inhibitory aspiny cells, and build_network_pc.py for excitatory spiny cells. For each cell, we need to change the ROTX,ROTZ,CELLID,POPNAME, accordingly. You can find these information in the "assets"->"morpholigies"=>human_celltypes_table.csv.
