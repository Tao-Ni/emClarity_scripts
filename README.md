# emClarity_scripts
scripts for emClarity

This script convert the emClarity database into pyTOM *em file, for Chimera Place Object plugin. 
The Place Object plugin can help visualize the localzation and orientation of particles in 3D, developed by Dr Qu Kun.

This script convert the emClarity dataset at the selected cycle to Dynamo convention first, then to pyTOM convention. Several functions from Dynamo are integrated to this script. For simplicity, only the
coordinate, cross-correlation coefficient and orientation information are converted. Other information are not converted.

Developed by Tao Ni, E-mail: taoni@hku.hk

To run the script:

(1) Have an emClarity database ready, an example database is provided: Emerge2C1allpca-150A.mat

(2) Fill the parameters in the script (For Example):

      subTomoMeta = load ('Emerge2C1allpca-150A.mat'); % load database.
      
      cycle_number = '009' ; % define which cycle the meta data shoudd be converted.
      
      stage_of_alignment = 'Avg_geometry'; % or 'Avg_geometry', or 'geometry', or 'RawAlign'
      
      TemplateSearch_bin = 4 ;% to generate model file at the binning of templateSearch.
      
      writeemClarityCsv = 0; % whether to write out csv file in emClarity format
      
      writepyTOMem = 1 ;% whether to write out pyTom *em file for visualization in Chimera
      
(3) Run the script in the MATLAB.
