# em2emClarity

This script aims to use the subregion em file to find the points within user-defined distance from the emClarity database.

Developed by Tao Ni, E-mail: taoni@hku.hk

NOTES:

1) if none of the  subtomogram in a subrgion is close to the points in the em file, the whole subregion will be removed (use rmfield command).

2) the subtomograms which do not satisfy the restraits, will be ignored (labelled as -9999 in column 26)

3) This script is usually used after manual inspection of the subregion em files and saved as <subregion>_binX_display.em
  
4) If the resulting emClarity database has problem for further process, one can write out the emClarity database to *csv files (for each subregion) and re-initialize the project database.

  
TO RUN:
  
1) Have the emClarity database ready, say emClarity_tutorial.mat

2) Have the pyTOM (*em) files ready, or IMOD model files ( this need to convert to txt file using model2point). Name convention of files:
  
     <subregion_name>_binX_display.em
    
     <subregion_name>_binX.txt
       
3) Put the emClarity database and subregion files in the same directory for MATLAB
       
4) Change the parameters in the em2emClarity.m script (for example):

       subTomoMeta = load('emClarity_tutorial.mat');
       
       output_file = 'emClarity_tutorial_clean.mat' ;
       
       template_search_bin = 4  % binning of the subregion em file.
       
       pixel_size = 1.34 % defined as pixel size in fixedStacks/
       
       cycle_number = '001' ; % define which cycle you want to expand
       
       stage_of_alignment = 'RawAlign' ;  % or 'Avg_geometry', 'RawAlign', 'geometry'
       
       use_modeltxt_file = 0;  % by default, use <subregion>_binX_display.em file
       
       min_distance_threshold = 300 % in angstrom. distance restraints
       
 5) Run the script using MATLAB (version after 2019).
      
