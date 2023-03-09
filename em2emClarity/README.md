# em2emClarity

This script aims to use the subregion em file to find the points within user-defined distance from the emClarity database.

Developed by Tao Ni, E-mail: taoni@hku.hk

NOTES:

1) if none of the  subtomogram in a subrgion is close to the points in the em file, the whole subregion will be removed (use rmfield command).

2) the subtomograms which do not satisfy the restraits, will be ignored (labelled as -9999 in column 26)

3) This script is usually used after manual inspection of the subregion em files and saved as <subregion>_binX_display.em
  
4) If the resulting emClarity database has problem for further process, one can write out the emClarity database to *csv files (for each subregion) and re-initialize the project database.
