%%% This script convert the emClarity database into pyTOM *em file, for Chimera Place Object plugin. 
%%% The Place Object plugin can help visualize the localzation and
%%% orientation of particles in 3D, developed by Dr Qu Kun.

%%% This script convert the emClarity dataset at the selected cycle to
%%% Dynamo convention first, then to pyTOM convention. Several functions
%%% from Dynamo are integrated to this script. For simplicity, only the
%%% coordinate, cross-correlation coefficient and orientation information
%%% are converted. Other information are not converted.

%%% by Tao Ni. Email: taoni@hku.hk


%%% clear the current workspace
clear;

%%% load emClarity database and define the cycle and stage you want work on %%%
subTomoMeta = load ('Emerge2C1allpca-150A.mat'); % load database.
cycle_number = '009' ; % define which cycle the meta data shoudd be converted.
stage_of_alignment = 'Avg_geometry'; % or 'Avg_geometry', or 'geometry', or 'RawAlign'
TemplateSearch_bin = 1 ;% to generate model file at the binning of templateSearch.
writeemClarityCsv = 0; % whether to write out csv file in emClarity format
writepyTOMem = 1 ;% whether to write out pyTom *em file for visualization in Chimera


%%%%%%%%%%%% No need to change after this line   %%%%%%%%%%%%%%%%
str_cycle_number = strcat('cycle',cycle_number);
meta = subTomoMeta.subTomoMeta.(str_cycle_number).(stage_of_alignment);
fn = fieldnames(meta);

for k=1:numel(fn)
    disp(fn{k});
    csv = meta.(fn{k});
    tomogram_IDX=k;
    csv_removed_outlier = csv(csv(:,26) ~= -9999,:);
    if (writeemClarityCsv)
        csvfilename = strcat(fn{k}, '_bin', num2str(TemplateSearch_bin),'.csv');
        coordinate_file = strcat(fn{k}, '_bin', num2str(TemplateSearch_bin),'.txt');
        tomogram_coordinates =  csv_removed_outlier(:,11:13) / TemplateSearch_bin;
        writematrix(csv_removed_outlier,csvfilename,'delimiter','tab');
        writematrix(tomogram_coordinates, coordinate_file,'delimiter','tab');
    end
    
    if (writepyTOMem)
        csv = csv_emClarity2dynamo(csv_removed_outlier,TemplateSearch_bin, tomogram_IDX);
        [motl, wedge]=dynamo__table2motl(csv);  % convert table to pyTOM format    
        %%% save *em into disk
        pyTOM_file= strcat(fn{k},'_bin',num2str(TemplateSearch_bin),'.em');
        dynamo_write_emfile(motl,pyTOM_file);
    end
end


%%%% emClarity_csv2dynamo_table
function dynamo_table = csv_emClarity2dynamo(emClarity_csv_original,tomogram_bin, tomogram_IDX)

    emClarity_csv = emClarity_csv_original(emClarity_csv_original(:,26) ~= -9999,:);  % only take the good ones
    number_of_subtomograms = size(emClarity_csv,1)
    
    dynamo_table = zeros(number_of_subtomograms,35);
    %initialize a dynamo convention table format
    dynamo_convention_base = [1 1 1 0 0 0 0 0 0 0.5 0 0 1 -60 60 -60 60 0 0 1 1 0 0 0 0 0 0 0 0 0 0 1 0 1 0];
    for i = 1:number_of_subtomograms
        %%% for dynamo convention table
        dynamo_table(i,:) = dynamo_convention_base;
        dynamo_table(i,24:26) = emClarity_csv(i,11:13) / tomogram_bin;  % consider the binning of dynamo tomogram
        dynamo_table(i,20) = tomogram_IDX; % tomogram ID
        dynamo_table(i,1) = i; % subtomogram ID
        dynamo_table(i,10) = emClarity_csv(i,1); % cross correlation coefficient
        dynamo_table(i,22) = emClarity_csv(i,26); % class number
        rotation_matrix_emClarity = emClarity_csv(i,17:25);
        dynamo_rm = emClarity2dynamo_rm(rotation_matrix_emClarity);
        dynamo_rm = dynamo_rm' ;
        [tdrot,tilt,narot]=local_matrix2euler(dynamo_rm);
        dynamo_table(i,7:9) = [tdrot,tilt,narot]; % rotation matrix    
    end
end

function dynamo_rm = emClarity2dynamo_rm(rotation_matrix_emClarity)
%%% basically reshape function in MATLAB.
%%% dynamo_rm = reshape(rotation_matrix_emClarity, [3,3]);
%%% modified from Alister Burt.
    dynamo_rm(1,1) = rotation_matrix_emClarity(1);
    dynamo_rm(2,1) = rotation_matrix_emClarity(2);
    dynamo_rm(3,1) = rotation_matrix_emClarity(3);
    dynamo_rm(1,2) = rotation_matrix_emClarity(4);
    dynamo_rm(2,2) = rotation_matrix_emClarity(5);
    dynamo_rm(3,2) = rotation_matrix_emClarity(6);
    dynamo_rm(1,3) = rotation_matrix_emClarity(7);
    dynamo_rm(2,3) = rotation_matrix_emClarity(8);
    dynamo_rm(3,3) = rotation_matrix_emClarity(9);
end

function [tdrot,tilt,narot] = local_matrix2euler(rm);
    tol=1e-4;

%warning('indetermination in defining narot and tdrot: rotation about z');

% rm(3,3) ~= -1
    if abs(rm(3,3)-1)<tol;
        tilt=0;
        narot=atan2(rm(2,1),rm(1,1))*180/pi;
        tdrot=0;
        return
    end

% rm(3,3) ~= -1
    if abs(rm(3,3)+1)<tol;
        tdrot=0;
        tilt=180;
        narot=atan2(rm(2,1),rm(1,1))*180/pi;
    
        return
    end

 % General case 
    tdrot  = atan2(rm(3,1),rm(3,2));
    tilt   = acos(rm(3,3));
    narot  = atan2(rm(1,3),-rm(2,3));

    tilt=tilt*180/pi;
    narot=narot*180/pi;
    tdrot=tdrot*180/pi;
end


function my_volume = dynamo_read_emfile(filename)

% reads an EM file into the MATLAB workspace
%
% INPUT
%         filename:         string with the name of a file containing a
%                           volume or image in EM format
%
% OUTPUT
%         volume:           2D or 3D Matlab array
%
% SYNTAX:
%        volume=dynamo_read_emfile(filename);
%
% NOTE:
% Using the file format conventions of the TOM PACKAGE



 
 
% Author: Daniel Castano-Diez, April 2012 (daniel.castano@unibas.ch)
% Copyright (c) 2012 Daniel Castano-Diez and Henning Stahlberg
% Center for Cellular Imaging and Nano Analytics 
% Biocenter, University of Basel


%f=sprintf('Reading EM-file: %s',filename);disp(f);

fid = fopen(filename,'r','ieee-le');
if fid==-1;
    error(sprintf('Unable to open: %s',filename));
end;

% according to EM file standards
machine_type   = fread(fid,[1],'uint8');
noinfo         = fread(fid,[2],'uint8');
data_type      = fread(fid,[1],'uint8');
fclose(fid);


switch machine_type
    
    case 0
        fid = fopen(filename,'r','ieee-be');
    case 3
        fid = fopen(filename,'r','ieee-be');
    case 5
        fid = fopen(filename,'r','ieee-be');
    case 6
        fid = fopen(filename,'r','ieee-le');
    otherwise
        my_message=sprintf('Format error at file %s:  unable to identify endianness for apparent machine type %d',filename,machine_type);
        error(my_message);
end


% Reads the header to determine how many entries are contained in the
% density map
words_in_header=128;
header = fread(fid,words_in_header,'uint32');
nx = header(2);
ny = header(3);
nz = header(4);
number_entries=nx*ny*nz;


% Reads the density map according to the determined data type
switch data_type
    case 1
        my_volume = uint8(reshape(fread(fid,number_entries,'uint8'),nx, ny, nz));
        
    case 2
        my_volume = reshape(fread(fid,number_entries,'uint16'),nx, ny, nz);
        
    case 5
        my_volume = reshape(fread(fid,number_entries,'float32'),nx, ny, nz);
        
    otherwise
        my_message=sprintf('Format error at file %s:  unable to identify apparent data type %d',filename,data_type);
        error(my_message);
        
end

fclose(fid);

end


function dynamo_write_emfile(volume, filename)
% Writes a file following the .em format 
% (Martinsried's convention)
% 
% dynamo_write_emfile(volume, filename)
 
 
% Author: Daniel Castano-Diez, April 2012 (daniel.castano@unibas.ch)
% Copyright (c) 2012 Daniel Castano-Diez and Henning Stahlberg
% Center for Cellular Imaging and Nano Analytics 
% Biocenter, University of Basel
 

[sx,sy,sz] = size(volume);

% writes in little endianness
fid = fopen(filename,'w','ieee-le');

if fid==-1     
    error('Cannot open file %s for writing.',filename);
end;

% first words in header:
fwrite(fid,6,'char');
fwrite(fid,0,'char');
fwrite(fid,0,'char');

% determines the volume type
if(isa(volume,'uint8'))
    fwrite(fid,1,'char');
else
    fwrite(fid,5,'char');
end;

% writes dimensions of the cube
fwrite(fid,sx,'uint32');
fwrite(fid,sy,'uint32');
fwrite(fid,sz,'uint32');


% rest of header remains zero
header_rest=496;
for i = 1:header_rest
    fwrite(fid,0,'char');
end


if isa(volume,'uint8')
    % respects original uint8
    fwrite(fid,volume,'uint8');
else
    % converts into float
    fwrite(fid,volume,'float');
end

fclose(fid);
end



function [motl,wedgelist] = dynamo__table2motl(table)
% full conversion of AV3-style motl to Dynamo-style table
% 
% INPUT
%
%     table            
%
%
% OUTPUT
%
%     table:             as command line variable
%     wedgelist:         as command line variable

% Note: this function uses the following AV3 convention:
% 
%       1         : Cross-Correlation Coefficient
%       2         : x-coordinate in full tomogram
%       3         : y-coordinate in full tomogram
%       4         : particle number
%       5         : running number of tomogram - used for wedgelist
%       6         : index of feature in tomogram (optional)
%       8         : x-coordinate in full tomogram
%       9         : y-coordinate in full tomogram
%       10        : z-coordinate in full tomogram
%       11        : x-shift in subvolume - AFTER rotation of template
%       12        : y-shift in subvolume - AFTER rotation of template
%       13        : z-shift in subvolume - AFTER rotation of template
%     ( 14        : x-shift in subvolume - BEFORE rotation of template )
%     ( 15        : y-shift in subvolume - BEFORE rotation of template )
%     ( 16        : z-shift in subvolume - BEFORE rotation of template )
%       17        : Phi (in deg)
%       18        : Psi
%       19        : Theta 
%       20        : class no
%
%   (This paragraph has been written from the supporting help in TOM package, 2008)

 
 
% Author: Daniel Castano-Diez, April 2012 (daniel.castano@unibas.ch)
% Copyright (c) 2012 Daniel Castano-Diez and Henning Stahlberg
% Center for Cellular Imaging and Nano Analytics 
% Biocenter, University of Basel

%table =dynamo_file_checkout(table);
N=size(table,1);

motl=zeros(20,N);

%table = dynamo_table_sort(table);
motl(17,:) = -table(:,9)'; % phi   <- -narot;
motl(18,:) = -table(:,7)'; % psi   <- -tdrot
motl(19,:) = -table(:,8)'; % theta <- -tilt

% insert shifts
motl(11,:) = table(:,4);
motl(12,:) = table(:,5);
motl(13,:) = table(:,6);

% identities
motl(4,:) = table(:,1);
% cc
motl(1,:)=table(:,10);
% volume positions
motl(8,:)  = table(:,24);
motl(9,:)  = table(:,25);
motl(10,:) = table(:,26);

% column 20:   tomogram
motl(5,:)  = table(:,20);
% column 21:   tomogram area
motl(6,:)  = table(:,21);
% column 22:   class
motl(20,:) = table(:,22);

wedgelist = zeros(N,3);
wedgelist(:,1) = table(:,1);
wedgelist(:,2) = table(:,14);
wedgelist(:,3) = table(:,15);
end
