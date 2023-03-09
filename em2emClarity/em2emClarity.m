%%% This script aims to use the subregion em file to find the points within 
%%% user-defined distance from the emClarity database.
%%% 1) if none of the  subtomogram in a subrgion is close to the points in
%%% the em file, the whole subregion will be removed (use rmfield command).
%%% 2) the subtomograms which do not satisfy the restraits, will be ignored
%%% (labelled as -9999 in column 26)
%%% 3) This script is usually used after manual inspection of the subregion
%%% em files and saved as <subregion>_binX_display.em
%%% by Tao Ni

%%% clear the current workspace
clear;
disp('Make sure your subregion em file name convention follows:')
disp('  ')
disp('<subregion>_binX_display.em')

%%% load emClarity database and define which cycle you want to work on %%%

subTomoMeta = load('alpha12.mat');
output_file = 'alpha12-300.mat'
template_search_bin = 4  % binning of the subregion em file.
pixel_size = 1.34 % defined as pixel size in fixedStacks/
cycle_number = '001' ; % define which cycle you want to expand
stage_of_alignment = 'RawAlign' ;  % or 'Avg_geometry', 'RawAlign', 'geometry'
use_modeltxt_file = 0;  % by default, use <subregion>_binX_display.em file
min_distance_threshold = 300 % in angstrom. distance restraints

%%%%%%%%%%%% probably no need to change after this line   %%%%%%%%%%%%%%%%

str_cycle_number = strcat('cycle',cycle_number);
meta = subTomoMeta.subTomoMeta.(str_cycle_number).(stage_of_alignment);
tiltGeometry = subTomoMeta.subTomoMeta.tiltGeometry;
reconGeometry = subTomoMeta.subTomoMeta.reconGeometry;
tomoName = subTomoMeta.subTomoMeta.mapBackGeometry.tomoName;
ctfGroupSize = subTomoMeta.subTomoMeta.ctfGroupSize;
mapBackGeometry = subTomoMeta.subTomoMeta.mapBackGeometry;
fn = fieldnames(meta);
min_distance_threshold = min_distance_threshold / pixel_size;


%%% main loop to update csv for each tomogram
particle_count_before_clean=0;
particle_count_after_clean=0;

for k=1:numel(fn)
    subregion_name = strcat(fn{k});
    emfile = strcat(subregion_name,'_bin',num2str(template_search_bin),'_display.em')
    modelFile = strcat(subregion_name,'_bin',num2str(template_search_bin),'.txt');
    
    %subregion_name = erase(name, str)
    csv_before_update = meta.(subregion_name);  % load the subregion from emClarity database
    particle_count_before_clean = particle_count_before_clean + size(csv_before_update(csv_before_update(:,26)~=-9999),1)
    if exist(emfile)
        molt = dynamo_read_emfile(emfile);
        model_file_txt = molt(8:10,:)';
        tomogram_outlier_removed = find_close_neighbors(csv_before_update, ...
                              model_file_txt,template_search_bin, min_distance_threshold);
        particle_count_after_clean = particle_count_after_clean + size(tomogram_outlier_removed(tomogram_outlier_removed(:,26)~=-9999),1)                      
        meta.(subregion_name) = tomogram_outlier_removed;
    elseif exist(modelFile) & (use_modeltxt_file)
        model_file_txt = load(modelFile);
        tomogram_outlier_removed = find_close_neighbors(csv_before_update, ...
                              model_file_txt,template_search_bin, min_distance_threshold);
        particle_count_after_clean = particle_count_after_clean + size(tomogram_outlier_removed(tomogram_outlier_removed(:,26)~=-9999),1)                      
        meta.(subregion_name) = tomogram_outlier_removed;
    else
        display("emfile and modelFile do not exist");
        csv_before_update(:,26) = -9999;
        tomogram_outlier_removed = csv_before_update;
        subregion_name = strcat(fn{k})
        meta=rmfield(meta,(subregion_name));
        tiltGeometry = rmfield(tiltGeometry, (subregion_name));
        reconGeometry = rmfield(reconGeometry, (subregion_name));
        tomoName = rmfield(tomoName, (subregion_name));
        ctfGroupSize = rmfield(ctfGroupSize, (subregion_name));
        
        %%% remove the corresponding subregion line in coords, this is the
        %%% major thing to update in the metadata.
        tmp_id = split(subregion_name, '_');
        id = str2num(tmp_id{end});
        tilt = strjoin(tmp_id(1:size(tmp_id,1)-1), '_');
        mapBackGeometry.(tilt).coords(id,:) = 0;
        mapBackGeometry.(tilt).nTomos = mapBackGeometry.(tilt).nTomos -1;
    end
end


subTomoMeta.subTomoMeta.(str_cycle_number).(stage_of_alignment) = meta;
subTomoMeta.subTomoMeta.tiltGeometry = tiltGeometry;
subTomoMeta.subTomoMeta.reconGeometry = reconGeometry;
subTomoMeta.subTomoMeta.mapBackGeometry.tomoName = tomoName;
subTomoMeta.subTomoMeta.ctfGroupSize = ctfGroupSize;
subTomoMeta.subTomoMeta.mapBackGeometry = mapBackGeometry;



%%% remove the subregions in the mapBackGeometry if they are not
%%% present in the current stage of alignment.
subregions = fieldnames(meta);
tilt2keep = cell(numel(subregions),1);  %%%% assign a cell to store the names of subregion2remove
for k = 1 : numel(subregions)
    subregion_name = subregions{k};
    tmp = split(subregion_name,'_');
    tilt2keep{k} = strjoin(tmp(1:size(tmp,1)-1), '_');  %%% extract the tilt series base name
end
   % tilt2keep = tilt2keep(~cellfun(@isempty,subregion2remove));  %%% remove the empty cell
tilt2keep = unique(tilt2keep); 
tilt2keep_tomoName=[tilt2keep;'tomoName'];
tilt2remove = setdiff(fieldnames(subTomoMeta.subTomoMeta.mapBackGeometry), tilt2keep_tomoName);
subTomoMeta.subTomoMeta.mapBackGeometry = rmfield(subTomoMeta.subTomoMeta.mapBackGeometry, tilt2remove);

subTomoMeta = subTomoMeta.subTomoMeta;
save(output_file,'subTomoMeta','-v7.3');

particle_count_before_clean
particle_count_after_clean


%%% define a function to find the nearest particles of each tomogram csv file
function  tomogram_outlier_removed = find_close_neighbors(csv_before_update, ...
                              model_file_txt,template_search_bin, min_distance_threshold)
    number_of_subvolumes = size(csv_before_update,1);
    for i=1:number_of_subvolumes
        if csv_before_update(i,26) ~= -9999
            distance = zeros(size(model_file_txt,1),1);
            for k = 1:size(model_file_txt,1)
                distance(k) = norm(csv_before_update(i, 11:13) - model_file_txt(k,1:3) * template_search_bin);
                if min(distance) > min_distance_threshold
                    csv_before_update(i, 26) = -9999;
                end
            end
        end
    end
    tomogram_outlier_removed = csv_before_update;               
end



 
 
% Author: Daniel Castano-Diez, April 2012 (daniel.castano@unibas.ch)
% Copyright (c) 2012 Daniel Castano-Diez and Henning Stahlberg
% Center for Cellular Imaging and Nano Analytics 
% Biocenter, University of Basel
% 
%
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
