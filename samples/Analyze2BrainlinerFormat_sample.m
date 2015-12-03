% Sample script to create brainliner data (HDF5).
%
% This script convert ANALYZE data to brainliner data format
% for Hand Shape Decoding http://brainliner.jp/data/brainliner/Hand_shape_decoding_Rock_Papers/
% 
% Required:
%                  spm5 - http://www.fil.ion.ucl.ac.uk/spm/
%   ATR-DNI CommonTools - https://github.com/ATR-DNI/CommonTools/tree/master/MATLAB
%
% ANALYZE data and ROI mask data are available from 
% http://brainliner.jp/data/brainliner/Hand_shape_decoding_Rock_Papers/
% Download rawData.tgz from Supplementary Files of above site.
%
%
% Use Matlab >= 2012a
%
% by brainliner-admin 2015/2/3 
% Modified 2015/12/1
%

clear ;

%% Param
% set path to the SPM5
if isempty( which('spm') )
    addpath(''); % specify the path to the SPM5
end
% set path to the ATR-DNI CommonTools
if isempty( which('writeHDF5FromStruct') )
    addpath(genpath('')); % specify the path to the ATR-DNI CommonTools
end

% path to the root directory of the data
dataPath  = './';
% path to the Analyze files
rawDataPath = [ dataPath, 'analyze/' ];
% path to the ROI mask data
RoiMaskPath = [ dataPath, 'roi/'];

% outputfile name
outputHdf5File = 'RockPaperScissors.h5';

%% experimental protocol
% Condition labels for each run:
% 1:'rest', 2:'rock', 3:'scissors', 4:'paper' 
labelsRunBlocks  = {...    % {ltypes}{runs}[1 x blocksPerRun]
	[ 1  3  2  4  2  4  3  1 ];...      % run 1
	[ 1  3  4  2  2  3  4  1 ];...      % run 2
	[ 1  2  3  4  2  4  3  1 ];...      % ...
	[ 1  4  2  3  2  3  4  1 ];...
	[ 1  3  4  2  2  3  4  1 ];...
	[ 1  4  3  2  4  2  3  1 ];...
	[ 1  4  2  3  2  4  3  1 ];...
	[ 1  3  4  2  3  2  4  1 ];...
	[ 1  3  2  4  2  4  3  1 ];...
	[ 1  2  3  4  3  4  2  1 ]};
nRun   = length(labelsRunBlocks ) ;
nBlock = length(labelsRunBlocks{1}) ;
% Number of samples per label:
samplesPerLabel = 4;

conditionLabel = zeros( nRun*nBlock*samplesPerLabel, 1);
runNumber = zeros(size(conditionLabel));
for runIdx = 1:nRun
    idxRange = (runIdx-1)*samplesPerLabel*nBlock+1:runIdx*samplesPerLabel*nBlock;
    conditionLabel(idxRange) = reshape( repmat(labelsRunBlocks{runIdx}, [ samplesPerLabel, 1]), samplesPerLabel*nBlock, 1);    
    runNumber(idxRange) = runIdx ;
end

%% load data
% brain data
files  = dir([rawDataPath, '*.img']);
nFiles = length(files);
runNames = {'a','b','c','d','e','f','g','h','i','j'};

cnt  = 0;
xdim = 64;
ydim = 64;
zdim = 50;
spaceReso = 3; % [mm]
dataSize = [ xdim ydim zdim ]';
nRun = length(runNames);
fmriTC = zeros(nFiles,xdim,ydim,zdim); % fMRI time course

for fileIdx = 1:nFiles
    fileInfo=spm_vol(sprintf([rawDataPath, files(fileIdx).name]));
    [y xyzLocation]=spm_read_vols(fileInfo);
    nVoxel = size(xyzLocation, 2);
    minPosi = min(xyzLocation, [], 2);
    
    for voxIdx = 1:nVoxel
        pIdx = ( xyzLocation(:, voxIdx) - minPosi )/spaceReso + 1;
        fmriTC(fileIdx,pIdx(1),pIdx(2),pIdx(3)) = y(voxIdx);
    end

    fprintf([ num2str(fileIdx) '/', num2str(nFiles), ' files are loaded.\n']);
end

% roi
roiList = { 'CB_LHand', 'CB_RHand', 'M1_LHand', 'M1_RHand', 'SMA_LHand', 'SMA_RHand' };
nRoi    = length(roiList);
roiMask = cell(nRoi, 1);
stats   = zeros(xdim, ydim, zdim);
posiIdx = zeros(1,3);
for roiIdx = 1:nRoi
    roiData = load([ RoiMaskPath, 'VOX_', roiList{roiIdx}, '.mat' ]);
    roiPosi = roiData.( ['VOX_', roiList{roiIdx}] )(1:3,:) ;
    statsTemp  =  roiData.( ['VOX_', roiList{roiIdx}] )(4,:) ;
    roiMask{roiIdx} = zeros(xdim, ydim, zdim);
    
    for idx = 1:size(roiPosi,2)
        posiIdx(:) = (roiPosi(:,idx) + spaceReso/2)/spaceReso + dataSize/2 ;        
        roiMask{roiIdx}(posiIdx(1),posiIdx(2),posiIdx(3)) = 1;
        stats( posiIdx(1),posiIdx(2),posiIdx(3) ) = statsTemp(idx) ;
    end
    
end

% location
location = zeros(xdim, ydim, zdim, 3 );
reso = 3 ; % space resolution [mm]
for xIdx = 1:xdim
    x = reso*(xIdx - ( xdim/2 + 0.5)) ;
    for yIdx = 1:ydim
       y = reso*(yIdx - ( ydim/2 + 0.5)) ;
       for zIdx = 1:zdim
            z = reso*(zIdx - ( zdim/2 + 0.5)) ;
            
            location(xIdx,yIdx,zIdx,1) = x;
            location(xIdx,yIdx,zIdx,2) = y;
            location(xIdx,yIdx,zIdx,3) = z;            
        end
    end
end


%% create hdf5 file
%if exist(hdf5File,'file')
%    delete(hdf5File);
%end

% fileHeader
data.fileHeader.affiliation = 'ATR';
data.fileHeader.date = '2010/05/11';
data.fileHeader.description = 'The fMRI data was obtained from subject while he make finger-form of rock/paper/scissors.';
data.fileHeader.license = 'pddl';
data.fileHeader.schema = {'default:http://brainliner.jp/1.0'};
data.fileHeader.schemaDef = {''};
data.fileHeader.subjectLabel = 'SS100511';
data.fileHeader.subjectSex = 'male' ;
data.fileHeader.subjectSpecies = 'human' ;
data.fileHeader.title = 'Rock-paper-scissors Movement Task' ;

% group1 ( voxel data )
data.group1.props.preprocessing.overview = ...
    'Collected in ATR Japan. These fMRI data were realigned, coregistered, and resliced by SPM5.' ;
for roiIdx = 1:nRoi
    data.group1.props.roi.(roiList{roiIdx}) = roiMask{roiIdx} ;
end
data.group1.props.stats.tval = stats ;
data.group1.props.coordinateSystem = 'individual';
data.group1.props.description = 'fMRI timeline';
data.group1.props.hardware = 'fMRI';
data.group1.props.location = location ;

data.group1.props.samplingRate = 0.2 ;
data.group1.props.startTime = 0 ;
data.group1.props.title = 'fMRI experimental data' ;
data.group1.props.type = 'fMRI';
data.group1.props.voxelSize = [ 3 3 3 ] ;

data.group1.data = fmriTC ;

% group2 ( label )
data.group2.props.description = '1 = rest, 2 = rock, 3 = scissors, 4 = paper' ;
data.group2.props.title = 'behavior label ( shape of hand )' ;
data.group2.props.type = 'movement' ;
data.group2.props.samplingRate = 0.2 ;
data.group2.data = conditionLabel ;

labelList = {'rest', 'rock', 'scissors', 'paper'}; 
data.group2.props.dataNames = labelList(data.group2.data);

% group3 ( run )
data.group3.props.description = 'Run number' ;
data.group3.props.title = 'Experimental design (runs)' ;
data.group3.props.type = 'expDesign' ;
data.group3.props.samplingRate = 0.2 ;
data.group3.data = runNumber ;


%% saveFile
display('data file created.');
display('save to H5 file');
writeHDF5FromStruct(outputHdf5File, data);



