% SynD - Synapse Detector
% 
% NOTE: This version was just rewritten to use object oriented
% matlab, it should be equivalent to old version.
% Please report any problems.
% 
% Johannes Hjorth
% hjorth@kth.se
%
% old: j.j.j.hjorth@damtp.cam.ac.uk
% old: johannes.hjorth@cncr.vu.nl
%
% Sabine Schmitz
% sabine.mosch@uni.lu
%
% old: sabine.schmitz@cncr.vu.nl
%
% The latest version of the code is available from:
% http://www.johanneshjorth.se/SynD
% 
% You can also find the code at:
%
% https://github.com/Hjorthmedh/SynD
%
% Please cite our article:
%
% Schmitz SK, Hjorth JJJ, Joemai RMS, Wijntjes R, Eijgenraam S,  
% de Bruijn P, Georgiou C, de Jong APH, van Ooyen A, Verhage M, 
% Cornelisse LN, Toonen RF, Veldkamp W 
% Automated Analysis of Neuronal Morphology, Synapse Number and 
% Synaptic Recruitment, Journal of Neuroscience Methods 195 (2011) 185–193
% DOI: 10.1016/j.jneumeth.2010.12.011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code also uses:
% LSM file toolbox by Peter Li
% Tiffread by Francois Nedelec
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% !!! If class new changed, constructor name must be updated
   
classdef SynD_extract2017 < handle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Description of data structures:
%
%
% Data directly relating to the image
%
% The image is stored in data.image, which is a four dimensional matrix
% data.height x data.width x 3 x data.num. The size of each pixel is 
% specified in data.xyRes. The three channels of the image are identified
% by detection.morphChannel, detection.synChannel and detection.XChannel.
% The colours used to display the image are optionally remapped 
% (if detection.enableRemapChan is set), The remapping  is specified by 
% data.remapChan, here red=1, green=2 and blue=3.
%
% SynD can read in tiff-stacks containing multiple images, but it only
% processes the first image. The dispInfo.curImage is a legacy from when
% the program allowed the user to select which of the images in a stack
% to use. The functionality can be reimplemented, since all other functions
% uses dispInfo.curImage to extract the relevant image to use.
%
% The soma is defined by data.somaMask which is data.height x data.width
% in size, pixels belonging to the soma are set to 1, the others are 0.
% More than one soma is allowed in the mask. In addition there is a 
% data.somaMeasureMask that optionally defines circular regions of interest
% in the soma that can be specified by the user. These ROIs are derived from
% data.measurePoint which specifies the pixel index of the center of the
% ROIs.
%
% The neurites are specified in data.neuriteMask where a non-zero values
% means a pixel is part of a neurite. During editing of a neurite mask
% a pixel belonging to a neurite pixel connected to a soma is specified
% by a 1, and for a non-connected neurite the pixel is set to 2. After the
% neurite editing all neurite pixels must be connected to a soma. If this
% is not the case then dispInfo.needCleaning is set.
%
% The shortest distance to the soma for each pixel in the neurite mask
% is calculated and saved in data.distMask. The data.distMask does not
% include pixels in the neurite padding. The neurite skeleton is calculated 
% and stored as a binary mask in data.skeleton.
%
% There are several representations of the synapse. During editing the
% data.synapseMask is a binary mask similar to the data.somaMask.
% A synapse has to be within detection.neuritePadding micrometers of
% a neurite, and it has to contain at least one synapse center. 
% After editing the center of each synapse is located, and the index
% (one dimensional offset into the data.synapseMask) of each synapse
% center is stored in data.synapseCenter. In addition the pixel of
% each synapse is stored in a cell array data.synapsePixels, in the
% case where a synapse region has more than one center, synapse pixels
% are assigned to the closest center.
%
% There is an exclude mask to mark regions of the image that should not
% be included in the processing stored in data.excludeMask. Note though
% that the data.excludeMask is not currently restored if the user does
% a batch reanalyse. The data.excludeMask, data.somaMeasureMask, 
% data.neuriteMask, and data.synapseMask have the same dimensions as 
% data.somaMask.
%
% There are a number of measures calculated for the neuron. A general rule
% is that all values internally are stored in either SI-units, in a few
% cases they are stored in pixels.
%
% For display purposes the maximal intensity of the red, green and blue
% channels are stored in data.maxRed, data.maxGreen and data.maxBlue.
% The intensity histogram for the image is specified by 
% data.intensityHistogramEdges, data.morphHist, data.synHist and 
% data.XHist and detection.intensityBinSize. 
%
%
% Display information
%
% The legacy variable data.curImg specifies which image (if there are
% multiple, should be displayed and processed).
%
% In each stage only the relevant channels are shown by default, this
% is tracked by dispInfo.showRed, dispInfo.showGreen and dispInfo.showBlue
% together with dispInfo.showMask and dispInfo.showSkeleton. The user can 
% change the contrast of the image channels by clicking on the histograms, 
% the intensity scaling is stored in  dispInfo.scaleRed, dispInfo.scaleGreen 
% and dispInfo.scaleBlue.
%
% The colour of the different masks are specified by dispInfo.somaColor,
% dispInfo.neuriteColor, dispInfo.synapseColor and 
% dispInfo.measurePointColor.
%
% The current axis of the image are specified by dispInfo.axis, the variable
% is empty if no image is displayed.
%
% The program is divided into five steps, the GUI shown is tracked by
% dispInfo.state, which is either 'load', 'soma', 'neurite', 'synapse'
% or 'analyse'. As the user can jump between the GUI state, the program
% also tracks the total progress in dispInfo.stage (values 1-5).
%
% The position of the figure on the screen is remembered between sessions
% and stored in dispInfo.figPosition.
%
% Additional output can be displayed by setting dispInfo.verbose to true.
%
%
% Editing information
%
% These variables are only relevant to the current sessions and are not
% saved. The callback functions that allow editing of the different masks
% keep track of their state using editInfo.mode (0 = not editing, 1 = adding
% to neurite mask, -1 = subtracting from neurite mask, 2/-2 = adding/deleting
% from soma mask, 3/-3 = adding/deleting from synapse mask.
%
% The handle of the line currently drawn is stored in editInfo.line and
% the colour and width is specified by editInfo.color and editInfo.width.
% The coordinates are stored in editInfo.xLine and editInfo.yLine.
%
% An undo function is implemented, and remembers editInfo.maxUndo steps.
% The old states are stored in editInfo.undo.
%
%
% Automatic detection
%
% The soma is detected by thresholding the image at the intensity specified
% in detection.morphThreshold. If detection.singleSoma is set then only the 
% largest connected component is kept, then the soma mask is eroded using a 
% disk of radius detection.somaErodeRadius (specified  in pixels).
% The threshold will be automatically detected using cross entropy
% minimization if either detection.autoSomaThreshold is set or if the user
% clicks the "auto" button.
%
% Optional ROIs are placed in the soma. At most detection.nMeasurePoints 
% disks of radius detection.measurePointRadius pixels are added. If unable
% to add the maximal number of ROIs the program will warn the user.
%
% The default way of detecting neurites uses steerable filters, tracing 
% the neurites from the soma and outwards as long as the cost of adding a
% pixel is lower than detection.maxAddCost. This can be turned off in
% by setting detection.useSteerableFilters to false, the program then falls 
% back to a simpler thresholding detection for the neurites. The size
% of the steerable filters (sigma) are defined by detection.filterSize 
% (in meters). The internal variable detection.pixelQue tracks which pixels
% remain to be used as seeding points, this vector is grows and shrinks
% while the neurites are tracked.
%
% The neurite directions and the similiarty to a neurite (ie bright region
% surrounded by dark borders) for each pixel as calculated by the steerable 
% filter are stored in data.dirVect and data.rigidity. The maximal value
% of the latter is stored in data.lambdaMax. The relative contributions
% of the directions and rigidity to the pixel add cost can be adjusted by
% changing detection.connectCostLambda. Optionally the filter can be 
% elongated by setting detection.steerableFilterAlpha to a non-zero value.
%
% When detecting the neurite skeleton any branches smaller than 
% detection.minProtusionLength (meters) are removed.
%
% The synapse regions are detected by first Wiener filtering the image
% if detection.wiener is set to a pixel width value (NaN means off).
% If detection.backgroundFiltering is set then a smoothed image defined
% by detection.smoothingRadius is subtracted from the synapse channel.
% The image is then thresholded at detection.synThreshold number of
% standard deviations above noise level (calculated by using the pixels
% in the neurite mask). Alternatively the user can use cross-entropy
% to calculate the synapse threshold automatically. Any regions smaller 
% than detection.minSynapseSize are removed.
%
% Synapse center detection uses a single synapse template generated from
% synapse regions with only one unique maxima. The size of the template is
% detection.maxRadie x detection.maxRadie pixels and it is stored in 
% data.meanSynapse. The data in the morphology channel is deconvolved by
% the data.meanSynapse to find the synapse centers.
%
% The program internally tracks if it is in batch reanalyse mode using
% detection.reanalyseFlag. It also counts the number of times it has
% exported data detection.numSaved (this is reset if the config file is 
% overwritten).
%
%
% Derived measures:
%
% Soma measures
%
% The basic size properties of the soma(s) are stored in data.somaArea,
% data.somaMajorAxisLength and data.somaMinorAxisLength, which each are
% as long as there are number of somas.
%
% The channel intensities for the optional regions of interest in the soma
% are stored in data.somaMeasureMorph, data.somaMeasureSyn and 
% data.somaMeasureX.
%
%
% Neurite measures
%
% Sholl analysis is done both on the neurites and the synapses. The bins
% are specified in data.shollEdges (bin width defined by 
% detection.shollBinSize), and the histograms for the number of
% neurites and synapses in each bin are stored in data.shollDendHist and
% data.shollSynHist.
%
% The channel intensity gradients along the neurites are defiend by 
% data.gradientEdges, data.morphGradientMean, data.synGradientMean
% and data.XGradientMean. There are also standard deviations and
% standard errors of the mean (replace names with Std and SEM).
%
%
% Synapse measures
%
% The area in micrometers of a synapse is specified in data.synapseArea,
% and data.synapseDist stores the shortest distance from the synapse center, 
% following a neurite, to the soma. The synapses are ordered in the same way
% as in data.synapseCenter. The average intensities of the three image
% channels are stored in data.synapseIntensityMorphMean, 
% data.synapseIntensitySynMean and data.synapseIntensityXMean. There are
% also three corresponding variables with the standard error of the mean 
% for each synapse.
%
% The single synapse profile displayed during synapse detection is stored
% in data.meanSynapseMorph, data.meanSynapseSyn and data.meanSynapseX.
%
% The Sholl analysis for the neurites also applies for the synapses, and
% the number of synapses within each bin of data.shollEdges are stored in
% data.shollSynHist. The average intensities of these synapses are stored
% in data.shollIntMorphMean, data.shollIntSynMean, data.shollIntXMean
% as well as the standard deviations and standard error of the means in
% corresponding variables (replace Mean by Std or SEM in the names).
%
%
% Export information
%
% The original names of the images used are stored in a cell array
% data.fileName. The image directory is saved in data.loadPath, and the 
% export directory in data.exportPath. The export file names are specified
% data.exportXMLfile, data.exportSaveFile (specifies the MAT-file), 
% data.exportNeuriteMaskFile and data.exportSynapseMaskFile.
%
% The program knows what information to export by looking at 
% exportInfo.saveMask, exportInfo.saveSholl, 
% exportInfo.saveIntensityHistogram, exportInfo.saveMat, 
% exportInfo.saveSomaMeasure, exportInfo.saveSynapseProfile and
% exportInfo.saveGradient.
%
%
% For additional information see our article Schmitz, Hjorth et al 2011
% or visit www.johanneshjorth.se/SynD or software.incf.org/software/synd
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  properties
    
    % Location variables
    SynDscript = [];
    SynDscriptFull = [];
    SyndScriptPath = [];

    licenseOk = 1;
    
    configFile = [];
    
    enableDebug = 1;
    debugFlag = 0;
    

    referenceStr = ['Schmitz SK, Hjorth JJJ, ' ...
                    'Joemai RMS, Wijntjes R, ' ...
                    'Eijgenraam S,  de Bruijn P, Georgiou C, ' ...
                    'de Jong APH, van Ooyen A, Verhage M, ' ...
                    'Cornelisse LN, Toonen RF, Veldkamp W ', ...
                    'Automated Analysis of Neuronal Morphology, ' ...
                    'Synapse Number and Synaptic Recruitment, ' ...
                    'J Neurosci Methods 195 (2011) 185–193'];

    referenceDOI = '10.1016/j.neumeth.2010.12.011';
    
    % Data structure which contains image data and results

    data = struct('image', [], ...
                  'height', 0, ...
                  'width', 0, ...
                  'num', 0, ...
                  'maxRed', 0, ...
                  'maxGreen', 0, ...
                  'maxBlue', 0, ...
                  'xyRes', 0.20e-6, ... % meters 
                  'somaMask', [], ...
                  'neuriteMask', [], ...
                  'synapseMask', [], ...
                  'somaMeasureMask', [], ...
                  'includeMask', [], ...
                  'synapseCenter', [], ...
                  'synapseArea', [], ...
                  'synapseDist', [], ...
                  'synapsePixels', [], ... % should be curly brackets
                  'neuriteLength', NaN, ...
                  'distMask', [], ...
                  'synapseIntensityMorphMean', [], ... % morph-channel
                  'synapseIntensityMorphSEM', [], ...  % morph-channel
                  'synapseIntensitySynMean', [], ...   % syn-channel
                  'synapseIntensitySynSEM', [], ...    % syn-channel
                  'synapseIntensityXMean', [], ...     % X-channel
                  'synapseIntensityXSEM', [], ...      % X-channel
                  'meanSynapse', [], ... % Template used for synapse center detection
                  'meanSynapseMorph', [], ...
                  'meanSynapseSyn', [], ...
                  'meanSynapseX', [], ...    % Profile of X channel for synapse
                  'meanSynapseProfileDist', [], ...
                  'meanSynapseProfileMorph', [], ...
                  'meanSynapseProfileSyn', [], ...
                  'meanSynapseProfileX', [], ...
                  'fileName', [], ...
                  'skeleton', [], ...
                  'shollEdges', [NaN NaN NaN], ...
                  'shollDendHist', [], ...
                  'shollSynHist', [], ...
                  'shollIntMorphMean', [], ...
                  'shollIntMorphStd', [], ...
                  'shollIntMorphSEM', [], ...
                  'shollIntSynMean', [], ...
                  'shollIntSynStd', [], ...
                  'shollIntSynSEM', [], ...
                  'shollIntXMean', [], ...
                  'shollIntXStd', [], ...
                  'shollIntXSEM', [], ...
                  'somaMeasureMorph', NaN, ...
                  'somaMeasureSyn', NaN, ...
                  'somaMeasureX', NaN, ...
                  'somaArea', [], ...
                  'somaMajorAxisLength', [], ...
                  'somaMinorAxisLength', [], ...
                  'measurePoint', [], ...
                  'intensityHistogramEdges', [], ...
                  'morphHist', [], ...
                  'synHist', [], ...
                  'XHist', [], ...
                  'gradientEdges', [], ...   % Intensity gradient along neurites
                  'morphGradientMean', [], ...
                  'morphGradientStd', [], ...
                  'morphGradientSEM', [], ...
                  'synGradientMean', [], ...
                  'synGradientStd', [], ...
                  'synGradientSEM', [], ...
                  'XGradientMean', [], ...
                  'XGradientStd', [], ...
                  'XGradientSEM', [], ...
                  'exportPath', pwd, ...
                  'exportXMLfile', [], ...
                  'exportSaveFile', [], ...
                  'exportNeuriteMaskFile', [], ...
                  'exportSynapseMaskFile', [], ...
                  'loadPath', pwd, ...
                  'dirVect', [], ...
                  'rigidity', [], ...
                  'lambdaMax', 0, ...
                  'somaIntensityThreshold', NaN, ...
                  'synapseIntensityThreshold', NaN);
    
        
    % Information used to display the image, only relevant to session

    dispInfo = struct('curImg', 1, ...
                      'showRed', 1, ...
                      'showGreen', 1, ...
                      'showBlue', 1, ...
                      'showMask', 1, ...
                      'scaleRed', 1, ...
                      'scaleGreen', 1, ...
                      'scaleBlue', 1, ...
                      'axis', [], ...
                      'state', 'load', ... % Which GUI view is currently displayed
                      'stage', 1, ... % How far along are we in processing?
                      'somaColor', [1 1 1], ...
                      'neuriteColor', NaN, ... %[1 1 1]*0.7, ...
                      'synapseColor', NaN, ...
                      'measurePointColor', NaN, ...
                      'defaultMeasurePointColor', [1 1 0], ...
                      'showSkeleton', 0, ...
                      'showSynapseProfileDataPoints', 1, ...
                      'needCleaning', 0, ...
                      'figPosition', [], ...
                      'verbose', 0);   % Add additional debug plots
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    % Editing information, only relevant to session

    editInfo = struct('mode', 0, ...
                      'line', [], ...
                      'color', [1 1 1], ...
                      'width', 3, ...
                      'defaultWidth', 3, ...
                      'xLine', [], ...
                      'yLine', [], ...
                      'undo', struct('somaMask', [], ...
                                     'neuriteMask', [], ...
                                     'synapseMask', [], ...
                                     'includeMask', [], ...
                                     'synapseCenter', [], ...
                                     'measurePoint', [], ...
                                     'description', [], ...
                                     'state', [], ...
                                     'stage', 0), ...
                      'maxUndo', 10, ...
                      'fileTypeOrder', [1 2]); % LSM=1 or TIFF=2 default file type?


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    exportInfo = struct('saveMask', 1, ...
                        'saveSholl', 1, ...
                        'saveIntensityHistogram', 1, ...
                        'saveMat', 1, ...
                        'saveSomaMeasure', 0, ...
                        'saveSynapseProfile', 1, ...
                        'saveGradient', 1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameters used for detection

    detection = struct('morphChannel', 1, ...
                       'synChannel', 2, ...   
                       'XChannel', 3, ...     
                       'remapChan', [1 2 3], ... % Normally R=ch1, G=ch2, B=ch3
                       'enableRemapChan', 'off', ...
                       'singleSoma', 1, ...
                       'measurePointRadius', 5, ...
                       'nMeasurePoints', 10, ...
                       'wienerSize', 7, ... % NaN to disable
                       'morphThreshold', 200, ... % 0.03?
                       'somaErodeRadius', 15, ...
                       'autoSomaThreshold', 0, ... % Use cross-entropy minimization?
                       'useSteerableFilters', 1, ... % Should we use steerable filters
                       'filterSize', 1.8e-6, ... % in meters, corresponding to about 6 pixels in normal resolution, ...
                       'maxAddCost', 0.9, ... % Maximal cost if we still want to add pixel
                       'connectCostLambda', 0.7, ... % Weighting to lambda (vs eigenvectorcost)  Meijering 2004 uses elongated filters, set it to -1/3 to use theirs. If you change this, also consider decreasing the maxAddCost (e.g. 0.85)
                       'steerableFilterAlpha', 0, ... % -1/3
                       'minProtrusionLength', 2e-6, ... % Min length of thin neurites added
                       'pixelQue', [], ... 
                       'smoothingRadius', 25, ...
                       'synThreshold', 0.6, ...  % How many standard deviations above noise?
                       'backgroundFiltering', true, ...  % Subtract smoothed image?
                       'minSynapseSize', 0.35e-12, ... % corresponds to < 4 pixels, ...
                       'neuritePadding', 1e-6, ...  % How many meters outside neurite are synapses allowed to be 
                       'singleSynapseRadie', NaN, ... %0.5e-6, ... % micrometers, NaN gives guestimation of profile
                       'maxRadie', 10, ... % Pixel size of synapse template
                       'excludeSomaSynapses', true, ...
                       'numSaved', 0, ... % How many times has the user saved a neuron
                       'trimNeuriteSize', 12, ...  % Neurites shorter than 10 pixels do not count
                       'shollBinSize', 5e-6, ... % Bin size for Sholl analysis
                       'intensityBinSize', 100, ...
                       'reanalyseFlag', 0);
    
    
    bilder = struct('logo',[],'icon',[]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Where to find updates and version files

    urlBase = 'http://www.johanneshjorth.se/files/SynD/';
    verFile = 'version.txt';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    handles = struct('menuHelp', [], ...
                     'menuItemHelp', [], ...
                     'menuItemVersion', [], ...
                     'menuItemAbout', [], ...
                     'menuFile', [], ...
                     'menuItemLoad', [], ...
                     'menuItemBatchAnalyse', [], ...
                     'menuItemBatchReanalyse', [], ...
                     'menuItemImportMask', [], ...
                     'menuItemExport', [], ...
                     'menuItemRestart', [], ...
                     'menuEdit', [], ...
                     'menuItemUndo', [], ...
                     'menuItemExclude', [], ...
                     'menuItemClearExclude', [], ...
                     'menuItemMaxProjection', [], ...
                     'menuSettings', [], ...
                     'autoThreshold', [], ...
                     'menuItemSteerable', [], ...
                     'menuItemBackgroundFiltering', [], ...
                     'menuItemSteerableSettings', [], ...
                     'menuItemSynapseSettings', [], ...
                     'menuItemShollSettings', [], ...
                     'menuItemSavePath', [], ...
                     'menuDebug', [], ...
                     'menuItemDebug', [], ...
                     'measurePointMenu', [], ...
                     'fig', [], ...
                     'image', [], ...
                     'redHist', [], ...
                     'greenHist', [], ...
                     'blueHist', [], ...
                     'credits', [], ...
                     'next', [], ...
                     'prev', [], ...
                     'num', [], ...
                     'loadNeuron', [], ...
                     'nextStage', [], ...
                     'detectSoma', [], ...
                     'editSoma', [], ...
                     'addSoma', [], ...
                     'deleteSoma', [], ...
                     'addMeasurePoints', [], ...
                     'editMeasurePoints', [], ...
                     'somaMeasureRadiusLabel', [], ...
                     'somaMeasureRadius', [], ...
                     'measureLabel', [], ...
                     'detectNeurite', [], ...
                     'editNeurite', [], ...
                     'addThinNeurites', [], ...
                     'addNeurite', [], ...
                     'clean', [], ...
                     'killBlob', [], ...
                     'extendNeurites', [], ...
                     'showSkeleton', [], ...
                     'calculateLabel', [], ...
                     'skipSynapseDetection', [], ...
                     'exportSimpleAverages', [], ...
                     'detectSynapses', [], ...
                     'editSynapse', [], ...
                     'killSynBlob', [], ...
                     'singleSynapse', [], ...
                     'morphLabel', [], ...
                     'morphChannel', [], ...
                     'morphChannelNumber', [], ...
                     'synLabel', [], ...
                     'synChannel', [], ...
                     'synChannelNumber', [], ...
                     'XLabel', [], ...
                     'XChannel', [], ...
                     'XChannelNumber', [], ...
                     'remapLabel', [], ...
                     'remap', [], ...
                     'XYresLabel', [], ...
                     'XYres', [], ...
                     'red', [], ...
                     'green', [], ...
                     'blue', [], ...
                     'redLabel', [], ...
                     'greenLabel', [], ...
                     'blueLabel', [], ...
                     'mask', [], ...
                     'maskLabel', [], ...
                     'morphThreshLabel', [], ...
                     'morphThresh', [], ...
                     'guessSomaThreshold', [], ...
                     'somaErodeLabel', [], ...
                     'somaErodeRadius', [], ...
                     'singleSomaLabel', [], ...
                     'singleSoma', [], ...
                     'synThreshLabel', [], ...
                     'synThresh', [], ...
                     'neuritePaddingLabel', [], ...
                     'neuritePadding', [], ...
                     'minSynapseSizeLabel', [], ...
                     'minSynapseSize', [], ...
                     'growThreshLabel', [], ...
                     'growThresh', [], ...
                     'filterSizeLabel', [], ...
                     'filterSize', [], ...
                     'exportData', [], ...
                     'restartButton', [], ...
                     'saveMaskLabel', [], ...
                     'saveMask', [], ...
                     'saveShollLabel', [], ...
                     'saveSholl', [], ...
                     'saveIntensityHistLabel', [], ...
                     'saveIntensityHist', [], ...
                     'saveMatLabel', [], ...
                     'saveMat', [], ...
                     'saveSomaMeasureLabel', [], ...
                     'saveSomaMeasure', [], ...
                     'saveSynapseProfileLabel', [], ...
                     'saveSynapseProfile', [], ...
                     'saveGradientLabel', [], ...
                     'saveGradient', [], ...
                     'loadStage', [], ...
                     'somaStage', [], ...
                     'neuriteStage', [], ...
                     'synapseStage', [], ...
                     'analyseStage', [], ...
                     'allIcons', []);
    

    %%% Helper variables for showGUI
    
    loadGUI = [];
    somaGUI = [];
    neuriteGUI = [];
    synapseGUI = [];
    analyseGUI = [];
    allGUI = [];
    allIcons = [];
    allIconCallback = {};
    allAxes = [];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  end
  
  methods
    
    % Constructor

    function obj = SynD_extract2017()
						 
      format compact
      
      % For some reason if I have these as cell arrays to begin with
      % the struct becomes empty instead of having one element
      % initialised with default data
      data.synapsePixels = {}
      data.fileName = {}    
  
      % To make sure the config file is in the right place...
      obj.SynDscript = mfilename();
      obj.SynDscriptFull = mfilename('fullpath');
      obj.SyndScriptPath = [];
      
      % Check the licenses
      if(~license('checkout','matlab'))
        disp('MATLAB license missing.')
        obj.licenseOk = 0;
      end
      
      if(~license('checkout','image_toolbox'))
        disp('Image Toolbox license missing.')
        obj.licenseOk = 0;
      end
      
      if(~license('checkout','statistics_toolbox'))
        disp('Statistics Toolbox license missing.')
        obj.licenseOk = 0;
      end
      
      if(~obj.licenseOk)
        disp('SynD does not have access to all required functions')
        disp('Terminating now because the needed functions are not available.')
        return
      end
      
      obj.setSynDpath()
      
      obj.configFile = strcat(obj.SyndScriptPath,'SynD_config.mat');

      % Load saved config values from file
      obj.loadConfig();
            
      obj.setupGUI();

      
      
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function setupGUI(obj)
      
      obj.handles.fig = figure('Name','SynD - Synapse Detection', ...
                               'MenuBar','none', ...
                               'Toolbar','figure', ...
                               'Position', [50 50 1150 680], ...
                               'Resize','on', ...
                               'ResizeFcn', @obj.resizeFunction, ...
			       'Visible','off');

      obj.handles.image = axes('Units','Pixels', ...
                               'Position', [50 50 615 615]);
      
      obj.handles.redHist = axes('Units','Pixels', ...
                                 'Visible','off', ...
                                 'Position', [700 310 424 100]);
      
      obj.handles.greenHist = axes('Units','Pixels', ...
                                   'Visible','off', ...
                                   'Position', [700 180 424 100]);
      
      obj.handles.blueHist = axes('Units','Pixels', ...
                                  'Visible','off', ...
                                  'Position', [700 50 424 100]);
      
      obj.handles.credits = uicontrol('Style', 'text', ...
                                      'String', 'Johannes Hjorth and Sabine Schmitz, 2010', ...
                                      'HorizontalAlignment', 'right', ...
                                      'Foregroundcolor', 0.7*[1 1 1], ...
                                      'Backgroundcolor', get(gcf,'color'), ...
                                      'Position', [925 5 210 15], ...
                                      'Fontsize',8);
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % Load GUI
      %
      
      obj.handles.next = uicontrol('Style','pushbutton', ...
                                   'String','<html>Next</html>', ...
                                   'Interruptible','off', ...
                                   'Visible','off', ...
                                   'Position', [965 417 60 20], ...
                                   'Callback', @obj.nextImage);
      
      obj.handles.prev = uicontrol('Style','pushbutton', ...
                                   'String','<html>Prev</html>', ...
                                   'Interruptible','off', ...
                                   'Visible','off', ...
                                   'Position', [895 417 60 20], ... 
                                   'Callback', @obj.prevImage);
      
      obj.handles.num = uicontrol('Style','text', ...
                                  'String', 'No image', ...
                                  'Position', [1045 417 100 20], ... 
                                  'Visible','off', ...
                                  'BackgroundColor', get(gcf,'color'), ...
                                  'FontSize', 12, ...
                                  'HorizontalAlignment', 'left');
      
      obj.handles.loadNeuron = uicontrol('Style','pushbutton', ...
                                         'String','<html><u>L</u>oad image</html>', ...
                                     'Interruptible','off', ...
                                         'Position', [680 530 120 35], ...
                                         'Fontsize', 12, ...
                                         'BackgroundColor', [1 0.6 0], ...
                                         'Callback', @obj.loadNeuron);
  
      obj.handles.nextStage = uicontrol('Style','pushbutton', ...
                                        'String','<html><u>N</u>ext</html>', ...
                                        'TooltipString', ...
                                        ['<html>Continue to next stage<br>' ...
                          'of analysis when this<br>' ...
                          'one is done</html>'], ...
                                        'Interruptible','off', ...
                                        'Position', [820 530 120 35], ...
                                        'Fontsize', 12, ...
                                        'BackgroundColor', get(gcf,'color'), ...
                                        'Callback', @obj.nextStage);
      
      % Detect soma GUI
      %
      
      obj.handles.detectSoma = uicontrol('Style','pushbutton', ...
                                         'String','<html><u>D</u>etect soma</html>', ...
                                         'Interruptible','off', ...
                                         'Visible','off', ...
                                         'Position', [680 530 125 35], ...
                                         'Fontsize', 12, ...
                                         'BackgroundColor', [1 0.6 0], ...
                                         'TooltipString', ...
                                         ['<html>Left click and draw to add soma,<br>' ...
                          'right click and draw to remove soma</html>'], ...
                                         'Callback', @obj.detectSoma);
      
      obj.handles.editSoma = uicontrol('Style','pushbutton', ...
                                       'String','<html><u>E</u>dit soma</html>', ...
                                       'Interruptible','off', ...
                                       'Visible','off', ...
                                       'Position', [680 420 120 30], ...
                                       'Fontsize', 12, ...
                                       'Callback', @obj.editSoma);
      
      obj.handles.addSoma = uicontrol('Style','pushbutton', ...
                                      'String','<html><u>A</u>dd soma</html>', ...
                                      'Interruptible','off', ...
                                      'Visible','off', ...
                                      'Position', [820 420 120 30], ...
                                      'TooltipString', ...
                                      'Manually add ROI around soma.', ...
                                      'Fontsize', 12, ...
                                      'Callback', @obj.addSoma);
      
      obj.handles.deleteSoma = uicontrol('Style','pushbutton', ...
                                         'String','<html>Delete soma</html>', ...
                                         'Interruptible','off', ...
                                         'Visible','off', ...
                                         'Position', [960 420 120 30], ...
                                         'Fontsize', 12, ...
                                         'Callback', @obj.deleteSoma);
      
      
      obj.handles.addMeasurePoints = uicontrol('Style','pushbutton', ...
                                               'String','<html>New measu<u>r</u>e</html>', ...
                                               'TooltipString', ...
                                               ['<html>Add 10 random regions to the soma<br>', ...
                          'to measure protein levels.</html>'], ...
                                               'Interruptible','off', ...
                                               'Visible','off', ...
                                               'Position', [680 310 120 30], ...
                                               'Fontsize', 12, ...
                                               'Callback', @obj.setSomaMeasurePoints);
      
      obj.handles.editMeasurePoints = uicontrol('Style','pushbutton', ...
                                                'String','<html>Add p<u>o</u>int</html>', ...
                                                'Interruptible','off', ...
                                                'Visible','off', ...
                                                'TooltipString', ...
                                                '<html>Add regions in soma</html>', ...
                                                'Position', [820 310 120 30], ...
                                                'Fontsize', 12, ...
                                                'Callback', @obj.addSomaMeasurePoints);
      
      obj.handles.somaMeasureRadiusLabel = uicontrol('Style','text', ...
                                                     'String', 'Radius:', ...
                                                     'Visible','off', ...
                                                     'Position', [960 310 50 20], ...
                                                     'BackgroundColor', get(gcf,'color'), ...
                                                     'FontSize', 10, ...
                                                     'HorizontalAlignment', 'left');
      
      obj.handles.somaMeasureRadius = uicontrol('Style','edit', ...
                                                'String', ...
                                                num2str(obj.detection.measurePointRadius), ...
                                                'Visible','off', ...
                                                'Position', [1010 312 50 20], ... 
                                                'TooltipString', ...
                                                '<html>Radius of soma measure</html>', ...
                                                'Interruptible', 'off', ...
                                                'Callback', @obj.setSomaMeasureRadius);
      
      

      obj.handles.measureLabel = uicontrol('Style','text', ...
                                           'String', 'Measure protein levels in soma:', ...
                                           'Visible','off', ...
                                           'Position', [680 350 270 20], ... 
                                           'BackgroundColor', get(gcf,'color'), ...
                                           'FontSize', 12, ...
                                           'HorizontalAlignment', 'left');
      
      
      % Detect neurite GUI
      % 

      obj.handles.detectNeurite = uicontrol('Style','pushbutton', ...
                                            'String','<html><u>D</u>etect neurites</html>', ...
                                            'Interruptible','off', ...
                                            'Visible','off', ...
                                            'Position', [680 531 125 35], ...
                                            'Fontsize', 12, ...
                                            'BackgroundColor', [1 0.6 0], ...
                                            'Callback', @obj.detectNeurite);
      
      obj.handles.editNeurite = uicontrol('Style','pushbutton', ...
                                          'String','<html><u>E</u>dit neurite</html>', ...
                                          'Interruptible','off', ...
                                          'Visible','off', ...
                                          'Position', [680 420 120 30], ...
                                          'Fontsize', 12, ...
                                          'TooltipString', ...
                                          ['<html>Left click and draw adds neurite,<br>' ...
                          'right click and draw removes neurite</html>'], ...
                                          'Callback', @obj.editNeurite);

      obj.handles.addThinNeurites = uicontrol('Style','pushbutton', ...
                                              'String','<html>Add <u>t</u>hin</html>', ...
                                              'Interruptible','off', ...
                                              'Visible','off', ...
                                              'Position', [680 380 120 30], ...
                                              'Fontsize', 12, ...
                                              'TooltipString', ...
                                              ['<html>Extend detection of neurites, using<br>' ...
                          'steerable filter with half the size.</html>'], ...
                                              'Callback', @obj.addThinNeurites);
      

      obj.handles.addNeurite = uicontrol('Style','pushbutton', ...
                                         'String','<html><u>A</u>dd neurite</html>', ...
                                         'Interruptible','off', ...
                                         'Visible','off', ...
                                         'Position', [820 420 100 30], ...
                                         'Fontsize', 12, ...
                                         'Visible','off', ...
                                         'TooltipString', ...
                                         'Click on neurite to auto detect', ...
                                         'Callback', @obj.addNeurite);
      
      obj.handles.clean = uicontrol('Style','pushbutton', ...
                                    'String','<html><u>C</u>lean</html>', ...
                                    'Interruptible','off', ...
                                    'Visible','off', ...
                                    'Position', [940 420 100 30], ...
                                    'Fontsize', 12, ...
                                    'Visible','off', ...
                                    'TooltipString', ...
                                    'Removes all unconnected neurite parts', ...
                                    'Callback', @obj.cleanMask);
      
      obj.handles.killBlob = uicontrol('Style','pushbutton', ...
                                       'String','<html>Bl<u>o</u>b erase</html>', ...
                                       'Interruptible','off', ...
                                       'Visible','off', ...
                                       'Position', [940 380 100 30], ...
                                       'Fontsize', 12, ...
                                       'Visible','off', ...
                                       'TooltipString', ...
                                       'Removes one unconnected neurite part', ...
                                       'Callback', @obj.killBlob);
      
      obj.handles.extendNeurites = uicontrol('Style','pushbutton', ...
                                             'String','<html>Extend</html>', ...
                                             'Interruptible','off', ...
                                             'Visible','off', ...
                                             'Position', [820 380 100 30], ...
                                             'Fontsize', 12, ...
                                             'Visible','off', ...
                                             'TooltipString', ...
                                             'Looks for unconnected neurites', ...
                                             'Callback', @obj.extendNeuritesHandler);
      
      obj.handles.showSkeleton = uicontrol('Style','pushbutton', ...
                                           'String','<html><u>S</u>keleton</html>', ...
                                           'Interruptible','off', ...
                                           'Visible','off', ...
                                           'Position', [940 340 100 30], ...
                                           'Fontsize', 12, ...
                                           'Visible','off', ...
                                           'TooltipString', ...
                                           'Shows the neurite skeleton.', ...
                                           'Callback', @obj.showSkeleton);
      
      
      obj.handles.calculateLabel = uicontrol('Style', 'Text', ...
                                             'String', '', ...
                                             'Position', [680 320 250 30], ...
                                             'BackgroundColor', get(gcf,'color'), ...
                                             'Visible','off', ...
                                             'FontSize', 10, ...
                                             'HorizontalAlignment', 'left');
      
      
      obj.handles.skipSynapseDetection = uicontrol('Style','pushbutton', ...
                                                   'String', '<html>Skip synapses</html>', ...
                                                   'Interruptible','off', ...
                                                   'Visible','off', ...
                                                   'Position', [680 50 100 30], ...
                                                   'Fontsize', 7, ...
                                                   'ForegroundColor',0.3*[1 1 1], ...
                                                   'Visible','off', ...
                                                   'TooltipString', ...
                                                   [ '<html>Goes directly to export,<br>' ...
                          'if you do not need to<br>' ...
                          'detect synapses etc.</html>'], ...
                                                   'Callback', @obj.skipSynapseDetection);
      
      obj.handles.exportSimpleAverages = uicontrol('Style','pushbutton', ...
                                                   'String', '<html>Export averages</html>', ...
                                                   'Interruptible','off', ...
                                                   'Visible','off', ...
                                                   'Position', [800 50 100 30], ...
                                                   'Fontsize', 7, ...
                                                   'ForegroundColor',0.3*[1 1 1], ...
                                                   'Visible','off', ...
                                                   'TooltipString', ...
                                                   [ '<html>Export average intensity under<br>' ...
                          'soma and neurite masks.</html>'], ...
                                                   'Callback', @obj.exportOnlyAverages);
      
      
      % Synapse GUI
      %

      obj.handles.detectSynapses = uicontrol('Style','pushbutton', ...
                                             'String','<html><u>D</u>etect synapses</html>', ...
                                             'Interruptible','off', ...
                                             'Visible','off', ...
                                             'Position', [680 531 125 35], ...
                                             'Fontsize', 12, ...
                                             'BackgroundColor', [1 0.6 0], ...
                                             'Callback', @obj.detectSynapses);



      obj.handles.editSynapse = uicontrol('Style','pushbutton', ...
                                          'String','<html><u>E</u>dit synapses</html>', ...
                                          'Interruptible','off', ...
                                          'Visible','off', ...
                                          'Position', [680 420 120 30], ...
                                          'Fontsize', 12, ...
                                          'TooltipString', ...
                                          ['<html>Left click and draw adds synapse<br>' ...
                          'area, right click and draw removes<br>' ... 
                          'synapse area</html>'], ...
                                          'Callback', @obj.editSynapse);
      
      obj.handles.killSynBlob = uicontrol('Style','pushbutton', ...
                                          'String','<html>Bl<u>o</u>b erase</html>', ...
                                          'Interruptible','off', ...
                                          'Visible','off', ...
                                          'Position', [820 420 100 30], ...
                                          'Fontsize', 12, ...
                                          'Visible','off', ...
                                          'TooltipString', ...
                                          'Removes one unconnected neurite part', ...
                                          'Callback', @obj.killSynBlob);
      
      
      obj.handles.singleSynapse = axes('Units','Pixels', ...
                                       'visible','off', ...
                                       'Position', [720 100 400 300]);
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      obj.handles.morphLabel = uicontrol('Style','text', ...
                                         'String', 'Morph:', ...
                                         'Visible','off', ...
                                         'Position', [680 490 60 20], ... 
                                         'BackgroundColor', get(gcf,'color'), ...
                                         'FontSize', 12, ...
                                         'HorizontalAlignment', 'left');
      
      obj.handles.morphChannel = uicontrol('Style','popupmenu', ...
                                           'String', {'R','G','B'}, ...
                                           'Position', [740 495 60 15], ... 
                                           'FontSize', 10, ...
                                           'TooltipString', ...
                                           'Channel colour of morphology staining', ...
                                           'Value', obj.detection.remapChan(1), ...
                                           'Callback', @obj.setMorphChannel);
      
      obj.handles.morphChannelNumber = uicontrol('Style','popupmenu', ...
                                                 'String', {'1','2','3'}, ...
                                                 'Enable', obj.detection.enableRemapChan, ...
                                                 'Position', [740 465 60 15], ... 
                                                 'FontSize', 10, ...
                                                 'TooltipString', ...
                                                 'Channel with morphology staining', ...
                                                 'Value', ...
                                                 obj.detection.morphChannel, ...
                                                 'Callback', @obj.remapChannels);
      
      
      obj.handles.synLabel = uicontrol('Style','text', ...
                                       'String', 'Syn:', ...
                                       'Position', [810 490 40 20], ...
                                       'BackgroundColor', get(gcf,'color'), ...
                                       'FontSize', 12, ...
                                       'HorizontalAlignment', 'left');
      
      obj.handles.synChannel = uicontrol('Style','popupmenu', ...
                                         'String', {'R','G','B'}, ...
                                         'Position', [850 495 60 15], ...
                                         'FontSize', 10, ...
                                         'TooltipString', ...
                                         'Channel colour of synapse staining', ...
                                         'Value', obj.detection.remapChan(2), ...
                                         'Callback', @obj.setSynChannel);
      
      obj.handles.synChannelNumber = uicontrol('Style','popupmenu', ...
                                               'String', {'1','2','3'}, ...
                                               'Enable', obj.detection.enableRemapChan, ...
                                               'Position', [850 465 60 15], ...
                                               'FontSize', 10, ...
                                               'TooltipString', ...
                                               'Channel with synapse staining', ...
                                               'Value', ...
                                               obj.detection.synChannel, ...
                                               'Callback', @obj.remapChannels);
      
      obj.handles.XLabel = uicontrol('Style','text', ...
                                     'String', 'X:', ...
                                     'Position', [920 490 20 20], ... 
                                     'BackgroundColor', get(gcf,'color'), ...
                                     'FontSize', 12, ...
                                     'HorizontalAlignment', 'left');

      obj.handles.XChannel = uicontrol('Style','popupmenu', ...
                                       'String', {'R','G','B'}, ...
                                       'Position', [940 495 60 15], ... 
                                       'FontSize', 10, ...
                                       'TooltipString', ...
                                       'Channel to be analysed', ...
                                       'Value', obj.detection.remapChan(3), ...
                                       'Callback', @obj.setXChannel);
      
      obj.handles.XChannelNumber = uicontrol('Style','popupmenu', ...
                                             'String', {'1','2','3'}, ...
                                             'Enable', obj.detection.enableRemapChan, ...
                                             'Position', [940 465 60 15], ... 
                                             'FontSize', 10, ...
                                             'TooltipString', ...
                                             'Channel to be analysed', ...
                                             'Value', ...
                                             obj.detection.XChannel, ...
                                             'Callback', @obj.remapChannels);
      
      obj.handles.remapLabel = uicontrol('Style','text', ...
                                         'String', 'Remap:', ...
                                         'Visible','off', ...
                                         'Position', [1020 460 60 20], ... 
                                         'BackgroundColor', get(gcf,'color'), ...
                                         'FontSize', 12, ...
                                         'HorizontalAlignment', 'left');
      
      obj.handles.remap = uicontrol('Style', 'checkbox', ...
                                    'Position', [1080 460 20 20], ... 
                                    'BackgroundColor', get(gcf,'color'), ...
                                    'TooltipString', ...
                                    'Enable remapping of channels', ...
                                    'Value', strcmp(obj.detection.enableRemapChan,'on'), ...
                                    'Callback', @obj.toggleRemap);
      

      obj.handles.XYresLabel = uicontrol('Style','text', ...
                                         'String', 'XY res (micrometer) :', ...
                                         'Position', [680 415 133 20], ...
                                         'BackgroundColor', get(gcf,'color'), ...
                                         'FontSize', 10, ...
                                         'HorizontalAlignment', 'left');

      try
        resStr = num2str(obj.data.xyRes*1e6);
      catch e
        getReport(e)
        keyboard
      end
      
      obj.handles.XYres= uicontrol('Style','edit', ...
                                   'String', resStr, ...
                                   'Position', [813 417 60 20], ... 
                                   'TooltipString', ...
                                   '<html>Image resolution<br>in micrometers</html>', ...
                                   'Interruptible', 'off', ...
                                   'Callback', @obj.setXYres);
      


      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      obj.handles.red = uicontrol('Style', 'checkbox', ...
                                  'Position', [1004 530 20 20], ... 
                                  'BackgroundColor', get(gcf,'color'), ...
                                  'TooltipString', ...
                                  'Show red channel', ...
                                  'Value', 1, ...
                                  'Callback', @obj.toggleRed);
      
      obj.handles.green = uicontrol('Style', 'checkbox', ...
                                    'Position', [1054 530 20 20], ...
                                    'BackgroundColor', get(gcf,'color'), ...
                                    'TooltipString', ...
                                    'Show green channel', ...
                                    'Value', 1, ...
                                    'Callback', @obj.toggleGreen);
      
      obj.handles.blue = uicontrol('Style', 'checkbox', ...
                                   'Position', [1104 530 20 20], ... 
                                   'BackgroundColor', get(gcf,'color'), ...
                                   'TooltipString', ...
                                   'Show blue channel', ...
                                   'Value', 1, ...
                                   'Callback', @obj.toggleBlue);
      
      obj.handles.redLabel = uicontrol('Style','text', ...
                                       'String', 'R', ...
                                       'Position', [984 530 20 20], ... 
                                       'BackgroundColor', get(gcf,'color'), ...
                                       'FontSize', 12, ...
                                       'HorizontalAlignment', 'left');
      
      obj.handles.greenLabel = uicontrol('Style','text', ...
                                         'String', 'G', ...
                                         'Position', [1034 530 20 20], ... 
                                         'BackgroundColor', get(gcf,'color'), ...
                                         'FontSize', 12, ...
                                         'HorizontalAlignment', 'left');
      
      obj.handles.blueLabel = uicontrol('Style','text', ...
                                        'String', 'B', ...
                                        'Position', [1084 530 20 20], ... 
                                        'BackgroundColor', get(gcf,'color'), ...
                                        'FontSize', 12, ...
                                        'HorizontalAlignment', 'left');
      
      obj.handles.mask = uicontrol('Style', 'checkbox', ...
                                   'Position', [1104 500 20 20], ... 
                                   'Visible','off', ...
                                   'BackgroundColor', get(gcf,'color'), ...
                                   'TooltipString', ...
                                   'Show mask', ...
                                   'Value', 1, ...
                                   'Callback', @obj.toggleMask);
      
      obj.handles.maskLabel = uicontrol('Style','text', ...
                                        'String', 'M', ...
                                        'Visible','off', ...
                                        'Position', [1084 500 20 20], ... 
                                        'BackgroundColor', get(gcf,'color'), ...
                                        'FontSize', 12, ...
                                        'HorizontalAlignment', 'left');
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      
      obj.handles.morphThreshLabel = uicontrol('Style','text', ...
                                               'String', 'Morph thresh:', ...
                                               'Visible','off', ...
                                               'Position', [680 460 100 20], ...
                                               'BackgroundColor', get(gcf,'color'), ...
                                               'FontSize', 10, ...
                                               'HorizontalAlignment', 'left');
      
      obj.handles.morphThresh = uicontrol('Style','edit', ...
                                          'String', num2str(obj.detection.morphThreshold), ...
                                          'Visible','off', ...
                                          'Enable', obj.onOff(~obj.detection.autoSomaThreshold), ...
                                          'Position', [780 462 50 20], ... 
                                          'TooltipString', ...
                                          ['<html>Intensity threshold for morphology:<br>' ...
                          'Relative threshold (fraction of max intensity) if smaller than 1.<br>' ...
                          'Absolute threshold if larger than 1.</html>'], ...
                                          'Interruptible', 'off', ...
                                          'Callback', @obj.setMorphThresh);

      obj.handles.guessSomaThreshold = uicontrol('Style','pushbutton', ...
                                                 'String','<html>Auto</html>', ...
                                                 'Interruptible','off', ...
                                                 'Visible','off', ...
                                                 'Position', [850 462 50 20], ...
                                                 'TooltipString', ...
                                                 [ '<html>Find threshold by minimizing<br>' ...
                          'cross-entropy of image.</html>'], ...
                                                 'Fontsize', 12, ...
                                                 'BackgroundColor', get(gcf,'color'), ...
                                                 'Callback', @obj.autoSomaThreshold);

      obj.handles.somaErodeLabel = uicontrol('Style','text', ...
                                             'String', 'Soma erode:', ...
                                             'Visible','off', ...
                                             'Position', [680 490 100 20], ... 
                                             'BackgroundColor', get(gcf,'color'), ...
                                             'FontSize', 10, ...
                                             'HorizontalAlignment', 'left');
      
      obj.handles.somaErodeRadius = uicontrol('Style','edit', ...
                                              'String', num2str(obj.detection.somaErodeRadius), ...
                                              'TooltipString', ...
                                              ['<html>Higher value removes more<br>' ...
                          'of the thresholded image</html>'], ...
                                              'Visible','off', ...
                                              'Position', [780 492 50 20], ... 
                                              'Interruptible', 'off', ...
                                              'Callback', @obj.setSomaErodeRadius);
      
      obj.handles.singleSomaLabel = uicontrol('Style','text', ...
                                              'String', 'Single soma:', ...
                                              'Visible','off', ...
                                              'Position', [850 490 100 20], ... 
                                              'BackgroundColor', get(gcf,'color'), ...
                                              'FontSize', 10, ...
                                              'HorizontalAlignment', 'left');

      obj.handles.singleSoma = uicontrol('Style', 'checkbox', ...
                                         'Visible','off', ...
                                         'Position', [940 492 20 20], ... 
                                         'BackgroundColor', get(gcf,'color'), ...
                                         'TooltipString', ...
                                         'Show red channel', ...
                                         'Value', obj.detection.singleSoma, ...
                                         'Callback', @obj.toggleSingleSoma);


      obj.handles.synThreshLabel = uicontrol('Style','text', ...
                                             'String', 'Syn thresh:', ...
                                             'Visible','off', ...
                                             'Position', [680 490 80 20], ... 
                                             'BackgroundColor', get(gcf,'color'), ...
                                             'FontSize', 10, ...
                                             'HorizontalAlignment', 'left');

      obj.handles.synThresh = uicontrol('Style','edit', ...
                                        'String', num2str(obj.detection.synThreshold), ...
                                        'Visible','off', ...
                                        'Position', [760 492 50 20], ... 
                                        'TooltipString', ...
                                        ['<html>Number of standard deviations<br>' ...
                          'above noise level required to<br>' ...
                          'be considered a synapse.<br>' ...
                          'If NaN uses cross-entropy for <br>' ...
                          'threshold obj.detection.<br>' ...
                          'If negative, use absolute value as<br>' ...
                          'absolute intensity threshold.</html>'], ...
                                        'Interruptible', 'off', ...
                                        'Callback', @obj.setSynThresh);
      
      obj.handles.neuritePaddingLabel = uicontrol('Style','text', ...
                                                  'String', 'Neurite padding:', ...
                                                  'Visible','off', ...
                                                  'Position', [830 490 115 20], ... 
                                                  'BackgroundColor', get(gcf,'color'), ...
                                                  'FontSize', 10, ...
                                                  'HorizontalAlignment', 'left');
      
      obj.handles.neuritePadding = uicontrol('Style','edit', ...
                                             'String', num2str(obj.detection.neuritePadding*1e6), ...
                                             'Visible','off', ...
                                             'Position', [945 492 50 20], ... 
                                             'TooltipString', ...
                                             ['<html>How many micrometers outside a neurite<br>' ...
                          'synapses are allowed to extend</html>'], ...
                                             'Interruptible', 'off', ...
                                             'Callback', @obj.setNeuritePadding);
      
      obj.handles.minSynapseSizeLabel = uicontrol('Style','text', ...
                                                  'String', 'Min size:', ...
                                                  'Visible','off', ...
                                                  'Position', [830 460 115 20], ... 
                                                  'BackgroundColor', get(gcf,'color'), ...
                                                  'FontSize', 10, ...
                                                  'HorizontalAlignment', 'left');
      
      obj.handles.minSynapseSize = uicontrol('Style','edit', ...
                                             'String', num2str(obj.detection.minSynapseSize*1e12), ...
                                             'Visible','off', ...
                                             'Position', [945 462 50 20], ... 
                                             'TooltipString', ...
                                             ['<html>Minimum size in micrometers square<br>' ...
                          'for detected synapses.</html>'], ...
                                             'Interruptible', 'off', ...
                                             'Callback', @obj.setSynapseMinSize);
      
      
      obj.handles.growThreshLabel = uicontrol('Style','text', ...
                                              'String', 'Max cost:', ...
                                              'Visible','off', ...
                                              'Position', [680 490 80 20], ...
                                              'BackgroundColor', get(gcf,'color'), ...
                                              'FontSize', 10, ...
                                              'HorizontalAlignment', 'left');
      
      obj.handles.growThresh = uicontrol('Style','edit', ...
                                         'String', num2str(obj.detection.maxAddCost), ...
                                         'Position', [760 492 50 20], ...
                                         'Visible','off', ...
                                         'TooltipString', ...
                                         ['<html>Value between 0 and 1,<br>', ...
                          'higher value means more<br>', ...
                          'neurite is included</html>'], ...
                                         'Interruptible', 'off', ...
                                         'Callback', @obj.setGrowThresh);
      
      obj.handles.filterSizeLabel = uicontrol('Style','text', ...
                                              'String', 'Filter size:', ...
                                              'Visible','off', ...
                                              'Position', [680 460 80 20], ...
                                              'BackgroundColor', get(gcf,'color'), ...
                                              'FontSize', 10, ...
                                              'HorizontalAlignment', 'left');
      
      obj.handles.filterSize = uicontrol('Style','edit', ...
                                         'String', num2str(obj.detection.filterSize*1e6), ...
                                         'Position', [760 462 50 20], ...
                                         'Visible','off', ...
                                         'TooltipString', ...
                                         ['<html>Filter size, width in micrometers.<br>', ...
                          'Lower values detect thinner neurites.</html>'], ...
                                         'Interruptible', 'off', ...
                                         'Callback', @obj.setFilterSize);


      % Export buttons
      %

      obj.handles.exportData = uicontrol('Style','pushbutton', ...
                                         'String','<html><u>E</u>xport data</html>', ...
                                         'Interruptible','off', ...
                                         'Visible','off', ...
                                         'Position', [680 530 120 35], ...
                                         'Fontsize', 12, ...
                                         'BackgroundColor', [0.8 0.8 0.8], ...
                                         'Callback', @obj.exportData);
      
      obj.handles.restartButton = uicontrol('Style','pushbutton', ...
                                            'String','<html>Restart SynD</html>', ...
                                            'Interruptible','off', ...
                                            'Visible','off', ...
                                            'Position', [950 530 120 35], ...
                                            'Fontsize', 12, ...
                                            'BackgroundColor', [0.8 0.0 0.0], ...
                                            'Callback', @obj.restartSynD);
      
      
      obj.handles.saveMaskLabel = uicontrol('Style','text', ...
                                            'String', 'Save mask', ...
                                            'Visible','off', ...
                                            'Position', [680 480 120 20], ... 
                                            'BackgroundColor', get(gcf,'color'), ...
                                            'FontSize', 12, ...
                                            'HorizontalAlignment', 'left');

      obj.handles.saveMask = uicontrol('Style', 'checkbox', ...
                                       'Visible','off', ...
                                       'Position', [810 480 20 20], ... 
                                       'BackgroundColor', get(gcf,'color'), ...
                                       'Value', obj.exportInfo.saveMask);
      
      obj.handles.saveShollLabel = uicontrol('Style','text', ...
                                             'String', 'Save Sholl', ...
                                             'Visible','off', ...
                                             'Position', [680 450 120 20], ... 
                                             'BackgroundColor', get(gcf,'color'), ...
                                             'FontSize', 12, ...
                                             'HorizontalAlignment', 'left');
      
      obj.handles.saveSholl = uicontrol('Style', 'checkbox', ...
                                        'Visible','off', ...
                                        'Position', [810 450 20 20], ... 
                                        'BackgroundColor', get(gcf,'color'), ...
                                        'Value', obj.exportInfo.saveSholl);
      
      obj.handles.saveIntensityHistLabel = uicontrol('Style','text', ...
                                                     'String', 'Save histogram', ...
                                                     'Visible','off', ...
                                                     'Position', [680 420 120 20], ... 
                                                     'BackgroundColor', get(gcf,'color'), ...
                                                     'FontSize', 12, ...
                                                     'HorizontalAlignment', 'left');
      
      obj.handles.saveIntensityHist = uicontrol('Style', 'checkbox', ...
                                                'Visible','off', ...
                                                'Position', [810 420 20 20], ... 
                                                'BackgroundColor', get(gcf,'color'), ...
                                                'Value', obj.exportInfo.saveIntensityHistogram);
      
      obj.handles.saveMatLabel = uicontrol('Style','text', ...
                                           'String', 'Save mat-file', ...
                                           'Visible','off', ...
                                           'Position', [680 390 120 20], ... 
                                           'BackgroundColor', get(gcf,'color'), ...
                                           'FontSize', 12, ...
                                           'HorizontalAlignment', 'left');
      
      obj.handles.saveMat = uicontrol('Style', 'checkbox', ...
                                      'Visible','off', ...
                                      'Position', [810 390 20 20], ... 
                                      'BackgroundColor', get(gcf,'color'), ...
                                      'Value', obj.exportInfo.saveMat);

      obj.handles.saveSomaMeasureLabel = uicontrol('Style','text', ...
                                                   'String', 'Save soma measure', ...
                                                   'Visible','off', ...
                                                   'Position', [680 360 120 20], ... 
                                                   'BackgroundColor', get(gcf,'color'), ...
                                                   'FontSize', 12, ...
                                                   'HorizontalAlignment', 'left');
      
      obj.handles.saveSomaMeasure = uicontrol('Style', 'checkbox', ...
                                              'Visible','off', ...
                                              'Position', [810 360 20 20], ... 
                                              'BackgroundColor', get(gcf,'color'), ...
                                              'Value', obj.exportInfo.saveSomaMeasure);

      obj.handles.saveSynapseProfileLabel = uicontrol('Style','text', ...
                                                      'String', 'Save synapse profile', ...
                                                      'Visible','off', ...
                                                      'Position', [680 330 120 20], ... 
                                                      'BackgroundColor', get(gcf,'color'), ...
                                                      'FontSize', 12, ...
                                                      'HorizontalAlignment', 'left');
      

      obj.handles.saveSynapseProfile = uicontrol('Style', 'checkbox', ...
                                                 'Visible','off', ...
                                                 'Position', [810 330 20 20], ... 
                                                 'BackgroundColor', get(gcf,'color'), ...
                                                 'Value', obj.exportInfo.saveSynapseProfile);
      
      
      obj.handles.saveGradientLabel = uicontrol('Style','text', ...
                                                'String', 'Save neurite gradient', ...
                                                'Visible','off', ...
                                                'Position', [680 300 120 20], ... 
                                                'BackgroundColor', get(gcf,'color'), ...
                                                'FontSize', 12, ...
                                                'HorizontalAlignment', 'left');
      
      
      obj.handles.saveGradient = uicontrol('Style', 'checkbox', ...
                                           'Visible','off', ...
                                           'Position', [810 300 20 20], ... 
                                           'BackgroundColor', get(gcf,'color'), ...
                                           'Value', obj.exportInfo.saveGradient);
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

      %%% Icons for switching "stages", ie load, soma, neurites, synapses, analyse
      % Swedish lesson: bilder = pictures  

      
      % Load the images
      
      
      try
        obj.bilder.icon{1} = imread('bilder/icon_load.png');
        obj.bilder.icon{2} = imread('bilder/icon_soma.png');
        obj.bilder.icon{3} = imread('bilder/icon_neurite.png');
        obj.bilder.icon{4} = imread('bilder/icon_synapse.png');
        obj.bilder.icon{5} = imread('bilder/icon_export.png');
      catch
        close all
        disp('Icons not found, aborting.')
        disp('Please download from www.johanneshjorth.se')
        return
      end

      obj.handles.loadStage = uicontrol('Style','pushbutton', ...
                                        'Cdata', obj.bilder.icon{1}, ...
                                        'TooltipString', ...
                                        'Load image with neuron', ...
                                        'Interruptible','off', ...
                                        'Position', [680 581 84 84], ...
                                        'Callback', @obj.loadStage);

      obj.handles.somaStage = uicontrol('Style','pushbutton', ...
                                        'Cdata', obj.desaturateIcon(obj.bilder.icon{2}), ...
                                        'TooltipString', ...
                                        '<html>Detect soma and<br> edit soma mask</html>', ...
                                        'Interruptible','off', ...
                                        'Position', [770 581 84 84]);, ...


          obj.handles.neuriteStage = uicontrol('Style','pushbutton', ...
                                               'Cdata', obj.desaturateIcon(obj.bilder.icon{3}), ...
                                               'TooltipString', ...
                                               ['<html>Detect neurites and<br>' ...
                          'edit neurite mask</html>'], ...
                                               'Interruptible','off', ...
                                               'Position', [860 581 84 84]);

      obj.handles.synapseStage = uicontrol('Style','pushbutton', ...
                                           'Cdata', obj.desaturateIcon(obj.bilder.icon{4}), ...
                                           'TooltipString', ...
                                           ['<html>Detect synapses and<br>' ...
                          'edit synapse mask</html>'], ...
                                           'Interruptible','off', ...
                                           'Position', [950 581 84 84]);

      obj.handles.analyseStage = uicontrol('Style','pushbutton', ...
                                           'Cdata', obj.desaturateIcon(obj.bilder.icon{5}), ...
                                           'TooltipString', ...
                                           'Export data to xml', ...
                                           'Interruptible','off', ...
                                           'Position', [1040 581 84 84]);


      % Helper variables for GUI
      obj.handles.allIcons = [obj.handles.loadStage, ...
                          obj.handles.somaStage, ...
                          obj.handles.neuriteStage, ...
                          obj.handles.synapseStage, ...
                          obj.handles.analyseStage];



      obj.loadGUI = [ obj.handles.loadNeuron, ...
                      obj.handles.prev, obj.handles.next, ...
                      obj.handles.num, ...
                      obj.handles.next, obj.handles.prev, ...
                      obj.handles.morphLabel, ...
                      obj.handles.morphChannel, ...
                      obj.handles.morphChannelNumber, ...
                      obj.handles.synLabel, ...
                      obj.handles.synChannel, ...
                      obj.handles.synChannelNumber, ...
                      obj.handles.remapLabel, obj.handles.remap, ...
                      obj.handles.XLabel, ...
                      obj.handles.XChannel, ...
                      obj.handles.XChannelNumber, ...
                      obj.handles.redLabel, obj.handles.red, ...
                      obj.handles.greenLabel, obj.handles.green, ...
                      obj.handles.blueLabel, obj.handles.blue, ...
                      obj.handles.XYresLabel, obj.handles.XYres, ...
                      obj.handles.nextStage, ...
                      obj.handles.redHist, obj.handles.greenHist, obj.handles.blueHist, ...
                    ];
      
      obj.somaGUI = [ obj.handles.detectSoma, obj.handles.editSoma, ...
                      obj.handles.addSoma, obj.handles.deleteSoma, ...
                      obj.handles.addMeasurePoints, obj.handles.editMeasurePoints, ...
                      obj.handles.measureLabel, ...
                      obj.handles.morphThreshLabel, obj.handles.morphThresh, ...
                      obj.handles.guessSomaThreshold, ...
                      obj.handles.singleSomaLabel, obj.handles.singleSoma, ...
                      obj.handles.somaErodeLabel, obj.handles.somaErodeRadius, ...
                      obj.handles.somaMeasureRadiusLabel, ...
                      obj.handles.somaMeasureRadius, ...
                      obj.handles.redLabel, obj.handles.red, ...
                      obj.handles.greenLabel, obj.handles.green, ...
                      obj.handles.blueLabel, obj.handles.blue, ...
                      obj.handles.maskLabel, obj.handles.mask, ...
                      obj.handles.nextStage, ...
                    ];

      obj.neuriteGUI = [ obj.handles.detectNeurite, obj.handles.editNeurite, ...
                         obj.handles.addThinNeurites, ...
                         obj.handles.morphThreshLabel, obj.handles.morphThresh, ...
                         obj.handles.clean, obj.handles.killBlob, ... 
                         obj.handles.addNeurite, obj.handles.extendNeurites, ...
                         obj.handles.showSkeleton, ...
                         obj.handles.growThreshLabel, obj.handles.growThresh, ...
                         obj.handles.filterSizeLabel, obj.handles.filterSize, ...
                         obj.handles.redLabel, obj.handles.red, ...
                         obj.handles.greenLabel, obj.handles.green, ...
                         obj.handles.blueLabel, obj.handles.blue, ...
                         obj.handles.maskLabel, obj.handles.mask, ...
                         obj.handles.nextStage, ...
                         obj.handles.calculateLabel, ...
                         obj.handles.skipSynapseDetection, ...
                         obj.handles.exportSimpleAverages, ...
                       ];
      
      obj.synapseGUI = [ obj.handles.detectSynapses, obj.handles.editSynapse...
                         obj.handles.killSynBlob, ...
                         obj.handles.synThreshLabel, obj.handles.synThresh, ...
                         obj.handles.neuritePaddingLabel, obj.handles.neuritePadding, ...
                         obj.handles.minSynapseSizeLabel, obj.handles.minSynapseSize, ...
                         obj.handles.redLabel, obj.handles.red, ...
                         obj.handles.greenLabel, obj.handles.green, ...
                         obj.handles.blueLabel, obj.handles.blue, ...
                         obj.handles.maskLabel, obj.handles.mask, ...
                         obj.handles.nextStage, ...
                         obj.handles.singleSynapse, ...
                       ];
      
      obj.analyseGUI = [ obj.handles.exportData, ...
                         obj.handles.restartButton, ...
                         obj.handles.saveMaskLabel, obj.handles.saveMask, ...
                         obj.handles.saveShollLabel, obj.handles.saveSholl, ...
                         obj.handles.saveIntensityHistLabel, ...
                         obj.handles.saveIntensityHist, ...
                         obj.handles.saveMatLabel, obj.handles.saveMat, ...
                         obj.handles.saveSomaMeasureLabel, obj.handles.saveSomaMeasure, ...
                         obj.handles.saveSynapseProfileLabel, ...
                         obj.handles.saveSynapseProfile, ...
                         obj.handles.saveGradientLabel, ...
                         obj.handles.saveGradient, ...
                       ];
      
      
      obj.allGUI = union(union(union(obj.loadGUI, obj.somaGUI), obj.neuriteGUI), ...
                         union(obj.synapseGUI, obj.analyseGUI));
      
      obj.allIcons = [obj.handles.loadStage obj.handles.somaStage ...
                      obj.handles.neuriteStage obj.handles.synapseStage ...
                      obj.handles.analyseStage];
      
      obj.allIconCallback = {@obj.loadStage @obj.somaStage ...
                          @obj.neuriteStage @obj.synapseStage ...
                          @obj.analyseStage};
      
      obj.allAxes = [obj.handles.image, obj.handles.singleSynapse, ... 
                     obj.handles.redHist, obj.handles.greenHist, obj.handles.blueHist];
      
      
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      % Make the figure resizeable by using normalisation
      set([obj.allGUI obj.allIcons obj.handles.fig obj.allAxes, ...
           obj.handles.credits], ...
          'Units','Normalized')
      
      % Set the figure position to what was previously read in from config file
      if(~isempty(obj.dispInfo.figPosition))
        set(obj.handles.fig,'position', obj.dispInfo.figPosition);
      end

      % Add callback function for keypresses
      set(obj.handles.fig,'WindowKeyPressFcn',@obj.captureKeyPress);

      obj.showGUI('load');
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      obj.handles.menuFile = uimenu(obj.handles.fig,'Label','File');
      obj.handles.menuItemLoad = uimenu(obj.handles.menuFile, ...
                                        'Label','Load neuron', ...
                                        'Interruptible','off', ...
                                        'Callback', @obj.loadNeuron);

      obj.handles.menuItemBatchAnalyse = uimenu(obj.handles.menuFile, ...
                                         'Label','Batch analyse neurons', ...
                                         'Interruptible','off', ...
                                         'Callback', @obj.batchAnalyse);
      
      obj.handles.menuItemBatchReanalyse = uimenu(obj.handles.menuFile, ...
                                         'Label','Batch re-analyse neurons', ...
                                         'Interruptible','off', ...
                                         'Callback', @obj.batchReAnalyse);
      
      obj.handles.menuItemImportMask = uimenu(obj.handles.menuFile, ...
                                              'Label','Import old masks', ...
                                              'Interruptible', 'off', ...
                                              'Callback', @obj.importMask);
      
      obj.handles.menuItemExport = uimenu(obj.handles.menuFile, ...
                                          'Label','Export XML', ...
                                          'Interruptible', 'off', ...
                                          'Callback', @obj.exportData);
      
      obj.handles.menuItemRestart = uimenu(obj.handles.menuFile, ...
                                           'Label','Restart', ...
                                           'Interruptible', 'off', ...
                                           'Callback', @obj.restartSynD);
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      obj.handles.menuEdit = uimenu(obj.handles.fig,'Label','Edit');
      
      obj.handles.menuItemUndo = uimenu(obj.handles.menuEdit, ...
                                        'Label', '<html><u>U</u>ndo: No history</html>', ...
                                        'Interruptible', 'off', ...
                                        'Callback', @obj.undoLastAction);
      
      obj.handles.menuItemExclude = uimenu(obj.handles.menuEdit, ...
                                           'Label', 'Exclude region', ...
                                           'Interruptible', 'off', ...
                                           'Callback', @obj.excludeRegion);
      
      obj.handles.menuItemClearExclude = uimenu(obj.handles.menuEdit, ...
                                                'Label', 'Reset exclude region', ...
                                                'Interruptible', 'off', ...
                                                'Callback', @obj.clearExcludeRegions);
      
      obj.handles.menuItemMaxProjection = uimenu(obj.handles.menuEdit, ...
                                                 'Label', 'Max projection', ...
                                                 'Interruptible', 'off', ...
                                                 'Callback', @obj.maxProjection);
      
      %  obj.handles.menuItemCollapse = uimenu(obj.handles.menuEdit, ...
      %				   'Label', 'Collapse image stack', ...
      %				   'Interruptible', 'off', ...
      %				   'Callback', @obj.collapseImageStack);

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      obj.handles.menuSettings = uimenu(obj.handles.fig,'Label','Settings');

      obj.handles.autoThreshold = uimenu(obj.handles.menuSettings, ...
                                         'Label','<html>Auto soma threshold</html>', ...
                                         'Interruptible', 'off', ...
                                         'Checked', obj.onOff(obj.detection.autoSomaThreshold), ...
                                         'Callback', @obj.toggleAutoSomaThreshold);

      obj.handles.menuItemSteerable = uimenu(obj.handles.menuSettings, ...
                                             'Label','<html>Use steerable <u>f</u>ilters</html>', ...
                                             'Interruptible', 'off', ...
                                             'Checked','on', ...
                                             'Callback', @obj.toggleSteerableFilters);
      
      obj.handles.menuItemBackgroundFiltering = ...
          uimenu(obj.handles.menuSettings, ...
                 'Label','Synapse background removal', ...
                 'Interruptible', 'off', ...
                 'Checked','on', ...
                 'Callback', @obj.toggleBackgroundRemoval);
      
      obj.handles.menuItemSteerableSettings = ...
          uimenu(obj.handles.menuSettings, ...
                 'Label','Steerable filter settings', ...
                 'Interruptible', 'off', ...
                 'Callback', @obj.steerableFilterSettings);

      obj.handles.menuItemSynapseSettings = ...
          uimenu(obj.handles.menuSettings, ...
                 'Label','Synapse detection settings', ...
                 'Interruptible', 'off', ...
                 'Callback', @obj.synapseSettings);
      
      obj.handles.menuItemShollSettings = ...
          uimenu(obj.handles.menuSettings, ...
                 'Label','Sholl settings', ...
                 'Interruptible', 'off', ...
                 'Callback', @obj.shollSettings);
      
      obj.handles.menuItemSavePath = ...
          uimenu(obj.handles.menuSettings, ...
                 'Label','Save SynD path permanently', ...
                 'Interruptible', 'off', ...
                 'Callback', @obj.saveSynDpath);
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      obj.handles.menuHelp = uimenu(obj.handles.fig,'Label','Help');
      
      webHelp = 'web(''http://www.johanneshjorth.se/SynD/'',''-browser'')';

      obj.handles.menuItemHelp = uimenu(obj.handles.menuHelp, ...
                                        'Label','How to use SynD', ...
                                        'Interruptible', 'off', ...
                                        'Callback', webHelp);

      obj.handles.menuItemVersion = uimenu(obj.handles.menuHelp, ...
                                           'Label','Get latest version', ...
                                           'Interruptible', 'off', ...
                                           'Callback', @obj.getLatestVersion); 

      obj.handles.menuItemAbout = uimenu(obj.handles.menuHelp, ...
                                         'Label', 'About SynD', ...
                                         'Interruptible', 'off', ...
                                         'Callback', @obj.aboutSynD);			     


      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      if(obj.enableDebug)

        obj.handles.menuDebug = uimenu(obj.handles.fig,'Label','Debug');
        obj.handles.menuItemDebug = uimenu(obj.handles.menuDebug, ...
                                           'Label','Keyboard', ...
                                           'Interruptible', 'off', ...
                                           'Callback', @obj.runDebug);
        
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      obj.handles.measurePointMenu = uicontextmenu;
      
      uimenu(obj.handles.measurePointMenu, ...
             'Label', 'Move measure point',...
             'Interruptible', 'off', ...
             'Callback',@obj.moveMeasureCallback);
      
      uimenu(obj.handles.measurePointMenu, ...
             'Label', 'Delete measure point',...
             'Interruptible', 'off', ...
             'Callback',@obj.deleteMeasureCallback);
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % Some dirty tricks to get the zoom and pan tools
      % Set up toolbar, only keep zoom in, zoom out and pan
      obj.handles.toolbar = findall(obj.handles.fig,'Type','uitoolbar');
      oldChild = allchild(obj.handles.toolbar);
      
      for i=1:length(oldChild)
        tmp = get(oldChild(i),'Tag');
        
        switch(tmp)
          case 'Exploration.ZoomIn';   
          case 'Exploration.ZoomOut';
          case 'Exploration.Pan';
            % Do nothing, we want to keep these
          otherwise 
            delete(oldChild(i));         % Remove the rest
        end
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % Display SynD logo
      obj.showLogo();

      % Display the graphics  					       
      set(obj.handles.fig,'visible','on')
      
    end
    













 



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function aboutSynD(obj,source, event)

    SynDversion = obj.getVersion();

    fprintf('Running version %s of SynD.\n', SynDversion);

    if~(isempty(obj.bilder.logo))
      miniLogo = imresize(obj.bilder.logo, [1 1]*50);
    else
      miniLogo = [];
    end

    aboutMsg = sprintf(['You are running SynD version %s\n\n' ...
                        'Johannes Hjorth\n' ...
                        'Sabine Schmitz\n\n\n' ...
                        'Reference:\n%s\n\n' ...
                        'DOI: %s'], ...
                       SynDversion, obj.referenceStr, obj.referenceDOI);
    
    uiwait(msgbox(aboutMsg,'About SynD','Custom',miniLogo));
    
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function SynDversion = getVersion(obj)

    try
     fid = fopen(strcat(obj.SyndscriptPath,obj.verFile),'r');
      SynDversion = fgets(fid);
      if(SynDversion(end) == char(10))
        SynDversion = SynDversion(1:end-1);
      end
      fclose(fid);
    catch
      SynDversion = '(unknown)';
      disp('No version information available.')
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function newVersion = askLatestVersion(obj)

    try
      newVersion = urlread(strcat(obj.urlBase,obj.verFile));

      if(newVersion(end) == char(10))
        newVersion = newVersion(1:end-1);
      end
    catch
      disp('Unable to check latest version')
      newVersion = '(unavailable)';
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function getLatestVersion(obj,source, event)

    try
      curVersion = obj.getVersion();
      newVersion = obj.askLatestVersion();

      fprintf('Current version: %s\nLatest version: %s\n', ...
	      curVersion, newVersion);

      if(strcmp(newVersion,curVersion))
        uiwait(helpdlg(sprintf(['SynD version %s ' ...
                                'is the latest version.\n'],curVersion), ...
		       'SynD update'));
      else
        newFile = sprintf('SynD-revision-%s.tar.bz2',newVersion);

        [outFile, outPath] = uiputfile('SynD-revision*.tar.bz2', ...
				       'Save new version of SynD', ...
				       newFile);

        urlwrite(strcat(obj.urlBase,newFile),strcat(outPath,outFile));

        uiwait(helpdlg(sprintf(['Downloaded %s, ' ...
				'please extract SynD from the file.'], ...
			       outFile), ...
		       'SynD update'));

      end

    catch
      disp('Unable to update automatically. Please see the INCF SynD page.');
      web('http://software.incf.org/software/synd/','-browser');
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function runDebug(obj,source, event)
    disp('Type return to exit debug mode')
    keyboard
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function s = restartSynD(obj,source, event)
      userInput = questdlg('Do you want to restart SynD?', ...
			   'Restart imminent', ...
			   'Yes','No','Maybe','Yes');

      switch(userInput)
        case 'Yes'
          disp('Restarting SynD.')
          close(obj.handles.fig)
          s = SynD_extract2017();
        case 'No'
          disp('Restart aborted.')
          s = obj;
        case 'Maybe'
          if(rand(1) < 0.5)
            disp('Restarting SynD.')
            close(obj.handles.fig)
            s = SynD_extract2017();
          else
            disp('Not restarting SynD.')
            s = obj;
          end
      end
    end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function toggleRed(obj, source, event)
    % If this function was called without a source, then toggle red
    % otherwise read in status from the clicked object
    if(~exist('source') | source ~= obj.handles.red)
      set(obj.handles.red,'Value',~get(obj.handles.red,'Value'));
    end

    obj.dispInfo.showRed = get(obj.handles.red,'Value');
  
    obj.showImage();
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function toggleGreen(obj, source, event)
    if(~exist('source') | source ~= obj.handles.green)
      set(obj.handles.green,'Value',~get(obj.handles.green,'Value'));
    end

    obj.dispInfo.showGreen = get(obj.handles.green,'Value');
    obj.showImage();
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function toggleBlue(obj, source, event)
    if(~exist('source') | source ~= obj.handles.blue)
      set(obj.handles.blue,'Value',~get(obj.handles.blue,'Value'));
    end

    obj.dispInfo.showBlue = get(obj.handles.blue,'Value');
    obj.showImage();
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function toggleMask(obj, source, event)
    if(~exist('source') | source ~= obj.handles.mask)
      set(obj.handles.mask,'Value',~get(obj.handles.mask,'Value'));
    end

    obj.dispInfo.showMask = get(obj.handles.mask,'Value');
  
    obj.showImage();
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setRGBvisible(obj,redFlag,greenFlag,blueFlag)
    obj.dispInfo.showRed = redFlag;
    set(obj.handles.red,'Value',redFlag);
    obj.dispInfo.showGreen = greenFlag;
    set(obj.handles.green,'Value',greenFlag);
    obj.dispInfo.showBlue = blueFlag;
    set(obj.handles.blue,'Value',blueFlag);

    obj.showImage();
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setChannelVisible(obj,morphFlag,synFlag,XFlag)

    visibleFlags = [0 0 0];

    if(morphFlag)
      visibleFlags(obj.detection.remapChan(1)) = 1;
    end

    if(synFlag)
      visibleFlags(obj.detection.remapChan(2)) = 1;
    end

    if(XFlag)
      visibleFlags(obj.detection.remapChan(3)) = 1;
    end

    obj.setRGBvisible(visibleFlags(1),visibleFlags(2),visibleFlags(3));

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setXYres(obj,source, event)
    tmp = 1e-6*str2num(get(obj.handles.XYres,'String'));

    if(~isempty(tmp) & 0 < tmp & tmp < Inf)
      obj.data.xyRes = tmp;
      if(isempty(obj.data.image))
        uicontrol(obj.handles.loadNeuron);
      else
        uicontrol(obj.handles.nextStage);
      end
    else
      uiwait(warndlg(sprintf('Incorrect resolution %s, using %s micrometers instead,', ...
			     get(obj.handles.XYres,'String'), ...
			     num2str(1e6*obj.data.xyRes)), ...
		     'Input error','modal'));

      set(obj.handles.XYres, 'String', num2str(obj.data.xyRes*1e6))

    end

    fprintf('XY-resolution set to %d m\n', obj.data.xyRes)

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setMorphChannel(obj, source, event)
    obj.detection.remapChan(1) = get(obj.handles.morphChannel,'Value');

    if(strcmp(obj.detection.enableRemapChan,'off'))
      % Use default: R maps to 1, G maps to 2, B maps to 3
      set(obj.handles.morphChannelNumber,'Value',obj.detection.remapChan(1));
      obj.detection.morphChannel = obj.detection.remapChan(1);

    end

    % Tell the user to redo everything...
    obj.dispInfo.stage = 2;

    obj.calcMaxIntensity();

    % Reset scaling
    obj.dispInfo.scaleRed = 1;
    obj.dispInfo.scaleBlue = 1;
    obj.dispInfo.scaleGreen = 1;

    obj.showGUI(obj.dispInfo.state);
    obj.showImage();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setSynChannel(obj, source, event)
    obj.detection.remapChan(2) = get(obj.handles.synChannel,'Value');
  
    if(strcmp(obj.detection.enableRemapChan,'off'))
      % Use default: R maps to 1, G maps to 2, B maps to 3
      set(obj.handles.synChannelNumber,'Value',obj.detection.remapChan(2));
      obj.detection.synChannel = obj.detection.remapChan(2);

    end

    % Tell the user to redo everything...
    obj.dispInfo.stage = 2;

    obj.calcMaxIntensity();

    % Reset scaling
    obj.dispInfo.scaleRed = 1;
    obj.dispInfo.scaleBlue = 1;
    obj.dispInfo.scaleGreen = 1;

    obj.showGUI(obj.dispInfo.state);
    obj.showImage();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setXChannel(obj, source, event)
    
    obj.detection.remapChan(3) = get(obj.handles.XChannel,'Value');

    if(strcmp(obj.detection.enableRemapChan,'off'))
      % Use default: R maps to 1, G maps to 2, B maps to 3
      set(obj.handles.XChannelNumber,'Value',obj.detection.remapChan(3));
      obj.detection.XChannel = obj.detection.remapChan(3);
    end

    % Tell the user to redo everything...
    obj.dispInfo.stage = 2;

    obj.calcMaxIntensity();

    % Reset scaling
    obj.dispInfo.scaleRed = 1;
    obj.dispInfo.scaleBlue = 1;
    obj.dispInfo.scaleGreen = 1;

    obj.showGUI(obj.dispInfo.state);
    obj.showImage();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function remapChannels(obj, source, event)

    obj.detection.morphChannel = get(obj.handles.morphChannelNumber,'Value');
    obj.detection.synChannel = get(obj.handles.synChannelNumber,'Value');
    obj.detection.XChannel = get(obj.handles.XChannelNumber,'Value');

    % Tell the user to redo everything...
    obj.dispInfo.stage = 2;

    obj.calcMaxIntensity();

    % Reset scaling
    obj.dispInfo.scaleRed = 1;
    obj.dispInfo.scaleBlue = 1;
    obj.dispInfo.scaleGreen = 1;

    obj.showGUI(obj.dispInfo.state);
    obj.showImage();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function toggleRemap(obj,source, event)
    if(get(obj.handles.remap,'Value'))
      obj.detection.enableRemapChan = 'on';

    else
      obj.detection.enableRemapChan = 'off';

      set(obj.handles.morphChannelNumber,'Value',obj.detection.remapChan(1));
      obj.detection.morphChannel = obj.detection.remapChan(1);

      set(obj.handles.synChannelNumber,'Value',obj.detection.remapChan(2));
      obj.detection.synChannel = obj.detection.remapChan(2);

      set(obj.handles.XChannelNumber,'Value',obj.detection.remapChan(3));
      obj.detection.XChannel = obj.detection.remapChan(3);

      obj.showImage();

      obj.calcMaxIntensity();

      % Reset scaling
      obj.dispInfo.scaleRed = 1;
      obj.dispInfo.scaleBlue = 1;
      obj.dispInfo.scaleGreen = 1;

      % Tell the user to redo everything...
      obj.dispInfo.stage = 2;
      obj.showGUI(obj.dispInfo.state);

    end

    set(obj.handles.morphChannelNumber,'enable',obj.detection.enableRemapChan);
    set(obj.handles.synChannelNumber,'enable',obj.detection.enableRemapChan);
    set(obj.handles.XChannelNumber,'enable',obj.detection.enableRemapChan);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function toggleSingleSoma(obj, source, event)
    % If this function was called without a source, then toggle red
    % otherwise read in status from the clicked object
    if(~exist('source') | source ~= obj.handles.singleSoma)
      set(obj.handles.singleSoma,'Value',~get(obj.handles.singleSoma,'Value'));
    end

    obj.detection.singleSoma = get(obj.handles.singleSoma,'Value');

    obj.detectSoma();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setMorphThresh(obj, source, event)
    tmp = str2num(get(obj.handles.morphThresh,'String'));

    if(~isempty(tmp) & 0 < tmp & tmp < Inf)
      obj.detection.morphThreshold = tmp;
      disp(sprintf('Using morphology threshold %d', obj.detection.morphThreshold))
    
      if(exist('source') & source == obj.handles.morphThresh)
        obj.detectSoma();
      end
    else
      uiwait(warndlg(sprintf('Incorrect threshold %s, using %s instead.', ...
			     get(obj.handles.morphThresh,'String'), ...
			     num2str(obj.detection.morphThreshold)), ...
		     'Input error','modal'));

      set(obj.handles.morphThresh,'String', num2str(obj.detection.morphThreshold));
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setSomaErodeRadius(obj, source, event)
    tmp = str2num(get(obj.handles.somaErodeRadius,'String'));

    if(~isempty(tmp) & 0 < tmp & tmp < Inf)
      obj.detection.somaErodeRadius = tmp;
      disp(sprintf('Using soma erode radius %d', obj.detection.somaErodeRadius))

      % Only call detect soma if the user pressed enter
      % (because detectSoma also calls this function to make sure it
      %  has latest user input).
      if(exist('source') & source == obj.handles.somaErodeRadius)
        obj.detectSoma();
      end

    else
      uiwait(warndlg(sprintf('Incorrect radius %s, using %s instead.', ...
			     get(obj.handles.somaErodeRadius,'String'), ...
			     num2str(obj.detection.somaErodeRadius)), ...
		     'Input error','modal'));

      set(obj.handles.somaErodeRadius,'String', num2str(obj.detection.somaErodeRadius));
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setSomaMeasureRadius(obj, source, event)

    tmp = str2num(get(obj.handles.somaMeasureRadius,'String'));
    if(~isempty(tmp) & 0 < tmp & tmp < Inf)
      obj.detection.measurePointRadius = tmp;
    else
      uiwait(warndlg(sprintf('Incorrect radius %s, using %d.', ...
			     get(obj.handles.somaMeasureRadius,'String'), ...
                             obj.detection.measurePointRadius)));
      set(obj.handles.somaMeasureRadius,'String', ...
                        num2str(obj.detection.measurePointRadius));
    end

    obj.makeMeasureMask();

    obj.showImage();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setSynThresh(obj, source, event)
    tmp = str2num(get(obj.handles.synThresh,'String'));

    [obj.detection.synThreshold, modFlag] = ...
        obj.sanitiseInput(get(obj.handles.synThresh,'String'), ...
                          obj.detection.synThreshold, ...
                          -inf, inf, false, true);

    if(modFlag)
      uiwait(warndlg(sprintf('Incorrect threshold %s, using %s instead.', ...
                             get(obj.handles.synThresh,'String'), ...
                             num2str(obj.detection.synThreshold)), ...
                     'Input error','modal'));

      set(obj.handles.synThresh, 'String', num2str(obj.detection.synThreshold));
      
    else
      if(exist('source') & source == obj.handles.synThresh)
        if(isnan(obj.detection.synThreshold))
          disp('Using cross-entropy for synapse threshold')
        else
          fprintf('Using synapse threshold %d\n', obj.detection.synThreshold)
        end

        obj.detectSynapses();

      end

    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setNeuritePadding(obj, source, event)
    tmp = str2num(get(obj.handles.neuritePadding, 'String'))*1e-6;

    if(~isempty(tmp) & 0 <= tmp & tmp < Inf)
      obj.detection.neuritePadding = tmp;
      disp(sprintf('Using %d pixels neurite padding', ...
		   round(obj.detection.neuritePadding/obj.data.xyRes)))

      if(exist('source') & source == obj.handles.neuritePadding)
        obj.detectSynapses();
      end
    else
      uiwait(warndlg(sprintf('Incorrect neurite padding %s, using %s instead.', ...
			     get(obj.handles.neuritePadding,'String'), ...
			     num2str(obj.detection.neuritePadding*1e6)), ...
		     'Input error','modal'));

      set(obj.handles.neuritePadding, 'String', num2str(obj.detection.neuritePadding*1e6));

    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setSynapseMinSize(obj, source, event)

    tmp = str2num(get(obj.handles.minSynapseSize, 'String'))*1e-12;

    if(~isempty(tmp) & 0 <= tmp & tmp < Inf)

      obj.detection.minSynapseSize = tmp;

      if(exist('source') & source == obj.handles.minSynapseSize)
        detectSynapses();
      end

    else
      set(obj.handles.minSynapseSize,'String', ...
                        num2str(obj.detection.minSynapseSize*1e12));
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setGrowThresh(obj, source, event)
    tmp = str2num(get(obj.handles.growThresh,'String'));

    if(~isempty(tmp) & 0 < tmp & tmp <= 1)
      obj.detection.maxAddCost = tmp;
      disp(sprintf('Using maximal grow cost of %.4f', obj.detection.maxAddCost))

      if(exist('source') & source == obj.handles.growThresh)
        obj.detectNeurite();
      end
    else

      uiwait(warndlg(sprintf('Incorrect cost %s, using %s instead.', ...
			     get(obj.handles.growThresh,'String'), ...
			     num2str(obj.detection.maxAddCost)), ...
		     'Input error','modal'));

      set(obj.handles.growThresh, 'String', num2str(obj.detection.maxAddCost));
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setFilterSize(obj, source, event)

    tmp = str2num(get(obj.handles.filterSize,'String'));

    if(~isempty(tmp) & 0.1 < tmp & tmp <= 20)
      obj.detection.filterSize = tmp*1e-6;
      disp(sprintf('Using filter size of %.1f micrometer', ...
		   obj.detection.filterSize*1e6))

      if(exist('source') & source == obj.handles.filterSize)
        obj.detectNeurite();
      end
    else

      uiwait(warndlg(sprintf('Incorrect filter sizet %s, using %s micrometer instead.', ...
			     get(obj.handles.filterSize,'String'), ...
			     num2str(obj.detection.filterSize*1e6)), ...
		     'Input error','modal'));

      set(obj.handles.filterSize, 'String', num2str(obj.detection.filterSize*1e6));
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function toggleSteerableFilters(obj, source, event)
    if(obj.detection.useSteerableFilters)
      obj.detection.useSteerableFilters = 0;
      set(obj.handles.menuItemSteerable,'Checked','off')
    else
      obj.detection.useSteerableFilters = 1;
      set(obj.handles.menuItemSteerable,'Checked','on')
    end

    obj.showGUI(obj.dispInfo.state);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function toggleAutoSomaThreshold(obj, source, event)

    if(obj.detection.autoSomaThreshold)
      obj.detection.autoSomaThreshold = 0;
    else
      obj.detection.autoSomaThreshold = 1;
    end

    set(obj.handles.autoThreshold, 'Checked', obj.onOff(obj.detection.autoSomaThreshold))

    set(obj.handles.morphThresh, 'Enable', obj.onOff(~obj.detection.autoSomaThreshold));

    obj.showGUI(obj.dispInfo.state);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function toggleBackgroundRemoval(obj, source, event)
    
    if(obj.detection.backgroundFiltering)
      obj.detection.backgroundFiltering = false;
      set(obj.handles.menuItemBackgroundFiltering,'Checked','off');
    else
      obj.detection.backgroundFiltering = true;
      set(obj.handles.menuItemBackgroundFiltering,'Checked','on');
    end

    obj.showGUI(obj.dispInfo.state);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function status = onOff(obj, flag)
    if(flag)
      status = 'on';
    else
      status = 'off';
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function editSoma(obj, source, event)

    obj.editInfo.mode = 2;
    obj.editInfo.defaultWidth = 8;

    obj.stopEdit(); 
    set(obj.handles.fig,'windowbuttondownfcn', @obj.startDrawLine)

    obj.dispInfo.somaColor = [1 1 1];
    obj.dispInfo.neuriteColor = NaN;
    obj.dispInfo.synapseColor = NaN;
    obj.dispInfo.measurePointColor = NaN;
    obj.showImage();

    % The user is responsible for drawing a soma if they click this button
    obj.activateStage('neurites');

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function addSoma(obj, source, event)

    set(obj.handles.fig,'CurrentAxes',obj.handles.image)    
    obj.saveUndo('add soma');

    obj.stopEdit(); % Clears mouse handlers

    set(obj.somaGUI,'enable','off')
    set(obj.handles.allIcons,'enable','off')
    mask = roipoly();
    set(obj.somaGUI,'enable','on')
    set(obj.handles.allIcons,'enable','on')

    if(~isempty(mask))
      obj.data.somaMask = obj.data.somaMask | mask;
      obj.activateStage('neurites');
    end

    obj.showImage();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function deleteSoma(obj, source, event)
    
    disp('Select a soma to delete.')
    [x,y,button] = ginput(1);

    a = axis;

    while(a(1) <= x & x <= a(2) ...
          & a(3) <= y & y <= a(4) ...
          & 1 <= x & x <= obj.data.width ...
          & 1 <= y & y <= obj.data.height ...
          & 1 == button)

      obj.saveUndo('delete soma');

      obj.data.somaMask = obj.cleanBlob(obj.data.somaMask,round(x),round(y));

      obj.data.somaMeasureMask = [];

      % Update image
      obj.showImage();

      % Check if we should remove one more
      disp('Select soma to remove.')
      [x,y,button] = ginput(1);
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function mask = cleanBlob(obj,mask,x,y)

    idx = sub2ind(size(mask),y,x);

    CC = bwconncomp(mask);

    for i = 1:length(CC.PixelIdxList)
      % Remove the blob that the user clicked on
      if(ismember(idx,CC.PixelIdxList{i}))
        mask(CC.PixelIdxList{i}) = 0;
      end
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setSomaMeasurePoints(obj,source, event)

    obj.stopEdit();

    if(isempty(obj.data.somaMask))
      disp('No soma detected yet.')
      return
    end

    measureDisk = strel('disk',obj.detection.measurePointRadius);
    measureDiskPadded = strel('disk',obj.detection.measurePointRadius*2);

    obj.saveUndo('new soma measure');

    measurePoint = [];

    [yS,xS] = ind2sub(size(obj.data.somaMask),find(obj.data.somaMask));
    somaRangeX = min(xS)-1:max(xS)+1;
    somaRangeY = min(yS)-1:max(yS)+1;

    somaRangeX = min(max(1,somaRangeX),size(obj.data.somaMask,2));
    somaRangeY = min(max(1,somaRangeY),size(obj.data.somaMask,1));

    % Zoom in on the soma
    obj.dispInfo.axis = [min(xS)-1, max(xS)+1, min(yS)-1, max(yS)+1];
    axis(obj.dispInfo.axis);

    miniSomaMask = obj.data.somaMask(somaRangeY,somaRangeX);
    possiblePoints = imerode(miniSomaMask,measureDisk);
    possiblePointIdx = find(possiblePoints);  

    removeMask = zeros(size(possiblePoints));

    warnedUser = 0;

    for i = 1:obj.detection.nMeasurePoints
      if(~isempty(possiblePointIdx))
        newPoint = possiblePointIdx(1+floor(length(possiblePointIdx)*rand(1)));
        measurePoint = [measurePoint, newPoint];    

        removeMask(:) = 0;
        removeMask(measurePoint) = 1;
        removeMask = imdilate(removeMask,measureDiskPadded);
        possiblePoints = possiblePoints-removeMask > 0;
        possiblePointIdx = find(possiblePoints);
      elseif(~warnedUser)
        uiwait(warndlg(sprintf('Unable to add %d points to soma', ...
			       obj.detection.nMeasurePoints-i+1), ...  
		       'SynD: Not enough space in soma'));
        warnedUser = 1;
        break;
      end

    end

    % Remap back to big soma mask
    [yM,xM] = ind2sub(size(miniSomaMask),measurePoint);
    yM = yM + min(somaRangeY)-1;
    xM = xM + min(somaRangeX)-1;

    obj.data.measurePoint = sub2ind(size(obj.data.somaMask),yM,xM);

    obj.makeMeasureMask();

    % Display the mask

    obj.dispInfo.somaColor = NaN;
    obj.dispInfo.measurePointColor = obj.dispInfo.defaultMeasurePointColor;
    obj.dispInfo.neuriteColor = NaN;
    obj.dispInfo.synapseColor = NaN;

    obj.showImage();

    obj.exportInfo.saveSomaMeasure = 1;
    set(obj.handles.saveSomaMeasure,'Value',obj.exportInfo.saveSomaMeasure);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function makeMeasureMask(obj)

    measureDisk = strel('disk',obj.detection.measurePointRadius);

    obj.data.somaMeasureMask = zeros(size(obj.data.somaMask));
    obj.data.somaMeasureMask(obj.data.measurePoint) = 1;
    obj.data.somaMeasureMask = imdilate(obj.data.somaMeasureMask,measureDisk);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function moveMeasureCallback(obj, source, event)

    pointIdx = obj.findClosestMeasurePoint();

    obj.saveUndo('move soma measure');
    
    disp('Move soma measure.')
    [x,y,idx] = obj.getImagePoint();

    if(isempty(idx))
      % Nothing to do, abort
      return
    end

    obj.data.measurePoint(pointIdx) = idx;
    obj.makeMeasureMask();

    obj.dispInfo.somaColor = NaN;
    obj.dispInfo.measurePointColor = obj.dispInfo.defaultMeasurePointColor;
    obj.dispInfo.neuriteColor = NaN;
    obj.dispInfo.synapseColor = NaN;

    obj.showImage();
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function deleteMeasureCallback(obj, source, event)

    obj.saveUndo('delete soma measure');

    pointIdx = obj.findClosestMeasurePoint();
    obj.data.measurePoint(pointIdx) = [];
    obj.makeMeasureMask();

    obj.dispInfo.somaColor = NaN;
    obj.dispInfo.measurePointColor = obj.dispInfo.defaultMeasurePointColor;
    obj.dispInfo.neuriteColor = NaN;
    obj.dispInfo.synapseColor = NaN;

    obj.showImage();
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [x,y,idx] = getImagePoint(obj)
    [x,y] = ginput(1);
    x = round(x);
    y = round(y);

    idx = sub2ind(size(obj.data.somaMask),y,x);

    a = axis;
    if(x < a(1) | a(2) < x | y < a(3) | a(4) < y)
      % We are outside the image
      x = []; y = []; idx = [];
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function pointIdx = findClosestMeasurePoint(obj)

    tmpXY = get(obj.handles.image,'CurrentPoint');
    x = round(tmpXY(1,1));
    y = round(tmpXY(1,2));

    [yP,xP] = ind2sub(size(obj.data.somaMask),obj.data.measurePoint);
    Pdist = sqrt((xP-x).^2+(yP-y).^2);

    [minValue, pointIdx] = min(Pdist);
    % [~, pointIdx] = min(Pdist); % Only works in 2009b and later

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function addSomaMeasurePoints(obj, source, event)

    if(isempty(obj.data.somaMask))
      disp('No soma mask set yet.')
      return
    end

    obj.saveUndo('add soma measure');

    % Zoom in on the soma

    [yS,xS] = ind2sub(size(obj.data.somaMask),find(obj.data.somaMask));

    obj.dispInfo.axis = [min(xS)-1, max(xS)+1, min(yS)-1, max(yS)+1];
    axis(obj.dispInfo.axis);

    obj.dispInfo.somaColor = NaN;
    obj.dispInfo.measurePointColor = obj.dispInfo.defaultMeasurePointColor;
    obj.dispInfo.neuriteColor = NaN;
    obj.dispInfo.synapseColor = NaN;

    obj.showImage();

    % Ask user for new start point
    disp('Add some measure')
    [x,y,idx] = obj.getImagePoint();

    if(isempty(idx))
      % Nothing to do, abort
      return
    end

    obj.data.measurePoint(end+1) = idx;
    obj.makeMeasureMask();

    obj.showImage();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function editNeurite(obj, source, event)

    obj.editInfo.mode = 1;
    obj.editInfo.defaultWidth = 3;

    obj.stopEdit();
    set(obj.handles.fig,'windowbuttondownfcn', @obj.startDrawLine)

    obj.dispInfo.somaColor = [1 1 1];
    obj.dispInfo.neuriteColor = [1 1 1]*0.7;
    obj.dispInfo.synapseColor = NaN;
    obj.dispInfo.measurePointColor = NaN;

    obj.showImage();

    % The user is responsible for drawing neurites if they click this button
    obj.activateStage('synapses');

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function editSynapse(obj, source, event)

    obj.editInfo.mode = 3;
    obj.editInfo.defaultWidth = 2;

    obj.stopEdit();
    set(obj.handles.fig,'windowbuttondownfcn', @obj.startDrawLine)

    obj.dispInfo.somaColor = NaN;
    obj.dispInfo.neuriteColor = NaN;
    obj.dispInfo.synapseColor = [1 1 1];
    obj.dispInfo.measurePointColor = NaN;
    obj.showImage();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function stopEdit(obj)

    % obj.editInfo.mode = 0;
    set(obj.handles.fig,'windowbuttondownfcn', []);
    set(obj.handles.fig,'windowbuttonmotionfcn', []);
    set(obj.handles.fig,'windowbuttonupfcn', []);

    obj.showImage();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function killSynBlob(obj, source, event)

    disp('Select blob to remove.')
    [x,y,button] = ginput(1);

    a = axis;

    while(a(1) <= x & x <= a(2) ...
	  & a(3) <= y & y <= a(4) ...
	  & 1 <= x & x <= obj.data.width ...
	  & 1 <= y & y <= obj.data.height ...
	  & 1 == button)

      obj.saveUndo('synapse blob removal');

      obj.data.synapseMask = obj.cleanBlob(obj.data.synapseMask,round(x),round(y));

      % Update image
      obj.showImage();

      disp('Select blob to remove.')
      % Check if we should remove one more
      [x,y,button] = ginput(1);
    end

    obj.locateSynapseCenters();
    obj.filterSynapseMask(obj.data.synapseMask);

    % Update synapse properties
    obj.calculateSynapseProperties();
    obj.showImage();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function steerableFilterSettings(obj, source, event)

    prompt = { 'Sigma (in micrometers)', 'Max growth cost' };
    defaultVal = { num2str(obj.detection.filterSize*1e6), ...
                   num2str(obj.detection.maxAddCost) };
    dialogName = 'Steerable filter settings';
    numLines = 1;

    answers = inputdlg(prompt, dialogName, numLines, defaultVal);

    if(~isempty(answers))

      [obj.detection.filterSize, modFlag1] = ...
          obj.sanitiseInput(answers{1}, obj.detection.filterSize*1e6, ...
                            0.1, 20, false);

      % Everything internally is stored in SI units.
      obj.detection.filterSize = obj.detection.filterSize * 1e-6;

      [obj.detection.maxAddCost, modFlag2] = ...
          obj.sanitiseInput(answers{2}, obj.detection.maxAddCost, ...
                            0, 1, false);

      if(modFlag1 | modFlag2)
        warnMsg = sprintf(['Sigma must 1-20, max growth cost 0-1.']);
		   
        uiwait(warndlg(warnMsg, 'Input sanitation','modal'));
   
      end

      % Update GUI
      set(obj.handles.growThresh,'String', num2str(obj.detection.maxAddCost));
      set(obj.handles.filterSize,'String', num2str(obj.detection.filterSize*1e6));

      % Update the filters
      if(~isempty(obj.data.image))
        obj.steerableFilter();
      end

    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function synapseSettings(obj, source, event)

    prompt = { 'Threshold (in std, NaN for auto)', ...
               'Min synapse size (um^2)', ...
               'Radius (in um, NaN for auto)', ...
               'Neurite padding (um)', ...
               'Include soma synapses (default 0=false)' };
    defaultVal = { num2str(obj.detection.synThreshold), ...
                   num2str(obj.detection.minSynapseSize*1e12), ...
                   num2str(obj.detection.singleSynapseRadie*1e6), ...
                   num2str(obj.detection.neuritePadding*1e6), ...
                   num2str(~obj.detection.excludeSomaSynapses) };
    dialogName = 'Synapse detection settings';
    numLines = 1;

    answers = inputdlg(prompt, dialogName, numLines, defaultVal);

    if(~isempty(answers))

      % Allows NaN = auto synapse threshold
      [obj.detection.synThreshold, modFlag1] = ...
          obj.sanitiseInput(answers{1}, obj.detection.synThreshold, ...
                        -inf, inf, false, true);
      
      [tmpSize, modFlag2] = ...
          obj.sanitiseInput(answers{2}, obj.detection.minSynapseSize*1e12, ...
                        0, inf, false, false);
      obj.detection.minSynapseSize = tmpSize / 1e12;
      
      % Allows NaN = auto synapse radie
      [tmpRadie, modFlag3] = ...
          obj.sanitiseInput(answers{3}, obj.detection.singleSynapseRadie*1e6, ...
                        0, inf, false, true);
      obj.detection.singleSynapseRadie = tmpRadie*1e-6;
      
      [obj.detection.neuritePadding, modFlag4] = ...
          obj.sanitiseInput(answers{4}, obj.detection.neuritePadding * 1e6, ...
		      0, 20, false, false);

      obj.detection.neuritePadding = obj.detection.neuritePadding * 1e-6;

      % Changed it to include soma synapses in settings, more intuitive
      [includeSomaSynapses,modFlag5] = ...
          obj.sanitiseInput(answers{5}, ~obj.detection.excludeSomaSynapses, ...
		      0, 1, true, false);
      obj.detection.excludeSomaSynapses = ~includeSomaSynapses;

      if(modFlag1 | modFlag2 | modFlag3 | modFlag4 | modFlag5)
        warnMsg = sprintf(['Synapse threshold must be 1 or larger '...
			   '(measured in std of noise, if negative absolute value is used as intensity threshold). ' ...
			   'Min synapse size is in micrometer^2. ' ...
                           'Single synapse radius (micrometers), ' ...
                           'if set to NaN tries to estimate it. ' ...
                           'Neurite padding must be between 0 and 20 pixels.'...
                           ' 1 = include, 0 = exclude soma synapses.' ...
			   ]);
	   
        uiwait(warndlg(warnMsg, 'Input sanitation','modal'));
   
      else
        % Update GUI
        set(obj.handles.synThresh,'String', ...
	    num2str(obj.detection.synThreshold));
        set(obj.handles.neuritePadding, 'String', ...
	    num2str(obj.detection.neuritePadding*1e6));
        set(obj.handles.minSynapseSize,'String', ...
	    num2str(obj.detection.minSynapseSize*1e12));
 
      end

      if(nnz(obj.data.neuriteMask))
        obj.detectSynapses();
      end

    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function shollSettings(obj, source, event)

    prompt = { 'Bin size (in micrometers)' };
    defaultVal = { num2str(obj.detection.shollBinSize*1e6) };

    dialogName = 'Sholl analysis settings';
    numLines = 1;

    answers = inputdlg(prompt, dialogName, numLines, defaultVal);

    if(~isempty(answers))

      [shollBinSize, modFlag1] = ...
          obj.sanitiseInput(answers{1}, obj.detection.shollBinSize*1e6, ...
		      1, 1000, false);

      obj.detection.shollBinSize = shollBinSize*1e-6;

      if(modFlag1)
        warnMsg = 'Please specify a bin size between 1 and 1000.';
        uiwait(warndlg(warnMsg, 'Input sanitation','modal'));
      end

    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Copied from EvA_extractMorphology.m
  % (also available from johanneshjorth.se)
  %
  % This function is used to sanitise the input, we require that
  % it is a number between minVal and maxVal, and we can also
  % require it to be an integer if we want.

  function [newVal,modFlag] = sanitiseInput(obj, inputString, oldVal, ...
					    minVal, maxVal, ...
					    integerFlag, allowNANFlag)

    if(~exist('allowNANFlag'))
      allowNANFlag = false;
    end

    readVal = str2num(inputString);
    modFlag = 0;

    if(isempty(readVal))
      % We end up here if user types in letters instead of numbers
      newVal = oldVal;
      modFlag = 1;     
    else
      if(integerFlag)
        newVal = round(readVal);
      else
        newVal = readVal;
      end     
 
      if(allowNANFlag & isnan(newVal))
        % We got NaN, and we allow it
      else
        newVal = max(minVal, newVal);
        newVal = min(maxVal, newVal);

        % Mark that we changed the input when sanitising it
        if(readVal ~= newVal)
          modFlag = 1;
        end

      end

    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function startDrawLine(obj, source, event)

    set(obj.handles.fig,'CurrentAxes',obj.handles.image)    

    % Where did the user click, within the image?
    xy = get(obj.handles.image,'CurrentPoint');
    x = xy(1,1);
    y = xy(1,2);

    a = axis();    

    if(x < a(1) | a(2) < x ...
       | y < a(3) | a(4) < y)

     % We are outside of the picture, abort.
     disp('Clicked outside of image')
     return;

    end 

    % What button was pressed (left draw, right erase)
    switch(get(obj.handles.fig,'selectiontype'))
      case 'normal'
        obj.editInfo.mode = abs(obj.editInfo.mode); % Drawing
        obj.editInfo.color = [1 1 1]*0.8;
        obj.editInfo.width = obj.editInfo.defaultWidth;
      case 'alt'
        obj.editInfo.mode = -abs(obj.editInfo.mode); % Erasing
        obj.editInfo.color = [1 1 1]*0.2;
        obj.editInfo.width = obj.editInfo.defaultWidth*2;
      otherwise
        disp('Unknown key pressed')  
        disp(get(obj.handles.fig,'selectiontype'))
        return
    end

    obj.editInfo.line = line(x,y, ...
                             'color',obj.editInfo.color, ...
                             'linewidth', obj.editInfo.width);
    set(obj.handles.fig,'windowbuttonmotionfcn',@obj.drawLineCallback)
    set(obj.handles.fig,'windowbuttonupfcn', @obj.endDrawLine)
    
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function drawLineCallback(obj, source, event)

    xy = get(obj.handles.image,'CurrentPoint');
    x = xy(1,1);
    y = xy(1,2);

    a = axis();    

    if(x < a(1) | a(2) < x ...
       | y < a(3) | a(4) < y)
      % We are outside figure, end drawing
      obj.endDrawLine();
     
    end

    try
      xOld = get(obj.editInfo.line,'Xdata');
      yOld = get(obj.editInfo.line,'Ydata');
    catch
      disp('Unable to read points from mouse.')
      obj.endDrawLine();
      return
    end

    if(1 <= x & x <= obj.data.width ...
       & 1 <= y & y <= obj.data.height)
      obj.editInfo.xLine = [xOld,x];
      obj.editInfo.yLine = [yOld,y];

      set(obj.editInfo.line, ...
          'Xdata', obj.editInfo.xLine, ...
          'Ydata', obj.editInfo.yLine);
    else
      % We were outside image, this should not happen
      % but put it here just to be on safe side.
      obj.endDrawLine();
    end

    drawnow

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function endDrawLine(obj, source, event)

    disp('End draw line called')

    % Remove drawing callback (but keep button down function)
    set(obj.handles.fig,'windowbuttonmotionfcn', [])
    set(obj.handles.fig,'windowbuttonupfcn', [])

    if(isempty(obj.editInfo.xLine))
      % Nothing to draw
      return;
    end


    % Modify the mask
    tmpMask = zeros(size(obj.data.neuriteMask));
   
    tmpMask(round(obj.editInfo.yLine(1)),round(obj.editInfo.xLine(1))) = 1;

    % There should be a built in function to connect points to lines
    % and return a binaries mask. This is a bad way to do a line, it will 
    % add extra points. Better to use Bresenhams algorithm...

    for i = 2:length(obj.editInfo.xLine)
      nDots = 2*max(abs(obj.editInfo.xLine(i)-obj.editInfo.xLine(i-1)), ...
		    abs(obj.editInfo.yLine(i)-obj.editInfo.yLine(i-1)));

      xDots = round(linspace(obj.editInfo.xLine(i-1),obj.editInfo.xLine(i),nDots));
      yDots = round(linspace(obj.editInfo.yLine(i-1),obj.editInfo.yLine(i),nDots));
      iDots = sub2ind(size(tmpMask),yDots,xDots);

      tmpMask(iDots) = 1;
    end

    % Use appropriate width
    a = axis;
    zoomLevel = max((a(2)-a(1))/obj.data.width,(a(4)-a(3))/obj.data.height);
    brushSize = max(round(obj.editInfo.width*zoomLevel),1);
    tmpMask = imdilate(tmpMask,strel('disk',brushSize));

    switch(obj.editInfo.mode)
      case 1
        % Add the pixels to neurite mask
        obj.saveUndo('add neurite');
        obj.data.neuriteMask = double((obj.data.neuriteMask + tmpMask) > 0);
        obj.updateNeuriteMask();
      case -1
        % Remove the pixels from neurite mask
        obj.saveUndo('remove neurite');
        obj.data.neuriteMask = double(((obj.data.neuriteMask>0) - tmpMask) > 0);
        obj.updateNeuriteMask();
      case 2
        % Add the pixels to soma mask
        obj.saveUndo('add soma');
        obj.data.somaMask = double((obj.data.somaMask + tmpMask) > 0);
      case -2
        % Remove the pixels from soma mask
        obj.saveUndo('remove soma');
        obj.data.somaMask = double((obj.data.somaMask - tmpMask) > 0);
      case 3
        % Add the pixels to synapse mask
        obj.saveUndo('add synapse');
        obj.data.synapseMask = double((obj.data.synapseMask + tmpMask) > 0);

        obj.locateSynapseCenters();
        obj.filterSynapseMask(obj.data.synapseMask);

        % Update synapse properties
        obj.calculateSynapseProperties();

      case -3
        % Remove the pixels from synapse mask
        obj.saveUndo('remove synapse');
        obj.data.synapseMask = double((obj.data.synapseMask - tmpMask) > 0);

        obj.locateSynapseCenters();
        obj.filterSynapseMask(obj.data.synapseMask);

        % Update synapse properties
        obj.calculateSynapseProperties();

      otherwise
        % Do nothing
    end


    % Display the result
    obj.data.distMask = obj.makeDistMask(obj.data.somaMask,obj.data.neuriteMask);
    obj.showImage();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function autoSomaThreshold(obj,source, event)

    % Find threshold that minimizes cross-entropy
    grayImage = obj.data.image(:,:,obj.detection.morphChannel, obj.dispInfo.curImg);
    obj.detection.morphThreshold  = obj.minimizeCrossEntropy(grayImage);

    % Update GUI
    set(obj.handles.morphThresh,'String',obj.detection.morphThreshold);

    if(exist('source'))
      obj.detectSoma();
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  function threshold = minimizeCrossEntropy(obj, grayImage)


    % Implementing Brink and Pendock 1996
    % "Minimum Cross-Entropy Threshold Selection", eq 11.
    % doi:10.1016/0031-3203(95)00066-6   

    offset = 1; % To avoid division by 0

    minG = min(grayImage(:)+offset);
    maxG = max(grayImage(:)+offset);
  
    edges = transpose(minG:maxG);
    freqG = histc(grayImage(:),edges);

    threshold = NaN;
    minH = inf;

    % Golden section search
    phi = (1 + sqrt(5))/2;
    resphi = 2 - phi;

    G1 = minG;
    G2 = minG + (maxG-minG)*resphi;
    G3 = maxG;

    H1 = obj.crossEntropy(G1,edges,freqG);
    H2 = obj.crossEntropy(G2,edges,freqG);
    H3 = obj.crossEntropy(G3,edges,freqG);

    Hctr = 3;

    while(abs(G1 - G3) > 0.5)

      G4 = G2  + resphi*(G3 - G2);

      H4 = obj.crossEntropy(G4,edges,freqG);
      Hctr = Hctr + 1;

      if(H4 < H2)
        G1 = G2; H1 = H2;
        G2 = G4; H2 = H4;  
        % G3 unchanged
      else
        G3 = G1; H3 = H1;
        % G2 unchanged
        G1 = G4; H1 = H4;       
      end
  
    end

    threshold = round((G1 + G3)/2);

    % fprintf('Threshold %.0f minimizes cross-entropy.\n',threshold)
    % fprintf('Speedup is %.0f times\n', (maxG-minG+1)/Hctr);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Helper function to minimizeCrossEntropy

  function H = crossEntropy(obj, thresh, edges, freqG)

      % The two distributions are separated by a threshold
      % idx = find(edges == thresh);
      tmp = abs(edges - thresh);
      [value, idx] = min(tmp);
      % [~, idx] = min(tmp); % Only works in 2009b and later

      % How to deal with val0(1) = 0, that will give H=NaN
      val0 = edges(1:idx-1);
      freq0 = freqG(1:idx-1);

      val1 = edges(idx:end);
      freq1 = freqG(idx:end);
      
      % Calculate mean values for the two distributions
      mu0 = sum(val0.*freq0)/sum(freq0);
      mu1 = sum(val1.*freq1)/sum(freq1);
      
      % Calculate cross-entropy measure
      H = nansum(freq0.*(mu0.*log(mu0./val0) + val0.*log(val0./mu0))) ...
          + nansum(freq1.*(mu1.*log(mu1./val1) + val1.*log(val1./mu1)));

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % !!! Orphaned function

  function threshold = minimizeCrossEntropyOldSlow(obj,grayImage)


    % Implementing Brink and Pendock 1996
    % "Minimum Cross-Entropy Threshold Selection", eq 11.
    % doi:10.1016/0031-3203(95)00066-6   

    offset = 1; % To avoid division by 0

    minG = min(grayImage(:)+offset);
    maxG = max(grayImage(:)+offset);
  
    edges = transpose(minG:maxG);
    freqG = histc(grayImage(:),edges);

    threshold = NaN;
    minH = inf;

    allH = [];
    
    for thresh = minG:maxG
  
      % The two distributions are separated by a threshold
      idx = find(edges == thresh);

      % How to deal with val0(1) = 0, that will give H=NaN
      val0 = edges(1:idx-1);
      freq0 = freqG(1:idx-1);

      val1 = edges(idx:end);
      freq1 = freqG(idx:end);
      
      % Calculate mean values for the two distributions
      mu0 = sum(val0.*freq0)/sum(freq0);
      mu1 = sum(val1.*freq1)/sum(freq1);
      
      % Calculate cross-entropy measure
      H = nansum(freq0.*(mu0.*log(mu0./val0) + val0.*log(val0./mu0))) ...
          + nansum(freq1.*(mu1.*log(mu1./val1) + val1.*log(val1./mu1)));
    
      allH = [allH,H];
  
      if(H < minH)
        minH = H;
        threshold = thresh - offset;
      end
  
    end

    % fprintf('Threshold %d minimizes cross-entropy.\n',threshold)

    if(1)
      f = gcf;
      figure, plot(minG:maxG,allH)
      figure(f)
    end
  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function detectNeuritesSimple(obj)

    % Read in new settings from GUI using callback function
    obj.setGrowThresh();
    obj.setFilterSize();

    % Extract selected channel for morphology detection
    tmp = obj.data.image(:,:,obj.detection.morphChannel, obj.dispInfo.curImg);

    if(~isnan(obj.detection.wienerSize))
      tmp = wiener2(tmp,obj.detection.wienerSize*[1 1]);
    end

    if(isempty(obj.data.includeMask))
      obj.data.includeMask = ones(size(tmp));
    end

    if(obj.detection.morphThreshold < 1)
      % If between 0 and 1, use relative threshold
      mask = tmp > obj.detection.morphThreshold*max(tmp(:));
    else
      mask = tmp > obj.detection.morphThreshold;
    end

    mask = mask .* obj.data.includeMask;

    se = strel('disk',1);
    mask = imopen(mask,se);

    obj.data.neuriteMask = mask;

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function detectNeurite(obj,source, event)

    obj.saveUndo('detect neurites');

    if(obj.detection.useSteerableFilters)
      % Replace the threshold mask for the neurite with that from
      % the steerable filters
      obj.steerableFilter();
      obj.growNeuriteFromSoma();
    else
      obj.detectNeuritesSimple();
    end

    % Mark unconnected components in neurite mask
    obj.updateNeuriteMask();

    obj.showImage();

    obj.activateStage('synapses');

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function extendNeuritesHandler(obj, source, event)
    obj.saveUndo('extend neurites');

    if(nnz(obj.data.neuriteMask))
      % We require that there are some neurite parts detected
      obj.extendNeurites(5e-6/obj.data.xyRes);

    else
      disp('Detect neurites first!')
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % This function looks for neurite components which are a few pixels away

  function extendNeurites(obj,seedDist)

    distMask = obj.makeDistMask(obj.data.somaMask,obj.data.neuriteMask);

    % Locate neurite end points, set inf and NaN to 0 and find
    % largest remaining distances, those are our end points.
    tmpMask = distMask;
    tmpMask(distMask == inf) = 0;
    tmpMask(isnan(distMask)) = 0;
    [endY,endX] = find(imregionalmax(tmpMask));

    seedPoints = zeros(size(obj.data.neuriteMask));

    seedList = [];
    endList = [];

    for iEnd = 1:length(endX)

      if(endX(iEnd) < seedDist | endX(iEnd) > obj.data.width - seedDist ...
         | endY(iEnd) < seedDist | endY(iEnd) > obj.data.height - seedDist)
        % Too close to border, skip this point
        continue
      end

      % What is direction of the neurite at the end point?
      endDir = squeeze(obj.data.dirVect(endY(iEnd),endX(iEnd),:));

      % Which way is forward, +endDir or -endDir
      % Probe neurite by taking small step in one direction
      checkX = round(endX(iEnd) + endDir(1)*2);
      checkY = round(endY(iEnd) + endDir(2)*2);
 
      % Did we go towards neurite, or away from neurite?
      if(obj.data.neuriteMask(checkY,checkX))
        % Wrong way... try other way. 

        walkDir = [-endDir(1); -endDir(2)];
      else
        % It was right way, take big step in that direction
        walkDir = [endDir(1); endDir(2)];

      end

      % In addition to walking away from neurite, we also want
      % to check a few additional pixels in an arc +/- 15 degrees
      rotAngles = pi/180 * linspace(-15,15,5);

      for iRot = 1:length(rotAngles)

        theta = rotAngles(iRot);
        rotMat = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        rotWalkDir = rotMat*walkDir*seedDist;

        seedX = round(endX(iEnd) + rotWalkDir(1));
        seedY = round(endY(iEnd) + rotWalkDir(2));

        endList = [endList; endX(iEnd) endY(iEnd)];
        seedList = [seedList; seedX seedY];

        seedPoints(seedY,seedX) = 1;

      end

    end

    % Grow neurites from the new seed points
    obj.growNeurite(seedPoints);

    % See which neurites are large enough, and reconnect them

    tmpMask = zeros(size(obj.data.neuriteMask));

    for iSeed = 1:size(seedList,1)
      seedX = seedList(iSeed,1);
      seedY = seedList(iSeed,2);

      % Is the seed point part of a neurite?
      if(obj.data.neuriteMask(seedY,seedX))

         % Also check that the neurite is large enough
         % !!! Not implemented yet.

         endX = endList(iSeed,1);
         endY = endList(iSeed,2);

         % Connect seed and end point
         nDots = max(abs(endX-seedX),abs(endY-seedY))*2;

         xDots = round(linspace(endX,seedX,nDots));
         yDots = round(linspace(endY,seedY,nDots));

         for iDot = 1:length(xDots)
           tmpMask(yDots(iDot),xDots(iDot)) = 1;
         end
      end

    end

    tmpMask = imdilate(tmpMask, strel('disk',3));

    obj.data.neuriteMask = obj.data.neuriteMask + tmpMask;

    obj.updateNeuriteMask();

    obj.data.distMask = obj.makeDistMask(obj.data.somaMask,obj.data.neuriteMask);

    % Show the final image

    obj.showImage();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function updateNeuriteMask(obj)

    obj.data.neuriteMask = 2*(obj.data.neuriteMask > 0);
    CC = bwconncomp(obj.data.neuriteMask);
    somaIdx = find(obj.data.somaMask);

    for i = 1:length(CC.PixelIdxList)
      if(nnz(ismember(CC.PixelIdxList{i},somaIdx)))
        obj.data.neuriteMask(CC.PixelIdxList{i}) = 1;
      end
    end

    if(nnz(obj.data.neuriteMask == 2))
      obj.dispInfo.needCleaning = 1;
      set(obj.handles.clean,'BackgroundColor', [1 0.6 0]);
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function cleanMask(obj,source, event)

    obj.saveUndo('clean neurite mask');

    % This function removes unconnected components from the neurite mask
    % Those with distance set at inf.

    obj.data.neuriteMask(obj.data.distMask .* obj.data.neuriteMask == inf) = 0;
    obj.data.neuriteMask(isnan(obj.data.distMask .* obj.data.neuriteMask)) = 0;

    obj.dispInfo.needCleaning = 0;
    set(obj.handles.clean,'BackgroundColor', [0.8 0.8 0.8]);

    obj.showImage();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function killBlob(obj, source, event)

    disp('Select blob to remove.')
    [x,y,button] = ginput(1);

    a = axis;

    while(a(1) <= x & x <= a(2) ...
          & a(3) <= y & y <= a(4) ...
          & 1 <= x & x <= obj.data.width ...
          & 1 <= y & y <= obj.data.height ...
          & 1 == button)

      obj.saveUndo('blob removal');

      obj.data.neuriteMask = obj.cleanBlob(obj.data.neuriteMask,round(x),round(y));

      % Update image
      obj.showImage();

      disp('Select blob to remove.')
      % Check if we should remove one more
      [x,y,button] = ginput(1);
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function addNeurite(obj, source, event)

    neuriteSeed = zeros(size(obj.data.neuriteMask));

    disp('Add neurite.')
    [x,y,button] = ginput(1);
    
    a = axis();

    while(button == 1 ...
          & 1 <= x & x <= obj.data.width ...
          & 1 <= y & y <= obj.data.height ...
          & a(1) <= x & x <= a(2) ...
          & a(3) <= y & y <= a(4))

      neuriteSeed(round(y),round(x)) = 1;

      fprintf(['Detecting neurite at x=%d,y=%d ' ...
	       '(click outside, or click right button to stop)\n'], ...
	      round(x), round(y))

      obj.saveUndo(sprintf('add neurite (%d,%d)', round(x),round(y)))
      obj.growNeurite(neuriteSeed);

      disp('Add neurite.')
      [x,y,button] = ginput(1);

    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function addThinNeurites(obj, source, event)

    obj.saveUndo('add thin neurites');

    oldFilterSize = obj.detection.filterSize;
    obj.detection.filterSize = obj.detection.filterSize / 2;
    obj.steerableFilter();

    oldNeuriteMask = obj.data.neuriteMask > 0;
    obj.growNeurite(oldNeuriteMask - imerode(oldNeuriteMask,strel('disk',1)));
    
    addedNeuriteBits = (obj.data.neuriteMask>0) - oldNeuriteMask;

    CC = bwconncomp(addedNeuriteBits);
    newProps = regionprops(CC, ...
			   'PixelIdxList', ...
			   'Area', ...
			   'MajorAxisLength');

    % If they are shorter than 2 micrometers, remove.
    for i = 1:length(newProps)
      if(newProps(i).MajorAxisLength < obj.detection.minProtrusionLength/obj.data.xyRes)
        obj.data.neuriteMask(newProps(i).PixelIdxList) = 0;
      end
    end

    obj.updateNeuriteMask();
    obj.showImage();

    obj.detection.filterSize = oldFilterSize;
    obj.steerableFilter();

  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function distMask = makeDistMask(obj, somaMask, neuriteMask)

    % Create a padded version of the masks, this allows us to avoid 
    % range checks at edges
    somaMaskPadded = zeros(size(somaMask)+2);
    somaMaskPadded(2:end-1,2:end-1) = somaMask;
    neuriteMaskPadded = zeros(size(neuriteMask)+2);
    neuriteMaskPadded(2:end-1,2:end-1) = neuriteMask;

    distMaskPadded = inf*ones(size(somaMaskPadded));
    distMaskPadded(find(somaMaskPadded)) = 0;

    % Find all the pixels that we need to calculate the distance for
    [yNeurite,xNeurite] = ind2sub(size(neuriteMaskPadded), ...
				  find(neuriteMaskPadded));
    [ySoma,xSoma] = ind2sub(size(somaMaskPadded), ...
			    find(somaMaskPadded));

    % Seed the algorithm with the pixels belonging to the soma
    % Then we want to add all the immediate neighbours to those pixels 
    % and then keep adding their neighbours, etc, while keeping track of
    % the distance.
    xPrev = xSoma;
    yPrev = ySoma;
    allDone = 0;

    neighOfsX = [-1 0 1 -1 0 1 -1 0 1];
    neighOfsY = [1 1 1 0 0 0 -1 -1 -1];
    neighDist = [sqrt(2) 1 sqrt(2) 1 0 1 sqrt(2) 1 sqrt(2)]*obj.data.xyRes;
    
    idxNew = NaN;
    iter = 0;

    while(~isempty(idxNew))
      idxNew = [];
      iter = iter + 1;

      % Get the neighbours to the previously updated pixels
      for i = 1:length(xPrev)
        xN = xPrev(i) + neighOfsX;
        yN = yPrev(i) + neighOfsY;
        dN = distMaskPadded(yPrev(i),xPrev(i)) + neighDist;

        for j = 1:length(xN)
	  % Was the neighbour part of the neuron, and did we find
	  % a shorter path to this pixel than it had before?
	  if(neuriteMaskPadded(yN(j),xN(j)) > 0 ...
	     & distMaskPadded(yN(j),xN(j)) > dN(j))
            distMaskPadded(yN(j),xN(j)) = dN(j);

            % Keep track of the newest added
            idxNew = [idxNew, yN(j) + (xN(j)-1)*size(somaMaskPadded,1)];
          end
        end
      end

     idxNew = unique(idxNew);

     % fprintf('Iter: %d, num changed: %d\n', iter, length(idxNew))

     % Set the x and y values of the pixels previously updated
     [yPrev,xPrev] = ind2sub(size(somaMaskPadded),idxNew);

    end

    % Extract the original part of the mask
    distMask = distMaskPadded(2:end-1,2:end-1);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function synDist = calculateSynapseDist(obj)

    paddedDistMask = obj.makeDistMask(obj.data.somaMask, obj.padNeuriteMask());
				    
    synDist = paddedDistMask(obj.data.synapseCenter);

    if(nnz(synDist == inf))
      disp('This should not happen, talk to Johannes! Argh...')
      keyboard
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function detectSynapses(obj, source, event)

    % Get values from GUI
    obj.setSynThresh();
    obj.setNeuritePadding();
    obj.setSynapseMinSize();

    obj.saveUndo('detect synapses');

    % obj.dispInfo.axis = [];

    obj.detectSynapseCenters();

    obj.calculateSynapseProperties();

    obj.showImage();
    obj.activateStage('analyse');

    text(20,20,sprintf('Found %d synapses', length(obj.data.synapseCenter)), ...
         'color','white')

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function detectSynapseCenters(obj)

    % This calculates a preliminary threshold
    obj.makeSynapseTemplate(); 

    noiseMask = obj.thresholdSynapseChannel();

    % This one will exclude soma if obj.detection.excludeSomaSynapses = true
    paddedNeuriteMask = obj.padNeuriteMask();

    % Synapse centres are on neurites or its padding, but outside soma
    obj.data.synapseMask = noiseMask .* paddedNeuriteMask;

    obj.locateSynapseCenters();

    obj.filterSynapseMask(noiseMask);

    if(obj.debugFlag)
      keyboard
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function noiseMask = thresholdSynapseChannel(obj)

    synChan = obj.data.image(:,:,obj.detection.synChannel, obj.dispInfo.curImg);

    neuriteIdx = find(obj.data.neuriteMask);

    if(obj.detection.backgroundFiltering)
      synChan = ...
        max(0, synChan -imfilter(synChan, ...
				 fspecial('gaussian', ...
					  [3 3]*obj.detection.smoothingRadius, ...
					  obj.detection.smoothingRadius)));
    end

    synOnNeurite = synChan(neuriteIdx);

    % Threshold based detection of the synapse pixels
    if(~isnan(obj.detection.synThreshold))
      if(obj.detection.synThreshold < 0)
        disp(['Negative threshold specified, using the absolute value as absolute threshold'])
        synThresh = abs(obj.detection.synThreshold);
      else
        % Normal case
        synThresh = mean(synOnNeurite(:)) ...
            + obj.detection.synThreshold*std(synOnNeurite(:));
      end
      
    else
      % Auto detect threshold
      synThresh = obj.minimizeCrossEntropy(synOnNeurite);
      fprintf('Using cross-entropy. ')
    end

    % Save the actual threshold used
    obj.data.synapseIntensityThreshold = synThresh;

    fprintf('Absolute synapse intensity threshold used: %.1f\n', synThresh)

    % Note synChan might be filtered (see above)
    noiseMask = synChan > synThresh;

    % Remove isolated pixels by enforcing minimum size on synapses
    CCs = bwconncomp(noiseMask);

    for i = 1:length(CCs.PixelIdxList)
      if(length(CCs.PixelIdxList{i}) < obj.detection.minSynapseSize/obj.data.xyRes^2)
       noiseMask(CCs.PixelIdxList{i}) = 0;
      end
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function paddedNeuriteMask = padNeuriteMask(obj)

    if(obj.detection.excludeSomaSynapses) 
      % This is default
      paddedNeuriteMask = imdilate((obj.data.neuriteMask>0) - obj.data.somaMask > 0, ...
                                   strel('disk', ...
                                         round(obj.detection.neuritePadding ...
                                               /obj.data.xyRes)));

      paddedNeuriteMask = (paddedNeuriteMask - obj.data.somaMask) > 0;
    else
      %Include soma
      paddedNeuriteMask = imdilate(obj.data.neuriteMask>0, ...
				   strel('disk', ...
					 round(obj.detection.neuritePadding ...
					       /obj.data.xyRes)));


    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function makeSynapseTemplate(obj)

    noiseMask = obj.thresholdSynapseChannel();

    paddedNeuriteMask = obj.padNeuriteMask();

    morphChan = obj.data.image(:,:,obj.detection.morphChannel, obj.dispInfo.curImg);
    synChan = obj.data.image(:,:,obj.detection.synChannel, obj.dispInfo.curImg);
    XChan = obj.data.image(:,:,obj.detection.XChannel, obj.dispInfo.curImg);

    % Synapses centres are on neurites, but outside soma
    obj.data.synapseMask = noiseMask .* paddedNeuriteMask;

    putativeCenters = imregionalmax(synChan) .* noiseMask;
    putativeCenters = putativeCenters .* obj.data.synapseMask;
    putativeCentersIdx = find(putativeCenters);

    CCp = bwconncomp(putativeCenters);

    labeledSynapses = bwlabel(obj.data.synapseMask);

    % Look for synapses with just one clear maxima, 
    % and also only one center per synapse
    % used for auto-detecting synapse kernel
    centerMask = zeros(size(obj.data.neuriteMask));

    for i = 1:length(CCp.PixelIdxList)
      if((length(CCp.PixelIdxList{i}) == 1) ...
         & (1 == nnz(labeledSynapses(CCp.PixelIdxList{i}) ...
                     == labeledSynapses(putativeCentersIdx))))
        centerMask(CCp.PixelIdxList{i}) = 1;
      end
    end

    % Average those clear synapses together, to get a template single
    % single synapse (that we will then deconvolve image with to 
    % locate the putative synapse centres.

    [yC,xC] = ind2sub(size(centerMask),find(centerMask));

    meanSynapseMorph = zeros(obj.detection.maxRadie*2+1,obj.detection.maxRadie*2+1);
    meanSynapseSyn = zeros(obj.detection.maxRadie*2+1,obj.detection.maxRadie*2+1);
    meanSynapseX = zeros(obj.detection.maxRadie*2+1,obj.detection.maxRadie*2+1);

    synCtr = 0;
    for i = 1:length(xC)
      xReg = (xC(i)-obj.detection.maxRadie):(xC(i)+obj.detection.maxRadie);
      yReg = (yC(i)-obj.detection.maxRadie):(yC(i)+obj.detection.maxRadie);

      if(1 <= min(xReg) & max(xReg) <= obj.data.width ...
         & 1 <= min(yReg) & max(yReg) <= obj.data.height)
        meanSynapseMorph = meanSynapseMorph + morphChan(yReg,xReg);
        meanSynapseSyn = meanSynapseSyn + synChan(yReg,xReg);
        meanSynapseX = meanSynapseX + XChan(yReg,xReg);
        synCtr = synCtr + 1;
      end
    end

    meanSynapseMorph = meanSynapseMorph / synCtr;
    meanSynapseSyn = meanSynapseSyn / synCtr;
    meanSynapseX = meanSynapseX / synCtr;

    fprintf('Found %d single synapses with clear peaks (and %d blobs)\n', ...
	    synCtr, length(CCp.PixelIdxList)-synCtr)

    if(synCtr == 0)

      if(nnz(synChan) & isnan(obj.detection.singleSynapseRadie))
        % There is at least a non-zero pixel in synapse channel
        % warn the user then, that the automatic kernel failed.
        % But do not warn, if the user has specified their own
        % kernel size.

        uiwait(errordlg(['The program found no unique local maximas in ' ...
  		         'any of the synapses. Ouch! Guestimating the ' ...
                         'single synapse by as a gaussian with standard ' ...
                         'deviation 1 pixel.'], ...
		         'No single synapses detected!', 'modal'))
      end

      disp('No single synapses detected, we have a problem!')
      disp('Guessing a synapse kernel...')

      [y,x] = meshgrid(-obj.detection.maxRadie:obj.detection.maxRadie, ...
                       -obj.detection.maxRadie:obj.detection.maxRadie);
      meanSynapse = exp(-x.^2-y.^2);

      obj.data.meanSynapseMorph = [];
      obj.data.meanSynapseSyn = [];
      obj.data.meanSynapseX = [];
    else
      obj.data.meanSynapseMorph = meanSynapseMorph;
      obj.data.meanSynapseSyn = meanSynapseSyn;
      obj.data.meanSynapseX = meanSynapseX;

      % Use average synapse for auto-detected template
      meanSynapse = meanSynapseSyn;;
    end

    % Did the user specify their own synapse template, if so use that 
    % for deconvolution instead.
    if(~isnan(obj.detection.singleSynapseRadie))

      fprintf('Using template synapse, radie %.2d micrometers\n', ...
	      1e6*obj.detection.singleSynapseRadie)
      r = obj.detection.singleSynapseRadie/obj.data.xyRes;

      [y,x] = meshgrid(-obj.detection.maxRadie:obj.detection.maxRadie, ...
		       -obj.detection.maxRadie:obj.detection.maxRadie);
      meanSynapse = exp(-(x/r).^2-(y/r).^2);

    end

    % Template used for synapse detection
    obj.data.meanSynapse = meanSynapse / max(meanSynapse(:));

    % How does the intensity decrease with distance from synapse center
    % for single synapses.
    obj.calculateSynapseProfile();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function calculateSynapseProfile(obj)

    if(~isempty(obj.data.meanSynapseSyn))

      [y,x] = meshgrid(-obj.detection.maxRadie:obj.detection.maxRadie, ...
		       -obj.detection.maxRadie:obj.detection.maxRadie);

      % Distance to center for each pixel in template
      d = sqrt(x.^2 + y.^2)*obj.data.xyRes;

      binSize = obj.data.xyRes; %obj.detection.synapseProfileBinSize;
      distCenters = 0:binSize:min(3e-6,max(d(:)));

      valMorph = NaN*zeros(size(distCenters));  % Morphology staining
      valSyn   = NaN*zeros(size(distCenters));  % Synapse staining
      valX     = NaN*zeros(size(distCenters));  % X staining intensity

      for i = 1:length(distCenters)
        idx = find(distCenters(i)-binSize/2 <= d ...
		   & d < distCenters(i)+binSize/2);
        valMorph(i) = mean(obj.data.meanSynapseMorph(idx));
        valSyn(i) = mean(obj.data.meanSynapseSyn(idx));
        valX(i) = mean(obj.data.meanSynapseX(idx));
      end

      obj.data.meanSynapseProfileDist = distCenters;
      obj.data.meanSynapseProfileMorph = valMorph;
      obj.data.meanSynapseProfileSyn = valSyn;
      obj.data.meanSynapseProfileX = valX;

      plotCol = [1 0 0; 0 1 0; 0 0 1];

      set(obj.handles.fig,'CurrentAxes',obj.handles.singleSynapse)

      if(obj.dispInfo.showSynapseProfileDataPoints)
        plot(d(:)*1e6,obj.data.meanSynapseMorph(:), ...
             '.', 'color',plotCol(obj.detection.remapChan(1),:));
        hold on
        plot(d(:)*1e6,obj.data.meanSynapseSyn(:), ...
             '.', 'color',plotCol(obj.detection.remapChan(2),:));
        plot(d(:)*1e6,obj.data.meanSynapseX(:), ...
             '.', 'color',plotCol(obj.detection.remapChan(3),:));
      end

      plot(distCenters*1e6,valMorph,'-', ...
           'linewidth',3, 'color',plotCol(obj.detection.remapChan(1),:));
      hold on
      plot(distCenters*1e6,valSyn,'-', ...
           'linewidth',3, 'color',plotCol(obj.detection.remapChan(2),:));
      plot(distCenters*1e6,valX,'-', ...
           'linewidth',3, 'color',plotCol(obj.detection.remapChan(3),:));
      
      hold off

      title('Single synapse intensity profile')
      xlabel('Distance (\mum)')
      ylabel('Intensity')

    else
      % No synapses, make sure all variables are cleared
      obj.data.meanSynapseProfileDist = [];
      obj.data.meanSynapseProfileMorph = [];
      obj.data.meanSynapseProfileSyn = [];
      obj.data.meanSynapseProfileX = [];

      set(obj.handles.fig,'CurrentAxes',obj.handles.singleSynapse)
      cla
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function locateSynapseCenters(obj)

    % This function assumes the synapse mask and synapse template is already 
    % calculated it only detects the centers.

    synChan = obj.data.image(:,:,obj.detection.synChannel, obj.dispInfo.curImg);

    % Deconvolve synapse image with the mean single synapse
    % [qImg,rImg] = deconvreg(synChan,obj.data.meanSynapse);
    qImg = deconvlucy(synChan,obj.data.meanSynapse);
    % figure, imagesc(qImg)

    % Find the local maximas, that are above noise level in 
    % synapse channel, and within neurite (excluding soma).
    qMask = imregionalmax(qImg);
    qMask = qMask .* obj.data.synapseMask;
    %qMask = imregionalmax(qImg.*obj.data.synapseMask);


    % This just marks centres of synapses, not the entire synapse.
    % The export function assumes there are only centres masked right now,
    % but perhaps change it to obj.data.synapseMaskCentre
    obj.data.synapseCenter = find(qMask);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function filterSynapseMask(obj,mask)

    % Important, the synapse can go outside the neurite marker
    % so we want to keep those parts also, but still exclude soma!

    % For auto detect, mask = noiseMask are pixels above synapse noise level

    % Update: They have to be within obj.detection.neuritePadding
    % of the neurite

    paddedNeuriteMask = obj.padNeuriteMask();

    CC = bwconncomp(mask .* paddedNeuriteMask);

    obj.data.synapseMask = zeros(size(obj.data.synapseMask));

    neuriteIdx = find(obj.data.neuriteMask);

    for i = 1:length(CC.PixelIdxList)
      % Require that putative synapse regions contain at least one center
      % and that in addition to all being within the padded region around
      % a neurite, at least one pixel has to be directly on the neurite.
      if(nnz(ismember(obj.data.synapseCenter,CC.PixelIdxList{i})) ...
         & nnz(ismember(neuriteIdx,CC.PixelIdxList{i})))
        obj.data.synapseMask(CC.PixelIdxList{i}) = 1;
      else
        % fprintf('Orphaned pixels: %d\n', length(CC.PixelIdxList{i}))
      end
    end


  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function calculateSomaMeasures(obj)

    % Calculate some measures, radie, area etc

    somaProp = regionprops(obj.data.somaMask, ...
			   'Area', ...
			   'MajorAxisLength', ...
			   'MinorAxisLength');

    obj.data.somaArea = cat(1,somaProp.Area)*obj.data.xyRes^2*1e12;
    obj.data.somaMajorAxisLength = cat(1,somaProp.MajorAxisLength)*obj.data.xyRes*1e6;
    obj.data.somaMinorAxisLength = cat(1,somaProp.MinorAxisLength)*obj.data.xyRes*1e6;

    if(isempty(obj.data.somaMeasureMask))
      disp('No soma measure regions marked.')
      return
    end

    tmpMorph = obj.data.image(:,:,obj.detection.morphChannel, obj.dispInfo.curImg);
    obj.data.somaMeasureMorph = tmpMorph(find(obj.data.somaMeasureMask));

    tmpSyn = obj.data.image(:,:,obj.detection.synChannel, obj.dispInfo.curImg);
    obj.data.somaMeasureSyn = tmpSyn(find(obj.data.somaMeasureMask));

    tmpX = obj.data.image(:,:,obj.detection.XChannel, obj.dispInfo.curImg);
    obj.data.somaMeasureX = tmpX(find(obj.data.somaMeasureMask));


  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function calculateSynapseProperties(obj)

    disp('Calculating synapse properties.')

    % Remove centres that lack synapse mask pixels below them
    goodIdx = find(obj.data.synapseMask(obj.data.synapseCenter));
    obj.data.synapseCenter = obj.data.synapseCenter(goodIdx);

    % Determine which synapse center each synapse pixel belongs to
    CC = bwconncomp(obj.data.synapseMask);

    % Clear the old synapse pixels... 
    obj.data.synapsePixels = {};

    for i = 1:length(CC.PixelIdxList)
      % idx is the synapse center number, centreIdx is the location of
      % the centre in the synapseMask.
      idx = find(ismember(obj.data.synapseCenter,CC.PixelIdxList{i}));

      % Find the pixels belonging to each synapse centre
      % First calculate distance to all centres for every pixel
      % in synapse

      centerIdx = obj.data.synapseCenter(idx);
      tmpDist = zeros(length(CC.PixelIdxList{i}),length(centerIdx));

      for j = 1:length(centerIdx)
        [yPix,xPix] = ind2sub(size(obj.data.synapseMask),CC.PixelIdxList{i});
        [yCent,xCent] = ind2sub(size(obj.data.synapseMask),centerIdx(j));

        tmpDist(:,j) = sqrt((xPix-xCent).^2+(yPix-yCent).^2);
      end

      % The first column of sortIdx contains the indexes of the closest
      % synapses.
      [sortedDist, sortIdx] = sort(tmpDist,2);

      for j = 1:length(centerIdx)
        obj.data.synapsePixels{idx(j)} = ...
          CC.PixelIdxList{i}(find(sortIdx(:,1) == j));

      end

    end

    % Debug plot to verify that the synapse picture was correct
    if(0)
      tmpMaskR = zeros(size(obj.data.neuriteMask));
      tmpMaskG = zeros(size(obj.data.neuriteMask));
      tmpMaskB = zeros(size(obj.data.neuriteMask));

      for i = 1:length(obj.data.synapsePixels)
        tmpMaskR(obj.data.synapsePixels{i}) = rand(1);
        tmpMaskG(obj.data.synapsePixels{i}) = rand(1);
        tmpMaskB(obj.data.synapsePixels{i}) = rand(1);
      end

      % Mark centres
      tmpMaskR(obj.data.synapseCenter) = 1;
      tmpMaskG(obj.data.synapseCenter) = 1;
      tmpMaskB(obj.data.synapseCenter) = 1;

      tmpMask = zeros(size(obj.data.image));
      tmpMask(:,:,1) = tmpMaskR;
      tmpMask(:,:,2) = tmpMaskG;
      tmpMask(:,:,3) = tmpMaskB;

      f = gcf;
      figure, imagesc(tmpMask)
      figure(f);
    end

    % Store the distances of the synapses to the soma also
    obj.data.synapseDist = obj.calculateSynapseDist();

    % Calculate the average intensity in each synapse for each channel
    obj.data.synapseIntensityMorphMean = NaN*zeros(length(obj.data.synapseCenter),1);
    obj.data.synapseIntensityMorphSEM  = NaN*zeros(length(obj.data.synapseCenter),1);
    obj.data.synapseIntensitySynMean   = NaN*zeros(length(obj.data.synapseCenter),1);
    obj.data.synapseIntensitySynSEM    = NaN*zeros(length(obj.data.synapseCenter),1);
    obj.data.synapseIntensityXMean     = NaN*zeros(length(obj.data.synapseCenter),1);
    obj.data.synapseIntensityXSEM      = NaN*zeros(length(obj.data.synapseCenter),1);
    obj.data.synapseArea = NaN*zeros(size(obj.data.synapseCenter));

    tmpMorph = double(obj.data.image(:,:,obj.detection.morphChannel));
    tmpSyn = double(obj.data.image(:,:,obj.detection.synChannel));
    tmpX = double(obj.data.image(:,:,obj.detection.XChannel));

    for i = 1:length(obj.data.synapsePixels)

      nPixels = length(obj.data.synapsePixels{i});

      obj.data.synapseIntensityMorphMean(i) ...
        = mean(tmpMorph(obj.data.synapsePixels{i}));

      obj.data.synapseIntensityMorphSEM(i) ...
          = std(tmpMorph(obj.data.synapsePixels{i})) / sqrt(nPixels);

      obj.data.synapseIntensitySynMean(i) ...
          = mean(tmpSyn(obj.data.synapsePixels{i}));

      obj.data.synapseIntensitySynSEM(i) ...
          = std(tmpSyn(obj.data.synapsePixels{i})) / sqrt(nPixels);

      obj.data.synapseIntensityXMean(i) ...
        = mean(tmpX(obj.data.synapsePixels{i}));

      obj.data.synapseIntensityXSEM(i) ...
        = std(tmpX(obj.data.synapsePixels{i})) / sqrt(nPixels);

      obj.data.synapseArea(i) ...
        = length(obj.data.synapsePixels{i})*(obj.data.xyRes^2)*1e12; % Micrometers^2

    end

    % Calculate how many of synapse pixels are at max value
    sIdx = find(obj.data.synapseMask);
    scIdx = obj.data.synapseCenter;

    fprintf('Synapse pixels: %.2f%% at %d (Syn), %.2f%% at %d (X)\n', ...
	    100*nnz(tmpSyn(sIdx) == max(tmpSyn(sIdx)))/numel(sIdx), ...
	    max(tmpSyn(sIdx)), ...
	    100*nnz(tmpX(sIdx) == max(tmpX(sIdx)))/numel(sIdx), ...
	    max(tmpX(sIdx)))

    fprintf('Synapse centers: %.2f%% at %d (Syn), %.2f%% at %d (X)\n', ...
	    100*nnz(tmpSyn(scIdx) == max(tmpSyn(scIdx))) ...
	      /numel(scIdx), ...
	    max(tmpSyn(scIdx)), ...
	    100*nnz(tmpX(scIdx) == max(tmpX(scIdx)))/numel(scIdx), ...
	    max(tmpX(scIdx)))

    if(obj.dispInfo.verbose)
      figure
      subplot(2,2,1);
      hist(tmpSyn(sIdx),100)
      hold on, plot(mean(tmpSyn(sIdx)),0,'r*')
      xlabel('Syn')

      subplot(2,2,3);
      hist(tmpSyn(scIdx),100)
      hold on, plot(mean(tmpSyn(scIdx)),0,'r*')
      xlabel('Syn (cent)')

      subplot(2,2,2);
      hist(tmpX(sIdx),100)
      hold on, plot(mean(tmpX(sIdx)),0,'r*')
      xlabel('X')

      subplot(2,2,4);
      hist(tmpX(scIdx),100)
      hold on, plot(mean(tmpX(scIdx)),0,'r*')
      xlabel('X (cent)')
      figure(obj.handles.fig)
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setFigureTitle(obj)

    set(obj.handles.fig,'name', sprintf('SynD - Synapse Detection - %s', ...
				    obj.data.fileName{obj.dispInfo.curImg}));

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function nextImage(obj, source, event)
    obj.dispInfo.curImg = min(obj.dispInfo.curImg + 1,obj.data.num);
    obj.setFigureTitle();				    
    obj.showImage();
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function prevImage(obj, source, event)
    obj.dispInfo.curImg = max(obj.dispInfo.curImg - 1,1);
    obj.setFigureTitle();				    
    obj.showImage();
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function showLogo(obj)

    set(obj.handles.fig,'CurrentAxes',obj.handles.image);
    try
      if(isempty(obj.bilder.logo))
        obj.bilder.logo = imread('bilder/SynD-logo.jpg');
      end

      set(obj.handles.image,'units','pixels')
      pos = get(obj.handles.image,'position');
      set(obj.handles.image,'units','normalized')

      width = pos(3) - pos(1);
      height = pos(4) - pos(2);

      imagesc(imresize(obj.bilder.logo, [height width]));
    catch e
      getReport(e)
      disp('Unable to load logo.')
    end

    axis off

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function loadNeuron(obj, source, event)

    curPwd = pwd;
    try
      cd(obj.data.loadPath);
    catch
      fprintf('Unable to change to %s\n', obj.data.loadPath)
    end

    fileTypes = { '*.lsm', 'LSM-image'; ...
		  '*.tif;*.tiff', 'Tiff-image' };
    fileTypes = fileTypes(obj.editInfo.fileTypeOrder,:);

    [imageFile, imagePath, filterIndex] = ...
      uigetfile(fileTypes, ...
		'Select neuron image', ...
		'MultiSelect', 'off');

    cd(curPwd);

    if(~iscell(imageFile) & imageFile == 0)
      % User pressed cancel
      return;
    end

    obj.data.loadPath = imagePath;

    obj.dispInfo.stage = 1;

    % Clearing all the old data
    obj.clearData();
   
    switch(obj.editInfo.fileTypeOrder(filterIndex))
      case 1
        obj.loadLSMfile(imageFile,imagePath);
      case 2
        obj.loadTiffFiles(imageFile,imagePath);
    end

    obj.data.neuriteMask = zeros(size(obj.data.image,1),size(obj.data.image,2));
    obj.data.somaMask    = zeros(size(obj.data.image,1),size(obj.data.image,2));
    obj.data.synapseMask = zeros(size(obj.data.image,1),size(obj.data.image,2));

    % Make sure the used filter is default next time
    switch(obj.editInfo.fileTypeOrder(filterIndex))
      case 1
        obj.editInfo.fileTypeOrder = [1 2];
      case 2
        obj.editInfo.fileTypeOrder = [2 1]; 
    end
    
    % Update the title of the window
    obj.setFigureTitle();				    

    obj.activateStage('soma');
    obj.showGUI('load');
    obj.showImage();
    obj.clearUndo();

    axis on

    % Allow the user to zoom using the mouse wheel
    obj.startZoomWheelMode();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Load tiff can handle multi-select

  function loadTiffFiles(obj,tiffFile,tiffPath)

    if(~iscell(tiffFile))
      tiffFile = {tiffFile};
    end

    for i = 1:length(tiffFile)
      fprintf('Loading %s\n', tiffFile{i})

      fileName = [tiffPath tiffFile{i}];
      tiffInfo = imfinfo(fileName);

      obj.data.height = tiffInfo(1).Height;
      obj.data.width = tiffInfo(1).Width;

      if(obj.data.num == 0)
        obj.data.image = zeros(obj.data.height,obj.data.width,3,1);
      end

      for j = 1:length(tiffInfo)
        obj.data.num = obj.data.num+1;
        tmp = imread(fileName,j);
        switch(size(tmp,3))

          case 1
            oldTmp = tmp;
            tmp = zeros(obj.data.height,obj.data.width,3);
            tmp(:,:,1) = oldTmp;

          case 2
            oldTmp = tmp;
            tmp = zeros(obj.data.height,obj.data.width,3);
            tmp(:,:,1) = oldTmp(:,:,1);
            tmp(:,:,2) = oldTmp(:,:,2);

          case 3
            % No modification needed, use tmp as it is.

          otherwise
            uiwait(warndlg(['SynD only supports 1 to 3 colour channels, ' ...
                            'only using three first channels.'],  ...
                           'Image error'));

            tmp = tmp(:,:,1:3);
        end

        obj.data.image(:,:,:,obj.data.num) = tmp;
        obj.data.fileName{end+1} = tiffFile{i};
      end

      % Should we collapse image stack
      if(obj.singleChannelStack())
        obj.collapseImageStack();
      end

      obj.calcMaxIntensity();

      % Reset scaling
      obj.dispInfo.scaleRed = 1;
      obj.dispInfo.scaleBlue = 1;
      obj.dispInfo.scaleGreen = 1;

    end

    obj.dispInfo.axis = [];

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function calcMaxIntensity(obj)

    disp('Calculating max intensity')
		      
    img = obj.getRemappedImage();

    % Recalculate the new R,G,B max intensities, used for scaling
    tmpR = img(:,:,1);
    obj.data.maxRed = double(max([tmpR(:); 1]));

    tmpG = img(:,:,2);
    obj.data.maxGreen = double(max([tmpG(:); 1]));

    tmpB = img(:,:,3);
    obj.data.maxBlue = double(max([tmpB(:); 1]));

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Load LSM can not handle multi-select

  function loadLSMfile(obj,lsmFile,lsmPath)

    fprintf('Loading %s\n', lsmFile)

    lsmFileName = strcat(lsmPath,lsmFile);
    [lsmInfo,scanInfo,imInfo] = lsminfo(lsmFileName);
    imgStack = tiffread(lsmFileName);

    obj.data.height = imgStack(1).height;
    obj.data.width = imgStack(1).width;
    obj.data.image = zeros(obj.data.height,obj.data.width,3,1);

    for j = 1:length(imgStack)

      if(iscell(imgStack(j).data))
        for i = 1:length(imgStack(j).data)
          obj.data.image(:,:,i,j) = imgStack(j).data{i};
        end
      else

        % Only one channel present, put it in red
        obj.data.image(:,:,1,j) = imgStack(j).data;      
      end

      if(j > 1)
        % Make sure voxel size etc matches
        if(imgStack(1).lsm.VoxelSizeX ~= imgStack(j).lsm.VoxelSizeX)
          disp('Voxel size inconsistent in the image! You are screwed.')
        end

      end

      obj.data.fileName{j} = lsmFile;

    end

    obj.data.num = length(imgStack);


    obj.data.xyRes = imgStack(1).lsm.VoxelSizeX;

    if(imgStack(1).lsm.VoxelSizeX ~= imgStack(1).lsm.VoxelSizeY)
      % If you ever see this warning, let me know.
      disp('Warning: X and Y resolution differ, code assumes they are same.')
      disp(sprintf('Using %d m', obj.data.xyRes))
    end

    obj.calcMaxIntensity();

    % Reset scaling
    obj.dispInfo.scaleRed = 1;
    obj.dispInfo.scaleBlue = 1;
    obj.dispInfo.scaleGreen = 1;

    obj.dispInfo.axis = [];

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [selectedFiles,MATpath,LSMpath] = selectBatchReanalyseFiles(obj)

    curPwd = pwd;
    try
      cd(obj.data.loadPath);
    catch
      fprintf('Unable to change to %s\n', obj.data.loadPath)
    end

    [saveFile, savePath, filterIndex] = ...
      uigetfile({ '*-save.mat', 'Old save-file' }, ...
		'Select old save mat file', ...
		'MultiSelect', 'on');

    cd(curPwd);

    if(~iscell(saveFile) & saveFile == 0)
      % User pressed cancel
      selectedFiles = {};
      return;
    elseif(~iscell(saveFile))
      selectedFiles = {saveFile}; 
    else
      selectedFiles = saveFile; 
    end


    MATpath = savePath; % Where the old mat save files were stored
    LSMpath = uigetdir(savePath,'LSM/TIFF directory');

    if(isempty(LSMpath) | isempty(MATpath))
      selectedFiles = {};
      disp('User aborted reanalysis')
    else
      if(LSMpath(end) ~= '/')
        LSMpath(end+1) = '/';
      end
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [selectedFiles,fileType] = selectBatchFiles(obj)

    curPwd = pwd;
    try
      cd(obj.data.loadPath);
    catch
      fprintf('Unable to change to %s\n', obj.data.loadPath)
    end

    fileTypes = { '*.lsm', 'LSM-image'; ...
		  '*.tif*', 'Tiff-image' };
    fileTypes = fileTypes(obj.editInfo.fileTypeOrder,:);

    [imageFile, imagePath, filterIndex] = ...
      uigetfile(fileTypes, ...
		'Select neuron image', ...
		'MultiSelect', 'on');

    cd(curPwd);

    if(~iscell(imageFile) & imageFile == 0)
      % User pressed cancel
      selectedFiles = {};
      return;
    elseif(~iscell(imageFile))
      selectedFiles = {imageFile}; 
    else
      selectedFiles = imageFile; 
    end

    obj.data.loadPath = imagePath;
    fileType = obj.editInfo.fileTypeOrder(filterIndex);

    % Make sure the used filter is default next time
    switch(fileType)
      case 1
        obj.editInfo.fileTypeOrder = [1 2];
      case 2
        obj.editInfo.fileTypeOrder = [2 1]; 
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function guessBatchExportFiles(obj)

    obj.data.exportXMLfile = strcat(obj.data.fileName{obj.dispInfo.curImg}, ...
				'-export.xml');

    obj.data.exportSaveFile = strrep(obj.data.exportXMLfile,'-export.xml','-save.mat');

    obj.data.exportNeuriteMaskFile = ...
      strrep(obj.data.exportXMLfile,'-export.xml','-neuritemask.tiff');
  
    obj.data.exportSynapseMaskFile = ...
      strrep(obj.data.exportXMLfile,'-export.xml','-synapsemask.tiff');

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function sccFlag = singleChannelStack(obj)

    sccFlag = 0;

    if(nnz(obj.data.image(:,:,2,1)) == 0 ...
       & nnz(obj.data.image(:,:,3,1)) == 0 ...
       & size(obj.data.image,4) > 1)

      % This is a single channel image stack, collapse!
      sccFlag = 1;
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function collapseImageStack(obj,source, event)

    if(obj.data.num == 0)
      disp('No images loaded, unable to callapse image stack.')
      return
    end

    disp('Collapsing image stack.')
    disp(['Creating a new image with three channels out of the' ...
          ' first channel of three images'])

    newImage = zeros(obj.data.height,obj.data.width,3,1);

    for i = 1:min(3,obj.data.num)
      newImage(:,:,i,1) = obj.data.image(:,:,1,i);
    end

    obj.data.image = newImage;
    obj.data.num = 1;
    obj.dispInfo.curImg = 1;

    obj.showImage();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function maxProjection(obj, source, event)

    newImage = zeros(obj.data.height,obj.data.width,3,1);

    for i = 1:size(obj.data.image,4)
      newImage = newImage + obj.data.image(:,:,:,i);
    end

    obj.data.image = newImage;
    obj.data.num = 1;
    obj.dispInfo.curImg = 1;

    obj.calcMaxIntensity();
    obj.showImage();

    obj.stopEdit();
    obj.activateStage('soma');
    obj.showGUI('load');


  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function batchAnalyse(obj, source, event)

    %%% Select neurons to analyse
    [selectedFiles,fileType] = obj.selectBatchFiles();

    %%% Load the neurons and analyse them one at a time.
    for iFiles = 1:length(selectedFiles)

      % Clean out the old
      obj.clearData();

      fprintf('Analysing %s\n', selectedFiles{iFiles})

      % Load file
      switch(fileType)
        case 1
          % LSM image
          obj.loadLSMfile(selectedFiles{iFiles},obj.data.loadPath);    
        case 2
          % TIFF image
          obj.loadTiffFiles(selectedFiles{iFiles},obj.data.loadPath);
      end

      % Analyse
      obj.showGUI('soma');

      obj.detectSoma();
      obj.showMasks();
      obj.drawnow();

      obj.showGUI('neurites');
      obj.detectNeurite();
      obj.cleanMask();

      obj.showMasks();
      obj.drawnow();

      obj.addThinNeurites();
      obj.cleanMask();

      obj.showMasks();
      obj.drawnow();

      obj.showGUI('synapses'); 
      obj.detectSynapses();

      obj.showMasks();
      obj.drawnow();

      obj.showGUI('analyse'); 
      obj.guessBatchExportFiles();
      obj.exportData();

    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function fileType = getFileType(obj, filename)

    idx = find(filename == '.',1,'last');
    fileEnding = filename(idx:end);

    if(strcmpi(fileEnding,'.tif'))
      fileType = 2;
    elseif(strcmpi(fileEnding,'.tiff'))
      fileType = 2;
    elseif(strcmpi(fileEnding,'.lsm'))
      fileType = 1;
    else
      fprintf('Unable file ending: %s\n', fileEnding)
      fileType = 0;
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function clearData(obj)

    disp('Clearing obj.data.')

    obj.data.num = 0;
    obj.data.image = [];
    obj.data.fileName = {};

    obj.data.somaMask = [];
    obj.data.somaMeasureMask = [];

    obj.data.neuriteMask = [];
    obj.data.skeleton = [];
    obj.data.includeMask = [];

    obj.data.synapseMask = [];
    obj.data.synapseCenter = [];
    obj.data.synapseArea = [];
    obj.data.synapseDist = [];
    obj.data.synapsePixels = {};
    obj.data.neuriteLength = NaN;

    obj.data.synapseIntensityMorphMean = []; % morph-channel
    obj.data.synapseIntensityMorphSEM = [];  % morph-channel
    obj.data.synapseIntensitySynMean = [];   % syn-channel
    obj.data.synapseIntensitySynSEM = [];    % syn-channel
    obj.data.synapseIntensityXMean = [];     % X-channel
    obj.data.synapseIntensityXSEM = [];      % X-channel

    obj.data.meanSynapse = [];

    obj.data.meanSynapseMorph = [];
    obj.data.meanSynapseSyn = [];
    obj.data.meanSynapseX = [];    % Profile of X channel for synapse

    obj.data.meanSynapseProfileDist = [];
    obj.data.meanSynapseProfileMorph = [];
    obj.data.meanSynapseProfileSyn = [];
    obj.data.meanSynapseProfileX = [];

    obj.data.skeleton = [];
    obj.data.shollEdges = [NaN NaN NaN];
    obj.data.shollDendHist = [];
    obj.data.shollSynHist = [];

    obj.data.shollIntMorphMean = [];
    obj.data.shollIntMorphStd = [];
    obj.data.shollIntMorphSEM = [];

    obj.data.shollIntSynMean = [];
    obj.data.shollIntSynStd = [];
    obj.data.shollIntSynSEM = [];

    obj.data.shollIntXMean = [];
    obj.data.shollIntXStd = [];
    obj.data.shollIntXSEM = [];

    obj.data.somaMeasureMorph = NaN;
    obj.data.somaMeasureSyn = NaN;
    obj.data.somaMeasureX = NaN;

    obj.data.somaArea = [];
    obj.data.somaMajorAxisLength = [];
    obj.data.somaMinorAxisLength = [];

    obj.data.measurePoint = [];

    obj.data.intensityHistogramEdges = [];
    obj.data.morphHist = [];
    obj.data.synHist = [];
    obj.data.XHist = [];

    % Intensity gradient along neurites
    obj.data.gradientEdges = [];

    obj.data.morphGradientMean = [];
    obj.data.morphGradientStd = [];
    obj.data.morphGradientSEM = [];

    obj.data.synGradientMean = [];
    obj.data.synGradientStd = [];
    obj.data.synGradientSEM = [];

    obj.data.XGradientMean = [];
    obj.data.XGradientStd = [];
    obj.data.XGradientSEM = [];

    obj.data.distMask = [];

    obj.data.dirVect = [];
    obj.data.rigidity = [];
    obj.data.lambdaMax = 0;

    % These are absolute thresholds used for detection
    obj.data.somaIntensityThreshold = NaN;
    obj.data.synapseIntensityThreshold = NaN;

    % Clear to prevent overwrites...
    obj.data.exportXMLfile = [];
    obj.data.exportSaveFile = [];
    obj.data.exportNeuriteMaskFile = [];
    obj.data.exportSynapseMaskFile = [];

    obj.dispInfo.curImg = 1;

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function batchReAnalyse(obj, source, event)

    obj.detection.reanalyseFlag = 1;

    doSynapseDetection = NaN;

    [selectedFiles, MATpath,LSMpath] = obj.selectBatchReanalyseFiles();

    oldExportPath = obj.data.exportPath;
    obj.data.exportPath = MATpath;

    %%% Load the neurons and analyse them one at a time.
    for iFiles = 1:length(selectedFiles)

      % Clean out the old
      obj.clearData();

      try
        old = load(strcat(MATpath,selectedFiles{iFiles}));
        old = old.old;

        try
          obj.detection.morphChannel = old.morphChannel;
          obj.detection.synChannel = old.synChannel;
          obj.detection.XChannel = old.XChannel;
          obj.detection.remapChan = old.remapChan;
        catch
          disp('Old channel info. Assuming not remapped')
          old.remapChan = [obj.detection.morphChannel, ...
			   obj.detection.synChannel, ...
			   obj.detection.XChannel];
        end

        try
          obj.detection.excludeSomaSynapses = old.excludeSomaSynapses;
        catch
          disp('Excluding soma synapses.')
          obj.detection.excludeSomaSynapses = true;
        end

        % Restore old XY-res (in case of a tiff-file)
        obj.data.xyRes = old.xyRes;

        fileType = obj.getFileType(old.fileName{1}); 
        switch(fileType)
          case 1
            fprintf('Loading %s%s\n',LSMpath,old.fileName{1})
            obj.loadLSMfile(old.fileName{1},LSMpath);    

          case 2 
            obj.loadTiffFiles(old.fileName{1},LSMpath);

          otherwise
            fprintf('Unknown file type: %s', old.fileName{1})
            continue
        end

        fprintf('Re-analysing %s\n', selectedFiles{iFiles})

        try
          obj.data.neuriteMask = old.neuriteMask;
          obj.data.somaMask = old.somaMask;
        catch
          % Old file names in old versions of program...
          obj.data.neuriteMask = old.neuriteMaskIdx;
          obj.data.somaMask = old.somaMaskIdx;
        end

        % Load the old intensity threshold used
        try
          obj.data.somaIntensityThreshold = old.somaIntensityThreshold;
        catch
          obj.data.somaIntensityThreshold = old.obj.detection.morphThreshold;
        end

        obj.exportInfo.saveSholl = 1;
        obj.exportInfo.saveSynapseProfile = 1;
        obj.exportInfo.saveMat = 1;
        obj.exportInfo.saveIntensityHistogram = 1;
        obj.exportInfo.saveGradient = 1;

        try
          obj.data.somaMeasureMask = old.somaMeasureMask;
          obj.exportInfo.saveSomaMeasure = 1;
        catch
          disp('No soma measure marked, ignoring.')
          obj.data.somaMeasureMask = [];
        end

        obj.writeExportSettings();

        obj.showGUI('soma');
        obj.showMasks();
        obj.drawnow();

        obj.calculateSomaMeasures();
        obj.dispInfo.measurePointColor = obj.dispInfo.defaultMeasurePointColor;
        obj.showMasks();
        obj.drawnow();

        obj.showGUI('neurites');
        obj.updateNeuriteMask();
        obj.data.distMask = obj.makeDistMask(obj.data.somaMask,obj.data.neuriteMask);

        obj.showMasks();
        obj.drawnow();

        if(~nnz(old.synapseMask))
          if(isnan(doSynapseDetection))
            % First time, ask user if they want to do synapse detection        
            answer = questdlg('Should I detect synapses?', ...
			      'Missing synapses', ...
			      'Detect all','Skip all',...
			      'Skip all');

            switch(answer)
              case 'Detect all'
                doSynapseDetection = 1;
              case 'Skip all'
                doSynapseDetection = 0;
              otherwise
                doSynapseDetection = 0;
            end

          end

        end

        if(nnz(old.synapseMask) | doSynapseDetection)
          obj.showGUI('synapses'); 
          obj.detectSynapses();

          obj.showMasks();
          obj.drawnow();

        end

        obj.activateStage('analyse')
        obj.showGUI('analyse'); 
        obj.guessBatchExportFiles();

        % Modifying the save-file so we do not overwrite the old save file
        obj.data.exportSaveFile = ...
            strrep(obj.data.exportSaveFile,'-save.mat','-re-save.mat');

        obj.exportData();


      catch e
        getReport(e)
        fprintf('Failed to re-analyse %s\n', selectedFiles{iFiles})

      end

    end

    % Restore the old export path
    obj.data.exportPath = oldExportPath;

    % Exit re-analyse mode.
    obj.detection.reanalyseFlag = 0;

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function mask = largestComponent(obj, mask)

    CC = bwconncomp(mask);

    maxCC = 0;
    maxCCidx = 0;
    for i = 1:length(CC.PixelIdxList)
      if(length(CC.PixelIdxList{i}) > maxCC)
        maxCCidx = i;
        maxCC = length(CC.PixelIdxList{i});
      end
    end

    % Keep only largest connected component when doing soma detection
    mask(:) = 0;
    if(maxCCidx)
      mask(CC.PixelIdxList{maxCCidx}) = 1;
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function detectSoma(obj, source, event)

    % Should we auto-detect the threshold?
    if(obj.detection.autoSomaThreshold)
      obj.autoSomaThreshold();
    end

    % Save the actual threshold used
    obj.data.somaIntensityThreshold = obj.detection.morphThreshold;

    % Read in latest inputs from GUI (so user doesnt have to hit enter)
    obj.setMorphThresh();
    obj.setSomaErodeRadius();

    obj.saveUndo('detect soma');

    % This gives soma and neurites
    obj.detectNeuritesSimple(); 

    % Next we erode away the neurites, a bit further down...

    if(obj.detection.singleSoma)
      % If only one soma, use largest component.
      mask = obj.largestComponent(obj.data.neuriteMask);
    else
      mask = obj.data.neuriteMask;
    end

    if(~nnz(obj.data.neuriteMask))
      uiwait(errordlg(['No pixels in channel selected for morphology. ' ...
                       'Check that the morphology channel is selected, ' ...
                       'also try lowering the threshold and soma erode radius.'], ...
		      'Bad channel selected', 'modal'))
      return
    end
    
    seSoma = strel('disk',obj.detection.somaErodeRadius);

    obj.data.somaMask = imopen(mask,seSoma);

    if(nnz(obj.data.somaMask) == 0)
      % Display the neurite mask for the ellusive soma
      obj.dispInfo.neuriteColor = [1 1 1]*0.7;
      obj.showImage();

      uiwait(errordlg(['No soma detected, reduce soma erode radius and ' ...
		       'verify you have the right channel for morphology'], ...
		      'Overly aggressive erosion', 'modal'))

      % Set focus to erode radius
      uicontrol(obj.handles.somaErodeRadius);
    else
      obj.dispInfo.neuriteColor = NaN;
      obj.data.neuriteMask = zeros(size(obj.data.neuriteMask));
      obj.showImage();

      obj.activateStage('neurites');

      CC = bwconncomp(obj.data.somaMask);
      text(20,20,sprintf('Found %d soma', CC.NumObjects), ...
	   'color','white')

    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function showHist(obj)

    nHist = 30;

    img = obj.getRemappedImage();

    % Create histogram over red channel
    % The part of the histogram that is saturated is marked black
    set(obj.handles.fig,'CurrentAxes',obj.handles.redHist)    

    tmp = img(:,:,1);
    try
      redEdges = linspace(0,max(1,obj.data.maxRed),nHist);
    catch e
        getReport(e)
        keyboard
    end
    nRed = histc(tmp(:),redEdges);
    dCenter = diff(redEdges(1:2))/2;
    idx = find(redEdges*obj.dispInfo.scaleRed <= obj.data.maxRed);

    bar(redEdges+dCenter,nRed, ...
        'facecolor',[0 0 0], 'edgecolor',[0 0 0])
    hold on
    bar(redEdges(idx)+dCenter,nRed(idx), ...
        'facecolor',[1 0 0], 'edgecolor',[1 0 0])
    hold off
    set(gca,'YScale','log')
    axis off
    a = axis; a(2) = max([obj.data.maxRed obj.data.maxGreen obj.data.maxBlue 1]); axis(a);

    set(obj.handles.fig,'CurrentAxes',obj.handles.greenHist)    
    tmp = img(:,:,2);
    greenEdges = linspace(0,max(1,obj.data.maxGreen),nHist);
    nGreen = histc(tmp(:),greenEdges);
    dCenter = diff(greenEdges(1:2))/2;
    idx = find(greenEdges*obj.dispInfo.scaleGreen <= obj.data.maxGreen);

    bar(greenEdges+dCenter,nGreen, ...
        'facecolor',[0 0 0], 'edgecolor',[0 0 0])
    hold on
    bar(greenEdges(idx)+dCenter,nGreen(idx), ...
        'facecolor',[0 1 0], 'edgecolor',[0 1 0])
    hold off
    set(gca,'YScale','log')
    axis off
    a = axis; a(2) = max([obj.data.maxRed obj.data.maxGreen obj.data.maxBlue]); axis(a);

    set(obj.handles.fig,'CurrentAxes',obj.handles.blueHist)    
    tmp = img(:,:,3);
    blueEdges = linspace(0,max(1,obj.data.maxBlue),nHist);
    nBlue = histc(tmp(:),blueEdges);
    dCenter = diff(blueEdges(1:2))/2;
    idx = find(blueEdges*obj.dispInfo.scaleBlue <= obj.data.maxBlue);

    bar(blueEdges+dCenter,nBlue, ...
        'facecolor',[0 0 0], 'edgecolor',[0 0 0])
    hold on
    bar(blueEdges(idx)+dCenter,nBlue(idx), ...
        'facecolor',[0 0 1], 'edgecolor',[0 0 1])
    hold off
    set(gca,'YScale','log')
    axis off
    a = axis; a(2) = max([obj.data.maxRed obj.data.maxGreen obj.data.maxBlue]); axis(a);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function showSkeleton(obj, source, event)
 
    % Recalculate the skeleton
    obj.data.skeleton = bwmorph(obj.data.neuriteMask-obj.data.somaMask>0,'skel',inf);
    obj.trimSkeleton();

    obj.dispInfo.showSkeleton = 1;

    obj.showMasks();

    % obj.dispInfo.showSkeleton = 0;

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % altFig is normaly not specified, but can be used to draw 
  % in alternative figure, for example when making images

  function showMasks(obj, altFig)

    if(exist('altFig'))
      % If we want to draw the masks in an alternative figure
      figure(altFig)
    else 
      set(obj.handles.fig,'CurrentAxes',obj.handles.image)    
    end
    tmp = obj.getImg();

    tmpR = tmp(:,:,1);
    tmpG = tmp(:,:,2);
    tmpB = tmp(:,:,3);

    % Should we show neurites?
    if(~isnan(obj.dispInfo.neuriteColor) & obj.dispInfo.showMask)
      tmpR(find(obj.data.neuriteMask == 1)) = obj.dispInfo.neuriteColor(1);
      tmpG(find(obj.data.neuriteMask == 1)) = obj.dispInfo.neuriteColor(2);
      tmpB(find(obj.data.neuriteMask == 1)) = obj.dispInfo.neuriteColor(3);

      tmpR(find(obj.data.neuriteMask == 2)) = 0.4*obj.dispInfo.neuriteColor(1);
      tmpG(find(obj.data.neuriteMask == 2)) = 0.4*obj.dispInfo.neuriteColor(2);
      tmpB(find(obj.data.neuriteMask == 2)) = 0.4*obj.dispInfo.neuriteColor(3);

      if(obj.dispInfo.showSkeleton & ~isempty(obj.data.skeleton))
        tmpR(find(obj.data.skeleton)) = obj.dispInfo.neuriteColor(1);
        tmpG(find(obj.data.skeleton)) = 0.5*obj.dispInfo.neuriteColor(1);
        tmpB(find(obj.data.skeleton)) = 0.5*obj.dispInfo.neuriteColor(1);
      end

    end

    % Should we display soma?
    if(~isnan(obj.dispInfo.somaColor) & obj.dispInfo.showMask)
      tmpR(find(obj.data.somaMask)) = obj.dispInfo.somaColor(1);
      tmpG(find(obj.data.somaMask)) = obj.dispInfo.somaColor(2);
      tmpB(find(obj.data.somaMask)) = obj.dispInfo.somaColor(3);
    end 

    % Should we mark synapses
    if(~isnan(obj.dispInfo.synapseColor) & obj.dispInfo.showMask)
      % Mark synapses
      tmpR(find(obj.data.synapseMask)) = 0.6*obj.dispInfo.synapseColor(1);
      tmpG(find(obj.data.synapseMask)) = 0.6*obj.dispInfo.synapseColor(2);
      tmpB(find(obj.data.synapseMask)) = 0.6*obj.dispInfo.synapseColor(3);

      % Mark synapse centres
      tmpR(obj.data.synapseCenter) = obj.dispInfo.synapseColor(1);
      tmpG(obj.data.synapseCenter) = obj.dispInfo.synapseColor(2);
      tmpB(obj.data.synapseCenter) = obj.dispInfo.synapseColor(3);
    end

    if(~isnan(obj.dispInfo.measurePointColor) & obj.dispInfo.showMask)
      % Mark soma measure points
      tmpR(find(obj.data.somaMeasureMask)) = obj.dispInfo.measurePointColor(1);
      tmpG(find(obj.data.somaMeasureMask)) = obj.dispInfo.measurePointColor(2);
      tmpB(find(obj.data.somaMeasureMask)) = obj.dispInfo.measurePointColor(3);
    end

    tmp(:,:,1) = tmpR;
    tmp(:,:,2) = tmpG;
    tmp(:,:,3) = tmpB;

    imagesc(tmp);
    axis equal

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Returns remapped but un-normalized image

  function img = getRemappedImage(obj)

     if(isempty(obj.data.image))
       img = zeros(0,0,3);
       return
     end

     origImg = squeeze(obj.data.image(:,:,:,obj.dispInfo.curImg));

     % Red here refers to the red in the original image, ie channel 1
     % Green refers to channel 2, and Blue to channel 3

     % However, to confuse things obj.dispInfo.showRed refers
     % to the red colour showed to the user.

     img = zeros(obj.data.height,obj.data.width,3);

     img(:,:,obj.detection.remapChan(1)) = origImg(:,:,obj.detection.morphChannel);
     img(:,:,obj.detection.remapChan(2)) = origImg(:,:,obj.detection.synChannel) ...
		                     + img(:,:,obj.detection.remapChan(2));
     img(:,:,obj.detection.remapChan(3)) = origImg(:,:,obj.detection.XChannel) ...
		                     + img(:,:,obj.detection.remapChan(3));

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Returns the properly scaled image

  function img = getImg(obj)

    img = obj.getRemappedImage();

    if(obj.data.maxRed > 0 & obj.dispInfo.showRed)
      img(:,:,1) = min(img(:,:,1) / double(obj.data.maxRed) * obj.dispInfo.scaleRed,1);
    else
      img(:,:,1) = 0;
    end

    if(obj.data.maxGreen > 0 & obj.dispInfo.showGreen)
      img(:,:,2) = min(img(:,:,2) / double(obj.data.maxGreen) * obj.dispInfo.scaleGreen,1);
    else
      img(:,:,2) = 0;
    end

    if(obj.data.maxBlue > 0 & obj.dispInfo.showBlue)
      img(:,:,3) = min(img(:,:,3) / double(obj.data.maxBlue) * obj.dispInfo.scaleBlue,1);
    else
      img(:,:,3) = 0;
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function showImage(obj)

    if(isempty(obj.data.image))
      % Nothing loaded yet, so nothing to display
      return;
    end

    set(obj.handles.fig,'CurrentAxes',obj.handles.image)    
    oldAxis = axis;

    switch(obj.dispInfo.state)
      case 'load'

        % Just display the raw image...
        tmp = obj.getImg();
        imagesc(tmp);
        axis equal

        set(obj.handles.num,'String', ...
                          sprintf('Image %d (%d)', obj.dispInfo.curImg, obj.data.num))
  
        obj.showHist();

        set(obj.handles.fig,'CurrentAxes',obj.handles.image)    

      case 'soma'

        obj.showMasks();

        if(obj.data.measurePoint)   
          hold on
          [yP,xP] = ind2sub(size(obj.data.somaMask),obj.data.measurePoint);
          contour(obj.data.somaMask,1,'color',[0.3 0.3 0.3], ...
                  'linewidth',2)
          p = plot(xP,yP,'ok');
          set(p,'UIContextMenu',obj.handles.measurePointMenu);

          hold off
        end

      case 'neurites'

        obj.showMasks();

      case 'synapses'

        obj.showMasks();
        hold on
        contour(obj.data.neuriteMask,1,'color',[0.3 0.3 0.3])
        hold off

      case 'analyse'

        % Just display the raw image...
        set(obj.handles.fig,'CurrentAxes',obj.handles.image)    
        tmp = obj.getImg();
        imagesc(tmp);
        axis equal

    end

    if(isempty(obj.dispInfo.axis))
      obj.dispInfo.axis = axis;
    else
      obj.dispInfo.axis = oldAxis;
      set(obj.handles.fig,'CurrentAxes',obj.handles.image)    
      axis(oldAxis);
    end


    drawnow

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function saveUndo(obj, description)
  
    if(isempty(obj.data.somaMask))
      % There is no soma detected, no point adding undo
      return;
    end

    clear tmp
    tmp.somaMask = obj.data.somaMask;
    tmp.neuriteMask = obj.data.neuriteMask;
    tmp.synapseMask = obj.data.synapseMask;
    tmp.includeMask = obj.data.includeMask;
    tmp.synapseCenter = obj.data.synapseCenter;
    tmp.measurePoint = obj.data.measurePoint;

    tmp.description = description;
    tmp.state = obj.dispInfo.state;
    tmp.stage = obj.dispInfo.stage;


    obj.editInfo.undo(end+1) = tmp;

    set(obj.handles.menuItemUndo, 'Label', ...
                      sprintf('<html><u>U</u>ndo %s</html>', obj.editInfo.undo(end).description))

    if(length(obj.editInfo.undo) > obj.editInfo.maxUndo)
      % disp('Max undo length reached, removing old history')
      obj.editInfo.undo(1) = [];
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function clearUndo(obj)

    obj.editInfo.undo = struct('somaMask', [], ...
                               'neuriteMask', [], ...
                               'synapseMask', [], ...
                               'includeMask', [], ...
                               'synapseCenter', [], ...
                               'measurePoint', [], ...
                               'description', [], ...
                               'state', [], ...
                               'stage', 0);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function undoLastAction(obj, source, event)

    if(length(obj.editInfo.undo) > 0)

      obj.data.somaMask = obj.editInfo.undo(end).somaMask;
      obj.data.neuriteMask = obj.editInfo.undo(end).neuriteMask;
      obj.data.synapseMask = obj.editInfo.undo(end).synapseMask;
      obj.data.includeMask = obj.editInfo.undo(end).includeMask;
      obj.data.synapseCenter = obj.editInfo.undo(end).synapseCenter;
      obj.data.measurePoint = obj.editInfo.undo(end).measurePoint;
      obj.makeMeasureMask();

      obj.dispInfo.state = obj.editInfo.undo(end).state;
      obj.dispInfo.stage = obj.editInfo.undo(end).stage;

      obj.editInfo.undo(end) = [];

      obj.showGUI(obj.dispInfo.state);

    end

    if(length(obj.editInfo.undo) > 0)
      set(obj.handles.menuItemUndo, 'Label', ...
                        sprintf('<html><u>U</u>ndo %s</html>', obj.editInfo.undo(end).description))
    else
      set(obj.handles.menuItemUndo,'Label','<html><u>U</u>ndo: No history</html>')
    end

    obj.showImage();
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function excludeRegion(obj, source, event)

    obj.saveUndo('add exclude region');

    obj.stopEdit();

    obj.excludeMask = roipoly();

    obj.data.includeMask = obj.data.includeMask & ~excludeMask;

    set(obj.handles.fig,'CurrentAxes',obj.handles.image)    
    imagesc(~obj.data.includeMask);
    pause(2);
    obj.detectSoma();
    obj.dispInfo.axis = []; % Reset view
    obj.showImage();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function clearExcludeRegions(obj, source, event)

    obj.saveUndo('reset exclude region');
    obj.data.includeMask = ones(size(obj.data.neuriteMask));

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function steerableFilter(obj)

    % Parts of this function provided by Matthew Down
    % as well as GaussianDerivatives3D and nms.

    img = obj.data.image(:,:,obj.detection.morphChannel,obj.dispInfo.curImg);

    [g0x g0y g0z g1x g1y g1z g2x g2y g2z] = ...
        obj.GaussianDerivatives3D(obj.detection.filterSize/obj.data.xyRes);

    N = 4/(obj.detection.filterSize/obj.data.xyRes)^4;
    Rxx = -N*conv2(conv2(img, g2x, 'same'), g0y, 'same');
    Rxy = -N*conv2(conv2(img, g1x, 'same'), g1y, 'same');
    Ryy = -N*conv2(conv2(img, g0x, 'same'), g2y, 'same');

    Rg = N*conv2(conv2(img, g0x, 'same'), g0y, 'same');

    % We use inverted 2nd derivative of gaussian,
    % Theory behind Meijering 2003 uses non-inverted, hence minus
    % since eigenProperties based on Meijering

    if(obj.detection.steerableFilterAlpha)
      alpha = obj.detection.steerableFilterAlpha;

      [eigenValues, eigenVectors, neuriteDir, ridgeMask] = ...
        obj.eigenProperties(Rxx+alpha*Ryy,Ryy+alpha*Rxx,Rxy*(1-alpha)); 
    else
      % Default case, no elongated filters
      [eigenValues, eigenVectors, neuriteDir, ridgeMask] = ...
        obj.eigenProperties(Rxx,Ryy,Rxy); 
    end

    obj.data.dirVect = squeeze(eigenVectors(:,:,1,:));
    obj.data.rigidity = squeeze(eigenValues(:,:,2));
    obj.data.lambdaMax = max(obj.data.rigidity(:));

    clear Rxx Rxy Ryy Rg

    % Debug figures
    if(0)
      % This is the angle that gives the maximal response from the filter
      Gmax = 0.5*atan2(2*Rxy,(Ryy-Rxx)); 

      % This is the value of the maximal response
      Rmax = Rxx/2+Ryy/2+(Rxx/2-Ryy/2).*cos(2*Gmax)-Rxy.*sin(2*Gmax); 

      % Dividing the sterable filter response by the gaussian respone. 
      % This calculation is optional but gives better results often
      % Rmax = Rmax./Rg; 

      % Negative y-axis due to image 0,0 being top left
      GX = cos(Gmax); 
      GY = -sin(Gmax); 

      figure, imagesc(Gmax)
      figure, imagesc(ridgeMask)
      figure, imagesc(neuriteDir)
      colorbar

      [xC,yC] = meshgrid(1:obj.data.width,1:obj.data.height);
      xV = squeeze(eigenVectors(:,:,1,1));
      yV = squeeze(eigenVectors(:,:,1,2));
      maskIdx = find(ridgeMask);

      hold on
      quiver(xC(maskIdx),yC(maskIdx), ...
	     xV(maskIdx),yV(maskIdx), ...
	     'color',[1 1 1]);

      %keyboard
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [g0x g0y g0z g1x g1y g1z g2x g2y g2z] = GaussianDerivatives3D(obj,sigma)

    W = 2*sqrt(2)*sigma;
    R = -W:W;
    k = 1/sigma^2;

    vecsize = size(R,2);

    g0 = exp(-k*R.^2);   
    g0y = reshape(g0, vecsize, 1, 1); 
    g0x = reshape(g0, 1, vecsize, 1); 
    g0z = reshape(g0, 1, 1, vecsize);

    g1 = -2*k*R.*exp(-k*R.^2); % Added minus sign compared to Meijering
    g1y = reshape(g1, vecsize, 1, 1); 
    g1x = reshape(g1, 1, vecsize, 1); 
    g1z = reshape(g1, 1, 1, vecsize);

    g2 = 2*k*(2*k*R.^2-1).*exp(-k*R.^2);
    g2y = reshape(g2, vecsize, 1, 1); 
    g2x = reshape(g2, 1, vecsize, 1); 
    g2z = reshape(g2, 1, 1, vecsize);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Calculate eigenvalues and eigenvectors for the filters

  function [eigenValues, eigenVectors, neuriteDir, ridgeMask] = ...
    eigenProperties(obj, Rxx,Ryy,Rxy)

    eigenValues = zeros(size(Rxx,1),size(Rxx,2),2);
    eigenVectors = zeros(size(Rxx,1),size(Rxx,2),2,2);
    neuriteDir = zeros(size(Rxx));

    if(1)
      % Homemade function to calculate eigenvalues and eigenvectors
      [eigenValues, eigenVectors] = obj.calcEigen2x2(Rxx,Rxy,Rxy,Ryy);
      neuriteDir = -atan(eigenVectors(:,:,1,2)./eigenVectors(:,:,1,1));
    else
      % Using built in functions, slower due to doing them one at a time
      for i = 1:size(Rxx,1)
        for j = 1:size(Rxx,2)
  	  H = [Rxx(i,j), Rxy(i,j); Rxy(i,j), Ryy(i,j)];
          [V,D] = eig(H);
          eigenValues(i,j,:) = [D(1,1),D(2,2)];
          eigenVectors(i,j,1,:) = V(:,1);
          eigenVectors(i,j,2,:) = V(:,2);

          neuriteDir(i,j) = -atan(eigenVectors(i,j,1,2) ...
				  /eigenVectors(i,j,1,1));
        end
      end    
    end

    ridgeMask = max(eigenValues,[],3) == max(abs(eigenValues),[],3);
    % ridgeMask = (eigenValues(:,:,2) < 0) % !!! This one should work, hmmm
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function cc = connectCost(obj, pointA, pointB)

    % Taken from Meijering et al 2003
    % Note that our filter is inverted.

    % We use > 0 instead of < 0, since our filter is minus that one
    % used in Meijering
    if(obj.data.rigidity(pointB(2),pointB(1)) > 0)
      ccLambda = 1 - obj.data.rigidity(pointB(2),pointB(1)) / obj.data.lambdaMax;
    else
      ccLambda = 1;
    end

    AB = pointB-pointA;
    dirAB = AB / norm(AB);

    tmpDirA = [obj.data.dirVect(pointA(2),pointA(1),1), ...
	        obj.data.dirVect(pointA(2),pointA(1),2)];

    wA = tmpDirA / norm(tmpDirA);
 
    tmp = abs(sum(dirAB.*wA));
    ccV = 0.5*(sqrt(1-tmp) + sqrt(1+tmp));

    cc = obj.detection.connectCostLambda * ccLambda ...
          + (1-obj.detection.connectCostLambda)*ccV;

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function growNeuriteFromSoma(obj)

    % Start by seeding neurite mask with soma mask
    obj.data.neuriteMask = obj.data.somaMask;

    obj.saveUndo('detect neurites')
    obj.growNeurite(obj.data.somaMask);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function growNeurite(obj, neuriteSeed)

    %profile clear
    %profile on

    set(obj.handles.calculateLabel, 'String', 'Locating neurites...');
    drawnow

    tic

    somaIdx = find(neuriteSeed);

    obj.detection.pixelQue = somaIdx;

    loopCtr = 0;

    while(length(obj.detection.pixelQue) > 0)      
      
      % Optimized this line by doing calculation explicit, see below
      %[yCenter,xCenter] = ind2sub(size(obj.data.neuriteMask), ...
      %				  obj.detection.pixelQue(1));

      xCenter = floor((obj.detection.pixelQue(1)-1) / size(obj.data.neuriteMask,1))+1;
      yCenter = mod(obj.detection.pixelQue(1)-1,size(obj.data.neuriteMask,1))+1;

      % Find valid neighours
      [neighX, neighY] = obj.processNeighbourhood(xCenter,yCenter);

      % Make sure they are not already qued

      % Optimized this line by calculating manually...
      % neighIdx = sub2ind(size(obj.data.neuriteMask),neighY,neighX);
      neighIdx = neighY+(neighX-1)*size(obj.data.neuriteMask,1);

      obj.detection.pixelQue = [obj.detection.pixelQue(2:end); neighIdx];

      % Mark them in the mask
      obj.data.neuriteMask(neighIdx) = 1;      

      loopCtr = loopCtr + 1;

      if(mod(loopCtr,5000) == 0)
        % figure, imagesc(obj.data.neuriteMask), axis equal, drawnow
        fprintf('%d iterations\n', loopCtr)
        set(obj.handles.calculateLabel, ...
            'String', sprintf('Locating neurites... (%d)', loopCtr));
        drawnow
      end
    end
 
    % Remove spurious pixels around the neurite
    se = strel('disk',1);
    obj.data.neuriteMask = imopen(obj.data.neuriteMask,se);

    fprintf('%d iterations, done.\n', loopCtr)

    obj.updateNeuriteMask();

    obj.data.distMask = obj.makeDistMask(obj.data.somaMask,obj.data.neuriteMask);

    toc

    set(obj.handles.calculateLabel, 'String', '');

    obj.showImage();

    %profview

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [neighX, neighY] = processNeighbourhood(obj, xCenter, yCenter)

    nHoodX = transpose([1 1  1 0  0 -1 -1 -1] + xCenter);
    nHoodY = transpose([1 0 -1 1 -1  1  0 -1] + yCenter);

    if(xCenter == 1 | xCenter == obj.data.width ...
       | yCenter == 1 | yCenter == obj.data.height)

      % Only if we are at border do this more time consuming check

      validIdx = find(1 <= nHoodX & nHoodX <= obj.data.width ...
		      & 1 <= nHoodY & nHoodY <= obj.data.height);

      nHoodX = nHoodX(validIdx);
      nHoodY = nHoodY(validIdx);
    end

    nHoodMask = zeros(size(nHoodX));

    for iN = 1:length(nHoodX)

      % My idea was to check if connectCost was lower than a certain value
      % then if that is the case, the neighbourhood pixels are added
      % and returned. The function that calls this function then adds these
      % pixels to obj.detection.pixelQue (if they are not already part of
      % neurite mask.

     if(obj.data.neuriteMask(nHoodY(iN),nHoodX(iN)) == 0 ...
         & obj.connectCost([xCenter yCenter], [nHoodX(iN) nHoodY(iN)]) ...
        < obj.detection.maxAddCost)
	 

        % disp('Adding pixel')
        nHoodMask(iN) = 1;    
      % else
      %  c = obj.connectCost([xCenter yCenter], [nHoodX(iN) nHoodY(iN)])
      %  n = obj.data.neuriteMask(nHoodY(iN),nHoodX(iN))
      end
 
    end

    neighX = nHoodX(find(nHoodMask));
    neighY = nHoodY(find(nHoodMask));

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [eigenVal, eigenVec] = calcEigen2x2(obj, a11,a12,a21,a22)

    eigenVal = zeros(size(a11,1),size(a11,2),2);
    eigenVec = zeros(size(a11,1),size(a11,2),2,2);

    tmp1 = (a11+a22)/2;
    tmp2 = sqrt(tmp1.^2 + (a12.*a21-a11.*a22));

    eigenVal(:,:,1) = tmp1-tmp2;
    eigenVal(:,:,2) = tmp1+tmp2;

    eigenVec(:,:,1,1) = 1;
    eigenVec(:,:,1,2) = (eigenVal(:,:,1)-a11)./a12;

    eigenVec(:,:,2,1) = 1;
    eigenVec(:,:,2,2) = a21./(eigenVal(:,:,2)-a22);

    % Normalise
    tmp = sqrt(eigenVec(:,:,1,1).^2+eigenVec(:,:,1,2).^2);
    eigenVec(:,:,1,1) = eigenVec(:,:,1,1)./tmp;
    eigenVec(:,:,1,2) = eigenVec(:,:,1,2)./tmp;

    tmp = sqrt(eigenVec(:,:,2,1).^2+eigenVec(:,:,2,2).^2);
    eigenVec(:,:,2,1) = eigenVec(:,:,2,1)./tmp;
    eigenVec(:,:,2,2) = eigenVec(:,:,2,2)./tmp;

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function calculateNeuriteLength(obj)

    skelIdx = find(obj.data.skeleton);

    tmpSkel = obj.data.skeleton;

    distTotal = 0;

    for i = 1:length(skelIdx)
      [yC,xC] = ind2sub(size(obj.data.skeleton),skelIdx(i));

      % Clipping so we do not get outside image
      yNeigh = min(max(yC-1:yC+1,1),size(obj.data.skeleton,1));
      xNeigh = min(max(xC-1:xC+1,1),size(obj.data.skeleton,2));

      for yN = yNeigh
        for xN = xNeigh
          if(tmpSkel(yN,xN))
	    distTotal = distTotal + sqrt((yC-yN)^2+(xC-xN)^2);
          end
        end
      end

      % Remove the pixel
      tmpSkel(yC,xC) = 0;

    end

    obj.data.neuriteLength = distTotal * obj.data.xyRes;

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function trimSkeleton(obj)

    branchPoints = bwmorph(obj.data.skeleton,'branchpoints');    
    endPoints = bwmorph(obj.data.skeleton,'endpoints');
    endPointIdx = find(endPoints);

    % We want to remove all the mini-branches smaller than the trim size.
    % To find them we remove the branch points, then calculate the locate
    % all the connected components smaller than the trim size, and if
    % they only have one branch point around them (imdilate 1), 
    % then we remove it.

    trimMask = double(obj.data.skeleton - bwmorph(branchPoints,'dilate') > 0);

    mCC = bwconncomp(trimMask);

    for i = 1:length(mCC.PixelIdxList)
      if(length(mCC.PixelIdxList{i}) < obj.detection.trimNeuriteSize)
        % Does it contain an endpoint?
        if(nnz(ismember(mCC.PixelIdxList{i},endPointIdx)))
          [y,x] = ind2sub(size(obj.data.skeleton),mCC.PixelIdxList{i});
          yAll = max(1,min([y+1;y+1;y+1;y;y;y;y-1;y-1;y-1],obj.data.height));
          xAll = max(1,min([x+1;x;x-1;x+1;x;x-1;x+1;x;x-1],obj.data.width));
          trimIdx = sub2ind(size(obj.data.skeleton),yAll,xAll);

          %fprintf('Trimming away %d pixels.\n', length(mCC.PixelIdxList{i}))

          % Yes it does, remove pixel from skeleton
          %obj.data.skeleton(mCC.PixelIdxList{i}) = 0;
          obj.data.skeleton(trimIdx) = 0;
        end    
      end
    end   

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Calculates the sholl analysis for a skeleton of the neurites
  function shollAnalysis(obj)

    obj.data.skeleton = bwmorph(obj.data.neuriteMask-obj.data.somaMask>0,'skel',inf);

    % This give the centroid of all pixels included in all somas
    somaProps = regionprops(obj.data.somaMask,'centroid');

    % Prune away tiny pertrusions, mostly artifacts
    obj.trimSkeleton();

    CC = bwconncomp(obj.data.somaMask);

    if(CC.NumObjects ~= 1)
      fprintf('Detected %d somas, aborting Sholl analysis.\n', ...
	      length(somaProps))
      obj.exportInfo.saveSholl = 0;

      % Save histogram with NaN
      maxDist = sqrt(obj.data.width^2+obj.data.height^2)*obj.data.xyRes;
      obj.data.shollEdges = 0:obj.detection.shollBinSize:maxDist;
      obj.data.shollDendHist = NaN*zeros(1,length(obj.data.shollEdges)-1);

      obj.data.shollIntMorphMean = NaN*zeros(1,length(obj.data.shollEdges)-1);
      obj.data.shollIntMorphStd = NaN*zeros(1,length(obj.data.shollEdges)-1);
      obj.data.shollIntMorphSEM = NaN*zeros(1,length(obj.data.shollEdges)-1);
      obj.data.shollIntSynMean = NaN*zeros(1,length(obj.data.shollEdges)-1);
      obj.data.shollIntSynStd = NaN*zeros(1,length(obj.data.shollEdges)-1);
      obj.data.shollIntSynSEM = NaN*zeros(1,length(obj.data.shollEdges)-1);
      obj.data.shollIntXMean = NaN*zeros(1,length(obj.data.shollEdges)-1);
      obj.data.shollIntXStd = NaN*zeros(1,length(obj.data.shollEdges)-1);
      obj.data.shollIntXSEM = NaN*zeros(1,length(obj.data.shollEdges)-1);

      obj.writeExportSettings();

      obj.data.skeleton = bwmorph(obj.data.neuriteMask-obj.data.somaMask>0,'skel',inf);
      obj.trimSkeleton();

      obj.calculateNeuriteLength();

      return
    else
      obj.exportInfo.saveSholl = 1;
      obj.writeExportSettings();
    end

    % This must be modified to handle multiple somas
    somaCenter = somaProps.Centroid;

    [y,x] = find(obj.data.somaMask);
    somaRadius = max(sqrt((x-somaCenter(1)).^2 + (y-somaCenter(2)).^2)) ...
		     * obj.data.xyRes;

    [y,x] = find(obj.data.skeleton);
    distSkeleton = sqrt((x-somaCenter(1)).^2 + (y-somaCenter(2)).^2) ...
			     * obj.data.xyRes;

    distToSoma = zeros(size(obj.data.neuriteMask));
    for i = 1:length(x)
      % We want to the soma radius to be defined as distance 0
      distToSoma(y(i),x(i)) = distSkeleton(i) - somaRadius;
    end

    maxDist = sqrt(obj.data.width^2+obj.data.height^2)*obj.data.xyRes;
    %obj.data.shollEdges = (round(somaRadius/5e-6)*5e-6+5e-6)...
    %			 :obj.detection.shollBinSize:maxDist;

    obj.data.shollEdges = 0:obj.detection.shollBinSize:maxDist;
    obj.data.shollDendHist = zeros(1,length(obj.data.shollEdges)-1);
    shollHalfWidth = sqrt(2)*obj.data.xyRes/2;

    % Soma radius defined as distance 0
    [ySyn,xSyn] = ind2sub(size(obj.data.synapseMask),obj.data.synapseCenter);
    synDist = sqrt((xSyn-somaCenter(1)).^2+(ySyn-somaCenter(2)).^2) ...
                * obj.data.xyRes - somaRadius;

    obj.data.shollIntMorphMean = zeros(1,length(obj.data.shollEdges)-1);
    obj.data.shollIntMorphStd = zeros(1,length(obj.data.shollEdges)-1);
    obj.data.shollIntMorphSEM = zeros(1,length(obj.data.shollEdges)-1);
    obj.data.shollIntSynMean = zeros(1,length(obj.data.shollEdges)-1);
    obj.data.shollIntSynStd = zeros(1,length(obj.data.shollEdges)-1);
    obj.data.shollIntSynSEM = zeros(1,length(obj.data.shollEdges)-1);
    obj.data.shollIntXMean = zeros(1,length(obj.data.shollEdges)-1);
    obj.data.shollIntXStd = zeros(1,length(obj.data.shollEdges)-1);
    obj.data.shollIntXSEM = zeros(1,length(obj.data.shollEdges)-1);


    for i = 1:length(obj.data.shollEdges)-1
      % Neurite Sholl analysis
      BW = (obj.data.shollEdges(i) <= distToSoma) ...
	    & (distToSoma < obj.data.shollEdges(i+1));
      CC = bwconncomp(BW);
      obj.data.shollDendHist(i) = length(CC.PixelIdxList);

      % Synapse 
      synIdx = find(obj.data.shollEdges(i) <= synDist ...
		    & synDist < obj.data.shollEdges(i+1));

      if(~isempty(synIdx))
        obj.data.shollIntMorphMean(i) = mean(obj.data.synapseIntensityMorphMean(synIdx));
        obj.data.shollIntMorphStd(i) = std(obj.data.synapseIntensityMorphMean(synIdx));
        obj.data.shollIntMorphSEM(i) = obj.data.shollIntMorphStd(i) /sqrt(length(synIdx));

        obj.data.shollIntSynMean(i) = mean(obj.data.synapseIntensitySynMean(synIdx));
        obj.data.shollIntSynStd(i) = std(obj.data.synapseIntensitySynMean(synIdx));
        obj.data.shollIntSynSEM(i) = obj.data.shollIntSynStd(i) /sqrt(length(synIdx));

        obj.data.shollIntXMean(i) = mean(obj.data.synapseIntensityXMean(synIdx));
        obj.data.shollIntXStd(i) = std(obj.data.synapseIntensityXMean(synIdx));
        obj.data.shollIntXSEM(i) = obj.data.shollIntXStd(i) /sqrt(length(synIdx));
      else
        % NaN fields will be left empty when exporting to XML
        obj.data.shollIntMorphMean(i) = NaN;
        obj.data.shollIntMorphStd(i) = NaN;
        obj.data.shollIntMorphSEM(i) = NaN;

        obj.data.shollIntSynMean(i) = NaN;
        obj.data.shollIntSynStd(i) = NaN;
        obj.data.shollIntSynSEM(i) = NaN;

        obj.data.shollIntXMean(i) = NaN;
        obj.data.shollIntXStd(i) = NaN;
        obj.data.shollIntXSEM(i) = NaN;
       end

    end

    % obj.data.synapseDist contains the arc length
    % synDist is the shortest distance from soma to synapse...

    obj.data.shollSynHist = histc(synDist,obj.data.shollEdges);
    obj.data.shollSynHist = obj.data.shollSynHist(1:end-1);

    obj.calculateNeuriteLength();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function calculateIntensityGradients(obj)

    if(isempty(obj.data.distMask))
      obj.data.distMask = obj.makeDistMask(obj.data.somaMask,obj.data.neuriteMask);
    end

    tmp = obj.data.distMask(:);
    tmp(find(tmp == inf)) = NaN;
    maxDist = max(tmp(:));

    obj.data.gradientEdges = ...
      0:obj.detection.shollBinSize:(maxDist+obj.detection.shollBinSize);

    obj.data.morphGradientMean = NaN*zeros(size(obj.data.gradientEdges));
    obj.data.morphGradientStd = NaN*zeros(size(obj.data.gradientEdges));
    obj.data.morphGradientSEM = NaN*zeros(size(obj.data.gradientEdges));

    obj.data.synGradientMean = NaN*zeros(size(obj.data.gradientEdges));
    obj.data.synGradientStd = NaN*zeros(size(obj.data.gradientEdges));
    obj.data.synGradientSEM = NaN*zeros(size(obj.data.gradientEdges));

    obj.data.XGradientMean = NaN*zeros(size(obj.data.gradientEdges));
    obj.data.XGradientStd = NaN*zeros(size(obj.data.gradientEdges));
    obj.data.XGradientSEM = NaN*zeros(size(obj.data.gradientEdges));

    tmpMorph = obj.data.image(:,:,obj.detection.morphChannel, obj.dispInfo.curImg);
    tmpSyn = obj.data.image(:,:,obj.detection.synChannel, obj.dispInfo.curImg);
    tmpX = obj.data.image(:,:,obj.detection.XChannel, obj.dispInfo.curImg);

    for i = 1:length(obj.data.gradientEdges)
      if(i == 1)
        idx = find(obj.data.distMask == 0);
      else
        idx = find(obj.data.gradientEdges(i-1) < obj.data.distMask ...
		   & obj.data.distMask <= obj.data.gradientEdges(i));
      end

      obj.data.morphGradientMean(i) = mean(tmpMorph(idx));
      obj.data.morphGradientStd(i) = std(tmpMorph(idx));
      obj.data.morphGradientSEM(i) = mean(tmpMorph(idx))/sqrt(length(idx));

      obj.data.synGradientMean(i) = mean(tmpSyn(idx));
      obj.data.synGradientStd(i) = std(tmpSyn(idx));
      obj.data.synGradientSEM(i) = std(tmpSyn(idx))/sqrt(length(idx));

      obj.data.XGradientMean(i) = mean(tmpX(idx));
      obj.data.XGradientStd(i) = std(tmpX(idx));
      obj.data.XGradientSEM(i) = std(tmpX(idx))/sqrt(length(idx));

    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function intensityHistograms(obj)

    maxInt = max([obj.data.maxRed,obj.data.maxGreen,obj.data.maxBlue]);

    obj.data.intensityHistogramEdges = ...
      0:obj.detection.intensityBinSize:(maxInt+obj.detection.intensityBinSize-1);

    obj.data.morphHist = histc(obj.data.synapseIntensityMorphMean, ...
                           obj.data.intensityHistogramEdges);

    obj.data.synHist = histc(obj.data.synapseIntensitySynMean, ...
                         obj.data.intensityHistogramEdges);

    obj.data.XHist = histc(obj.data.synapseIntensityXMean, ...
		       obj.data.intensityHistogramEdges);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function readExportSettings(obj)

    obj.exportInfo.saveMask = get(obj.handles.saveMask,'Value');
    obj.exportInfo.saveSholl = get(obj.handles.saveSholl,'Value');
    obj.exportInfo.saveIntensityHistogram = get(obj.handles.saveIntensityHist,'Value');
    obj.exportInfo.saveMat = get(obj.handles.saveMat,'Value');
    obj.exportInfo.saveSomaMeasure = get(obj.handles.saveSomaMeasure,'Value');
    obj.exportInfo.saveSynapseProfile = get(obj.handles.saveSynapseProfile,'Value');
    obj.exportInfo.saveGradient = get(obj.handles.saveGradient,'Value');

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function writeExportSettings(obj)

    set(obj.handles.saveMask,'Value',obj.exportInfo.saveMask);
    set(obj.handles.saveSholl,'Value',obj.exportInfo.saveSholl);
    set(obj.handles.saveIntensityHist,'Value',obj.exportInfo.saveIntensityHistogram);
    set(obj.handles.saveMat,'Value',obj.exportInfo.saveMat);
    set(obj.handles.saveSomaMeasure,'Value',obj.exportInfo.saveSomaMeasure);
    set(obj.handles.saveSynapseProfile,'Value',obj.exportInfo.saveSynapseProfile);
    set(obj.handles.saveGradient,'Value',obj.exportInfo.saveGradient);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function exportData(obj, source, event)

    if(exist('source'))
      % The export function was called from a button, ask for save file
      obj.setExportFiles();
    end

    if(isempty(obj.data.exportXMLfile))
       % User aborted save.
      return
    end

    if(~obj.detection.reanalyseFlag)
      progressBar = waitbar(0,sprintf('Exporting %s...', ...
				      strrep(obj.data.fileName{obj.dispInfo.curImg}, ...
					     '_','\_')));
    else
      % Do not show a progress bar if reanalysing, to prevent focus stealing
      progressBar = [];
    end

    barFig = gcf;
    nSteps = 6;

    obj.readExportSettings();

    % Calculate obj.data...
    if(obj.exportInfo.saveSholl)
      obj.shollAnalysis();
      obj.exportWaitBar(1/nSteps, progressBar);
    else
      % Sholl normally calculates skeleton, but if we do not do that
      % lets do it here.
      
      obj.data.skeleton = bwmorph(obj.data.neuriteMask-obj.data.somaMask>0,'skel',inf);
      obj.trimSkeleton();

      obj.calculateNeuriteLength();

    end

    if(obj.exportInfo.saveIntensityHistogram)
      obj.intensityHistograms();
      obj.exportWaitBar(2/nSteps, progressBar);
    end

    obj.calculateSomaMeasures();
    obj.exportWaitBar(3/nSteps, progressBar);

    if(obj.exportInfo.saveGradient)
      obj.calculateIntensityGradients();
    end

    % exportMasks();

    obj.exportToXML();
    obj.exportWaitBar(4/nSteps, progressBar);

    % Also save the data in matlab readable format, just in case
    if(obj.exportInfo.saveMat)
      obj.saveData();
    end
    obj.exportWaitBar(5/nSteps, progressBar);

    if(obj.exportInfo.saveMask & ~isempty(progressBar))
      % Only need to return focus if we had a progressbar
      figure(obj.handles.fig)
      obj.saveMask();
    end

    if(~isempty(progressBar))
      % Only required if we have a progress bar
      figure(barFig)			    
    end
    obj.exportWaitBar(6/nSteps, progressBar);

    % Increase the saved counter
    obj.detection.numSaved = obj.detection.numSaved + 1;

    obj.saveConfig();

    delete(progressBar);
    % uiwait(helpdlg('All done. Tack för idag!', 'SynD'));

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setExportFiles(obj)

    exportFile = strcat(obj.data.fileName{obj.dispInfo.curImg}, ...
			'-export.xml');

    curPwd = pwd;
    try
      cd(obj.data.exportPath)
    catch
      fprintf('Unable to change to %s\n', obj.data.exportPath)
    end

    [obj.data.exportXMLfile, xmlExportPath] = ...
	  uiputfile('*-export.xml', ...
		    'Export synapse data', ...
		    exportFile);

    cd(curPwd);

    if(isempty(obj.data.exportXMLfile) | obj.data.exportXMLfile == 0)
      obj.data.exportXMLfile = [];
      disp('Export aborted.')
      return
    end

    obj.data.exportPath = xmlExportPath;

    if(strfind(obj.data.exportXMLfile,'-export.xml'))
      obj.data.exportSaveFile = ...
          strrep(obj.data.exportXMLfile,'-export.xml','-save.mat');
      obj.data.exportNeuriteMaskFile = ...
          strrep(obj.data.exportXMLfile,'-export.xml','-neuritemask.tiff');
      obj.data.exportSynapseMaskFile = ...
          strrep(obj.data.exportXMLfile,'-export.xml','-synapsemask.tiff');
    else
      obj.data.exportSaveFile = ...
          strcat(obj.data.exportXMLfile,'-save.mat');
      obj.data.exportNeuriteMaskFile = ...
          strcat(obj.data.exportXMLfile,'-neuritemask.tiff');
      obj.data.exportSynapseMaskFile = ...
        strcat(obj.data.exportXMLfile,'-synapsemask.tiff');
    end


  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function saveMask(obj)

    %%% Save neurite and soma masks

    % Get morphology image with true morphology colour
    obj.setChannelVisible(1,0,0); 
    morphImg = obj.getImg();

    % Show neurite mask with synapse colour
    morphImg(:,:,obj.detection.synChannel) = obj.data.neuriteMask;

    % Show soma with X-channel colour
    morphImg(:,:,obj.detection.XChannel) = obj.data.somaMask*0.5;

    imwrite(morphImg,strcat(obj.data.exportPath,obj.data.exportNeuriteMaskFile), ...
	    'tif', 'compression','lzw');

    %%% Save synapse picture

    obj.setChannelVisible(0,1,0); 
    synImg = obj.getImg();

    % Use morph channel to mark synapses
    synImg(:,:,obj.detection.morphChannel) = 0.5*obj.data.synapseMask;

    % Use X-channel to mark synapse centres
    tmp = zeros(obj.data.height,obj.data.width);
    tmp(obj.data.synapseCenter) = 1;
    synImg(:,:,obj.detection.XChannel) = tmp;

    imwrite(synImg,strcat(obj.data.exportPath,obj.data.exportSynapseMaskFile), ...
	    'tif', 'compression','lzw');


    obj.setChannelVisible(1,1,1); 

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function exportWaitBar(obj, progress,progressBar)
    if(~isempty(progressBar))
      waitbar(progress,progressBar);
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function saveData(obj)

    clear old

    old.image = obj.data.image;

    old.xyRes = obj.data.xyRes;
    old.fileName = obj.data.fileName;
    old.num = obj.data.num;

    old.somaMask = obj.data.somaMask;
    old.neuriteMask = obj.data.neuriteMask;
    old.synapseMask = obj.data.synapseMask;

    old.synapseCenter = obj.data.synapseCenter;
    old.synapseArea = obj.data.synapseArea;
    old.synapseDist = obj.data.synapseDist;

    old.meanSynapseMorph = obj.data.meanSynapseMorph;
    old.meanSynapseSyn = obj.data.meanSynapseSyn;
    old.meanSynapseX = obj.data.meanSynapseX;
    old.maxRadie = obj.detection.maxRadie;
    old.excludeSomaSynapses = obj.detection.excludeSomaSynapses;

    old.meanSynapseProfileDist = obj.data.meanSynapseProfileDist;
    old.meanSynapseProfileMorph = obj.data.meanSynapseProfileMorph;
    old.meanSynapseProfileSyn = obj.data.meanSynapseProfileSyn;
    old.meanSynapseProfileX = obj.data.meanSynapseProfileX;

    old.measurePoint = obj.data.measurePoint;
    old.somaMeasureMask = obj.data.somaMeasureMask;

    old.synapseIntensityMorphMean = obj.data.synapseIntensityMorphMean;
    old.synapseIntensityMorphSEM = obj.data.synapseIntensityMorphSEM;
    old.synapseIntensitySynMean = obj.data.synapseIntensitySynMean;
    old.synapseIntensitySynSEM = obj.data.synapseIntensitySynSEM;
    old.synapseIntensityXMean = obj.data.synapseIntensityXMean;
    old.synapseIntensityXSEM = obj.data.synapseIntensityXSEM;

    old.morphChannel = obj.detection.morphChannel;
    old.synChannel = obj.detection.synChannel;
    old.XChannel = obj.detection.XChannel;
    old.remapChan = obj.detection.remapChan;
    old.filterSize = obj.detection.filterSize;

    old.totDendLength = obj.data.neuriteLength;
    old.intensityHistogramEdges = obj.data.intensityHistogramEdges;
    old.morphHist = obj.data.morphHist;
    old.synHist = obj.data.synHist;
    old.XHist = obj.data.XHist;

    old.somaMeasureMorph = obj.data.somaMeasureMorph;
    old.somaMeasureSyn = obj.data.somaMeasureSyn;
    old.somaMeasureX = obj.data.somaMeasureX;

    old.somaArea = obj.data.somaArea;
    old.somaMajorAxisLength = obj.data.somaMajorAxisLength;
    old.somaMinorAxisLength = obj.data.somaMinorAxisLength;

    old.shollEdges = obj.data.shollEdges;
    old.shollDendHist = obj.data.shollDendHist;
    old.shollSynHist = obj.data.shollSynHist;

    old.shollIntMorphMean = obj.data.shollIntMorphMean;
    old.shollIntMorphStd = obj.data.shollIntMorphStd;
    old.shollIntMorphSEM = obj.data.shollIntMorphSEM;
    old.shollIntSynMean = obj.data.shollIntSynMean;
    old.shollIntSynStd = obj.data.shollIntSynStd;
    old.shollIntSynSEM = obj.data.shollIntSynSEM;
    old.shollIntXMean = obj.data.shollIntXMean;
    old.shollIntXStd = obj.data.shollIntXStd;
    old.shollIntXSEM = obj.data.shollIntXSEM;

    old.gradientEdges = obj.data.gradientEdges;

    old.morphGradientMean = obj.data.morphGradientMean;
    old.morphGradientStd = obj.data.morphGradientStd; 
    old.morphGradientSEM = obj.data.morphGradientSEM;

    old.synGradientMean = obj.data.synGradientMean;
    old.synGradientStd = obj.data.synGradientStd;
    old.synGradientSEM = obj.data.synGradientSEM;

    old.XGradientMean = obj.data.XGradientMean;
    old.XGradientStd = obj.data.XGradientStd;
    old.XGradientSEM = obj.data.XGradientSEM;

    old.dataImg = obj.dispInfo.curImg;
    old.detection = obj.detection;

    old.somaIntensityThreshold = obj.data.somaIntensityThreshold;
    old.synapseIntensityThreshold = obj.data.synapseIntensityThreshold;

    % Save version information also
    old.version = obj.getVersion();

    if(obj.data.exportSaveFile)
      fprintf('Saving data to %s\n', obj.data.exportSaveFile);
      save(strcat(obj.data.exportPath,obj.data.exportSaveFile),'old')
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function importMask(obj, source, event)

    if(isempty(obj.data.image))
      
      disp('Unable to import mask, you must load corresponding lsm file first')
      beep

      % uiwait(warndlg('You must load the corresponding lsm file first.', ...
      %		     'Unable to import masks','modal'));

      return
    end

    curDir = pwd;

    try
      cd(obj.data.exportPath)
    catch
      fprintf('Unable to change to %s\n', obj.data.exportPath)
    end

    % Ask the user which file to load
    [dataFile,dataPath] = uigetfile('*-save.mat','Select old save file');

    obj.data.exportPath = dataPath;
    cd(curDir);

    old = load(strcat(dataPath,dataFile));
    try
      old = old.old;
    catch
      fprintf('Unable to load file %s!\n', dataFile)
      return
    end


    if(~strcmp(old.fileName{1},obj.data.fileName{1}))
      uiwait(warndlg(sprintf('Image file loaded: %s\nMask file loaded: %s', ...
			     obj.data.fileName{1},old.fileName{1}), ...
		     'SynD : Mask import file mismatch'))

      obj.data.somaIntensityThreshold = NaN;
    else

      % Load the old intensity threshold used
      try
        obj.data.somaIntensityThreshold = old.somaIntensityThreshold;
      catch
        obj.data.somaIntensityThreshold = old.obj.detection.morphThreshold;
      end

    end

    obj.saveUndo('loading masks')

    try
      obj.data.neuriteMask = old.neuriteMask;
      obj.data.somaMask = old.somaMask;
    catch
      % Old format had wrong name, revision 162...
      obj.data.neuriteMask = old.neuriteMaskIdx;
      obj.data.somaMask = old.somaMaskIdx;     
    end

    try
      obj.data.somaMeasureMask = old.somaMeasureMask;
      obj.exportInfo.saveSomaMeasure = 1;
      obj.writeExportSettings();
    catch
      disp('No soma measure marked, ignoring.')
      obj.data.somaMeasureMask = [];
    end

    if(obj.data.height ~= size(obj.data.neuriteMask,1) ...
       | obj.data.width ~= size(obj.data.neuriteMask,2))
      disp('Neurite and soma mask size does not match current image')

      % Clear the masks...
      obj.data.neuriteMask = [];
      obj.data.somaMask = [];
      activateStage('soma');
      return
    end

    if(~strcmp(obj.data.fileName{1},old.fileName{1}))
      fprintf('Image file name: %s\nOld mask file name: %s\n', ...
	      obj.data.fileName{1}, old.fileName{1})
      beep
    end

    obj.updateNeuriteMask();

    obj.data.distMask = obj.makeDistMask(obj.data.somaMask,obj.data.neuriteMask);

    obj.activateStage('synapses');

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function exportMasks(obj)
    % !!! To be continued
    disp('exportMasks not implemented yet')
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % This exports to xml, readable by excel

  function exportToXML(obj)

    if(isempty(obj.data.exportXMLfile))
      return
    end

    CC = bwconncomp(obj.data.somaMask);
    nSoma = CC.NumObjects; 

    if(nSoma > 1 & obj.detection.singleSoma ...
       & ~obj.detection.reanalyseFlag)
      % If we are reanalysing we are reusing an old soma mask, so the user
      % should know about this by now...

      uiwait(warndlg(sprintf('Single soma selected, but found %d somas.', ...
			     nSoma), ...
		     'Soma inconsistency error','modal'));
    end


    % Header info, to make excel feel comfortable with our obj.data...
    docNode = obj.makeXMLheader();

    % Create first sheet -- File info

    columnName = {'File name','Pixel size [um]', ...
		  'Morphology channel', ...
		  'Synapse channel', ...
		  'X channel', ...
                  'Number of somas', ...
		  'Soma synapses'};

    if(length(obj.data.fileName) > 1)
      disp('Warning, export only obj.handles one input file');
    end

    if(obj.detection.excludeSomaSynapses)
      somaSynapsesStatus = 'excluded';
    else
      somaSynapsesStatus = 'included';
    end

    columnData = {{obj.data.fileName{obj.dispInfo.curImg}}, ...
		  [obj.data.xyRes*1e6], ...
		  [obj.detection.morphChannel], ...
		  [obj.detection.synChannel], ...
		  [obj.detection.XChannel], ...
                  [nSoma], ...
		  {somaSynapsesStatus}};

    obj.makeXMLsheet(docNode,'File info', ...
                     columnName, columnData);


    % Write second sheet of data, synapse info

    columnName = {'Synapse number', ...
		  'X (micrometers)', ...
		  'Y (micrometers)', ...
		  'Area (micrometers^2)', ...
		  'Morph intensity (mean)', ...
		  'Morph intensity (SEM)', ...
		  'Synapse intensity (mean)', ...
		  'Synapse intensity (SEM)', ...
		  'X intensity (mean)', ...
		  'X intensity (SEM)', ...
		  'Syn/Morph', ... 
		  'X/Morph', ...
		  'X/Syn' };

    [y,x] = ind2sub(size(obj.data.synapseMask),obj.data.synapseCenter);

    % To avoid division by zero we discard those data points
    okMidx = find(obj.data.synapseIntensityMorphMean ~= 0);
    okSidx = find(obj.data.synapseIntensitySynMean ~= 0);

    columnData = {1:length(obj.data.synapseCenter), ...
		  x, y, ...
		  obj.data.synapseArea, ...
		  obj.data.synapseIntensityMorphMean, ...
		  obj.data.synapseIntensityMorphSEM, ...
		  obj.data.synapseIntensitySynMean, ...
		  obj.data.synapseIntensitySynSEM, ...
		  obj.data.synapseIntensityXMean, ...
		  obj.data.synapseIntensityXSEM, ...
		  obj.data.synapseIntensitySynMean(okMidx) ...
		   ./obj.data.synapseIntensityMorphMean(okMidx), ...
		  obj.data.synapseIntensityXMean(okMidx) ...
		   ./obj.data.synapseIntensityMorphMean(okMidx), ...
		  obj.data.synapseIntensityXMean(okSidx) ...
		   ./obj.data.synapseIntensitySynMean(okSidx) };


    obj.makeXMLsheet(docNode,'Raw synapse data', ...
                     columnName, columnData);
    
    columnName = { 'Number of synapses', ...
		   'Total dendritic length (micrometer)', ...
		   'Synapses/micrometer', ...
		   'Mean area (micrometer^2)', ...
		   'Mean morph channel', ...
		   'Mean synapse channel', ...
		   'Mean X channel', ...
		   'Mean syn/morph', ...
		   'Mean X/morph', ...
		   'Mean X/syn' };

    nSynapses = length(obj.data.synapseCenter);
    totLength = obj.data.neuriteLength*1e6;


    columnData = { nSynapses, ...
		   totLength, ...
		   nSynapses/totLength, ...
		   mean(obj.data.synapseArea), ...
		   mean(obj.data.synapseIntensityMorphMean), ...
		   mean(obj.data.synapseIntensitySynMean), ...  
		   mean(obj.data.synapseIntensityXMean), ...  
		   nanmean(obj.data.synapseIntensitySynMean(okMidx) ...
			   ./obj.data.synapseIntensityMorphMean(okMidx)), ...
		   nanmean(obj.data.synapseIntensityXMean(okMidx) ...
			   ./obj.data.synapseIntensityMorphMean(okMidx)), ...
		   nanmean(obj.data.synapseIntensityXMean(okSidx) ...
			   ./obj.data.synapseIntensitySynMean(okSidx)) };

    obj.makeXMLsheet(docNode,'Synapse averages', ...
                     columnName, columnData);

    if(obj.exportInfo.saveIntensityHistogram)

      % Write histogram with channel intensities

      columnName = { 'From', 'To', ...
		     'Morph', 'Syn', 'X', ...
		     'Cum morph', 'Cum syn', 'Cum X' };
				      

      columnData = { obj.data.intensityHistogramEdges(1:end-1), ...
		     obj.data.intensityHistogramEdges(2:end), ...
		     obj.data.morphHist(1:end-1), ...
		     obj.data.synHist(1:end-1), ...
		     obj.data.XHist(1:end-1), ...
		     cumsum(obj.data.morphHist(1:end-1)), ...
		     cumsum(obj.data.synHist(1:end-1)), ...
		     cumsum(obj.data.XHist(1:end-1)) };

      obj.makeXMLsheet(docNode,'Intensity histograms', ...
                       columnName, columnData);
    end


    if(obj.exportInfo.saveSholl)
      % Add Sholl analysis
      columnName = { 'From (mu)', 'To (mu)', ...
		     'Dend count', ...
		     'Synapse count', ...
		     'Morph (mean)', ...
		     'Morph (std)', ...
		     'Morph (SEM)', ...
		     'Syn (mean)', ...
		     'Syn (std)', ...
		     'Syn (SEM)', ...
		     'X (mean)', ...
		     'X (std)', ...
		     'X (SEM)', ...
		     'Mean (syn/morph)', ...
		     'Mean (X/morph)', ...
		     'Mean (X/syn)', ...
                   };

      columnData = { obj.data.shollEdges(1:end-1)*1e6, obj.data.shollEdges(2:end)*1e6, ...
		     obj.data.shollDendHist, ...
		     obj.data.shollSynHist, ...
		     obj.data.shollIntMorphMean, ...
		     obj.data.shollIntMorphStd, ...
		     obj.data.shollIntMorphSEM, ...
		     obj.data.shollIntSynMean, ...
		     obj.data.shollIntSynStd, ...
		     obj.data.shollIntSynSEM, ...
		     obj.data.shollIntXMean, ...
		     obj.data.shollIntXStd, ...
		     obj.data.shollIntXSEM, ...
		     obj.data.shollIntSynMean./obj.data.shollIntMorphMean, ...
		     obj.data.shollIntXMean./obj.data.shollIntMorphMean, ...
		     obj.data.shollIntXMean./obj.data.shollIntSynMean, ...
                   };

      obj.makeXMLsheet(docNode,'Sholl analysis', ...
                       columnName, columnData);

    end
		   

    if(obj.exportInfo.saveSomaMeasure)
      columnName = { 'Soma morph (mean)', ...
                     'Soma morph (std)', ...
                     'Soma morph (SEM)', ...
                     'Soma syn (mean)', ...
                     'Soma syn (std)', ...
                     'Soma syn (SEM)', ...
                     'Soma X (mean)', ...
                     'Soma X (std)', ...
                     'Soma X (SEM)' };
      
       columnData = { mean(obj.data.somaMeasureMorph), ...
                      std(obj.data.somaMeasureMorph), ...
                      std(obj.data.somaMeasureMorph) ...
                      /sqrt(length(obj.data.somaMeasureMorph)), ...	 
                      mean(obj.data.somaMeasureSyn), ...
                      std(obj.data.somaMeasureSyn), ...		    
                      std(obj.data.somaMeasureSyn) ...
                      / sqrt(length(obj.data.somaMeasureSyn)), ...       
                      mean(obj.data.somaMeasureX), ...
                      std(obj.data.somaMeasureX), ...		    
                      std(obj.data.somaMeasureX) ...
                      / sqrt(length(obj.data.somaMeasureX)), ...
                    };
       
    else
      columnName = {};
      columnData = {};
    end

    columnName2 = {
		     'Soma Area (mean)', ...
		     'Soma Area (std)', ...
		     'Soma Area (SEM)', ...
		     'Soma Major Axis (mean)', ...
		     'Soma Major Axis (std)', ...
		     'Soma Major Axis (SEM)', ...
		     'Soma Minor Axis (mean)', ...
		     'Soma Minor Axis (std)', ...
		     'Soma Minor Axis (SEM)', ...
                   };

   columnData2 = {
		  mean(obj.data.somaArea), ...
		  std(obj.data.somaArea), ...
		  std(obj.data.somaArea) ...
		    / sqrt(length(obj.data.somaArea)), ...
		  mean(obj.data.somaMajorAxisLength), ...
		  std(obj.data.somaMajorAxisLength), ...
		  std(obj.data.somaMajorAxisLength) ...
  		    / sqrt(length(obj.data.somaMajorAxisLength)), ...
		  mean(obj.data.somaMinorAxisLength), ...
		  std(obj.data.somaMinorAxisLength), ...
		  std(obj.data.somaMinorAxisLength) ...
  		    / sqrt(length(obj.data.somaMinorAxisLength)), ...
                 };

    % Concatenate cell arrays
    columnName = {columnName{:},columnName2{:}};
    columnData = {columnData{:},columnData2{:}}; 

    obj.makeXMLsheet(docNode,'Soma measures', ...
                     columnName, columnData);



    if(obj.exportInfo.saveSynapseProfile)

      columnName = {'Distance (micrometer)', ...
		    'Morphology', ...
		    'Synapse', ...
		    'X', ...
		    'Syn/Morph', ...
		    'X/Morph', ...
		    'X/Syn' };

      columnData = { obj.data.meanSynapseProfileDist*1e6, ...
		     obj.data.meanSynapseProfileMorph, ...
		     obj.data.meanSynapseProfileSyn, ...
		     obj.data.meanSynapseProfileX, ...
		     obj.data.meanSynapseProfileSyn ...
		     ./obj.data.meanSynapseProfileMorph, ...
		     obj.data.meanSynapseProfileX ...
		     ./obj.data.meanSynapseProfileMorph, ...
		     obj.data.meanSynapseProfileX ...
		     ./obj.data.meanSynapseProfileSyn, ...
                   };


      obj.makeXMLsheet(docNode,'Synapse profile', ...
                       columnName, columnData);
    end

    if(obj.exportInfo.saveGradient)
      columnName = { 'From (mu)','To (mu)', ...
		     'Morph (mean)', 'Morph (std)', 'Morph (SEM)', ...
		     'Syn (mean)', 'Syn (std)', 'Syn (SEM)', ...
		     'X (mean)', 'X (std)', 'X (SEM)', ...	     
		     'Syn/Morph (mean)', ...
		     'X/Moph (mean)', ...
                   };

      columnData = { [0,obj.data.gradientEdges(1:end-1)]*1e6, ...
		     obj.data.gradientEdges*1e6,...
		     obj.data.morphGradientMean, ...
		     obj.data.morphGradientStd, ...
		     obj.data.morphGradientSEM, ...
		     obj.data.synGradientMean, ...
		     obj.data.synGradientStd, ...
		     obj.data.synGradientSEM, ...
		     obj.data.XGradientMean, ...
		     obj.data.XGradientStd, ...
		     obj.data.XGradientSEM, ...
		     obj.data.synGradientMean./obj.data.morphGradientMean, ...
		     obj.data.XGradientMean./obj.data.morphGradientMean, ...
                   };

      obj.makeXMLsheet(docNode,'Neurite gradients', ...
                       columnName, columnData);

    end

    % Write all to disk

    fprintf('Exporting data to %s\n',obj.data.exportXMLfile);
    xmlwrite(strcat(obj.data.exportPath,obj.data.exportXMLfile),docNode);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function docNode = makeXMLheader(obj)

    docNode = com.mathworks.xml.XMLUtils.createDocument('Workbook');
    docRootNode = docNode.getDocumentElement();

    docRootNode.setAttribute('xmlns','urn:schemas-microsoft-com:office:spreadsheet');
    docRootNode.setAttribute('xmlns:o','urn:schemas-microsoft-com:office:office');
    docRootNode.setAttribute('xmlns:x','urn:schemas-microsoft-com:office:excel');
    docRootNode.setAttribute('xmlns:ss','urn:schemas-microsoft-com:office:spreadsheet');

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Helper function to write a sheet of data
  % columnData is a cell array, where each element is either a vector 
  % or a cell array. Data in vectors get written as numbers, and data
  % in cell arrays get written as string to the xml file.


  function makeXMLsheet(obj,docNode, ...
                        sheetName, ...
                        columnHeader, ...
                        columnData)

    docRootNode = docNode.getDocumentElement();

    docSheet = docNode.createElement('Worksheet');
    docRootNode.appendChild(docSheet);
    docSheet.setAttribute('ss:Name',sheetName);

    docTable = docNode.createElement('Table');
    docSheet.appendChild(docTable);

    for i = 1:length(columnHeader)
      docTable.appendChild(docNode.createElement('Column'));
    end

    docHeader = docNode.createElement('Row');
    docTable.appendChild(docHeader);

   try

    for i = 1:length(columnHeader)

      docCell = docNode.createElement('Cell');
      docHeader.appendChild(docCell);

      docData = docNode.createElement('Data');
      docData.setAttribute('ss:Type','String');
      % docObj.Data.setAttribute('ss:Type','Number');
      docCell.appendChild(docData);

      docData.appendChild(docNode.createTextNode(columnHeader{i}));

    end

    maxRows = 0;

    for j = 1:length(columnData)
      maxRows = max(maxRows, length(columnData{j}));

      if(length(columnData{j}) ~= numel(columnData{j}))
        disp('Warning, XML write called with matrix, this is bad! Report it!')
      end

    end

    for j = 1:maxRows
      docRow = docNode.createElement('Row');
      docTable.appendChild(docRow);

      for i = 1:length(columnData)
        docCell = docNode.createElement('Cell');
        docRow.appendChild(docCell);

        docData = docNode.createElement('Data');
        docCell.appendChild(docData);

        if(length(columnData{i}) < j)
          % This column has less rows of data than other columns
          %disp('Short column, leaving additional rows empty.')

          docData.setAttribute('ss:Type','String');
          docData.appendChild(docNode.createTextNode(''));

          continue
        end

        if(isa(columnData{i},'double'))
          % We have numerical data
          if(columnData{i}(j) < inf)
            docData.setAttribute('ss:Type','Number');
            docData.appendChild(docNode.createTextNode(num2str(columnData{i}(j))));
          elseif(isnan(columnData{i}(j)))
            % Leave empty if NaN
            docData.setAttribute('ss:Type','String');
            docData.appendChild(docNode.createTextNode(''));
          else
	    % Excel 2003 cant handle INF, lets give it a HUGE number!
            docData.setAttribute('ss:Type','Number');
            docData.appendChild(docNode.createTextNode('1e300'));
          end
        else
          % We have string data
          docData.setAttribute('ss:Type','String');
          docData.appendChild(docNode.createTextNode(columnData{i}{j}));
        end
      end
    end

   % Temp debug, 
   catch exception
     getReport(exception)
     disp('If you see this, talk to Johannes, he can save your data!')
     keyboard
   end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function skipSynapseDetection(obj, source, event)

    obj.activateStage('analyse');
    obj.stopEdit();

    if(max(max(obj.data.distMask .* obj.data.neuriteMask)) == inf)
      uiwait(errordlg(['Manually edit the neurite mask to connect the dark ' ...
                       'grey neurites to the soma or use the clean button ' ...
                       'to automatically remove unconnected components!'], ...
 	               'Unconnected components', 'modal'))
      return             
    end

    obj.dispInfo.axis = []; % Reset view
    obj.showGUI('analyse');

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function exportOnlyAverages(obj, source, event)

    exportAverageFile = strcat(obj.data.fileName{obj.dispInfo.curImg}, ...
			       '-onlyAverage.csv');

    curPwd = pwd;
    try
      cd(obj.data.exportPath)
    catch
      fprintf('Unable to change to %s\n', obj.data.exportPath)
    end

    [exportAverageFile, savePath] = ...
	  uiputfile('*-onlyAverage.csv', ...
		    'Export mask averages', ...
		    exportAverageFile);

    cd(curPwd);

    if(isempty(exportAverageFile) | exportAverageFile == 0)
      disp('Export aborted.')
      return
    end

    obj.data.exportPath = savePath;

    tmpMorph = obj.data.image(:,:,obj.detection.morphChannel, obj.dispInfo.curImg);
    tmpSyn = obj.data.image(:,:,obj.detection.synChannel, obj.dispInfo.curImg);
    tmpX = obj.data.image(:,:,obj.detection.XChannel, obj.dispInfo.curImg);

    fid = fopen(strcat(obj.data.exportPath,exportAverageFile),'w');
    fprintf(fid, ['File,Soma area,Soma M,Soma M (SEM),Soma S,Soma S (SEM),' ...
                  'Soma X,Soma X (SEM),Neurite area,Neurite M,' ...
                  'Neurite M (SEM),Neurite S,Neurite S (SEM),Neurite X,' ...
                  'Neurite X (SEM)\n']);

    fprintf(fid,'%s,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n', ...
            obj.data.fileName{obj.dispInfo.curImg}, ...
            nnz(obj.data.somaMask)*(obj.data.xyRes*1e6)^2, ...
            mean(tmpMorph(find(obj.data.somaMask))), ...
            std(tmpMorph(find(obj.data.somaMask))) ...
            /sqrt(nnz(obj.data.somaMask)), ...
            mean(tmpSyn(find(obj.data.somaMask))), ...
            std(tmpSyn(find(obj.data.somaMask))) ...
            /sqrt(nnz(obj.data.somaMask)), ...
            mean(tmpX(find(obj.data.somaMask))), ...
            std(tmpX(find(obj.data.somaMask))) ...
            /sqrt(nnz(obj.data.somaMask)), ...
            nnz(obj.data.neuriteMask)*(obj.data.xyRes*1e6)^2, ...
            mean(tmpMorph(find(obj.data.neuriteMask))), ...
            std(tmpMorph(find(obj.data.neuriteMask))) ...
            /sqrt(nnz(obj.data.neuriteMask)), ...
            mean(tmpSyn(find(obj.data.neuriteMask))), ...
            std(tmpSyn(find(obj.data.neuriteMask))) ...
            /sqrt(nnz(obj.data.neuriteMask)), ...
            mean(tmpX(find(obj.data.neuriteMask))), ...
            std(tmpX(find(obj.data.neuriteMask))) ...
            /sqrt(nnz(obj.data.neuriteMask)));
    
    fclose(fid);

    obj.saveConfig();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function nextStage(obj,source, event)

    % Make sure the next stage is active before advancing

    % Turn of any motion obj.handles
    set(obj.handles.fig,'windowbuttonmotionfcn', [])

    switch(obj.dispInfo.state)
      case 'load'
        obj.stopEdit();

        if(obj.dispInfo.stage >= 2)
          obj.dispInfo.axis = []; % Reset view
          obj.showGUI('soma');
          if(~nnz(obj.data.somaMask))
            % No pixels in the soma mask, try and detect it
            obj.detectSoma();
          end
        end

      case 'soma'
        obj.stopEdit();

        if(obj.dispInfo.stage >= 3)
          if(nnz(obj.data.somaMask))      
            obj.dispInfo.axis = []; % Reset view  
            obj.showGUI('neurites');
            if(~nnz(obj.data.neuriteMask))
              obj.detectNeurite();
            end
          else
            % There are no pixels in the soma mask...
            uiwait(errordlg('You must detect or draw a soma!', ...
			    'No soma detected', 'modal'))
            return
          end
        end

      case 'neurites'
        obj.stopEdit();

        if(obj.dispInfo.stage >= 4)
          % Verify that there are no unconnected neurites...
	      
          if(max(max(obj.data.distMask .* obj.data.neuriteMask)) == inf)
            uiwait(errordlg(['Manually edit the neurite mask to connect the dark ' ...
                             'grey neurites to the soma or use the clean button ' ...
                             'to automatically remove unconnected components!'], ...
			    'Unconnected components', 'modal'))
            return             
          end

          obj.dispInfo.axis = []; % Reset view
          obj.showGUI('synapses');
          if(~nnz(obj.data.synapseMask))
            obj.detectSynapses();
          end
        end

      case 'synapses'
        obj.stopEdit();

        if(obj.dispInfo.stage >= 5)
          obj.dispInfo.axis = []; % Reset view
          obj.showGUI('analyse');
        end

      case 'analyse'
        obj.stopEdit();

        obj.dispInfo.axis = []; % Reset view
        % We are already at last stage, do nothing
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function activeIcon(obj, iconNumber)

    for i = 1:length(iconNumber)
      set(obj.handles.allIcons(iconNumber(i)),'backgroundcolor',[1 0 0]);
    end

    restOfIcons = setdiff(1:length(obj.handles.allIcons),iconNumber);
    for i = 1:length(restOfIcons)
      set(obj.handles.allIcons(restOfIcons(i)),'backgroundcolor',get(gcf,'color'));
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function showGUI(obj,GUIstate)

    makeVisible = [];
    makeInvisible = [];

    % Update icons so they are correct based on current stage
    % of analysis.
    obj.setIcons();

    switch(GUIstate)
      case 'load'
        % The GUI shown when we are loading image and selecting
        % which frames in a stack to display.
        % Also used to select which channel are which

        makeVisible = obj.loadGUI;
        makeInvisible = setdiff(obj.allGUI,obj.loadGUI);

        if(obj.data.num <= 1)
          multiImgButtons = [obj.handles.next,obj.handles.prev,obj.handles.num];
          makeVisible = setdiff(makeVisible,multiImgButtons);
          makeInvisible = union(makeInvisible, multiImgButtons);
        end

        obj.dispInfo.state = 'load';  
        obj.setChannelVisible(1,1,1);
        obj.nextStageButtonStatus(1);

        % Used for handling the histograms
        set(obj.handles.fig,'WindowButtonDownFcn', @obj.mouseHandler);

        set(obj.handles.XYres,'String', num2str(obj.data.xyRes*1e6));

        obj.activeIcon(1);
        uicontrol(obj.handles.loadNeuron);

      case 'soma'
        % The GUI shown when we are locating the soma

        if(strcmp(obj.dispInfo.state,'load'))
          % Moving from load? If so remove histogram listener
          set(obj.handles.fig,'WindowButtonDownFcn', []);
        end

        makeVisible = obj.somaGUI;
        makeInvisible = setdiff(obj.allGUI,obj.somaGUI);
        obj.dispInfo.state = 'soma';  
        obj.setChannelVisible(1,0,0);
        obj.nextStageButtonStatus(2);

        obj.dispInfo.somaColor = [1 1 1];
        obj.dispInfo.neuriteColor = NaN;
        obj.dispInfo.synapseColor = NaN;
        obj.dispInfo.measurePointColor = NaN;

        obj.activeIcon(2);
        % uicontrol(obj.handles.detectSoma);

      case 'neurites'
        % The GUI to locate neurites and edit neurite mask

        if(strcmp(obj.dispInfo.state,'load'))
          % Moving from load? If so remove histogram listener
          set(obj.handles.fig,'WindowButtonDownFcn', []);
        end

        makeVisible = obj.neuriteGUI;
        makeInvisible = setdiff(obj.allGUI,obj.neuriteGUI);

        % Alternate which GUI is shown depending on settings
        if(obj.detection.useSteerableFilters)
          cfgObj.Handles = [obj.handles.morphThresh, obj.handles.morphThreshLabel];
        else
          cfgObj.Handles = [obj.handles.growThresh, ...
                        obj.handles.growThreshLabel, ...
                        obj.handles.filterSize, ...
                        obj.handles.filterSizeLabel, ...
                        obj.handles.addNeurite];
        end

        makeVisible = setdiff(makeVisible,cfgObj.Handles);
        makeInvisible = union(makeInvisible,cfgObj.Handles);

        obj.dispInfo.state = 'neurites';  
        obj.setChannelVisible(1,0,0);
        obj.nextStageButtonStatus(3);

        obj.dispInfo.somaColor = [1 1 1];
        obj.dispInfo.neuriteColor = [1 1 1]*0.7;
        obj.dispInfo.synapseColor = NaN;
        obj.dispInfo.measurePointColor = NaN;

        obj.activeIcon(3);
        % uicontrol(obj.handles.detectNeurite);

      case 'synapses'
        % The GUI for synapse detection

        if(strcmp(obj.dispInfo.state,'load'))
          % Moving from load? If so remove histogram listener
          set(obj.handles.fig,'WindowButtonDownFcn', []);
        end

        makeVisible = obj.synapseGUI;
        makeInvisible = setdiff(obj.allGUI,obj.synapseGUI);

        obj.dispInfo.state = 'synapses';  
        obj.setChannelVisible(0,1,0);
        obj.nextStageButtonStatus(4);

        obj.dispInfo.somaColor = NaN;
        obj.dispInfo.neuriteColor = NaN;
        obj.dispInfo.synapseColor = [1 1 1];
        obj.dispInfo.measurePointColor = NaN;

        obj.activeIcon(4);
        % uicontrol(obj.handles.detectSynapses);

      case 'analyse'
        % Analyse and export GUI

        % Reset viewpoint
        obj.dispInfo.axis = [];

        if(strcmp(obj.dispInfo.state,'load'))
          % Moving from load? If so remove histogram listener
          set(obj.handles.fig,'WindowButtonDownFcn', []);
        end

        makeVisible = obj.analyseGUI;
        makeInvisible = setdiff(obj.allGUI,obj.analyseGUI);
        obj.dispInfo.state = 'analyse';  
        obj.setChannelVisible(1,1,1);

        obj.activeIcon(5);
        uicontrol(obj.handles.exportData);

      otherwise
        fprintf('Unkown GUI setup: %s\n', GUIstate)
    end

    set(makeInvisible,'visible','off');
    for i = 1:length(makeInvisible)
      children = get(makeInvisible(i),'children');
      set(children,'visible','off')
    end

    set(makeVisible,'visible','on');

    for i = 1:length(makeVisible)
      children = get(makeVisible(i),'children');
      set(children,'visible','on')
    end

    % Hide axes of the histograms, always.
    set([obj.handles.redHist, obj.handles.greenHist, obj.handles.blueHist],'visible','off');

    obj.showImage();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function mouseHandler(obj, source, event)

     if(isempty(obj.data.image))
       % Nothing loaded yet, ignore
       return
     end

    [insideRed, xRed, yRed] = obj.checkInside(obj.handles.redHist);
    [insideGreen, xGreen, yGreen] = obj.checkInside(obj.handles.greenHist);
    [insideBlue, xBlue, yBlue] = obj.checkInside(obj.handles.blueHist);

    if(insideRed)
      obj.dispInfo.scaleRed = obj.data.maxRed / xRed;
      fprintf('Setting red scale to %d\n', obj.dispInfo.scaleRed)
    elseif(insideGreen)
      obj.dispInfo.scaleGreen = obj.data.maxGreen / xGreen;
      fprintf('Setting green scale to %d\n', obj.dispInfo.scaleGreen)
    elseif(insideBlue)
      obj.dispInfo.scaleBlue = obj.data.maxBlue / xBlue;
      fprintf('Setting blue scale to %d\n', obj.dispInfo.scaleBlue)
    end
 
    obj.showImage();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Check that we clicked inside axis of figure

  function [insideFlag, x, y] = checkInside(obj, testHandle)

    tmpXY = get(testHandle,'CurrentPoint');
    x = tmpXY(1,1); 
    y = tmpXY(1,2);

    set(obj.handles.fig,'CurrentAxes',testHandle)
    a = axis();

    if(a(1) <= x & x <= a(2) & a(3) <= y & y <= a(4))
      insideFlag = 1;
    else
      insideFlag = 0;
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function startZoomWheelMode(obj,source,event)

    % disp('Zoom mode active, click to center, scroll wheel to zoom')
    set(obj.handles.fig,'WindowScrollWheelFcn', @obj.zoomHandler);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function zoomHandler(obj, source, event)

    set(obj.handles.fig,'CurrentAxes',obj.handles.image)

    a = axis();

    xCenter = (a(1)+a(2))/2;
    yCenter = (a(3)+a(4))/2;

    aWidth = a(2)-a(1);
    aHeight = a(4)-a(3);

    % Where did the user click, within the image?
    xy = get(obj.handles.image,'CurrentPoint');
    xCenter = xy(1,1);
    yCenter = xy(1,2);

    if(xCenter < a(1) | a(2) < xCenter ...
       | yCenter < a(3) | a(4) < yCenter)

      disp('Outside of figure, ignoring')
      % set(obj.handles.fig,'WindowScrollWheelFcn', '');

      return
    end

    aHeight = aHeight * (1.25)^event.VerticalScrollCount;
    aWidth = aWidth * (1.25)^event.VerticalScrollCount;

    xCenter = min(max(aWidth/2+1,xCenter),obj.data.width-aWidth/2);
    yCenter = min(max(aHeight/2+1,yCenter),obj.data.height-aHeight/2);

    if(aWidth > obj.data.width | aHeight > obj.data.height)
      a = [1 obj.data.width 1 obj.data.height];
    else
      a(1) = xCenter-aWidth/2;
      a(2) = xCenter+aWidth/2;          
      a(3) = yCenter-aHeight/2;
      a(4) = yCenter+aHeight/2;
    end

    axis(a);
 
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % !!! Orphaned function

  % Extra function to calculate synapse threshold sensitivity

  function [nSynapse, synArea, synStd] = synapseThresholdSensitivity(obj,threshRange)

    oldThresh = obj.detection.synThreshold;

    nSynapse = NaN*ones(size(threshRange));
    synArea = NaN*ones(size(threshRange));
    synStd = NaN*ones(size(threshRange));

    for iT = 1:length(threshRange)
      fprintf('Thresh: %d\n', threshRange(iT))
      synStd(iT) = threshRange(iT);

      obj.detection.synThreshold = threshRange(iT);
      obj.detectSynapseCenters();
      obj.calculateSynapseProperties();
      nSynapse(iT) = length(obj.data.synapseCenter);
      synArea(iT) = nnz(obj.data.synapseMask);
    end

    obj.detection.synThreshold = oldThresh;

    figure
    plot(synStd,nSynapse,'k-')
    xlabel('Synapse threshold (std)')
    ylabel('Synapse count')

    figure
    plot(synStd,synArea,'k-')
    xlabel('Synapse threshold (std)')
    ylabel('Total synapse area')

    [transpose(synStd),transpose(nSynapse), transpose(synArea)]

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function resizeFunction(obj, source, event)

    % Resize the bitmaps for the icons
    obj.setIcons();

    if(isempty(obj.data.image))
      % disp('Resizing and redisplaying logo.')
      obj.showLogo();
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setIcons(obj)

    % We want the icons to fit the push buttons

    % Find out how big the button is
    set(obj.handles.loadStage,'units','pixels')
    buttonPos = get(obj.handles.loadStage,'position');
    buttonSize = buttonPos(3:4)-4;
    set(obj.handles.loadStage,'units','normalized')

    % Make sure to desaturete the icons that should be
    maxStage = 5;
    stageFlag = (1:maxStage) <= obj.dispInfo.stage;

    for i = 1:length(obj.bilder.icon)
      obj.bilder.iconResize{i} = imresize(obj.bilder.icon{i},buttonSize([2 1]));
    end

    for i = 1:length(stageFlag)
      if(stageFlag(i))
        set(obj.allIcons(i), ...
            'Callback',obj.allIconCallback{i}, ...
            'Cdata', obj.bilder.iconResize{i});
				      
      else
        set(obj.allIcons(i), ...
            'Callback',[], ...
            'Cdata', obj.desaturateIcon(obj.bilder.iconResize{i}));
      end
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function loadStage(obj, source, event)
    obj.dispInfo.axis = []; % Reset view
    obj.stopEdit();
    obj.showGUI('load');
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function somaStage(obj, source, event)
    obj.dispInfo.axis = []; % Reset view
    obj.stopEdit();
    obj.showGUI('soma');
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function neuriteStage(obj, source, event)
    obj.dispInfo.axis = []; % Reset view
    obj.stopEdit();
    obj.showGUI('neurites');
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function synapseStage(obj, source, event)
    obj.dispInfo.axis = []; % Reset view
    obj.stopEdit();
    obj.showGUI('synapses');
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function analyseStage(obj, source, event)
    obj.dispInfo.axis = []; % Reset view
    obj.stopEdit();
    obj.showGUI('analyse');
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function activateStage(obj, stage)

    switch(stage)
      case 'load'
        obj.dispInfo.stage = 1;

      case 'soma'
        set(obj.handles.loadNeuron, 'BackgroundColor', get(gcf,'color'));
        obj.dispInfo.stage = 2;

      case 'neurites'
        set(obj.handles.detectSoma, 'BackgroundColor', get(gcf,'color'));
        obj.dispInfo.stage = 3;

      case 'synapses'
        set(obj.handles.detectNeurite, 'BackgroundColor', get(gcf,'color'));
        obj.dispInfo.stage = 4;

      case 'analyse'
        set(obj.handles.detectSynapses, 'BackgroundColor', get(gcf,'color'));
        obj.dispInfo.stage = 5;

    end

    obj.showGUI(obj.dispInfo.state);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function bwIcon = desaturateIcon(obj, colourIcon)

    bwIcon = zeros(size(colourIcon));

    for i = 1:3
      bwIcon(:,:,i) = mean(colourIcon,3)/255.0;
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function loadConfig(obj)

    if(exist(obj.configFile))
      fprintf('Loading configuration from %s\n', obj.configFile)
      cfg = load(obj.configFile);

      try
        obj.detection.morphChannel = cfg.old.morphChannel;
        obj.detection.synChannel = cfg.old.synChannel;
        obj.detection.XChannel = cfg.old.XChannel;
        obj.detection.remapChan = cfg.old.remapChan;
        obj.detection.enableRemapChan = cfg.old.enableRemapChan;

        obj.detection.filterSize = cfg.old.filterSize;
        obj.detection.singleSoma = cfg.old.singleSoma;

        obj.data.xyRes = cfg.old.xyRes;

        obj.detection.somaErodeRadius = cfg.old.somaErodeRadius;
        obj.detection.morphThreshold = cfg.old.morphThreshold;
        obj.detection.autoSomaThreshold = cfg.old.autoSomaThreshold;

        obj.detection.measurePointRadius = cfg.old.measurePointRadius;
        obj.detection.nMeasurePoints = cfg.old.nMeasurePoints;

        obj.detection.synThreshold = cfg.old.synThreshold;
        obj.detection.minSynapseSize = cfg.old.minSynapseSize;
        obj.detection.singleSynapseRadie = cfg.old.singleSynapseRadie;
        obj.detection.neuritePadding = cfg.old.neuritePadding;
        obj.detection.maxRadie = cfg.old.maxRadie;
        obj.detection.excludeSomaSynapses = cfg.old.excludeSomaSynapses;

        obj.detection.shollBinSize = cfg.old.shollBinSize;

        obj.exportInfo.saveMask = cfg.old.saveMask;

        % obj.exportInfo.saveSholl = cfg.old.saveSholl;
        obj.exportInfo.saveSholl = 1; % Always set it to save

        obj.exportInfo.saveIntensityHistogram = cfg.old.saveIntHist;
        obj.exportInfo.saveMat = cfg.old.saveMat;

        obj.data.loadPath = cfg.old.loadPath;
        obj.data.exportPath = cfg.old.exportPath;

        obj.editInfo.fileTypeOrder = cfg.old.fileTypeOrder;

        obj.dispInfo.figPosition = cfg.old.figPosition;

        obj.detection.numSaved = cfg.old.numSaved;

      catch exception
        getReport(exception)
        disp('Config file belonged to old version.')
      end
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function saveSynDpath(obj,source, event)
    disp('Saving current SynD path permanently')
    savepath
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setSynDpath(obj)

    obj.SyndScriptPath = obj.SynDscriptFull(1:(max(strfind(obj.SynDscriptFull, ...
						 obj.SynDscript))-1));

    SynDpath = { obj.SyndScriptPath, ...
                 strcat(obj.SyndScriptPath,'LSM/lsm'), ...
                 strcat(obj.SyndScriptPath,'LSM/cstruct'), ...
                 strcat(obj.SyndScriptPath,'bilder') };

    fprintf('Set SynD work path to %s\n', SynDpath{1})

    for i = 1:length(SynDpath)
      addpath(SynDpath{i});
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function saveConfig(obj)

    fprintf('Saving old configuration to %s\n', obj.configFile)

    old.morphChannel = obj.detection.morphChannel;
    old.synChannel = obj.detection.synChannel;
    old.XChannel = obj.detection.XChannel;
    old.remapChan = obj.detection.remapChan;
    old.enableRemapChan = obj.detection.enableRemapChan;

    old.filterSize = obj.detection.filterSize;
    old.singleSoma = obj.detection.singleSoma;

    old.xyRes = obj.data.xyRes;

    old.somaErodeRadius = obj.detection.somaErodeRadius;
    old.morphThreshold = obj.detection.morphThreshold;
    old.autoSomaThreshold = obj.detection.autoSomaThreshold;

    old.measurePointRadius = obj.detection.measurePointRadius;
    old.nMeasurePoints = obj.detection.nMeasurePoints;

    old.synThreshold = obj.detection.synThreshold;
    old.minSynapseSize = obj.detection.minSynapseSize;
    old.singleSynapseRadie = obj.detection.singleSynapseRadie;
    old.neuritePadding = obj.detection.neuritePadding;
    old.maxRadie = obj.detection.maxRadie;
    old.excludeSomaSynapses = obj.detection.excludeSomaSynapses;

    old.shollBinSize = obj.detection.shollBinSize;

    old.saveMask = obj.exportInfo.saveMask;
    old.saveSholl = obj.exportInfo.saveSholl;
    old.saveIntHist = obj.exportInfo.saveIntensityHistogram;
    old.saveMat = obj.exportInfo.saveMat;
    old.loadPath = obj.data.loadPath;
    old.exportPath = obj.data.exportPath;

    old.fileTypeOrder = obj.editInfo.fileTypeOrder;

    old.figPosition = get(obj.handles.fig,'position');

    old.numSaved = obj.detection.numSaved;

    try 
      save(obj.configFile, 'old');
    catch
      fprintf('Unable to save config file %s.\n', obj.configFile)
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function nextStageButtonStatus(obj,dispStage)
    if(dispStage < obj.dispInfo.stage)
      % We can advance to next stage
      set(obj.handles.nextStage, 'BackgroundColor', [1 0.6 0]);

      % Set focus to next button
      uicontrol(obj.handles.nextStage);
    else
      % Still not finished with this stage
      set(obj.handles.nextStage, 'BackgroundColor', [0.8 0.8 0.8]);
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function captureKeyPress(obj,source, event)

    try
      % Was it a keypress?
      event.Key;
    catch
      % ginput messes with scroll wheel handling, fixing it!
      obj.startZoomWheelMode();
      obj.zoomHandler(source, event);
      return
    end

    % Shared keys for all
    switch(event.Key)
      case 'r'
        obj.toggleRed();
      case 'g'
        obj.toggleGreen();
      case 'b'
        obj.toggleBlue();
    end

    switch(event.Character)
      case '+'

        if(isempty(obj.data.image))
          return
        end

        % Zoom in
        set(obj.handles.fig,'CurrentAxes',obj.handles.image)    
        disp('Zoom mode.')
        [xCenter,yCenter,buttonPressed] = ginput(1);

        a = axis;

        if(xCenter < a(1) | a(2) < xCenter ...
           | yCenter < a(3) | a(4) < yCenter)
          return
        end

        while(1 <= buttonPressed & buttonPressed <= 3)
           buttonPressed

           a = axis;

          if(buttonPressed == 1)
            % Zoom in
            aWidth = (a(2)-a(1))/2;
            aHeight = (a(4)-a(3))/2;
          else
            aWidth = (a(2)-a(1))*2;
            aHeight = (a(4)-a(3))*2;
          end

          xCenter = min(max(aWidth/2+1,xCenter),obj.data.width-aWidth/2);
          yCenter = min(max(aHeight/2+1,yCenter),obj.data.height-aHeight/2);

          if(aWidth > obj.data.width | aHeight > obj.data.height)
            a = [1 obj.data.width 1 obj.data.height];
          else
            a(1) = xCenter-aWidth/2;
            a(2) = xCenter+aWidth/2;          
            a(3) = yCenter-aHeight/2;
            a(4) = yCenter+aHeight/2;
          end
          axis(a);

          fprintf('Zoom mode.')
          [xCenter,yCenter,buttonPressed] = ginput(1);

          if(xCenter < a(1) | a(2) < xCenter ...
             | yCenter < a(3) | a(4) < yCenter)
            return
          end
         end

      case '='
        if(isempty(obj.data.image))
          return
        end

        % Zoom out
        set(obj.handles.fig,'CurrentAxes',obj.handles.image)    
        a = axis;
        xCenter = mean(a(1:2));
        yCenter = mean(a(3:4));
        aWidth = min((a(2)-a(1))*2,obj.data.width);
        aHeight = min((a(4)-a(3))*2,obj.data.height);

        xCenter = min(max(aWidth/2+1,xCenter),obj.data.width-aWidth/2);
        yCenter = min(max(aHeight/2+1,yCenter),obj.data.height-aHeight/2);

        a(1) = xCenter-aWidth/2;
        a(2) = xCenter+aWidth/2;          
        a(3) = yCenter-aHeight/2;
        a(4) = yCenter+aHeight/2;
        axis(a);

    end

    % Separate keys depending on state
    switch(obj.dispInfo.state)
      case 'load'
        switch(event.Key)
          case 'l'
            obj.loadNeuron();
          case 'n'
            obj.nextStage();
          case 'x'
            % Focus to the XY-res
	    uicontrol(obj.handles.XYres)
        end
      case 'soma'
        switch(event.Key)
          case 'a'
            obj.addSoma();
          case 'd'
            obj.detectSoma();
          case 'n'
            obj.nextStage();
          case 'e'
            obj.editSoma();
          case 'm'
            obj.toggleMask();
          case 's'
            % Change focus
            uicontrol(obj.handles.somaErodeRadius);
          case 't'
            uicontrol(obj.handles.morphThresh);
          case 'u'
            obj.undoLastAction();
          case 'r'
            obj.setSomaMeasurePoints();
          case 'o'
            obj.addSomaMeasurePoints();
        end
     
      case 'neurites'
        switch(event.Key)
          case 'd'
            obj.detectNeurite();
          case 'n'
            obj.nextStage();
          case 'e'
            obj.editNeurite();
          case 'c'
            obj.cleanMask();
          case 'a'
            if(obj.detection.useSteerableFilters)
              obj.addNeurite();
            end
          case 'o'
            obj.killBlob();
          case 'm'
            obj.toggleMask();
          case 'f'
            if(obj.detection.useSteerableFilters)
              uicontrol(obj.handles.growThresh);
            else
              uicontrol(obj.handles.morphThresh);
            end
          case 't'
            obj.addThinNeurites();
          case 'u'
            obj.undoLastAction();
          case 's'
            obj.showSkeleton();
        end

      case 'synapses'
        switch(event.Key)
          case 'd'
            obj.detectSynapses();
          case 'n'
            obj.nextStage();
          case 'e'
            obj.editSynapse();
          case 'm'
            obj.toggleMask();
          case 't'
            uicontrol(obj.handles.synThresh);
          case 'u'
            obj.undoLastAction();
        end

      case 'analyse'
        switch(event.Key)
          case 'e'
            obj.exportData();
        end

    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
      
end
