% SynD - Synapse Detector
% 
% --- THIS IS THE OLD VERSION, SynD_extractNEW.m is updated to use
% object oriented matlab
% 
% Johannes Hjorth
% j.j.j.hjorth@damtp.cam.ac.uk
% johannes.hjorth@cncr.vu.nl
% hjorth@kth.se
%
% Sabine Schmitz
% sabine.schmitz@cncr.vu.nl
%
% The latest version of the code is available from:
% http://www.johanneshjorth.se/SynD
% 
% You can also find the code at:
%
% http://software.incf.org/software/synd
% http://www.cncr.nl/resources
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

function SynD_extract()

  disp('--> This version of SynD might not work on new versions of matlab.')
  disp('--> The latest version is: SynD_extract2017.m')
  disp(' ')
   
  % To make sure the config file is in the right place...
  SynDscript = mfilename();
  SynDscriptFull = mfilename('fullpath');
  SyndScriptPath = [];

  % close all
  format compact

  licenseOk = 1;

  if(~license('checkout','matlab'))
    disp('MATLAB license missing.')
    licenseOk = 0;
  end

  if(~license('checkout','image_toolbox'))
    disp('Image Toolbox license missing.')
    licenseOk = 0;
  end
  
  if(~license('checkout','statistics_toolbox'))
    disp('Statistics Toolbox license missing.')
    licenseOk = 0;
  end

  if(~licenseOk)
    disp('SynD does not have access to all required functions')
    disp('Terminating now because the needed functions are not available.')
    return
  end

  setSynDpath()

  configFile = strcat(SyndScriptPath,'SynD_config.mat');

  enableDebug = 1;

  referenceStr = ['Schmitz SK, Hjorth JJJ, ' ...
                  'Joemai RMS, Wijntjes R, ' ...
                  'Eijgenraam S,  de Bruijn P, Georgiou C, ' ...
                  'de Jong APH, van Ooyen A, Verhage M, ' ...
                  'Cornelisse LN, Toonen RF, Veldkamp W ', ...
                  'Automated Analysis of Neuronal Morphology, ' ...
                  'Synapse Number and Synaptic Recruitment, ' ...
                  'J Neurosci Methods 195 (2011) 185–193'];

  referenceDOI = '10.1016/j.neumeth.2010.12.011';

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

  % Data structure which contains image data and results

  data.image = [];
  data.height = 0;
  data.width = 0;
  data.num = 0;
  data.maxRed = 0;
  data.maxGreen = 0;
  data.maxBlue = 0;
  data.xyRes = 0.20e-6; % meters 
  data.somaMask = [];
  data.neuriteMask = [];
  data.synapseMask = [];
  data.somaMeasureMask = [];
  data.includeMask = [];
  data.synapseCenter = [];
  data.synapseArea = [];
  data.synapseDist = [];
  data.synapsePixels = {};
  data.neuriteLength = NaN;
  data.distMask = [];

  data.synapseIntensityMorphMean = []; % morph-channel
  data.synapseIntensityMorphSEM = [];  % morph-channel
  data.synapseIntensitySynMean = [];   % syn-channel
  data.synapseIntensitySynSEM = [];    % syn-channel
  data.synapseIntensityXMean = [];     % X-channel
  data.synapseIntensityXSEM = [];      % X-channel

  data.meanSynapse = []; % Template used for synapse center detection

  data.meanSynapseMorph = [];
  data.meanSynapseSyn = [];
  data.meanSynapseX = [];    % Profile of X channel for synapse

  data.meanSynapseProfileDist = [];
  data.meanSynapseProfileMorph = [];
  data.meanSynapseProfileSyn = [];
  data.meanSynapseProfileX = [];

  data.fileName = {};

  data.skeleton = [];
  data.shollEdges = [NaN NaN NaN];
  data.shollDendHist = [];
  data.shollSynHist = [];

  data.shollIntMorphMean = [];
  data.shollIntMorphStd = [];
  data.shollIntMorphSEM = [];

  data.shollIntSynMean = [];
  data.shollIntSynStd = [];
  data.shollIntSynSEM = [];

  data.shollIntXMean = [];
  data.shollIntXStd = [];
  data.shollIntXSEM = [];

  data.somaMeasureMorph = NaN;
  data.somaMeasureSyn = NaN;
  data.somaMeasureX = NaN;

  data.somaArea = [];
  data.somaMajorAxisLength = [];
  data.somaMinorAxisLength = [];

  data.measurePoint = [];

  data.intensityHistogramEdges = [];
  data.morphHist = [];
  data.synHist = [];
  data.XHist = [];

  % Intensity gradient along neurites
  data.gradientEdges = [];

  data.morphGradientMean = [];
  data.morphGradientStd = [];
  data.morphGradientSEM = [];

  data.synGradientMean = [];
  data.synGradientStd = [];
  data.synGradientSEM = [];

  data.XGradientMean = [];
  data.XGradientStd = [];
  data.XGradientSEM = [];

  data.exportPath = pwd;
  data.exportXMLfile = [];
  data.exportSaveFile = [];
  data.exportNeuriteMaskFile = [];
  data.exportSynapseMaskFile = [];

  data.loadPath = pwd;

  data.dirVect = [];
  data.rigidity = [];
  data.lambdaMax = 0;

  data.somaIntensityThreshold = NaN;
  data.synapseIntensityThreshold = NaN;

  % Information used to display the image, only relevant to session

  dispInfo.curImg = 1;
  dispInfo.showRed = 1;
  dispInfo.showGreen = 1;
  dispInfo.showBlue = 1;
  dispInfo.showMask = 1;

  dispInfo.scaleRed = 1;
  dispInfo.scaleGreen = 1;
  dispInfo.scaleBlue = 1;

  dispInfo.axis = [];
  dispInfo.state = 'load'; % Which GUI view is currently displayed
  dispInfo.stage = 1; % How far along are we in processing?

  dispInfo.somaColor = [1 1 1];
  dispInfo.neuriteColor = NaN; %[1 1 1]*0.7;
  dispInfo.synapseColor = NaN;
  dispInfo.measurePointColor = NaN;
  dispInfo.defaultMeasurePointColor = [1 1 0];

  dispInfo.showSkeleton = 0;
  dispInfo.showSynapseProfileDataPoints = 1;

  dispInfo.needCleaning = 0;

  dispInfo.figPosition = [];

  % Add additional debug plots
  dispInfo.verbose = 0;

  % Editing information, only relevant to session

  editInfo.mode = 0;
  editInfo.line = [];
  editInfo.color = [1 1 1];
  editInfo.width = 3;
  editInfo.defaultWidth = 3;
  editInfo.xLine = [];
  editInfo.yLine = [];
  editInfo.undo = struct('somaMask', [], ...
			 'neuriteMask', [], ...
			 'synapseMask', [], ...
			 'includeMask', [], ...
			 'synapseCenter', [], ...
			 'measurePoint', [], ...
			 'description', [], ...
			 'state', [], ...
			 'stage', 0);

  editInfo.maxUndo = 10;

  editInfo.fileTypeOrder = [1 2]; % LSM=1 or TIFF=2 default file type?


  exportInfo.saveMask = 1;
  exportInfo.saveSholl = 1;
  exportInfo.saveIntensityHistogram = 1;
  exportInfo.saveMat= 1;
  exportInfo.saveSomaMeasure = 0;
  exportInfo.saveSynapseProfile = 1;
  exportInfo.saveGradient = 1;

  % Parameters used for detection

  detection.morphChannel = 1; 
  detection.synChannel = 2;   
  detection.XChannel = 3;     

  detection.remapChan = [1 2 3]; % Normally R=ch1, G=ch2, B=ch3
  detection.enableRemapChan = 'off';

  detection.singleSoma = 1;

  detection.measurePointRadius = 5;
  detection.nMeasurePoints = 10;

  detection.wienerSize = 7; % NaN to disable
  detection.morphThreshold = 200; % 0.03?
  detection.somaErodeRadius = 15;
  detection.autoSomaThreshold = 0; % Use cross-entropy minimization?

  detection.useSteerableFilters = 1; % Should we use steerable filters
  detection.filterSize = 1.8e-6; % in meters, corresponding to about 6 pixels in normal resolution;
  detection.maxAddCost = 0.9; % Maximal cost if we still want to add pixel
  detection.connectCostLambda = 0.7; % Weighting to lambda (vs eigenvectorcost)

  % Meijering 2004 uses elongated filters, set it to -1/3 to use theirs
  % If you change this, also consider decreasing the maxAddCost (e.g. 0.85)
  detection.steerableFilterAlpha = 0; % -1/3; 

  detection.minProtrusionLength = 2e-6; % Min length of thin neurites added
  detection.pixelQue = []; 

  detection.smoothingRadius = 25;

  detection.synThreshold = 0.6;  % How many standard deviations above noise?
  detection.backgroundFiltering = true;  % Subtract smoothed image?
  detection.minSynapseSize = 0.35e-12; % corresponds to < 4 pixels;
  detection.neuritePadding = 1e-6;  % How many meters outside neurite
                                    % are synapses allowed to be 
  detection.singleSynapseRadie = NaN; %0.5e-6; % micrometers, NaN gives guestimation of profile
  detection.maxRadie = 10; % Pixel size of synapse template
  detection.excludeSomaSynapses = true;

  detection.numSaved = 0; % How many times has the user saved a neuron

  % Neurites shorter than 10 pixels do not count
  detection.trimNeuriteSize = 12;

  % Bin size for Sholl analysis
  detection.shollBinSize = 5e-6;

  detection.intensityBinSize = 100;

  % detection.synapseProfileBinSize = 0.5e-6; % we use data.xyRes instead

  detection.reanalyseFlag = 0;

  % Load saved config values from file
  loadConfig();

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  handles = struct([]);
  
  handles(1).fig = figure('Name','SynD - Synapse Detection', ...
		       'MenuBar','none', ...
		       'Toolbar','figure', ...
		       'Position', [50 50 1150 680], ...
		       'Resize','on', ...
		       'ResizeFcn',@resizeFunction);

  handles.image = axes('Units','Pixels', ...
		       'Position', [50 50 615 615]);

  handles.redHist = axes('Units','Pixels', ...
			 'Visible','off', ...
			 'Position', [700 310 424 100]);

  handles.greenHist = axes('Units','Pixels', ...
			   'Visible','off', ...
			   'Position', [700 180 424 100]);

  handles.blueHist = axes('Units','Pixels', ...
			  'Visible','off', ...
			  'Position', [700 50 424 100]);

  handles.credits = uicontrol('Style', 'text', ...
			      'String', 'Johannes Hjorth and Sabine Schmitz, 2010', ...
			      'HorizontalAlignment', 'right', ...
			      'Foregroundcolor', 0.7*[1 1 1], ...
			      'Backgroundcolor', get(gcf,'color'), ...
			      'Position', [925 5 210 15], ...
			      'Fontsize',8);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Load GUI
  %

  handles.next = uicontrol('Style','pushbutton', ...
			   'String','<html>Next</html>', ...
			   'Interruptible','off', ...
			   'Visible','off', ...
			   'Position', [965 417 60 20], ...
			   'Callback', @nextImage);

  handles.prev = uicontrol('Style','pushbutton', ...
			   'String','<html>Prev</html>', ...
			   'Interruptible','off', ...
			   'Visible','off', ...
			   'Position', [895 417 60 20], ... 
			   'Callback', @prevImage);

  handles.num = uicontrol('Style','text', ...
			  'String', 'No image', ...
			  'Position', [1045 417 100 20], ... 
			  'Visible','off', ...
			  'BackgroundColor', get(gcf,'color'), ...
			  'FontSize', 12, ...
			  'HorizontalAlignment', 'left');

  handles.loadNeuron = uicontrol('Style','pushbutton', ...
				 'String','<html><u>L</u>oad image</html>', ...
				 'Interruptible','off', ...
				 'Position', [680 530 120 35], ...
				 'Fontsize', 12, ...
				 'BackgroundColor', [1 0.6 0], ...
				 'Callback', @loadNeuron);

  handles.nextStage = uicontrol('Style','pushbutton', ...
				'String','<html><u>N</u>ext</html>', ...
				'TooltipString', ...
				['<html>Continue to next stage<br>' ...
				 'of analysis when this<br>' ...
			         'one is done</html>'], ...
				'Interruptible','off', ...
				'Position', [820 530 120 35], ...
				'Fontsize', 12, ...
				'BackgroundColor', [0.8 0.8 0.8], ...
				'Callback', @nextStage);


  % Detect soma GUI
  %

  handles.detectSoma = uicontrol('Style','pushbutton', ...
				 'String','<html><u>D</u>etect soma</html>', ...
				 'Interruptible','off', ...
				 'Visible','off', ...
				 'Position', [680 530 125 35], ...
				 'Fontsize', 12, ...
				 'BackgroundColor', [1 0.6 0], ...
				 'TooltipString', ...
                                 ['<html>Left click and draw to add soma,<br>' ...
				  'right click and draw to remove soma</html>'], ...
				 'Callback', @detectSoma);

  handles.editSoma = uicontrol('Style','pushbutton', ...
			       'String','<html><u>E</u>dit soma</html>', ...
			       'Interruptible','off', ...
			       'Visible','off', ...
			       'Position', [680 420 120 30], ...
			       'Fontsize', 12, ...
			       'Callback', @editSoma);

  handles.addSoma = uicontrol('Style','pushbutton', ...
			      'String','<html><u>A</u>dd soma</html>', ...
			      'Interruptible','off', ...
			      'Visible','off', ...
			      'Position', [820 420 120 30], ...
			      'TooltipString', ...
			      'Manually add ROI around soma.', ...
			      'Fontsize', 12, ...
			      'Callback', @addSoma);

  handles.deleteSoma = uicontrol('Style','pushbutton', ...
				 'String','<html>Delete soma</html>', ...
				 'Interruptible','off', ...
				 'Visible','off', ...
				 'Position', [960 420 120 30], ...
				 'Fontsize', 12, ...
				 'Callback', @deleteSoma);


  handles.addMeasurePoints = uicontrol('Style','pushbutton', ...
				       'String','<html>New measu<u>r</u>e</html>', ...
				       'TooltipString', ...
				       ['<html>Add 10 random regions to the soma<br>', ...
										     'to measure protein levels.</html>'], ...
				       'Interruptible','off', ...
				       'Visible','off', ...
				       'Position', [680 310 120 30], ...
				       'Fontsize', 12, ...
				       'Callback', @setSomaMeasurePoints);

  handles.editMeasurePoints = uicontrol('Style','pushbutton', ...
					'String','<html>Add p<u>o</u>int</html>', ...
					'Interruptible','off', ...
					'Visible','off', ...
					'TooltipString', ...
				        '<html>Add regions in soma</html>', ...
					'Position', [820 310 120 30], ...
					'Fontsize', 12, ...
					'Callback', @addSomaMeasurePoints);

  handles.somaMeasureRadiusLabel = uicontrol('Style','text', ...
					     'String', 'Radius:', ...
					     'Visible','off', ...
					     'Position', [960 310 50 20], ...
					     'BackgroundColor', get(gcf,'color'), ...
					     'FontSize', 10, ...
					     'HorizontalAlignment', 'left');

  handles.somaMeasureRadius = uicontrol('Style','edit', ...
					'String', ...
					num2str(detection.measurePointRadius), ...
					'Visible','off', ...
					'Position', [1010 312 50 20], ... 
					'TooltipString', ...
					'<html>Radius of soma measure</html>', ...
					'Interruptible', 'off', ...
					'Callback', @setSomaMeasureRadius);



  handles.measureLabel = uicontrol('Style','text', ...
				   'String', 'Measure protein levels in soma:', ...
				   'Visible','off', ...
				   'Position', [680 350 270 20], ... 
				   'BackgroundColor', get(gcf,'color'), ...
				   'FontSize', 12, ...
				   'HorizontalAlignment', 'left');


  % Detect neurite GUI
  % 

  handles.detectNeurite = uicontrol('Style','pushbutton', ...
				    'String','<html><u>D</u>etect neurites</html>', ...
				    'Interruptible','off', ...
				    'Visible','off', ...
				    'Position', [680 531 125 35], ...
				    'Fontsize', 12, ...
				    'BackgroundColor', [1 0.6 0], ...
				    'Callback', @detectNeurite);

  handles.editNeurite = uicontrol('Style','pushbutton', ...
				  'String','<html><u>E</u>dit neurite</html>', ...
				  'Interruptible','off', ...
				  'Visible','off', ...
				  'Position', [680 420 120 30], ...
				  'Fontsize', 12, ...
				  'TooltipString', ...
                                  ['<html>Left click and draw adds neurite,<br>' ...
				   'right click and draw removes neurite</html>'], ...
				  'Callback', @editNeurite);

  handles.addThinNeurites = uicontrol('Style','pushbutton', ...
				      'String','<html>Add <u>t</u>hin</html>', ...
				      'Interruptible','off', ...
				      'Visible','off', ...
				      'Position', [680 380 120 30], ...
				      'Fontsize', 12, ...
				      'TooltipString', ...
				      ['<html>Extend detection of neurites, using<br>' ...
				             'steerable filter with half the size.</html>'], ...
				      'Callback', @addThinNeurites);


  handles.addNeurite = uicontrol('Style','pushbutton', ...
				 'String','<html><u>A</u>dd neurite</html>', ...
				 'Interruptible','off', ...
				 'Visible','off', ...
				 'Position', [820 420 100 30], ...
				 'Fontsize', 12, ...
				 'Visible','off', ...
				 'TooltipString', ...
				 'Click on neurite to auto detect', ...
				 'Callback', @addNeurite);

  handles.clean = uicontrol('Style','pushbutton', ...
			    'String','<html><u>C</u>lean</html>', ...
			    'Interruptible','off', ...
			    'Visible','off', ...
			    'Position', [940 420 100 30], ...
			    'Fontsize', 12, ...
			    'Visible','off', ...
			    'TooltipString', ...
			    'Removes all unconnected neurite parts', ...
			    'Callback', @cleanMask);

  handles.killBlob = uicontrol('Style','pushbutton', ...
			       'String','<html>Bl<u>o</u>b erase</html>', ...
			       'Interruptible','off', ...
			       'Visible','off', ...
			       'Position', [940 380 100 30], ...
			       'Fontsize', 12, ...
			       'Visible','off', ...
			       'TooltipString', ...
			       'Removes one unconnected neurite part', ...
			       'Callback', @killBlob);

  handles.extendNeurites = uicontrol('Style','pushbutton', ...
				     'String','<html>Extend</html>', ...
				     'Interruptible','off', ...
				     'Visible','off', ...
				     'Position', [820 380 100 30], ...
				     'Fontsize', 12, ...
				     'Visible','off', ...
				     'TooltipString', ...
				     'Looks for unconnected neurites', ...
				     'Callback', @extendNeuritesHandler);

  handles.showSkeleton = uicontrol('Style','pushbutton', ...
				   'String','<html><u>S</u>keleton</html>', ...
			       'Interruptible','off', ...
			       'Visible','off', ...
			       'Position', [940 340 100 30], ...
			       'Fontsize', 12, ...
			       'Visible','off', ...
			       'TooltipString', ...
			       'Shows the neurite skeleton.', ...
			       'Callback', @showSkeleton);


  handles.calculateLabel = uicontrol('Style', 'Text', ...
				     'String', '', ...
				     'Position', [680 320 250 30], ...
				     'BackgroundColor', get(gcf,'color'), ...
				     'Visible','off', ...
				     'FontSize', 10, ...
				     'HorizontalAlignment', 'left');


  handles.skipSynapseDetection = uicontrol('Style','pushbutton', ...
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
					'Callback', @skipSynapseDetection);
			
  handles.exportSimpleAverages = uicontrol('Style','pushbutton', ...
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
					'Callback', @exportOnlyAverages);
		

  % Synapse GUI
  %

  handles.detectSynapses = uicontrol('Style','pushbutton', ...
				     'String','<html><u>D</u>etect synapses</html>', ...
				     'Interruptible','off', ...
				     'Visible','off', ...
				     'Position', [680 531 125 35], ...
				     'Fontsize', 12, ...
				     'BackgroundColor', [1 0.6 0], ...
				     'Callback', @detectSynapses);



  handles.editSynapse = uicontrol('Style','pushbutton', ...
				  'String','<html><u>E</u>dit synapses</html>', ...
				  'Interruptible','off', ...
				  'Visible','off', ...
				  'Position', [680 420 120 30], ...
				  'Fontsize', 12, ...
				  'TooltipString', ...
				  ['<html>Left click and draw adds synapse<br>' ...
				   'area, right click and draw removes<br>' ... 
				   'synapse area</html>'], ...
				  'Callback', @editSynapse);

  handles.killSynBlob = uicontrol('Style','pushbutton', ...
				  'String','<html>Bl<u>o</u>b erase</html>', ...
				  'Interruptible','off', ...
				  'Visible','off', ...
				  'Position', [820 420 100 30], ...
				  'Fontsize', 12, ...
				  'Visible','off', ...
				  'TooltipString', ...
				  'Removes one unconnected neurite part', ...
				  'Callback', @killSynBlob);


  handles.singleSynapse = axes('Units','Pixels', ...
			       'visible','off', ...
			       'Position', [720 100 400 300]);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  handles.morphLabel = uicontrol('Style','text', ...
				 'String', 'Morph:', ...
				 'Visible','off', ...
				 'Position', [680 490 60 20], ... 
				 'BackgroundColor', get(gcf,'color'), ...
				 'FontSize', 12, ...
				 'HorizontalAlignment', 'left');

  handles.morphChannel = uicontrol('Style','popupmenu', ...
				   'String', {'R','G','B'}, ...
				   'Position', [740 495 60 15], ... 
				   'FontSize', 10, ...
				   'TooltipString', ...
				   'Channel colour of morphology staining', ...
				   'Value', detection.remapChan(1), ...
				   'Callback', @setMorphChannel);

  handles.morphChannelNumber = uicontrol('Style','popupmenu', ...
					 'String', {'1','2','3'}, ...
					 'Enable', detection.enableRemapChan, ...
					 'Position', [740 465 60 15], ... 
					 'FontSize', 10, ...
					 'TooltipString', ...
					 'Channel with morphology staining', ...
					 'Value', ...
					 detection.morphChannel, ...
					 'Callback', @remapChannels);


  handles.synLabel = uicontrol('Style','text', ...
			       'String', 'Syn:', ...
			       'Position', [810 490 40 20], ...
			       'BackgroundColor', get(gcf,'color'), ...
			       'FontSize', 12, ...
			       'HorizontalAlignment', 'left');

  handles.synChannel = uicontrol('Style','popupmenu', ...
				 'String', {'R','G','B'}, ...
				 'Position', [850 495 60 15], ...
				 'FontSize', 10, ...
				 'TooltipString', ...
				 'Channel colour of synapse staining', ...
				 'Value', detection.remapChan(2), ...
				 'Callback', @setSynChannel);

  handles.synChannelNumber = uicontrol('Style','popupmenu', ...
				       'String', {'1','2','3'}, ...
				       'Enable', detection.enableRemapChan, ...
				       'Position', [850 465 60 15], ...
				       'FontSize', 10, ...
				       'TooltipString', ...
				       'Channel with synapse staining', ...
				       'Value', ...
				       detection.synChannel, ...
				       'Callback', @remapChannels);

  handles.XLabel = uicontrol('Style','text', ...
			     'String', 'X:', ...
			     'Position', [920 490 20 20], ... 
			     'BackgroundColor', get(gcf,'color'), ...
			     'FontSize', 12, ...
			     'HorizontalAlignment', 'left');

  handles.XChannel = uicontrol('Style','popupmenu', ...
			       'String', {'R','G','B'}, ...
			       'Position', [940 495 60 15], ... 
			       'FontSize', 10, ...
			       'TooltipString', ...
			       'Channel to be analysed', ...
			       'Value', detection.remapChan(3), ...
			       'Callback', @setXChannel);

  handles.XChannelNumber = uicontrol('Style','popupmenu', ...
				     'String', {'1','2','3'}, ...
				     'Enable', detection.enableRemapChan, ...
				     'Position', [940 465 60 15], ... 
				     'FontSize', 10, ...
				     'TooltipString', ...
				     'Channel to be analysed', ...
				     'Value', ...
				     detection.XChannel, ...
				     'Callback', @remapChannels);

  handles.remapLabel = uicontrol('Style','text', ...
				 'String', 'Remap:', ...
				 'Visible','off', ...
				 'Position', [1020 460 60 20], ... 
				 'BackgroundColor', get(gcf,'color'), ...
				 'FontSize', 12, ...
				 'HorizontalAlignment', 'left');

  handles.remap = uicontrol('Style', 'checkbox', ...
			    'Position', [1080 460 20 20], ... 
			    'BackgroundColor', get(gcf,'color'), ...
			    'TooltipString', ...
			    'Enable remapping of channels', ...
			    'Value', strcmp(detection.enableRemapChan,'on'), ...
			    'Callback', @toggleRemap);


  handles.XYresLabel = uicontrol('Style','text', ...
				 'String', 'XY res (micrometer) :', ...
				 'Position', [680 415 133 20], ...
				 'BackgroundColor', get(gcf,'color'), ...
				 'FontSize', 10, ...
				 'HorizontalAlignment', 'left');

  handles.XYres= uicontrol('Style','edit', ...
			   'String', num2str(data.xyRes*1e6), ...
			   'Position', [813 417 60 20], ... 
			   'TooltipString', ...
			   '<html>Image resolution<br>in micrometers</html>', ...
			   'Interruptible', 'off', ...
			   'Callback', @setXYres);



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  handles.red = uicontrol('Style', 'checkbox', ...
			  'Position', [1004 530 20 20], ... 
			  'BackgroundColor', get(gcf,'color'), ...
			  'TooltipString', ...
			  'Show red channel', ...
			  'Value', 1, ...
			  'Callback', @toggleRed);

  handles.green = uicontrol('Style', 'checkbox', ...
			    'Position', [1054 530 20 20], ...
			    'BackgroundColor', get(gcf,'color'), ...
			    'TooltipString', ...
			    'Show green channel', ...
			    'Value', 1, ...
			    'Callback', @toggleGreen);

  handles.blue = uicontrol('Style', 'checkbox', ...
			   'Position', [1104 530 20 20], ... 
			   'BackgroundColor', get(gcf,'color'), ...
			   'TooltipString', ...
			   'Show blue channel', ...
			   'Value', 1, ...
			   'Callback', @toggleBlue);

  handles.redLabel = uicontrol('Style','text', ...
			       'String', 'R', ...
			       'Position', [984 530 20 20], ... 
			       'BackgroundColor', get(gcf,'color'), ...
			       'FontSize', 12, ...
			       'HorizontalAlignment', 'left');

  handles.greenLabel = uicontrol('Style','text', ...
				 'String', 'G', ...
				 'Position', [1034 530 20 20], ... 
				 'BackgroundColor', get(gcf,'color'), ...
				 'FontSize', 12, ...
				 'HorizontalAlignment', 'left');

  handles.blueLabel = uicontrol('Style','text', ...
				'String', 'B', ...
				'Position', [1084 530 20 20], ... 
				'BackgroundColor', get(gcf,'color'), ...
				'FontSize', 12, ...
				'HorizontalAlignment', 'left');

  handles.mask = uicontrol('Style', 'checkbox', ...
			   'Position', [1104 500 20 20], ... 
			   'Visible','off', ...
			   'BackgroundColor', get(gcf,'color'), ...
			   'TooltipString', ...
			   'Show mask', ...
			   'Value', 1, ...
			   'Callback', @toggleMask);

  handles.maskLabel = uicontrol('Style','text', ...
				'String', 'M', ...
				'Visible','off', ...
				'Position', [1084 500 20 20], ... 
				'BackgroundColor', get(gcf,'color'), ...
				'FontSize', 12, ...
				'HorizontalAlignment', 'left');

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
  handles.morphThreshLabel = uicontrol('Style','text', ...
				       'String', 'Morph thresh:', ...
				       'Visible','off', ...
				       'Position', [680 460 100 20], ...
				       'BackgroundColor', get(gcf,'color'), ...
				       'FontSize', 10, ...
				       'HorizontalAlignment', 'left');

  handles.morphThresh = uicontrol('Style','edit', ...
				  'String', num2str(detection.morphThreshold), ...
				  'Visible','off', ...
				  'Enable', onOff(~detection.autoSomaThreshold), ...
				  'Position', [780 462 50 20], ... 
				  'TooltipString', ...
				  ['<html>Intensity threshold for morphology:<br>' ...
				   'Relative threshold (fraction of max intensity) if smaller than 1.<br>' ...
                                   'Absolute threshold if larger than 1.</html>'], ...
				  'Interruptible', 'off', ...
				  'Callback', @setMorphThresh);

  handles.guessSomaThreshold = uicontrol('Style','pushbutton', ...
					 'String','<html>Auto</html>', ...
					 'Interruptible','off', ...
					 'Visible','off', ...
					 'Position', [850 462 50 20], ...
					 'TooltipString', ...
					 [ '<html>Find threshold by minimizing<br>' ...
					   'cross-entropy of image.</html>'], ...
					 'Fontsize', 12, ...
					 'BackgroundColor', get(gcf,'color'), ...
					 'Callback', @autoSomaThreshold);

  handles.somaErodeLabel = uicontrol('Style','text', ...
				     'String', 'Soma erode:', ...
				     'Visible','off', ...
				     'Position', [680 490 100 20], ... 
				     'BackgroundColor', get(gcf,'color'), ...
				     'FontSize', 10, ...
				     'HorizontalAlignment', 'left');

  handles.somaErodeRadius = uicontrol('Style','edit', ...
				      'String', num2str(detection.somaErodeRadius), ...
				      'TooltipString', ...
				      ['<html>Higher value removes more<br>' ...
				       'of the thresholded image</html>'], ...
				      'Visible','off', ...
				      'Position', [780 492 50 20], ... 
				      'Interruptible', 'off', ...
				      'Callback', @setSomaErodeRadius);

  handles.singleSomaLabel = uicontrol('Style','text', ...
				      'String', 'Single soma:', ...
				      'Visible','off', ...
				      'Position', [850 490 100 20], ... 
				      'BackgroundColor', get(gcf,'color'), ...
				      'FontSize', 10, ...
				      'HorizontalAlignment', 'left');

  handles.singleSoma = uicontrol('Style', 'checkbox', ...
				 'Visible','off', ...
				 'Position', [940 492 20 20], ... 
				 'BackgroundColor', get(gcf,'color'), ...
				 'TooltipString', ...
				 'Show red channel', ...
				 'Value', detection.singleSoma, ...
				 'Callback', @toggleSingleSoma);


  handles.synThreshLabel = uicontrol('Style','text', ...
				     'String', 'Syn thresh:', ...
				     'Visible','off', ...
				     'Position', [680 490 80 20], ... 
				     'BackgroundColor', get(gcf,'color'), ...
				     'FontSize', 10, ...
				     'HorizontalAlignment', 'left');

  handles.synThresh = uicontrol('Style','edit', ...
				'String', num2str(detection.synThreshold), ...
				'Visible','off', ...
				'Position', [760 492 50 20], ... 
				'TooltipString', ...
                                ['<html>Number of standard deviations<br>' ...
                                 'above noise level required to<br>' ...
                                 'be considered a synapse.<br>' ...
                                 'If NaN uses cross-entropy for <br>' ...
                                 'threshold detection.<br>' ...
                                 'If negative, use absolute value as<br>' ...
                                 'absolute intensity threshold.</html>'], ...
				'Interruptible', 'off', ...
				'Callback', @setSynThresh);

  handles.neuritePaddingLabel = uicontrol('Style','text', ...
					  'String', 'Neurite padding:', ...
					  'Visible','off', ...
					  'Position', [830 490 115 20], ... 
					  'BackgroundColor', get(gcf,'color'), ...
					  'FontSize', 10, ...
					  'HorizontalAlignment', 'left');

  handles.neuritePadding = uicontrol('Style','edit', ...
				     'String', num2str(detection.neuritePadding*1e6), ...
				     'Visible','off', ...
				     'Position', [945 492 50 20], ... 
				     'TooltipString', ...
				     ['<html>How many micrometers outside a neurite<br>' ...
				      'synapses are allowed to extend</html>'], ...
				     'Interruptible', 'off', ...
				     'Callback', @setNeuritePadding);

  handles.minSynapseSizeLabel = uicontrol('Style','text', ...
					  'String', 'Min size:', ...
					  'Visible','off', ...
					  'Position', [830 460 115 20], ... 
					  'BackgroundColor', get(gcf,'color'), ...
					  'FontSize', 10, ...
					  'HorizontalAlignment', 'left');

  handles.minSynapseSize = uicontrol('Style','edit', ...
				     'String', num2str(detection.minSynapseSize*1e12), ...
				     'Visible','off', ...
				     'Position', [945 462 50 20], ... 
				     'TooltipString', ...
				     ['<html>Minimum size in micrometers square<br>' ...
				      'for detected synapses.</html>'], ...
				     'Interruptible', 'off', ...
				     'Callback', @setSynapseMinSize);


  handles.growThreshLabel = uicontrol('Style','text', ...
				      'String', 'Max cost:', ...
				      'Visible','off', ...
				      'Position', [680 490 80 20], ...
				      'BackgroundColor', get(gcf,'color'), ...
				      'FontSize', 10, ...
				      'HorizontalAlignment', 'left');

  handles.growThresh = uicontrol('Style','edit', ...
				 'String', num2str(detection.maxAddCost), ...
				 'Position', [760 492 50 20], ...
				 'Visible','off', ...
				 'TooltipString', ...
				 ['<html>Value between 0 and 1,<br>', ...
                                  'higher value means more<br>', ...
			          'neurite is included</html>'], ...
				 'Interruptible', 'off', ...
				 'Callback', @setGrowThresh);

  handles.filterSizeLabel = uicontrol('Style','text', ...
				      'String', 'Filter size:', ...
				      'Visible','off', ...
				      'Position', [680 460 80 20], ...
				      'BackgroundColor', get(gcf,'color'), ...
				      'FontSize', 10, ...
				      'HorizontalAlignment', 'left');

  handles.filterSize = uicontrol('Style','edit', ...
				 'String', num2str(detection.filterSize*1e6), ...
				 'Position', [760 462 50 20], ...
				 'Visible','off', ...
				 'TooltipString', ...
				 ['<html>Filter size, width in micrometers.<br>', ...
                                  'Lower values detect thinner neurites.</html>'], ...
				 'Interruptible', 'off', ...
				 'Callback', @setFilterSize);


  % Export buttons
  %

  handles.exportData = uicontrol('Style','pushbutton', ...
				 'String','<html><u>E</u>xport data</html>', ...
				 'Interruptible','off', ...
				 'Visible','off', ...
				 'Position', [680 530 120 35], ...
				 'Fontsize', 12, ...
				 'BackgroundColor', [0.8 0.8 0.8], ...
				 'Callback', @exportData);

  handles.restartButton = uicontrol('Style','pushbutton', ...
				    'String','<html>Restart SynD</html>', ...
				    'Interruptible','off', ...
				    'Visible','off', ...
				    'Position', [950 530 120 35], ...
				    'Fontsize', 12, ...
				    'BackgroundColor', [0.8 0.0 0.0], ...
				    'Callback', @restartSynD);


  handles.saveMaskLabel = uicontrol('Style','text', ...
				    'String', 'Save mask', ...
				    'Visible','off', ...
				    'Position', [680 480 120 20], ... 
				    'BackgroundColor', get(gcf,'color'), ...
				    'FontSize', 12, ...
				    'HorizontalAlignment', 'left');

  handles.saveMask = uicontrol('Style', 'checkbox', ...
			       'Visible','off', ...
			       'Position', [810 480 20 20], ... 
			       'BackgroundColor', get(gcf,'color'), ...
			       'Value', exportInfo.saveMask);

  handles.saveShollLabel = uicontrol('Style','text', ...
				     'String', 'Save Sholl', ...
				     'Visible','off', ...
				     'Position', [680 450 120 20], ... 
				     'BackgroundColor', get(gcf,'color'), ...
				     'FontSize', 12, ...
				     'HorizontalAlignment', 'left');

  handles.saveSholl = uicontrol('Style', 'checkbox', ...
				'Visible','off', ...
				'Position', [810 450 20 20], ... 
				'BackgroundColor', get(gcf,'color'), ...
				'Value', exportInfo.saveSholl);

  handles.saveIntensityHistLabel = uicontrol('Style','text', ...
					     'String', 'Save histogram', ...
					     'Visible','off', ...
					     'Position', [680 420 120 20], ... 
					     'BackgroundColor', get(gcf,'color'), ...
					     'FontSize', 12, ...
					     'HorizontalAlignment', 'left');

  handles.saveIntensityHist = uicontrol('Style', 'checkbox', ...
					'Visible','off', ...
					'Position', [810 420 20 20], ... 
					'BackgroundColor', get(gcf,'color'), ...
					'Value', exportInfo.saveIntensityHistogram);

  handles.saveMatLabel = uicontrol('Style','text', ...
				   'String', 'Save mat-file', ...
				   'Visible','off', ...
				   'Position', [680 390 120 20], ... 
				   'BackgroundColor', get(gcf,'color'), ...
				   'FontSize', 12, ...
				   'HorizontalAlignment', 'left');

  handles.saveMat = uicontrol('Style', 'checkbox', ...
			      'Visible','off', ...
			      'Position', [810 390 20 20], ... 
			      'BackgroundColor', get(gcf,'color'), ...
			      'Value', exportInfo.saveMat);

  handles.saveSomaMeasureLabel = uicontrol('Style','text', ...
					   'String', 'Save soma measure', ...
					   'Visible','off', ...
					   'Position', [680 360 120 20], ... 
					   'BackgroundColor', get(gcf,'color'), ...
					   'FontSize', 12, ...
					   'HorizontalAlignment', 'left');

  handles.saveSomaMeasure = uicontrol('Style', 'checkbox', ...
				      'Visible','off', ...
				      'Position', [810 360 20 20], ... 
				      'BackgroundColor', get(gcf,'color'), ...
				      'Value', exportInfo.saveSomaMeasure);

  handles.saveSynapseProfileLabel = uicontrol('Style','text', ...
					      'String', 'Save synapse profile', ...
					      'Visible','off', ...
					      'Position', [680 330 120 20], ... 
					      'BackgroundColor', get(gcf,'color'), ...
					      'FontSize', 12, ...
					      'HorizontalAlignment', 'left');


  handles.saveSynapseProfile = uicontrol('Style', 'checkbox', ...
					 'Visible','off', ...
					 'Position', [810 330 20 20], ... 
					 'BackgroundColor', get(gcf,'color'), ...
					 'Value', exportInfo.saveSynapseProfile);


  handles.saveGradientLabel = uicontrol('Style','text', ...
					'String', 'Save neurite gradient', ...
					'Visible','off', ...
					'Position', [680 300 120 20], ... 
					'BackgroundColor', get(gcf,'color'), ...
					'FontSize', 12, ...
					'HorizontalAlignment', 'left');


  handles.saveGradient = uicontrol('Style', 'checkbox', ...
				   'Visible','off', ...
				   'Position', [810 300 20 20], ... 
				   'BackgroundColor', get(gcf,'color'), ...
				   'Value', exportInfo.saveGradient);



  %%% Icons for switching "stages", ie load, soma, neurites, synapses, analyse
  % Swedish lesson: bilder = pictures  

  bilder.logo = [];

  try
    bilder.icon{1} = imread('bilder/icon_load.png');
    bilder.icon{2} = imread('bilder/icon_soma.png');
    bilder.icon{3} = imread('bilder/icon_neurite.png');
    bilder.icon{4} = imread('bilder/icon_synapse.png');
    bilder.icon{5} = imread('bilder/icon_export.png');
  catch
    close all
    disp('Icons not found, aborting.')
    disp('Please download from www.johanneshjorth.se')
    return
  end
  
  handles(1).loadStage = uicontrol('Style','pushbutton', ...
                                'Cdata', bilder.icon{1}, ...
                                'TooltipString', ...
                                'Load image with neuron', ...
                                'Interruptible','off', ...
                                'Position', [680 581 84 84], ...
                                'Callback', @loadStage);

  handles.somaStage = uicontrol('Style','pushbutton', ...
				'Cdata', desaturateIcon(bilder.icon{2}), ...
				'TooltipString', ...
				'<html>Detect soma and<br> edit soma mask</html>', ...
				'Interruptible','off', ...
				'Position', [770 581 84 84]);, ...


  handles.neuriteStage = uicontrol('Style','pushbutton', ...
				   'Cdata', desaturateIcon(bilder.icon{3}), ...
				   'TooltipString', ...
				   ['<html>Detect neurites and<br>' ...
                                    'edit neurite mask</html>'], ...
				   'Interruptible','off', ...
				   'Position', [860 581 84 84]);

  handles.synapseStage = uicontrol('Style','pushbutton', ...
				   'Cdata', desaturateIcon(bilder.icon{4}), ...
				   'TooltipString', ...
				   ['<html>Detect synapses and<br>' ...
				    'edit synapse mask</html>'], ...
				   'Interruptible','off', ...
				   'Position', [950 581 84 84]);

  handles.analyseStage = uicontrol('Style','pushbutton', ...
				   'Cdata', desaturateIcon(bilder.icon{5}), ...
				   'TooltipString', ...
				   'Export data to xml', ...
				   'Interruptible','off', ...
				   'Position', [1040 581 84 84]);


  handles.allIcons = [handles.loadStage, ...
 	 	      handles.somaStage, ...
		      handles.neuriteStage, ...
		      handles.synapseStage, ...
		      handles.analyseStage];

  %%% Helper variables for showGUI

  loadGUI = [ handles.loadNeuron, ...
              handles.prev, handles.next, ...
	      handles.num, ...
	      handles.next, handles.prev, ...
	      handles.morphLabel, ...
	      handles.morphChannel, ...
	      handles.morphChannelNumber, ...
	      handles.synLabel, ...
              handles.synChannel, ...
              handles.synChannelNumber, ...
	      handles.remapLabel, handles.remap, ...
	      handles.XLabel, ...
              handles.XChannel, ...
              handles.XChannelNumber, ...
	      handles.redLabel, handles.red, ...
	      handles.greenLabel, handles.green, ...
	      handles.blueLabel, handles.blue, ...
 	      handles.XYresLabel, handles.XYres, ...
  	      handles.nextStage, ...
 	      handles.redHist, handles.greenHist, handles.blueHist, ...
            ];

  somaGUI = [ handles.detectSoma, handles.editSoma, ...
 	      handles.addSoma, handles.deleteSoma, ...
 	      handles.addMeasurePoints, handles.editMeasurePoints, ...
	      handles.measureLabel, ...
	      handles.morphThreshLabel, handles.morphThresh, ...
	      handles.guessSomaThreshold, ...
 	      handles.singleSomaLabel, handles.singleSoma, ...
	      handles.somaErodeLabel, handles.somaErodeRadius, ...
 	      handles.somaMeasureRadiusLabel, ...
	      handles.somaMeasureRadius, ...
	      handles.redLabel, handles.red, ...
	      handles.greenLabel, handles.green, ...
	      handles.blueLabel, handles.blue, ...
	      handles.maskLabel, handles.mask, ...
  	      handles.nextStage, ...
            ];

  neuriteGUI = [ handles.detectNeurite, handles.editNeurite, ...
		 handles.addThinNeurites, ...
		 handles.morphThreshLabel, handles.morphThresh, ...
  		 handles.clean, handles.killBlob, ... 
		 handles.addNeurite, handles.extendNeurites, ...
		 handles.showSkeleton, ...
		 handles.growThreshLabel, handles.growThresh, ...
	 	 handles.filterSizeLabel, handles.filterSize, ...
		 handles.redLabel, handles.red, ...
 	         handles.greenLabel, handles.green, ...
	         handles.blueLabel, handles.blue, ...
	 	 handles.maskLabel, handles.mask, ...
  	 	 handles.nextStage, ...
 	 	 handles.calculateLabel, ...
		 handles.skipSynapseDetection, ...
 		 handles.exportSimpleAverages, ...
               ];

  synapseGUI = [ handles.detectSynapses, handles.editSynapse...
		 handles.killSynBlob, ...
		 handles.synThreshLabel, handles.synThresh, ...
		 handles.neuritePaddingLabel, handles.neuritePadding, ...
		 handles.minSynapseSizeLabel, handles.minSynapseSize, ...
  	         handles.redLabel, handles.red, ...
	         handles.greenLabel, handles.green, ...
	         handles.blueLabel, handles.blue, ...
	 	 handles.maskLabel, handles.mask, ...
 		 handles.nextStage, ...
	 	 handles.singleSynapse, ...
               ];

  analyseGUI = [ handles.exportData, ...
 		 handles.restartButton, ...
 	 	 handles.saveMaskLabel, handles.saveMask, ...
		 handles.saveShollLabel, handles.saveSholl, ...
  		 handles.saveIntensityHistLabel, ...
		 handles.saveIntensityHist, ...
		 handles.saveMatLabel, handles.saveMat, ...
		 handles.saveSomaMeasureLabel, handles.saveSomaMeasure, ...
 	 	 handles.saveSynapseProfileLabel, ...
		 handles.saveSynapseProfile, ...
		 handles.saveGradientLabel, ...
		 handles.saveGradient, ...
               ];


  allGUI = union(union(union(loadGUI, somaGUI), neuriteGUI), ...
		 union(synapseGUI, analyseGUI));

  allIcons = [handles.loadStage handles.somaStage ...
	      handles.neuriteStage handles.synapseStage ...
	      handles.analyseStage];

  allIconCallback = {@loadStage @somaStage ...
		     @neuriteStage @synapseStage ...
		     @analyseStage};

  allAxes = [handles.image, handles.singleSynapse, ... 
 	     handles.redHist, handles.greenHist, handles.blueHist];

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Make the figure resizeable by using normalisation
  set([allGUI allIcons handles.fig allAxes, ...
       handles.credits], ...
      'Units','Normalized')

  % Set the figure position to what was previously read in from config file
  if(~isempty(dispInfo.figPosition))
    set(handles.fig,'position', dispInfo.figPosition);
  end

  % Add callback function for keypresses
  set(handles.fig,'WindowKeyPressFcn',@captureKeyPress);

  showGUI('load');

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  handles.menuFile = uimenu(handles.fig,'Label','File');
  handles.menuItemLoad = uimenu(handles.menuFile, ...
				'Label','Load neuron', ...
				'Interruptible','off', ...
				'Callback', @loadNeuron);

  handles.menuItemBatch = uimenu(handles.menuFile, ...
				 'Label','Batch analyse neurons', ...
				 'Interruptible','off', ...
				 'Callback', @batchAnalyse);

  handles.menuItemBatch = uimenu(handles.menuFile, ...
				 'Label','Batch re-analyse neurons', ...
				 'Interruptible','off', ...
				 'Callback', @batchReAnalyse);

  handles.menuItemImportMask = uimenu(handles.menuFile, ...
				      'Label','Import old masks', ...
				      'Interruptible', 'off', ...
				      'Callback', @importMask);

  handles.menuItemExport = uimenu(handles.menuFile, ...
				  'Label','Export XML', ...
				  'Interruptible', 'off', ...
				  'Callback', @exportData);

  handles.menuItemRestart = uimenu(handles.menuFile, ...
				   'Label','Restart', ...
				   'Interruptible', 'off', ...
				   'Callback', @restartSynD);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  handles.menuEdit = uimenu(handles.fig,'Label','Edit');

  handles.menuItemUndo = uimenu(handles.menuEdit, ...
				'Label', '<html><u>U</u>ndo: No history</html>', ...
				'Interruptible', 'off', ...
				'Callback', @undoLastAction);

  handles.menuItemExclude = uimenu(handles.menuEdit, ...
				   'Label', 'Exclude region', ...
				   'Interruptible', 'off', ...
				   'Callback', @excludeRegion);

  handles.menuItemClearExclude = uimenu(handles.menuEdit, ...
					'Label', 'Reset exclude region', ...
					'Interruptible', 'off', ...
					'Callback', @clearExcludeRegions);

  handles.menuItemMaxProjection = uimenu(handles.menuEdit, ...
					 'Label', 'Max projection', ...
					 'Interruptible', 'off', ...
					 'Callback', @maxProjection);

  

%  handles.menuItemCollapse = uimenu(handles.menuEdit, ...
%				   'Label', 'Collapse image stack', ...
%				   'Interruptible', 'off', ...
%				   'Callback', @collapseImageStack);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  handles.menuSettings = uimenu(handles.fig,'Label','Settings');

  handles.autoThreshold = uimenu(handles.menuSettings, ...
				 'Label','<html>Auto soma threshold</html>', ...
				 'Interruptible', 'off', ...
				 'Checked', onOff(detection.autoSomaThreshold), ...
				 'Callback', @toggleAutoSomaThreshold);


  handles.menuItemSteerable = uimenu(handles.menuSettings, ...
				     'Label','<html>Use steerable <u>f</u>ilters</html>', ...
				     'Interruptible', 'off', ...
				     'Checked','on', ...
				     'Callback', @toggleSteerableFilters);

  handles.menuItemBackgroundFiltering = ...
    uimenu(handles.menuSettings, ...
	   'Label','Synapse background removal', ...
	   'Interruptible', 'off', ...
	   'Checked','on', ...
	   'Callback', @toggleBackgroundRemoval);

  handles.menuItemSteerableSettings = ...
    uimenu(handles.menuSettings, ...
	   'Label','Steerable filter settings', ...
	   'Interruptible', 'off', ...
	   'Callback', @steerableFilterSettings);


  handles.menuItemSynapseSettings = ...
    uimenu(handles.menuSettings, ...
	   'Label','Synapse detection settings', ...
	   'Interruptible', 'off', ...
	   'Callback', @synapseSettings);

  handles.menuItemShollSettings = ...
    uimenu(handles.menuSettings, ...
	   'Label','Sholl settings', ...
	   'Interruptible', 'off', ...
	   'Callback', @shollSettings);

  handles.menuItemSavePath = ...
    uimenu(handles.menuSettings, ...
	   'Label','Save SynD path permanently', ...
	   'Interruptible', 'off', ...
	   'Callback', @saveSynDpath);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Where to find updates and version files

  urlBase = 'http://www.johanneshjorth.se/files/SynD/';
  verFile = 'version.txt';

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  handles.menuHelp = uimenu(handles.fig,'Label','Help');

  webHelp = 'web(''http://www.johanneshjorth.se/SynD/'',''-browser'')';

  handles.menuItemHelp = uimenu(handles.menuHelp, ...
				'Label','How to use SynD', ...
				'Interruptible', 'off', ...
				'Callback', webHelp);

  handles.menuItemHelp = uimenu(handles.menuHelp, ...
				'Label','Get latest version', ...
				'Interruptible', 'off', ...
				'Callback', @getLatestVersion); 

  handles.menuAbout = uimenu(handles.menuHelp, ...
			     'Label', 'About SynD', ...
			     'Interruptible', 'off', ...
			     'Callback', @aboutSynD);			     


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if(enableDebug)

    handles.menuDebug = uimenu(handles.fig,'Label','Debug');
    handles.menuItemDebug = uimenu(handles.menuDebug, ...
				   'Label','Keyboard', ...
				   'Interruptible', 'off', ...
				   'Callback', @runDebug);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  handles.measurePointMenu = uicontextmenu;

  uimenu(handles.measurePointMenu, ...
	 'Label', 'Move measure point',...
	 'Interruptible', 'off', ...
	 'Callback',@moveMeasureCallback);

  uimenu(handles.measurePointMenu, ...
	 'Label', 'Delete measure point',...
	 'Interruptible', 'off', ...
	 'Callback',@deleteMeasureCallback);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Some dirty tricks to get the zoom and pan tools
  % Set up toolbar, only keep zoom in, zoom out and pan
  handles.toolbar = findall(handles.fig,'Type','uitoolbar');
  oldChild = allchild(handles.toolbar);

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
  showLogo();

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function aboutSynD(source, event)

    SynDversion = getVersion();

    fprintf('Running version %s of SynD.\n', SynDversion);

    if~(isempty(bilder.logo))
      miniLogo = imresize(bilder.logo, [1 1]*50);
    else
      miniLogo = [];
    end

    aboutMsg = sprintf(['You are running SynD version %s\n\n' ...
 	  	        'Johannes Hjorth\n' ...
		        'Sabine Schmitz\n\n\n' ...
                        'Reference:\n%s\n\n' ...
                        'DOI: %s'], ...
		       SynDversion, referenceStr, referenceDOI);

    uiwait(msgbox(aboutMsg,'About SynD','Custom',miniLogo));

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function SynDversion = getVersion()

    try
     fid = fopen(strcat(SyndScriptPath,verFile),'r');
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

  function newVersion = askLatestVersion()

    try
      newVersion = urlread(strcat(urlBase,verFile));

      if(newVersion(end) == char(10))
        newVersion = newVersion(1:end-1);
      end
    catch
      disp('Unable to check latest version')
      newVersion = '(unavailable)';
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function getLatestVersion(source, event)

    try
      curVersion = getVersion();
      newVersion = askLatestVersion();

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

        urlwrite(strcat(urlBase,newFile),strcat(outPath,outFile));

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

  function runDebug(source, event)
    disp('Type return to exit debug mode')
    keyboard
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function restartSynD(source, event)
      userInput = questdlg('Do you want to restart SynD?', ...
			   'Restart imminent', ...
			   'Yes','No','Maybe','Yes');

      switch(userInput)
        case 'Yes'
          disp('Restarting SynD.')
          close(handles.fig)
	  SynD_extract();
        case 'No'
 	  disp('Restart aborted.')
        case 'Maybe'
	if(rand(1) < 0.5)
          disp('Restarting SynD.')
          close(handles.fig)
	  SynD_extract();
        else
	  disp('Not restarting SynD.')
        end
      end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function toggleRed(source, event)
    % If this function was called without a source, then toggle red
    % otherwise read in status from the clicked object
    if(~exist('source') | source ~= handles.red)
      set(handles.red,'Value',~get(handles.red,'Value'));
    end

    dispInfo.showRed = get(handles.red,'Value');
  
    showImage();
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function toggleGreen(source, event)
    if(~exist('source') | source ~= handles.green)
      set(handles.green,'Value',~get(handles.green,'Value'));
    end

    dispInfo.showGreen = get(handles.green,'Value');
    showImage();
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function toggleBlue(source, event)
    if(~exist('source') | source ~= handles.blue)
      set(handles.blue,'Value',~get(handles.blue,'Value'));
    end

    dispInfo.showBlue = get(handles.blue,'Value');
    showImage();
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function toggleMask(source, event)
    if(~exist('source') | source ~= handles.mask)
      set(handles.mask,'Value',~get(handles.mask,'Value'));
    end

    dispInfo.showMask = get(handles.mask,'Value');
  
    showImage();
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setRGBvisible(redFlag,greenFlag,blueFlag)
    dispInfo.showRed = redFlag;
    set(handles.red,'Value',redFlag);
    dispInfo.showGreen = greenFlag;
    set(handles.green,'Value',greenFlag);
    dispInfo.showBlue = blueFlag;
    set(handles.blue,'Value',blueFlag);

    showImage();
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setChannelVisible(morphFlag,synFlag,XFlag)

    visibleFlags = [0 0 0];

    if(morphFlag)
      visibleFlags(detection.remapChan(1)) = 1;
    end

    if(synFlag)
      visibleFlags(detection.remapChan(2)) = 1;
    end

    if(XFlag)
      visibleFlags(detection.remapChan(3)) = 1;
    end

    setRGBvisible(visibleFlags(1),visibleFlags(2),visibleFlags(3));

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setXYres(source, event)
    tmp = 1e-6*str2num(get(handles.XYres,'String'));

    if(~isempty(tmp) & 0 < tmp & tmp < Inf)
      data.xyRes = tmp;
      if(isempty(data.image))
        uicontrol(handles.loadNeuron);
      else
	uicontrol(handles.nextStage);
      end
    else
      uiwait(warndlg(sprintf('Incorrect resolution %s, using %s micrometers instead,', ...
			     get(handles.XYres,'String'), ...
			     num2str(1e6*data.xyRes)), ...
		     'Input error','modal'));

      set(handles.XYres, 'String', num2str(data.xyRes*1e6))

    end

    fprintf('XY-resolution set to %d m\n', data.xyRes)

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setMorphChannel(source, event)
    detection.remapChan(1) = get(handles.morphChannel,'Value');

    if(strcmp(detection.enableRemapChan,'off'))
      % Use default: R maps to 1, G maps to 2, B maps to 3
      set(handles.morphChannelNumber,'Value',detection.remapChan(1));
      detection.morphChannel = detection.remapChan(1);

    end

    % Tell the user to redo everything...
    dispInfo.stage = 2;

    calcMaxIntensity();

    % Reset scaling
    dispInfo.scaleRed = 1;
    dispInfo.scaleBlue = 1;
    dispInfo.scaleGreen = 1;

    showGUI(dispInfo.state);
    showImage();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setSynChannel(source, event)
    detection.remapChan(2) = get(handles.synChannel,'Value');
  
    if(strcmp(detection.enableRemapChan,'off'))
      % Use default: R maps to 1, G maps to 2, B maps to 3
      set(handles.synChannelNumber,'Value',detection.remapChan(2));
      detection.synChannel = detection.remapChan(2);

    end

    % Tell the user to redo everything...
    dispInfo.stage = 2;

    calcMaxIntensity();

    % Reset scaling
    dispInfo.scaleRed = 1;
    dispInfo.scaleBlue = 1;
    dispInfo.scaleGreen = 1;

    showGUI(dispInfo.state);
    showImage();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setXChannel(source, event)
    detection.remapChan(3) = get(handles.XChannel,'Value');

    if(strcmp(detection.enableRemapChan,'off'))
      % Use default: R maps to 1, G maps to 2, B maps to 3
      set(handles.XChannelNumber,'Value',detection.remapChan(3));
      detection.XChannel = detection.remapChan(3);

    end

    % Tell the user to redo everything...
    dispInfo.stage = 2;

    calcMaxIntensity();

    % Reset scaling
    dispInfo.scaleRed = 1;
    dispInfo.scaleBlue = 1;
    dispInfo.scaleGreen = 1;

    showGUI(dispInfo.state);
    showImage();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function remapChannels(source, event)

    detection.morphChannel = get(handles.morphChannelNumber,'Value');
    detection.synChannel = get(handles.synChannelNumber,'Value');
    detection.XChannel = get(handles.XChannelNumber,'Value');

    % Tell the user to redo everything...
    dispInfo.stage = 2;

    calcMaxIntensity();

    % Reset scaling
    dispInfo.scaleRed = 1;
    dispInfo.scaleBlue = 1;
    dispInfo.scaleGreen = 1;

    showGUI(dispInfo.state);
    showImage();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function toggleRemap(source, event)
    if(get(handles.remap,'Value'))
      detection.enableRemapChan = 'on';

    else
      detection.enableRemapChan = 'off';

      set(handles.morphChannelNumber,'Value',detection.remapChan(1));
      detection.morphChannel = detection.remapChan(1);

      set(handles.synChannelNumber,'Value',detection.remapChan(2));
      detection.synChannel = detection.remapChan(2);

      set(handles.XChannelNumber,'Value',detection.remapChan(3));
      detection.XChannel = detection.remapChan(3);

      showImage();

      calcMaxIntensity();

      % Reset scaling
      dispInfo.scaleRed = 1;
      dispInfo.scaleBlue = 1;
      dispInfo.scaleGreen = 1;

      % Tell the user to redo everything...
      dispInfo.stage = 2;
      showGUI(dispInfo.state);

    end

    set(handles.morphChannelNumber,'enable',detection.enableRemapChan);
    set(handles.synChannelNumber,'enable',detection.enableRemapChan);
    set(handles.XChannelNumber,'enable',detection.enableRemapChan);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function toggleSingleSoma(source, event)
    % If this function was called without a source, then toggle red
    % otherwise read in status from the clicked object
    if(~exist('source') | source ~= handles.singleSoma)
      set(handles.singleSoma,'Value',~get(handles.singleSoma,'Value'));
    end

    detection.singleSoma = get(handles.singleSoma,'Value');

    detectSoma();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setMorphThresh(source, event)
    tmp = str2num(get(handles.morphThresh,'String'));

    if(~isempty(tmp) & 0 < tmp & tmp < Inf)
      detection.morphThreshold = tmp;
      disp(sprintf('Using morphology threshold %d', detection.morphThreshold))
    
      if(exist('source') & source == handles.morphThresh)
        detectSoma();
      end
    else
      uiwait(warndlg(sprintf('Incorrect threshold %s, using %s instead.', ...
			     get(handles.morphThresh,'String'), ...
			     num2str(detection.morphThreshold)), ...
		     'Input error','modal'));

      set(handles.morphThresh,'String', num2str(detection.morphThreshold));
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setSomaErodeRadius(source, event)
    tmp = str2num(get(handles.somaErodeRadius,'String'));

    if(~isempty(tmp) & 0 < tmp & tmp < Inf)
      detection.somaErodeRadius = tmp;
      disp(sprintf('Using soma erode radius %d', detection.somaErodeRadius))

      % Only call detect soma if the user pressed enter
      % (because detectSoma also calls this function to make sure it
      %  has latest user input).
      if(exist('source') & source == handles.somaErodeRadius)
        detectSoma();
      end

    else
      uiwait(warndlg(sprintf('Incorrect radius %s, using %s instead.', ...
			     get(handles.somaErodeRadius,'String'), ...
			     num2str(detection.somaErodeRadius)), ...
		     'Input error','modal'));

      set(handles.somaErodeRadius,'String', num2str(detection.somaErodeRadius));
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setSomaMeasureRadius(source, event)
    tmp = str2num(get(handles.somaMeasureRadius,'String'));
    if(~isempty(tmp) & 0 < tmp & tmp < Inf)
      detection.measurePointRadius = tmp;
    else
      uiwait(warndlg(sprintf('Incorrect radius %s, using %d.', ...
			     get(handles.somaMeasureRadius,'String'), ...
			     detection.measurePointRadius)));
    set(handles.somaMeasureRadius,'String', ...
	num2str(detection.measurePointRadius));
    end

    makeMeasureMask();

    showImage();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setSynThresh(source, event)
    tmp = str2num(get(handles.synThresh,'String'));

    [detection.synThreshold, modFlag] = ...
      sanitiseInput(get(handles.synThresh,'String'), ...
		   detection.synThreshold, ...
		   -inf, inf, false, true);

    if(modFlag)
      uiwait(warndlg(sprintf('Incorrect threshold %s, using %s instead.', ...
			     get(handles.synThresh,'String'), ...
			     num2str(detection.synThreshold)), ...
		     'Input error','modal'));

      set(handles.synThresh, 'String', num2str(detection.synThreshold));

    else
      if(exist('source') & source == handles.synThresh)
        if(isnan(detection.synThreshold))
          disp('Using cross-entropy for synapse threshold')
        else
          fprintf('Using synapse threshold %d\n', detection.synThreshold)
        end

        detectSynapses();

      end

    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setNeuritePadding(source, event)
    tmp = str2num(get(handles.neuritePadding, 'String'))*1e-6;

    if(~isempty(tmp) & 0 <= tmp & tmp < Inf)
      detection.neuritePadding = tmp;
      disp(sprintf('Using %d pixels neurite padding', ...
		   round(detection.neuritePadding/data.xyRes)))

      if(exist('source') & source == handles.neuritePadding)
        detectSynapses();
      end
    else
      uiwait(warndlg(sprintf('Incorrect neurite padding %s, using %s instead.', ...
			     get(handles.neuritePadding,'String'), ...
			     num2str(detection.neuritePadding*1e6)), ...
		     'Input error','modal'));

      set(handles.neuritePadding, 'String', num2str(detection.neuritePadding*1e6));

    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setSynapseMinSize(source, event)

    tmp = str2num(get(handles.minSynapseSize, 'String'))*1e-12;

    if(~isempty(tmp) & 0 <= tmp & tmp < Inf)

      detection.minSynapseSize = tmp;

      if(exist('source') & source == handles.minSynapseSize)
        detectSynapses();
      end

    else
      set(handles.minSynapseSize,'String', ...
	  num2str(detection.minSynapseSize*1e12));
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setGrowThresh(source, event)
    tmp = str2num(get(handles.growThresh,'String'));

    if(~isempty(tmp) & 0 < tmp & tmp <= 1)
      detection.maxAddCost = tmp;
      disp(sprintf('Using maximal grow cost of %.4f', detection.maxAddCost))

      if(exist('source') & source == handles.growThresh)
        detectNeurite();
      end
    else

      uiwait(warndlg(sprintf('Incorrect cost %s, using %s instead.', ...
			     get(handles.growThresh,'String'), ...
			     num2str(detection.maxAddCost)), ...
		     'Input error','modal'));

      set(handles.growThresh, 'String', num2str(detection.maxAddCost));
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setFilterSize(source, event)

    tmp = str2num(get(handles.filterSize,'String'));

    if(~isempty(tmp) & 0.1 < tmp & tmp <= 20)
      detection.filterSize = tmp*1e-6;
      disp(sprintf('Using filter size of %.1f micrometer', ...
		   detection.filterSize*1e6))

      if(exist('source') & source == handles.filterSize)
        detectNeurite();
      end
    else

      uiwait(warndlg(sprintf('Incorrect filter sizet %s, using %s micrometer instead.', ...
			     get(handles.filterSize,'String'), ...
			     num2str(detection.filterSize*1e6)), ...
		     'Input error','modal'));

      set(handles.filterSize, 'String', num2str(detection.filterSize*1e6));
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function toggleSteerableFilters(source, event)
    if(detection.useSteerableFilters)
      detection.useSteerableFilters = 0;
      set(handles.menuItemSteerable,'Checked','off')
    else
      detection.useSteerableFilters = 1;
      set(handles.menuItemSteerable,'Checked','on')
    end

    showGUI(dispInfo.state);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function toggleAutoSomaThreshold(source, event)

    if(detection.autoSomaThreshold)
      detection.autoSomaThreshold = 0;
    else
      detection.autoSomaThreshold = 1;
    end

    set(handles.autoThreshold, 'Checked', onOff(detection.autoSomaThreshold))

    set(handles.morphThresh, 'Enable', onOff(~detection.autoSomaThreshold));

    showGUI(dispInfo.state);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function toggleBackgroundRemoval(source, event)
    if(detection.backgroundFiltering)
      detection.backgroundFiltering = false;
      set(handles.menuItemBackgroundFiltering,'Checked','off');
    else
      detection.backgroundFiltering = true;
      set(handles.menuItemBackgroundFiltering,'Checked','on');
    end

    showGUI(dispInfo.state);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function status = onOff(flag)
    if(flag)
      status = 'on';
    else
      status = 'off';
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function editSoma(source, event)

    editInfo.mode = 2;
    editInfo.defaultWidth = 8;

    stopEdit(); 
    set(handles.fig,'windowbuttondownfcn', @startDrawLine)

    dispInfo.somaColor = [1 1 1];
    dispInfo.neuriteColor = NaN;
    dispInfo.synapseColor = NaN;
    dispInfo.measurePointColor = NaN;
    showImage();

    % The user is responsible for drawing a soma if they click this button
    activateStage('neurites');

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function addSoma(source, event)

    set(handles.fig,'CurrentAxes',handles.image)    
    saveUndo('add soma');

    stopEdit(); % Clears mouse handlers

    set(somaGUI,'enable','off')
    set(handles.allIcons,'enable','off')
    mask = roipoly();
    set(somaGUI,'enable','on')
    set(handles.allIcons,'enable','on')

    if(~isempty(mask))
      data.somaMask = data.somaMask | mask;
      activateStage('neurites');
    end

    showImage();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function deleteSoma(source, event)
    
    disp('Select a soma to delete.')
    [x,y,button] = ginput(1);

    a = axis;

    while(a(1) <= x & x <= a(2) ...
	  & a(3) <= y & y <= a(4) ...
	  & 1 <= x & x <= data.width ...
	  & 1 <= y & y <= data.height ...
	  & 1 == button)

      saveUndo('delete soma');

      data.somaMask = cleanBlob(data.somaMask,round(x),round(y));

      data.somaMeasureMask = [];

      % Update image
      showImage();

      % Check if we should remove one more
      disp('Select soma to remove.')
      [x,y,button] = ginput(1);
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function mask = cleanBlob(mask,x,y)

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

  function setSomaMeasurePoints(source, event)

    stopEdit();

    if(isempty(data.somaMask))
      disp('No soma detected yet.')
      return
    end

    measureDisk = strel('disk',detection.measurePointRadius);
    measureDiskPadded = strel('disk',detection.measurePointRadius*2);

    saveUndo('new soma measure');

    measurePoint = [];

    [yS,xS] = ind2sub(size(data.somaMask),find(data.somaMask));
    somaRangeX = min(xS)-1:max(xS)+1;
    somaRangeY = min(yS)-1:max(yS)+1;

    somaRangeX = min(max(1,somaRangeX),size(data.somaMask,2));
    somaRangeY = min(max(1,somaRangeY),size(data.somaMask,1));

    % Zoom in on the soma
    dispInfo.axis = [min(xS)-1, max(xS)+1, min(yS)-1, max(yS)+1];
    axis(dispInfo.axis);

    miniSomaMask = data.somaMask(somaRangeY,somaRangeX);
    possiblePoints = imerode(miniSomaMask,measureDisk);
    possiblePointIdx = find(possiblePoints);  

    removeMask = zeros(size(possiblePoints));

    warnedUser = 0;

    for i = 1:detection.nMeasurePoints
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
			       detection.nMeasurePoints-i+1), ...  
		       'SynD: Not enough space in soma'));
        warnedUser = 1;
        break;
      end

    end

    % Remap back to big soma mask
    [yM,xM] = ind2sub(size(miniSomaMask),measurePoint);
    yM = yM + min(somaRangeY)-1;
    xM = xM + min(somaRangeX)-1;

    data.measurePoint = sub2ind(size(data.somaMask),yM,xM);

    makeMeasureMask();

    % Display the mask

    dispInfo.somaColor = NaN;
    dispInfo.measurePointColor = dispInfo.defaultMeasurePointColor;
    dispInfo.neuriteColor = NaN;
    dispInfo.synapseColor = NaN;

    showImage();

    exportInfo.saveSomaMeasure = 1;
    set(handles.saveSomaMeasure,'Value',exportInfo.saveSomaMeasure);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function makeMeasureMask()

    measureDisk = strel('disk',detection.measurePointRadius);

    data.somaMeasureMask = zeros(size(data.somaMask));
    data.somaMeasureMask(data.measurePoint) = 1;
    data.somaMeasureMask = imdilate(data.somaMeasureMask,measureDisk);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function moveMeasureCallback(source, event)

    pointIdx = findClosestMeasurePoint();

    saveUndo('move soma measure');
    
    disp('Move soma measure.')
    [x,y,idx] = getImagePoint();

    if(isempty(idx))
      % Nothing to do, abort
      return
    end

    data.measurePoint(pointIdx) = idx;
    makeMeasureMask();

    dispInfo.somaColor = NaN;
    dispInfo.measurePointColor = dispInfo.defaultMeasurePointColor;
    dispInfo.neuriteColor = NaN;
    dispInfo.synapseColor = NaN;

    showImage();
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function deleteMeasureCallback(source, event)

    saveUndo('delete soma measure');

    pointIdx = findClosestMeasurePoint();
    data.measurePoint(pointIdx) = [];
    makeMeasureMask();

    dispInfo.somaColor = NaN;
    dispInfo.measurePointColor = dispInfo.defaultMeasurePointColor;
    dispInfo.neuriteColor = NaN;
    dispInfo.synapseColor = NaN;

    showImage();
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [x,y,idx] = getImagePoint()
    [x,y] = ginput(1);
    x = round(x);
    y = round(y);

    idx = sub2ind(size(data.somaMask),y,x);

    a = axis;
    if(x < a(1) | a(2) < x | y < a(3) | a(4) < y)
      % We are outside the image
      x = []; y = []; idx = [];
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function pointIdx = findClosestMeasurePoint()

    tmpXY = get(handles.image,'CurrentPoint');
    x = round(tmpXY(1,1));
    y = round(tmpXY(1,2));

    [yP,xP] = ind2sub(size(data.somaMask),data.measurePoint);
    Pdist = sqrt((xP-x).^2+(yP-y).^2);

    [minValue, pointIdx] = min(Pdist);
    % [~, pointIdx] = min(Pdist); % Only works in 2009b and later

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function addSomaMeasurePoints(source, event)

    if(isempty(data.somaMask))
      disp('No soma mask set yet.')
      return
    end

    saveUndo('add soma measure');

    % Zoom in on the soma

    [yS,xS] = ind2sub(size(data.somaMask),find(data.somaMask));

    dispInfo.axis = [min(xS)-1, max(xS)+1, min(yS)-1, max(yS)+1];
    axis(dispInfo.axis);

    dispInfo.somaColor = NaN;
    dispInfo.measurePointColor = dispInfo.defaultMeasurePointColor;
    dispInfo.neuriteColor = NaN;
    dispInfo.synapseColor = NaN;

    showImage();

    % Ask user for new start point
    disp('Add some measure')
    [x,y,idx] = getImagePoint();

    if(isempty(idx))
      % Nothing to do, abort
      return
    end

    data.measurePoint(end+1) = idx;
    makeMeasureMask();

    showImage();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function editNeurite(source, event)

    editInfo.mode = 1;
    editInfo.defaultWidth = 3;

    stopEdit();
    set(handles.fig,'windowbuttondownfcn', @startDrawLine)

    dispInfo.somaColor = [1 1 1];
    dispInfo.neuriteColor = [1 1 1]*0.7;
    dispInfo.synapseColor = NaN;
    dispInfo.measurePointColor = NaN;

    showImage();

    % The user is responsible for drawing neurites if they click this button
    activateStage('synapses');

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function editSynapse(source, event)

    editInfo.mode = 3;
    editInfo.defaultWidth = 2;

    stopEdit();
    set(handles.fig,'windowbuttondownfcn', @startDrawLine)

    dispInfo.somaColor = NaN;
    dispInfo.neuriteColor = NaN;
    dispInfo.synapseColor = [1 1 1];
    dispInfo.measurePointColor = NaN;
    showImage();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function stopEdit()

    % editInfo.mode = 0;
    set(handles.fig,'windowbuttondownfcn', []);
    set(handles.fig,'windowbuttonmotionfcn', []);
    set(handles.fig,'windowbuttonupfcn', []);

    showImage();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function killSynBlob(source, event)

    disp('Select blob to remove.')
    [x,y,button] = ginput(1);

    a = axis;

    while(a(1) <= x & x <= a(2) ...
	  & a(3) <= y & y <= a(4) ...
	  & 1 <= x & x <= data.width ...
	  & 1 <= y & y <= data.height ...
	  & 1 == button)

      saveUndo('synapse blob removal');

      data.synapseMask = cleanBlob(data.synapseMask,round(x),round(y));

      % Update image
      showImage();

      disp('Select blob to remove.')
      % Check if we should remove one more
      [x,y,button] = ginput(1);
    end

    locateSynapseCenters();
    filterSynapseMask(data.synapseMask);

    % Update synapse properties
    calculateSynapseProperties();
    showImage();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function steerableFilterSettings(source, event)

    prompt = { 'Sigma (in micrometers)', 'Max growth cost' };
    defaultVal = { num2str(detection.filterSize*1e6), ...
		   num2str(detection.maxAddCost) };
    dialogName = 'Steerable filter settings';
    numLines = 1;

    answers = inputdlg(prompt, dialogName, numLines, defaultVal);

    if(~isempty(answers))

      [detection.filterSize, modFlag1] = ...
	sanitiseInput(answers{1}, detection.filterSize*1e6, ...
		      0.1, 20, false);

      % Everything internally is stored in SI units.
      detection.filterSize = detection.filterSize * 1e-6;

      [detection.maxAddCost, modFlag2] = ...
	sanitiseInput(answers{2}, detection.maxAddCost, ...
		      0, 1, false);

      if(modFlag1 | modFlag2)
        warnMsg = sprintf(['Sigma must 1-20, max growth cost 0-1.']);
		   
        uiwait(warndlg(warnMsg, 'Input sanitation','modal'));
   
      end

      % Update GUI
      set(handles.growThresh,'String', num2str(detection.maxAddCost));
      set(handles.filterSize,'String', num2str(detection.filterSize*1e6));

      % Update the filters
      if(~isempty(data.image))
        steerableFilter();
      end

    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function synapseSettings(source, event)

    prompt = { 'Threshold (in std, NaN for auto)', ...
	       'Min synapse size (um^2)', ...
	       'Radius (in um, NaN for auto)', ...
	       'Neurite padding (um)', ...
               'Include soma synapses (default 0=false)' };
    defaultVal = { num2str(detection.synThreshold), ...
		   num2str(detection.minSynapseSize*1e12), ...
		   num2str(detection.singleSynapseRadie*1e6), ...
		   num2str(detection.neuritePadding*1e6), ...
		   num2str(~detection.excludeSomaSynapses) };
    dialogName = 'Synapse detection settings';
    numLines = 1;

    answers = inputdlg(prompt, dialogName, numLines, defaultVal);

    if(~isempty(answers))

      % Allows NaN = auto synapse threshold
      [detection.synThreshold, modFlag1] = ...
          sanitiseInput(answers{1}, detection.synThreshold, ...
                        -inf, inf, false, true);

      [tmpSize, modFlag2] = ...
	sanitiseInput(answers{2}, detection.minSynapseSize*1e12, ...
		      0, inf, false, false);
      detection.minSynapseSize = tmpSize / 1e12;

      % Allows NaN = auto synapse radie
      [tmpRadie, modFlag3] = ...
        sanitiseInput(answers{3}, detection.singleSynapseRadie*1e6, ...
		      0, inf, false, true);
      detection.singleSynapseRadie = tmpRadie*1e-6;

      [detection.neuritePadding, modFlag4] = ...
	sanitiseInput(answers{4}, detection.neuritePadding * 1e6, ...
		      0, 20, false, false);

      detection.neuritePadding = detection.neuritePadding * 1e-6;

      % Changed it to include soma synapses in settings, more intuitive
      [includeSomaSynapses,modFlag5] = ...
	sanitiseInput(answers{5}, ~detection.excludeSomaSynapses, ...
		      0, 1, true, false);
      detection.excludeSomaSynapses = ~includeSomaSynapses;

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
        set(handles.synThresh,'String', ...
	    num2str(detection.synThreshold));
        set(handles.neuritePadding, 'String', ...
	    num2str(detection.neuritePadding*1e6));
        set(handles.minSynapseSize,'String', ...
	    num2str(detection.minSynapseSize*1e12));
 
      end

      if(nnz(data.neuriteMask))
        detectSynapses();
      end

    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function shollSettings(source, event)

    prompt = { 'Bin size (in micrometers)' };
    defaultVal = { num2str(detection.shollBinSize*1e6) };

    dialogName = 'Sholl analysis settings';
    numLines = 1;

    answers = inputdlg(prompt, dialogName, numLines, defaultVal);

    if(~isempty(answers))

      [shollBinSize, modFlag1] = ...
	sanitiseInput(answers{1}, detection.shollBinSize*1e6, ...
		      1, 1000, false);

      detection.shollBinSize = shollBinSize*1e-6;

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

  function [newVal,modFlag] = sanitiseInput(inputString, oldVal, ...
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

  function startDrawLine(source, event)

    set(handles.fig,'CurrentAxes',handles.image)    

    % Where did the user click, within the image?
    xy = get(handles.image,'CurrentPoint');
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
    switch(get(handles.fig,'selectiontype'))
      case 'normal'
        editInfo.mode = abs(editInfo.mode); % Drawing
	editInfo.color = [1 1 1]*0.8;
        editInfo.width = editInfo.defaultWidth;
      case 'alt'
        editInfo.mode = -abs(editInfo.mode); % Erasing
        editInfo.color = [1 1 1]*0.2;
        editInfo.width = editInfo.defaultWidth*2;
      otherwise
        disp('Unknown key pressed')  
        disp(get(handles.fig,'selectiontype'))
        return
    end

    editInfo.line = line(x,y, ...
			 'color',editInfo.color, ...
			 'linewidth', editInfo.width);
    set(handles.fig,'windowbuttonmotionfcn',@drawLineCallback)
    set(handles.fig,'windowbuttonupfcn', @endDrawLine)

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function drawLineCallback(source, event)

    xy = get(handles.image,'CurrentPoint');
    x = xy(1,1);
    y = xy(1,2);

    a = axis();    

    if(x < a(1) | a(2) < x ...
       | y < a(3) | a(4) < y)
      % We are outside figure, end drawing
      endDrawLine();
     
    end

    try
      xOld = get(editInfo.line,'Xdata');
      yOld = get(editInfo.line,'Ydata');
    catch
      disp('Unable to read points from mouse.')
      endDrawLine();
      return
    end

    if(1 <= x & x <= data.width ...
       & 1 <= y & y <= data.height)
      editInfo.xLine = [xOld,x];
      editInfo.yLine = [yOld,y];

      set(editInfo.line, ...
	  'Xdata', editInfo.xLine, ...
	  'Ydata', editInfo.yLine);
    else
      % We were outside image, this should not happen
      % but put it here just to be on safe side.
      endDrawLine();
    end

    drawnow

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function endDrawLine(source, event)

    disp('End draw line called')

    % Remove drawing callback (but keep button down function)
    set(handles.fig,'windowbuttonmotionfcn', [])
    set(handles.fig,'windowbuttonupfcn', [])

    if(isempty(editInfo.xLine))
      % Nothing to draw
      return;
    end


    % Modify the mask
    tmpMask = zeros(size(data.neuriteMask));
   
    tmpMask(round(editInfo.yLine(1)),round(editInfo.xLine(1))) = 1;

    % There should be a built in function to connect points to lines
    % and return a binaries mask. This is a bad way to do a line, it will 
    % add extra points. Better to use Bresenhams algorithm...

    for i = 2:length(editInfo.xLine)
      nDots = 2*max(abs(editInfo.xLine(i)-editInfo.xLine(i-1)), ...
		    abs(editInfo.yLine(i)-editInfo.yLine(i-1)));

      xDots = round(linspace(editInfo.xLine(i-1),editInfo.xLine(i),nDots));
      yDots = round(linspace(editInfo.yLine(i-1),editInfo.yLine(i),nDots));
      iDots = sub2ind(size(tmpMask),yDots,xDots);

      tmpMask(iDots) = 1;
    end

    % Use appropriate width
    a = axis;
    zoomLevel = max((a(2)-a(1))/data.width,(a(4)-a(3))/data.height);
    brushSize = max(round(editInfo.width*zoomLevel),1);
    tmpMask = imdilate(tmpMask,strel('disk',brushSize));

    switch(editInfo.mode)
      case 1
        % Add the pixels to neurite mask
        saveUndo('add neurite');
        data.neuriteMask = double((data.neuriteMask + tmpMask) > 0);
        updateNeuriteMask();
      case -1
        % Remove the pixels from neurite mask
        saveUndo('remove neurite');
        data.neuriteMask = double(((data.neuriteMask>0) - tmpMask) > 0);
        updateNeuriteMask();
      case 2
        % Add the pixels to soma mask
        saveUndo('add soma');
        data.somaMask = double((data.somaMask + tmpMask) > 0);
      case -2
        % Remove the pixels from soma mask
        saveUndo('remove soma');
        data.somaMask = double((data.somaMask - tmpMask) > 0);
      case 3
        % Add the pixels to synapse mask
        saveUndo('add synapse');
        data.synapseMask = double((data.synapseMask + tmpMask) > 0);

        locateSynapseCenters();
        filterSynapseMask(data.synapseMask);

        % Update synapse properties
	calculateSynapseProperties();

      case -3
        % Remove the pixels from synapse mask
        saveUndo('remove synapse');
        data.synapseMask = double((data.synapseMask - tmpMask) > 0);

        locateSynapseCenters();
        filterSynapseMask(data.synapseMask);

        % Update synapse properties
	calculateSynapseProperties();

      otherwise
        % Do nothing
    end


    % Display the result
    data.distMask = makeDistMask(data.somaMask,data.neuriteMask);
    showImage();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function autoSomaThreshold(source, event)

    % Find threshold that minimizes cross-entropy
    grayImage = data.image(:,:,detection.morphChannel, dispInfo.curImg);
    detection.morphThreshold  = minimizeCrossEntropy(grayImage);

    % Update GUI
    set(handles.morphThresh,'String',detection.morphThreshold);

    if(exist('source'))
      detectSoma();
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  function threshold = minimizeCrossEntropy(grayImage)


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

    H1 = crossEntropy(G1,edges,freqG);
    H2 = crossEntropy(G2,edges,freqG);
    H3 = crossEntropy(G3,edges,freqG);

    Hctr = 3;

    while(abs(G1 - G3) > 0.5)

      G4 = G2  + resphi*(G3 - G2);

      H4 = crossEntropy(G4,edges,freqG);
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

  function H = crossEntropy(thresh, edges, freqG)

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

  function threshold = minimizeCrossEntropyOldSlow(grayImage)


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

  function detectNeuritesSimple()

    % Read in new settings from GUI using callback function
    setGrowThresh();
    setFilterSize();

    % Extract selected channel for morphology detection
    tmp = data.image(:,:,detection.morphChannel, dispInfo.curImg);

    if(~isnan(detection.wienerSize))
      tmp = wiener2(tmp,detection.wienerSize*[1 1]);
    end

    if(isempty(data.includeMask))
      data.includeMask = ones(size(tmp));
    end

    if(detection.morphThreshold < 1)
      % If between 0 and 1, use relative threshold
      mask = tmp > detection.morphThreshold*max(tmp(:));
    else
      mask = tmp > detection.morphThreshold;
    end

    mask = mask .* data.includeMask;

    se = strel('disk',1);
    mask = imopen(mask,se);

    data.neuriteMask = mask;

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function detectNeurite(source, event)

    saveUndo('detect neurites');

    if(detection.useSteerableFilters)
      % Replace the threshold mask for the neurite with that from
      % the steerable filters
      steerableFilter();
      growNeuriteFromSoma();
    else
      detectNeuritesSimple();
    end

    % Mark unconnected components in neurite mask
    updateNeuriteMask();

    showImage();

    activateStage('synapses');

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function extendNeuritesHandler(source, event)
    saveUndo('extend neurites');

    if(nnz(data.neuriteMask))
      % We require that there are some neurite parts detected
      extendNeurites(5e-6/data.xyRes);

    else
      disp('Detect neurites first!')
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % This function looks for neurite components which are a few pixels away

  function extendNeurites(seedDist)

    distMask = makeDistMask(data.somaMask,data.neuriteMask);

    % Locate neurite end points, set inf and NaN to 0 and find
    % largest remaining distances, those are our end points.
    tmpMask = distMask;
    tmpMask(distMask == inf) = 0;
    tmpMask(isnan(distMask)) = 0;
    [endY,endX] = find(imregionalmax(tmpMask));

    seedPoints = zeros(size(data.neuriteMask));

    seedList = [];
    endList = [];

    for iEnd = 1:length(endX)

      if(endX(iEnd) < seedDist | endX(iEnd) > data.width - seedDist ...
	 | endY(iEnd) < seedDist | endY(iEnd) > data.height - seedDist)
        % Too close to border, skip this point
        continue
      end

      % What is direction of the neurite at the end point?
      endDir = squeeze(data.dirVect(endY(iEnd),endX(iEnd),:));

      % Which way is forward, +endDir or -endDir
      % Probe neurite by taking small step in one direction
      checkX = round(endX(iEnd) + endDir(1)*2);
      checkY = round(endY(iEnd) + endDir(2)*2);
 
      % Did we go towards neurite, or away from neurite?
      if(data.neuriteMask(checkY,checkX))
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
    growNeurite(seedPoints);

    % See which neurites are large enough, and reconnect them

    tmpMask = zeros(size(data.neuriteMask));

    for iSeed = 1:size(seedList,1)
      seedX = seedList(iSeed,1);
      seedY = seedList(iSeed,2);

      % Is the seed point part of a neurite?
      if(data.neuriteMask(seedY,seedX))

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

    data.neuriteMask = data.neuriteMask + tmpMask;

    updateNeuriteMask();

    data.distMask = makeDistMask(data.somaMask,data.neuriteMask);

    % Show the final image

    showImage();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function updateNeuriteMask()

    data.neuriteMask = 2*(data.neuriteMask > 0);
    CC = bwconncomp(data.neuriteMask);
    somaIdx = find(data.somaMask);

    for i = 1:length(CC.PixelIdxList)
      if(nnz(ismember(CC.PixelIdxList{i},somaIdx)))
	data.neuriteMask(CC.PixelIdxList{i}) = 1;
      end
    end

    if(nnz(data.neuriteMask == 2))
      dispInfo.needCleaning = 1;
      set(handles.clean,'BackgroundColor', [1 0.6 0]);
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function cleanMask(source, event)

    saveUndo('clean neurite mask');

    % This function removes unconnected components from the neurite mask
    % Those with distance set at inf.

    data.neuriteMask(data.distMask .* data.neuriteMask == inf) = 0;
    data.neuriteMask(isnan(data.distMask .* data.neuriteMask)) = 0;

    dispInfo.needCleaning = 0;
    set(handles.clean,'BackgroundColor', [0.8 0.8 0.8]);

    showImage();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function killBlob(source, event)

    disp('Select blob to remove.')
    [x,y,button] = ginput(1);

    a = axis;

    while(a(1) <= x & x <= a(2) ...
	  & a(3) <= y & y <= a(4) ...
	  & 1 <= x & x <= data.width ...
	  & 1 <= y & y <= data.height ...
	  & 1 == button)

      saveUndo('blob removal');

      data.neuriteMask = cleanBlob(data.neuriteMask,round(x),round(y));

      % Update image
      showImage();

      disp('Select blob to remove.')
      % Check if we should remove one more
      [x,y,button] = ginput(1);
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function addNeurite(source, event)

    neuriteSeed = zeros(size(data.neuriteMask));

    disp('Add neurite.')
    [x,y,button] = ginput(1);
    
    a = axis();

    while(button == 1 ...
	  & 1 <= x & x <= data.width ...
	  & 1 <= y & y <= data.height ...
	  & a(1) <= x & x <= a(2) ...
	  & a(3) <= y & y <= a(4))

      neuriteSeed(round(y),round(x)) = 1;

      fprintf(['Detecting neurite at x=%d,y=%d ' ...
	       '(click outside, or click right button to stop)\n'], ...
	      round(x), round(y))

      saveUndo(sprintf('add neurite (%d,%d)', round(x),round(y)))
      growNeurite(neuriteSeed);

      disp('Add neurite.')
      [x,y,button] = ginput(1);

    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function addThinNeurites(source, event)

    saveUndo('add thin neurites');

    oldFilterSize = detection.filterSize;
    detection.filterSize = detection.filterSize / 2;
    steerableFilter();

    oldNeuriteMask = data.neuriteMask > 0;
    growNeurite(oldNeuriteMask - imerode(oldNeuriteMask,strel('disk',1)));
    
    addedNeuriteBits = (data.neuriteMask>0) - oldNeuriteMask;

    CC = bwconncomp(addedNeuriteBits);
    newProps = regionprops(CC, ...
			   'PixelIdxList', ...
			   'Area', ...
			   'MajorAxisLength');

    % If they are shorter than 2 micrometers, remove.
    for i = 1:length(newProps)
      if(newProps(i).MajorAxisLength < detection.minProtrusionLength/data.xyRes)
	data.neuriteMask(newProps(i).PixelIdxList) = 0;
      end
    end

    updateNeuriteMask();
    showImage();

    detection.filterSize = oldFilterSize;
    steerableFilter();

  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function distMask = makeDistMask(somaMask, neuriteMask)

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
    neighDist = [sqrt(2) 1 sqrt(2) 1 0 1 sqrt(2) 1 sqrt(2)]*data.xyRes;
    
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

  function synDist = calculateSynapseDist()

    paddedDistMask = makeDistMask(data.somaMask, padNeuriteMask());
				    
    synDist = paddedDistMask(data.synapseCenter);

    if(nnz(synDist == inf))
      disp('This should not happen, talk to Johannes! Argh...')
      keyboard
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function detectSynapses(source, event)

    % Get values from GUI
    setSynThresh();
    setNeuritePadding();
    setSynapseMinSize();

    saveUndo('detect synapses');

    % dispInfo.axis = [];

    detectSynapseCenters();

    calculateSynapseProperties();

    showImage();
    activateStage('analyse');

    text(20,20,sprintf('Found %d synapses', length(data.synapseCenter)), ...
	 'color','white')

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  debugFlag = 0;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function detectSynapseCenters()

    % This calculates a preliminary threshold
    makeSynapseTemplate(); 

    noiseMask = thresholdSynapseChannel();

    % This one will exclude soma if detection.excludeSomaSynapses = true
    paddedNeuriteMask = padNeuriteMask();

    % Synapse centres are on neurites or its padding, but outside soma
    data.synapseMask = noiseMask .* paddedNeuriteMask;

    locateSynapseCenters();

    filterSynapseMask(noiseMask);

    if(debugFlag)
      keyboard
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function noiseMask = thresholdSynapseChannel()

    synChan = data.image(:,:,detection.synChannel, dispInfo.curImg);

    neuriteIdx = find(data.neuriteMask);

    if(detection.backgroundFiltering)
      synChan = ...
        max(0, synChan -imfilter(synChan, ...
				 fspecial('gaussian', ...
					  [3 3]*detection.smoothingRadius, ...
					  detection.smoothingRadius)));
    end

    synOnNeurite = synChan(neuriteIdx);

    % Threshold based detection of the synapse pixels
    if(~isnan(detection.synThreshold))
      if(detection.synThreshold < 0)
        disp(['Negative threshold specified, using the absolute value as absolute threshold'])
        synThresh = abs(detection.synThreshold);
      else
        % Normal case
        synThresh = mean(synOnNeurite(:)) ...
            + detection.synThreshold*std(synOnNeurite(:));
      end
      
    else
      % Auto detect threshold
      synThresh = minimizeCrossEntropy(synOnNeurite);
      fprintf('Using cross-entropy. ')
    end

    % Save the actual threshold used
    data.synapseIntensityThreshold = synThresh;

    fprintf('Absolute synapse intensity threshold used: %.1f\n', synThresh)

    % Note synChan might be filtered (see above)
    noiseMask = synChan > synThresh;

    % Remove isolated pixels by enforcing minimum size on synapses
    CCs = bwconncomp(noiseMask);

    for i = 1:length(CCs.PixelIdxList)
      if(length(CCs.PixelIdxList{i}) < detection.minSynapseSize/data.xyRes^2)
       noiseMask(CCs.PixelIdxList{i}) = 0;
      end
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function paddedNeuriteMask = padNeuriteMask()

    if(detection.excludeSomaSynapses) 
      % This is default
      paddedNeuriteMask = imdilate((data.neuriteMask>0) - data.somaMask > 0, ...
				   strel('disk', ...
					 round(detection.neuritePadding ...
					       /data.xyRes)));

      paddedNeuriteMask = (paddedNeuriteMask - data.somaMask) > 0;
    else
      %Include soma
      paddedNeuriteMask = imdilate(data.neuriteMask>0, ...
				   strel('disk', ...
					 round(detection.neuritePadding ...
					       /data.xyRes)));


    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function makeSynapseTemplate()

    noiseMask = thresholdSynapseChannel();

    paddedNeuriteMask = padNeuriteMask();

    morphChan = data.image(:,:,detection.morphChannel, dispInfo.curImg);
    synChan = data.image(:,:,detection.synChannel, dispInfo.curImg);
    XChan = data.image(:,:,detection.XChannel, dispInfo.curImg);

    % Synapses centres are on neurites, but outside soma
    data.synapseMask = noiseMask .* paddedNeuriteMask;

    putativeCenters = imregionalmax(synChan) .* noiseMask;
    putativeCenters = putativeCenters .* data.synapseMask;
    putativeCentersIdx = find(putativeCenters);

    CCp = bwconncomp(putativeCenters);

    labeledSynapses = bwlabel(data.synapseMask);

    % Look for synapses with just one clear maxima, 
    % and also only one center per synapse
    % used for auto-detecting synapse kernel
    centerMask = zeros(size(data.neuriteMask));

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

    meanSynapseMorph = zeros(detection.maxRadie*2+1,detection.maxRadie*2+1);
    meanSynapseSyn = zeros(detection.maxRadie*2+1,detection.maxRadie*2+1);
    meanSynapseX = zeros(detection.maxRadie*2+1,detection.maxRadie*2+1);

    synCtr = 0;
    for i = 1:length(xC)
      xReg = (xC(i)-detection.maxRadie):(xC(i)+detection.maxRadie);
      yReg = (yC(i)-detection.maxRadie):(yC(i)+detection.maxRadie);

      if(1 <= min(xReg) & max(xReg) <= data.width ...
         & 1 <= min(yReg) & max(yReg) <= data.height)
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

      if(nnz(synChan) & isnan(detection.singleSynapseRadie))
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

      [y,x] = meshgrid(-detection.maxRadie:detection.maxRadie, ...
		       -detection.maxRadie:detection.maxRadie);
      meanSynapse = exp(-x.^2-y.^2);

      data.meanSynapseMorph = [];
      data.meanSynapseSyn = [];
      data.meanSynapseX = [];
    else
      data.meanSynapseMorph = meanSynapseMorph;
      data.meanSynapseSyn = meanSynapseSyn;
      data.meanSynapseX = meanSynapseX;

      % Use average synapse for auto-detected template
      meanSynapse = meanSynapseSyn;;
    end

    % Did the user specify their own synapse template, if so use that 
    % for deconvolution instead.
    if(~isnan(detection.singleSynapseRadie))

      fprintf('Using template synapse, radie %.2d micrometers\n', ...
	      1e6*detection.singleSynapseRadie)
      r = detection.singleSynapseRadie/data.xyRes;

      [y,x] = meshgrid(-detection.maxRadie:detection.maxRadie, ...
		       -detection.maxRadie:detection.maxRadie);
      meanSynapse = exp(-(x/r).^2-(y/r).^2);

    end

    % Template used for synapse detection
    data.meanSynapse = meanSynapse / max(meanSynapse(:));

    % How does the intensity decrease with distance from synapse center
    % for single synapses.
    calculateSynapseProfile();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function calculateSynapseProfile()

    if(~isempty(data.meanSynapseSyn))

      [y,x] = meshgrid(-detection.maxRadie:detection.maxRadie, ...
		       -detection.maxRadie:detection.maxRadie);

      % Distance to center for each pixel in template
      d = sqrt(x.^2 + y.^2)*data.xyRes;

      binSize = data.xyRes; %detection.synapseProfileBinSize;
      distCenters = 0:binSize:min(3e-6,max(d(:)));

      valMorph = NaN*zeros(size(distCenters));  % Morphology staining
      valSyn   = NaN*zeros(size(distCenters));  % Synapse staining
      valX     = NaN*zeros(size(distCenters));  % X staining intensity

      for i = 1:length(distCenters)
        idx = find(distCenters(i)-binSize/2 <= d ...
		   & d < distCenters(i)+binSize/2);
        valMorph(i) = mean(data.meanSynapseMorph(idx));
        valSyn(i) = mean(data.meanSynapseSyn(idx));
        valX(i) = mean(data.meanSynapseX(idx));
      end

      data.meanSynapseProfileDist = distCenters;
      data.meanSynapseProfileMorph = valMorph;
      data.meanSynapseProfileSyn = valSyn;
      data.meanSynapseProfileX = valX;

      plotCol = [1 0 0; 0 1 0; 0 0 1];

      set(handles.fig,'CurrentAxes',handles.singleSynapse)

      if(dispInfo.showSynapseProfileDataPoints)
	plot(d(:)*1e6,data.meanSynapseMorph(:), ...
	     '.', 'color',plotCol(detection.remapChan(1),:));
        hold on
	plot(d(:)*1e6,data.meanSynapseSyn(:), ...
	     '.', 'color',plotCol(detection.remapChan(2),:));
	plot(d(:)*1e6,data.meanSynapseX(:), ...
	     '.', 'color',plotCol(detection.remapChan(3),:));
      end

      plot(distCenters*1e6,valMorph,'-', ...
	   'linewidth',3, 'color',plotCol(detection.remapChan(1),:));
      hold on
      plot(distCenters*1e6,valSyn,'-', ...
	   'linewidth',3, 'color',plotCol(detection.remapChan(2),:));
      plot(distCenters*1e6,valX,'-', ...
	   'linewidth',3, 'color',plotCol(detection.remapChan(3),:));

      hold off

      title('Single synapse intensity profile')
      xlabel('Distance (\mum)')
      ylabel('Intensity')

    else
      % No synapses, make sure all variables are cleared
      data.meanSynapseProfileDist = [];
      data.meanSynapseProfileMorph = [];
      data.meanSynapseProfileSyn = [];
      data.meanSynapseProfileX = [];

      set(handles.fig,'CurrentAxes',handles.singleSynapse)
      cla
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function locateSynapseCenters()

    % This function assumes the synapse mask and synapse template is already 
    % calculated it only detects the centers.

    synChan = data.image(:,:,detection.synChannel, dispInfo.curImg);

    % Deconvolve synapse image with the mean single synapse
    % [qImg,rImg] = deconvreg(synChan,data.meanSynapse);
    qImg = deconvlucy(synChan,data.meanSynapse);
    % figure, imagesc(qImg)

    % Find the local maximas, that are above noise level in 
    % synapse channel, and within neurite (excluding soma).
    qMask = imregionalmax(qImg);
    qMask = qMask .* data.synapseMask;
    %qMask = imregionalmax(qImg.*data.synapseMask);


    % This just marks centres of synapses, not the entire synapse.
    % The export function assumes there are only centres masked right now,
    % but perhaps change it to data.synapseMaskCentre
    data.synapseCenter = find(qMask);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function filterSynapseMask(mask)

    % Important, the synapse can go outside the neurite marker
    % so we want to keep those parts also, but still exclude soma!

    % For auto detect, mask = noiseMask are pixels above synapse noise level

    % Update: They have to be within detection.neuritePadding
    % of the neurite

    paddedNeuriteMask = padNeuriteMask();

    CC = bwconncomp(mask .* paddedNeuriteMask);

    data.synapseMask = zeros(size(data.synapseMask));

    neuriteIdx = find(data.neuriteMask);

    for i = 1:length(CC.PixelIdxList)
      % Require that putative synapse regions contain at least one center
      % and that in addition to all being within the padded region around
      % a neurite, at least one pixel has to be directly on the neurite.
      if(nnz(ismember(data.synapseCenter,CC.PixelIdxList{i})) ...
	 & nnz(ismember(neuriteIdx,CC.PixelIdxList{i})))
        data.synapseMask(CC.PixelIdxList{i}) = 1;
      else
	% fprintf('Orphaned pixels: %d\n', length(CC.PixelIdxList{i}))
      end
    end


  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function calculateSomaMeasures()

    % Calculate some measures, radie, area etc

    somaProp = regionprops(data.somaMask, ...
			   'Area', ...
			   'MajorAxisLength', ...
			   'MinorAxisLength');

    data.somaArea = cat(1,somaProp.Area)*data.xyRes^2*1e12;
    data.somaMajorAxisLength = cat(1,somaProp.MajorAxisLength)*data.xyRes*1e6;
    data.somaMinorAxisLength = cat(1,somaProp.MinorAxisLength)*data.xyRes*1e6;

    if(isempty(data.somaMeasureMask))
      disp('No soma measure regions marked.')
      return
    end

    tmpMorph = data.image(:,:,detection.morphChannel, dispInfo.curImg);
    data.somaMeasureMorph = tmpMorph(find(data.somaMeasureMask));

    tmpSyn = data.image(:,:,detection.synChannel, dispInfo.curImg);
    data.somaMeasureSyn = tmpSyn(find(data.somaMeasureMask));

    tmpX = data.image(:,:,detection.XChannel, dispInfo.curImg);
    data.somaMeasureX = tmpX(find(data.somaMeasureMask));


  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function calculateSynapseProperties()

    disp('Calculating synapse properties.')

    % Remove centres that lack synapse mask pixels below them
    goodIdx = find(data.synapseMask(data.synapseCenter));
    data.synapseCenter = data.synapseCenter(goodIdx);

    % Determine which synapse center each synapse pixel belongs to
    CC = bwconncomp(data.synapseMask);

    % Clear the old synapse pixels... 
    data.synapsePixels = {};

    for i = 1:length(CC.PixelIdxList)
      % idx is the synapse center number, centreIdx is the location of
      % the centre in the synapseMask.
      idx = find(ismember(data.synapseCenter,CC.PixelIdxList{i}));

      % Find the pixels belonging to each synapse centre
      % First calculate distance to all centres for every pixel
      % in synapse

      centerIdx = data.synapseCenter(idx);
      tmpDist = zeros(length(CC.PixelIdxList{i}),length(centerIdx));

      for j = 1:length(centerIdx)
        [yPix,xPix] = ind2sub(size(data.synapseMask),CC.PixelIdxList{i});
        [yCent,xCent] = ind2sub(size(data.synapseMask),centerIdx(j));

        tmpDist(:,j) = sqrt((xPix-xCent).^2+(yPix-yCent).^2);
      end

      % The first column of sortIdx contains the indexes of the closest
      % synapses.
      [sortedDist, sortIdx] = sort(tmpDist,2);

      for j = 1:length(centerIdx)
        data.synapsePixels{idx(j)} = ...
          CC.PixelIdxList{i}(find(sortIdx(:,1) == j));

      end

    end

    % Debug plot to verify that the synapse picture was correct
    if(0)
      tmpMaskR = zeros(size(data.neuriteMask));
      tmpMaskG = zeros(size(data.neuriteMask));
      tmpMaskB = zeros(size(data.neuriteMask));

      for i = 1:length(data.synapsePixels)
        tmpMaskR(data.synapsePixels{i}) = rand(1);
        tmpMaskG(data.synapsePixels{i}) = rand(1);
        tmpMaskB(data.synapsePixels{i}) = rand(1);
      end

      % Mark centres
      tmpMaskR(data.synapseCenter) = 1;
      tmpMaskG(data.synapseCenter) = 1;
      tmpMaskB(data.synapseCenter) = 1;

      tmpMask = zeros(size(data.image));
      tmpMask(:,:,1) = tmpMaskR;
      tmpMask(:,:,2) = tmpMaskG;
      tmpMask(:,:,3) = tmpMaskB;

      f = gcf;
      figure, imagesc(tmpMask)
      figure(f);
    end

    % Store the distances of the synapses to the soma also
    data.synapseDist = calculateSynapseDist();

    % Calculate the average intensity in each synapse for each channel
    data.synapseIntensityMorphMean = NaN*zeros(length(data.synapseCenter),1);
    data.synapseIntensityMorphSEM  = NaN*zeros(length(data.synapseCenter),1);
    data.synapseIntensitySynMean   = NaN*zeros(length(data.synapseCenter),1);
    data.synapseIntensitySynSEM    = NaN*zeros(length(data.synapseCenter),1);
    data.synapseIntensityXMean     = NaN*zeros(length(data.synapseCenter),1);
    data.synapseIntensityXSEM      = NaN*zeros(length(data.synapseCenter),1);
    data.synapseArea = NaN*zeros(size(data.synapseCenter));

    tmpMorph = double(data.image(:,:,detection.morphChannel));
    tmpSyn = double(data.image(:,:,detection.synChannel));
    tmpX = double(data.image(:,:,detection.XChannel));

    for i = 1:length(data.synapsePixels)

      nPixels = length(data.synapsePixels{i});

      data.synapseIntensityMorphMean(i) ...
        = mean(tmpMorph(data.synapsePixels{i}));

      data.synapseIntensityMorphSEM(i) ...
	= std(tmpMorph(data.synapsePixels{i})) / sqrt(nPixels);

      data.synapseIntensitySynMean(i) ...
	= mean(tmpSyn(data.synapsePixels{i}));

      data.synapseIntensitySynSEM(i) ...
        = std(tmpSyn(data.synapsePixels{i})) / sqrt(nPixels);

      data.synapseIntensityXMean(i) ...
        = mean(tmpX(data.synapsePixels{i}));

      data.synapseIntensityXSEM(i) ...
        = std(tmpX(data.synapsePixels{i})) / sqrt(nPixels);

      data.synapseArea(i) ...
        = length(data.synapsePixels{i})*(data.xyRes^2)*1e12; % Micrometers^2

    end

    % Calculate how many of synapse pixels are at max value
    sIdx = find(data.synapseMask);
    scIdx = data.synapseCenter;

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

    if(dispInfo.verbose)
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
      figure(handles.fig)
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setFigureTitle()

    set(handles.fig,'name', sprintf('SynD - Synapse Detection - %s', ...
				    data.fileName{dispInfo.curImg}));

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function nextImage(source, event)
    dispInfo.curImg = min(dispInfo.curImg + 1,data.num);
    setFigureTitle();				    
    showImage();
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function prevImage(source, event)
    dispInfo.curImg = max(dispInfo.curImg - 1,1);
    setFigureTitle();				    
    showImage();
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function showLogo()

    set(handles.fig,'CurrentAxes',handles.image);
    try
      if(isempty(bilder.logo))
        bilder.logo = imread('bilder/SynD-logo.jpg');
      end

      set(handles.image,'units','pixels')
      pos = get(handles.image,'position');
      set(handles.image,'units','normalized')

      width = pos(3) - pos(1);
      height = pos(4) - pos(2);

      imagesc(imresize(bilder.logo, [height width]));
    catch e
      getReport(e)
      disp('Unable to load logo.')
    end

    axis off

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function loadNeuron(source, event)

    curPwd = pwd;
    try
      cd(data.loadPath);
    catch
      fprintf('Unable to change to %s\n', data.loadPath)
    end

    fileTypes = { '*.lsm', 'LSM-image'; ...
		  '*.tif;*.tiff', 'Tiff-image' };
    fileTypes = fileTypes(editInfo.fileTypeOrder,:);

    [imageFile, imagePath, filterIndex] = ...
      uigetfile(fileTypes, ...
		'Select neuron image', ...
		'MultiSelect', 'off');

    cd(curPwd);

    if(~iscell(imageFile) & imageFile == 0)
      % User pressed cancel
      return;
    end

    data.loadPath = imagePath;

    dispInfo.stage = 1;

    % Clearing all the old data
    clearData();
   
    switch(editInfo.fileTypeOrder(filterIndex))
      case 1
	loadLSMfile(imageFile,imagePath);
      case 2
        loadTiffFiles(imageFile,imagePath);
    end

    data.neuriteMask = zeros(size(data.image,1),size(data.image,2));
    data.somaMask    = zeros(size(data.image,1),size(data.image,2));
    data.synapseMask = zeros(size(data.image,1),size(data.image,2));

    % Make sure the used filter is default next time
    switch(editInfo.fileTypeOrder(filterIndex))
      case 1
	editInfo.fileTypeOrder = [1 2];
      case 2
	editInfo.fileTypeOrder = [2 1]; 
    end
    
    % Update the title of the window
    setFigureTitle();				    

    activateStage('soma');
    showGUI('load');
    showImage();
    clearUndo();

    axis on

    % Allow the user to zoom using the mouse wheel
    startZoomWheelMode();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Load tiff can handle multi-select

  function loadTiffFiles(tiffFile,tiffPath)

    if(~iscell(tiffFile))
      tiffFile = {tiffFile};
    end

    for i = 1:length(tiffFile)
      fprintf('Loading %s\n', tiffFile{i})

      fileName = [tiffPath tiffFile{i}];
      tiffInfo = imfinfo(fileName);

      data.height = tiffInfo(1).Height;
      data.width = tiffInfo(1).Width;

      if(data.num == 0)
        data.image = zeros(data.height,data.width,3,1);
      end

      for j = 1:length(tiffInfo)
	data.num = data.num+1;
        tmp = imread(fileName,j);
        switch(size(tmp,3))

          case 1
	    oldTmp = tmp;
            tmp = zeros(data.height,data.width,3);
	    tmp(:,:,1) = oldTmp;

          case 2
	    oldTmp = tmp;
            tmp = zeros(data.height,data.width,3);
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

        data.image(:,:,:,data.num) = tmp;
        data.fileName{end+1} = tiffFile{i};
      end

      % Should we collapse image stack
      if(singleChannelStack())
	collapseImageStack();
      end

      calcMaxIntensity();

      % Reset scaling
      dispInfo.scaleRed = 1;
      dispInfo.scaleBlue = 1;
      dispInfo.scaleGreen = 1;

    end

    dispInfo.axis = [];

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function calcMaxIntensity()

    img = getRemappedImage();

    % Recalculate the new R,G,B max intensities, used for scaling
    tmpR = img(:,:,1);
    data.maxRed = double(max([tmpR(:); 1]));

    tmpG = img(:,:,2);
    data.maxGreen = double(max([tmpG(:); 1]));

    tmpB = img(:,:,3);
    data.maxBlue = double(max([tmpB(:); 1]));

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Load LSM can not handle multi-select

  function loadLSMfile(lsmFile,lsmPath)

    fprintf('Loading %s\n', lsmFile)

    lsmFileName = strcat(lsmPath,lsmFile);
    [lsmInfo,scanInfo,imInfo] = lsminfo(lsmFileName);
    imgStack = tiffread(lsmFileName);

    data.height = imgStack(1).height;
    data.width = imgStack(1).width;
    data.image = zeros(data.height,data.width,3,1);

    for j = 1:length(imgStack)

      if(iscell(imgStack(j).data))
        for i = 1:length(imgStack(j).data)
          data.image(:,:,i,j) = imgStack(j).data{i};
        end
      else

        % Only one channel present, put it in red
        data.image(:,:,1,j) = imgStack(j).data;      
      end

      if(j > 1)
        % Make sure voxel size etc matches
        if(imgStack(1).lsm.VoxelSizeX ~= imgStack(j).lsm.VoxelSizeX)
          disp('Voxel size inconsistent in the image! You are screwed.')
        end

      end

      data.fileName{j} = lsmFile;

    end

    data.num = length(imgStack);


    data.xyRes = imgStack(1).lsm.VoxelSizeX;

    if(imgStack(1).lsm.VoxelSizeX ~= imgStack(1).lsm.VoxelSizeY)
      % If you ever see this warning, let me know.
      disp('Warning: X and Y resolution differ, code assumes they are same.')
      disp(sprintf('Using %d m', data.xyRes))
    end

    calcMaxIntensity();

    % Reset scaling
    dispInfo.scaleRed = 1;
    dispInfo.scaleBlue = 1;
    dispInfo.scaleGreen = 1;

    dispInfo.axis = [];

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [selectedFiles,MATpath,LSMpath] = selectBatchReanalyseFiles()

    curPwd = pwd;
    try
      cd(data.loadPath);
    catch
      fprintf('Unable to change to %s\n', data.loadPath)
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

  function [selectedFiles,fileType] = selectBatchFiles()

    curPwd = pwd;
    try
      cd(data.loadPath);
    catch
      fprintf('Unable to change to %s\n', data.loadPath)
    end

    fileTypes = { '*.lsm', 'LSM-image'; ...
		  '*.tif*', 'Tiff-image' };
    fileTypes = fileTypes(editInfo.fileTypeOrder,:);

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

    data.loadPath = imagePath;
    fileType = editInfo.fileTypeOrder(filterIndex);

    % Make sure the used filter is default next time
    switch(fileType)
      case 1
	editInfo.fileTypeOrder = [1 2];
      case 2
	editInfo.fileTypeOrder = [2 1]; 
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function guessBatchExportFiles()

    data.exportXMLfile = strcat(data.fileName{dispInfo.curImg}, ...
				'-export.xml');

    data.exportSaveFile = strrep(data.exportXMLfile,'-export.xml','-save.mat');

    data.exportNeuriteMaskFile = ...
      strrep(data.exportXMLfile,'-export.xml','-neuritemask.tiff');
  
    data.exportSynapseMaskFile = ...
      strrep(data.exportXMLfile,'-export.xml','-synapsemask.tiff');

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function sccFlag = singleChannelStack()

    sccFlag = 0;

    if(nnz(data.image(:,:,2,1)) == 0 ...
       & nnz(data.image(:,:,3,1)) == 0 ...
       & size(data.image,4) > 1)

      % This is a single channel image stack, collapse!
      sccFlag = 1;
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function collapseImageStack(source, event)

    if(data.num == 0)
      disp('No images loaded, unable to callapse image stack.')
      return
    end

    disp('Collapsing image stack.')
    disp(['Creating a new image with three channels out of the' ...
          ' first channel of three images'])

    newImage = zeros(data.height,data.width,3,1);

    for i = 1:min(3,data.num)
      newImage(:,:,i,1) = data.image(:,:,1,i);
    end

    data.image = newImage;
    data.num = 1;
    dispInfo.curImg = 1;

    showImage();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function maxProjection(source, event)

    newImage = zeros(data.height,data.width,3,1);

    for i = 1:size(data.image,4)
      newImage = newImage + data.image(:,:,:,i);
    end

    data.image = newImage;
    data.num = 1;
    dispInfo.curImg = 1;

    calcMaxIntensity();
    showImage();

    stopEdit();
    activateStage('soma');
    showGUI('load');


  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function batchAnalyse(source, event)

    %%% Select neurons to analyse
    [selectedFiles,fileType] = selectBatchFiles();

    %%% Load the neurons and analyse them one at a time.
    for iFiles = 1:length(selectedFiles)

      % Clean out the old
      clearData();

      fprintf('Analysing %s\n', selectedFiles{iFiles})

      % Load file
      switch(fileType)
        case 1
          % LSM image
          loadLSMfile(selectedFiles{iFiles},data.loadPath);    
        case 2
          % TIFF image
          loadTiffFiles(selectedFiles{iFiles},data.loadPath);
      end

      % Analyse
      showGUI('soma');

      detectSoma();
      showMasks();
      drawnow();

      showGUI('neurites');
      detectNeurite();
      cleanMask();

      showMasks();
      drawnow();

      addThinNeurites();
      cleanMask();

      showMasks();
      drawnow();

      showGUI('synapses'); 
      detectSynapses();

      showMasks();
      drawnow();

      showGUI('analyse'); 
      guessBatchExportFiles();
      exportData();

    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function fileType = getFileType(filename)

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

  function clearData()

    disp('Clearing data.')

    data.num = 0;
    data.image = [];
    data.fileName = {};

    data.somaMask = [];
    data.somaMeasureMask = [];

    data.neuriteMask = [];
    data.skeleton = [];
    data.includeMask = [];

    data.synapseMask = [];
    data.synapseCenter = [];
    data.synapseArea = [];
    data.synapseDist = [];
    data.synapsePixels = {};
    data.neuriteLength = NaN;

    data.synapseIntensityMorphMean = []; % morph-channel
    data.synapseIntensityMorphSEM = [];  % morph-channel
    data.synapseIntensitySynMean = [];   % syn-channel
    data.synapseIntensitySynSEM = [];    % syn-channel
    data.synapseIntensityXMean = [];     % X-channel
    data.synapseIntensityXSEM = [];      % X-channel

    data.meanSynapse = [];

    data.meanSynapseMorph = [];
    data.meanSynapseSyn = [];
    data.meanSynapseX = [];    % Profile of X channel for synapse

    data.meanSynapseProfileDist = [];
    data.meanSynapseProfileMorph = [];
    data.meanSynapseProfileSyn = [];
    data.meanSynapseProfileX = [];

    data.skeleton = [];
    data.shollEdges = [NaN NaN NaN];
    data.shollDendHist = [];
    data.shollSynHist = [];

    data.shollIntMorphMean = [];
    data.shollIntMorphStd = [];
    data.shollIntMorphSEM = [];

    data.shollIntSynMean = [];
    data.shollIntSynStd = [];
    data.shollIntSynSEM = [];

    data.shollIntXMean = [];
    data.shollIntXStd = [];
    data.shollIntXSEM = [];

    data.somaMeasureMorph = NaN;
    data.somaMeasureSyn = NaN;
    data.somaMeasureX = NaN;

    data.somaArea = [];
    data.somaMajorAxisLength = [];
    data.somaMinorAxisLength = [];

    data.measurePoint = [];

    data.intensityHistogramEdges = [];
    data.morphHist = [];
    data.synHist = [];
    data.XHist = [];

    % Intensity gradient along neurites
    data.gradientEdges = [];

    data.morphGradientMean = [];
    data.morphGradientStd = [];
    data.morphGradientSEM = [];

    data.synGradientMean = [];
    data.synGradientStd = [];
    data.synGradientSEM = [];

    data.XGradientMean = [];
    data.XGradientStd = [];
    data.XGradientSEM = [];

    data.distMask = [];

    data.dirVect = [];
    data.rigidity = [];
    data.lambdaMax = 0;

    % These are absolute thresholds used for detection
    data.somaIntensityThreshold = NaN;
    data.synapseIntensityThreshold = NaN;

    % Clear to prevent overwrites...
    data.exportXMLfile = [];
    data.exportSaveFile = [];
    data.exportNeuriteMaskFile = [];
    data.exportSynapseMaskFile = [];

    dispInfo.curImg = 1;

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function batchReAnalyse(source, event)

    detection.reanalyseFlag = 1;

    doSynapseDetection = NaN;

    [selectedFiles, MATpath,LSMpath] = selectBatchReanalyseFiles();

    oldExportPath = data.exportPath;
    data.exportPath = MATpath;

    %%% Load the neurons and analyse them one at a time.
    for iFiles = 1:length(selectedFiles)

      % Clean out the old
      clearData();

      try
        old = load(strcat(MATpath,selectedFiles{iFiles}));
        old = old.old;

        try
          detection.morphChannel = old.morphChannel;
          detection.synChannel = old.synChannel;
          detection.XChannel = old.XChannel;
          detection.remapChan = old.remapChan;
        catch
          disp('Old channel info. Assuming not remapped')
          old.remapChan = [detection.morphChannel, ...
			   detection.synChannel, ...
			   detection.XChannel];
        end

        try
          detection.excludeSomaSynapses = old.excludeSomaSynapses;
        catch
          disp('Excluding soma synapses.')
	  detection.excludeSomaSynapses = true;
        end

        % Restore old XY-res (in case of a tiff-file)
        data.xyRes = old.xyRes;

        fileType = getFileType(old.fileName{1}); 
        switch(fileType)
          case 1
	    fprintf('Loading %s%s\n',LSMpath,old.fileName{1})
            loadLSMfile(old.fileName{1},LSMpath);    

          case 2 
            loadTiffFiles(old.fileName{1},LSMpath);

          otherwise
	    fprintf('Unknown file type: %s', old.fileName{1})
            continue
        end

        fprintf('Re-analysing %s\n', selectedFiles{iFiles})

        try
          data.neuriteMask = old.neuriteMask;
          data.somaMask = old.somaMask;
        catch
          % Old file names in old versions of program...
          data.neuriteMask = old.neuriteMaskIdx;
          data.somaMask = old.somaMaskIdx;
        end

        % Load the old intensity threshold used
        try
          data.somaIntensityThreshold = old.somaIntensityThreshold;
        catch
          data.somaIntensityThreshold = old.detection.morphThreshold;
        end

        exportInfo.saveSholl = 1;
        exportInfo.saveSynapseProfile = 1;
        exportInfo.saveMat = 1;
        exportInfo.saveIntensityHistogram = 1;
        exportInfo.saveGradient = 1;

        try
          data.somaMeasureMask = old.somaMeasureMask;
          exportInfo.saveSomaMeasure = 1;
        catch
          disp('No soma measure marked, ignoring.')
          data.somaMeasureMask = [];
        end

        writeExportSettings();

        showGUI('soma');
        showMasks();
        drawnow();

        calculateSomaMeasures();
        dispInfo.measurePointColor = dispInfo.defaultMeasurePointColor;
        showMasks();
        drawnow();

        showGUI('neurites');
        updateNeuriteMask();
        data.distMask = makeDistMask(data.somaMask,data.neuriteMask);

        showMasks();
        drawnow();

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
          showGUI('synapses'); 
          detectSynapses();

          showMasks();
          drawnow();

        end

        activateStage('analyse')
        showGUI('analyse'); 
        guessBatchExportFiles();

        % Modifying the save-file so we do not overwrite the old save file
        data.exportSaveFile = ...
	  strrep(data.exportSaveFile,'-save.mat','-re-save.mat');

        exportData();


      catch e
        getReport(e)
        fprintf('Failed to re-analyse %s\n', selectedFiles{iFiles})

      end

    end

    % Restore the old export path
    data.exportPath = oldExportPath;

    % Exit re-analyse mode.
    detection.reanalyseFlag = 0;

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function mask = largestComponent(mask)

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

  function detectSoma(source, event)

    % Should we auto-detect the threshold?
    if(detection.autoSomaThreshold)
      autoSomaThreshold();
    end

    % Save the actual threshold used
    data.somaIntensityThreshold = detection.morphThreshold;

    % Read in latest inputs from GUI (so user doesnt have to hit enter)
    setMorphThresh();
    setSomaErodeRadius();

    saveUndo('detect soma');

    % This gives soma and neurites
    detectNeuritesSimple(); 

    % Next we erode away the neurites, a bit further down...

    if(detection.singleSoma)
      % If only one soma, use largest component.
      mask = largestComponent(data.neuriteMask);
    else
      mask = data.neuriteMask;
    end

    if(~nnz(data.neuriteMask))
      uiwait(errordlg(['No pixels in channel selected for morphology. ' ...
                       'Check that the morphology channel is selected, ' ...
                       'also try lowering the threshold and soma erode radius.'], ...
		      'Bad channel selected', 'modal'))
      return
    end
    
    seSoma = strel('disk',detection.somaErodeRadius);

    data.somaMask = imopen(mask,seSoma);

    if(nnz(data.somaMask) == 0)
      % Display the neurite mask for the ellusive soma
      dispInfo.neuriteColor = [1 1 1]*0.7;
      showImage();

      uiwait(errordlg(['No soma detected, reduce soma erode radius and ' ...
		       'verify you have the right channel for morphology'], ...
		      'Overly aggressive erosion', 'modal'))

      % Set focus to erode radius
      uicontrol(handles.somaErodeRadius);
    else
      dispInfo.neuriteColor = NaN;
      data.neuriteMask = zeros(size(data.neuriteMask));
      showImage();

      activateStage('neurites');

      CC = bwconncomp(data.somaMask);
      text(20,20,sprintf('Found %d soma', CC.NumObjects), ...
	   'color','white')

    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function showHist()

    nHist = 30;

    img = getRemappedImage();

    % Create histogram over red channel
    % The part of the histogram that is saturated is marked black
    set(handles.fig,'CurrentAxes',handles.redHist)    

    tmp = img(:,:,1);
    redEdges = linspace(0,max(1,data.maxRed),nHist);
    nRed = histc(tmp(:),redEdges);
    dCenter = diff(redEdges(1:2))/2;
    idx = find(redEdges*dispInfo.scaleRed <= data.maxRed);

    bar(redEdges+dCenter,nRed, ...
	'facecolor',[0 0 0], 'edgecolor',[0 0 0])
    hold on
    bar(redEdges(idx)+dCenter,nRed(idx), ...
	'facecolor',[1 0 0], 'edgecolor',[1 0 0])
    hold off
    set(gca,'YScale','log')
    axis off
    a = axis; a(2) = max([data.maxRed data.maxGreen data.maxBlue 1]); axis(a);

    set(handles.fig,'CurrentAxes',handles.greenHist)    
    tmp = img(:,:,2);
    greenEdges = linspace(0,max(1,data.maxGreen),nHist);
    nGreen = histc(tmp(:),greenEdges);
    dCenter = diff(greenEdges(1:2))/2;
    idx = find(greenEdges*dispInfo.scaleGreen <= data.maxGreen);

    bar(greenEdges+dCenter,nGreen, ...
	'facecolor',[0 0 0], 'edgecolor',[0 0 0])
    hold on
    bar(greenEdges(idx)+dCenter,nGreen(idx), ...
	'facecolor',[0 1 0], 'edgecolor',[0 1 0])
    hold off
    set(gca,'YScale','log')
    axis off
    a = axis; a(2) = max([data.maxRed data.maxGreen data.maxBlue]); axis(a);

    set(handles.fig,'CurrentAxes',handles.blueHist)    
    tmp = img(:,:,3);
    blueEdges = linspace(0,max(1,data.maxBlue),nHist);
    nBlue = histc(tmp(:),blueEdges);
    dCenter = diff(blueEdges(1:2))/2;
    idx = find(blueEdges*dispInfo.scaleBlue <= data.maxBlue);

    bar(blueEdges+dCenter,nBlue, ...
	'facecolor',[0 0 0], 'edgecolor',[0 0 0])
    hold on
    bar(blueEdges(idx)+dCenter,nBlue(idx), ...
	'facecolor',[0 0 1], 'edgecolor',[0 0 1])
    hold off
    set(gca,'YScale','log')
    axis off
    a = axis; a(2) = max([data.maxRed data.maxGreen data.maxBlue]); axis(a);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function showSkeleton(source, event)
 
    % Recalculate the skeleton
    data.skeleton = bwmorph(data.neuriteMask-data.somaMask>0,'skel',inf);
    trimSkeleton();

    dispInfo.showSkeleton = 1;

    showMasks();

    % dispInfo.showSkeleton = 0;

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % altFig is normaly not specified, but can be used to draw 
  % in alternative figure, for example when making images

  function showMasks(altFig)

    if(exist('altFig'))
      % If we want to draw the masks in an alternative figure
      figure(altFig)
    else 
      set(handles.fig,'CurrentAxes',handles.image)    
    end
    tmp = getImg();

    tmpR = tmp(:,:,1);
    tmpG = tmp(:,:,2);
    tmpB = tmp(:,:,3);

    % Should we show neurites?
    if(~isnan(dispInfo.neuriteColor) & dispInfo.showMask)
      tmpR(find(data.neuriteMask == 1)) = dispInfo.neuriteColor(1);
      tmpG(find(data.neuriteMask == 1)) = dispInfo.neuriteColor(2);
      tmpB(find(data.neuriteMask == 1)) = dispInfo.neuriteColor(3);

      tmpR(find(data.neuriteMask == 2)) = 0.4*dispInfo.neuriteColor(1);
      tmpG(find(data.neuriteMask == 2)) = 0.4*dispInfo.neuriteColor(2);
      tmpB(find(data.neuriteMask == 2)) = 0.4*dispInfo.neuriteColor(3);

      if(dispInfo.showSkeleton & ~isempty(data.skeleton))
        tmpR(find(data.skeleton)) = dispInfo.neuriteColor(1);
        tmpG(find(data.skeleton)) = 0.5*dispInfo.neuriteColor(1);
        tmpB(find(data.skeleton)) = 0.5*dispInfo.neuriteColor(1);
      end

    end

    % Should we display soma?
    if(~isnan(dispInfo.somaColor) & dispInfo.showMask)
      tmpR(find(data.somaMask)) = dispInfo.somaColor(1);
      tmpG(find(data.somaMask)) = dispInfo.somaColor(2);
      tmpB(find(data.somaMask)) = dispInfo.somaColor(3);
    end 

    % Should we mark synapses
    if(~isnan(dispInfo.synapseColor) & dispInfo.showMask)
      % Mark synapses
      tmpR(find(data.synapseMask)) = 0.6*dispInfo.synapseColor(1);
      tmpG(find(data.synapseMask)) = 0.6*dispInfo.synapseColor(2);
      tmpB(find(data.synapseMask)) = 0.6*dispInfo.synapseColor(3);

      % Mark synapse centres
      tmpR(data.synapseCenter) = dispInfo.synapseColor(1);
      tmpG(data.synapseCenter) = dispInfo.synapseColor(2);
      tmpB(data.synapseCenter) = dispInfo.synapseColor(3);
    end

    if(~isnan(dispInfo.measurePointColor) & dispInfo.showMask)
      % Mark soma measure points
      tmpR(find(data.somaMeasureMask)) = dispInfo.measurePointColor(1);
      tmpG(find(data.somaMeasureMask)) = dispInfo.measurePointColor(2);
      tmpB(find(data.somaMeasureMask)) = dispInfo.measurePointColor(3);
    end

    tmp(:,:,1) = tmpR;
    tmp(:,:,2) = tmpG;
    tmp(:,:,3) = tmpB;

    imagesc(tmp);
    axis equal

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Returns remapped but un-normalized image

  function img = getRemappedImage()

     if(isempty(data.image))
       img = zeros(0,0,3);
       return
     end

     origImg = squeeze(data.image(:,:,:,dispInfo.curImg));

     % Red here refers to the red in the original image, ie channel 1
     % Green refers to channel 2, and Blue to channel 3

     % However, to confuse things dispInfo.showRed refers
     % to the red colour showed to the user.

     img = zeros(data.height,data.width,3);

     img(:,:,detection.remapChan(1)) = origImg(:,:,detection.morphChannel);
     img(:,:,detection.remapChan(2)) = origImg(:,:,detection.synChannel) ...
		                     + img(:,:,detection.remapChan(2));
     img(:,:,detection.remapChan(3)) = origImg(:,:,detection.XChannel) ...
		                     + img(:,:,detection.remapChan(3));

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Returns the properly scaled image

  function img = getImg()

    img = getRemappedImage();

    if(data.maxRed > 0 & dispInfo.showRed)
      img(:,:,1) = min(img(:,:,1) / double(data.maxRed) * dispInfo.scaleRed,1);
    else
      img(:,:,1) = 0;
    end

    if(data.maxGreen > 0 & dispInfo.showGreen)
      img(:,:,2) = min(img(:,:,2) / double(data.maxGreen) * dispInfo.scaleGreen,1);
    else
      img(:,:,2) = 0;
    end

    if(data.maxBlue > 0 & dispInfo.showBlue)
      img(:,:,3) = min(img(:,:,3) / double(data.maxBlue) * dispInfo.scaleBlue,1);
    else
      img(:,:,3) = 0;
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function showImage()

    if(isempty(data.image))
      % Nothing loaded yet, so nothing to display
      return;
    end

    set(handles.fig,'CurrentAxes',handles.image)    
    oldAxis = axis;

    switch(dispInfo.state)
      case 'load'

        % Just display the raw image...
        tmp = getImg();
        imagesc(tmp);
        axis equal

        set(handles.num,'String', ...
	    sprintf('Image %d (%d)', dispInfo.curImg, data.num))
  
        showHist();

        set(handles.fig,'CurrentAxes',handles.image)    

      case 'soma'

        showMasks();

        if(data.measurePoint)   
          hold on
          [yP,xP] = ind2sub(size(data.somaMask),data.measurePoint);
          contour(data.somaMask,1,'color',[0.3 0.3 0.3], ...
		  'linewidth',2)
          p = plot(xP,yP,'ok');
          set(p,'UIContextMenu',handles.measurePointMenu);

          hold off
        end

      case 'neurites'

        showMasks();

      case 'synapses'

        showMasks();
        hold on
        contour(data.neuriteMask,1,'color',[0.3 0.3 0.3])
        hold off

      case 'analyse'

        % Just display the raw image...
        set(handles.fig,'CurrentAxes',handles.image)    
        tmp = getImg();
        imagesc(tmp);
        axis equal

    end

    if(isempty(dispInfo.axis))
      dispInfo.axis = axis;
    else
      dispInfo.axis = oldAxis;
      set(handles.fig,'CurrentAxes',handles.image)    
      axis(oldAxis);
    end


    drawnow

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function saveUndo(description)
  
    if(isempty(data.somaMask))
      % There is no soma detected, no point adding undo
      return;
    end

    clear tmp
    tmp.somaMask = data.somaMask;
    tmp.neuriteMask = data.neuriteMask;
    tmp.synapseMask = data.synapseMask;
    tmp.includeMask = data.includeMask;
    tmp.synapseCenter = data.synapseCenter;
    tmp.measurePoint = data.measurePoint;

    tmp.description = description;
    tmp.state = dispInfo.state;
    tmp.stage = dispInfo.stage;


    editInfo.undo(end+1) = tmp;

    set(handles.menuItemUndo, 'Label', ...
	sprintf('<html><u>U</u>ndo %s</html>', editInfo.undo(end).description))

    if(length(editInfo.undo) > editInfo.maxUndo)
      % disp('Max undo length reached, removing old history')
      editInfo.undo(1) = [];
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function clearUndo()

    editInfo.undo = struct('somaMask', [], ...
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

  function undoLastAction(source, event)

    if(length(editInfo.undo) > 0)

      data.somaMask = editInfo.undo(end).somaMask;
      data.neuriteMask = editInfo.undo(end).neuriteMask;
      data.synapseMask = editInfo.undo(end).synapseMask;
      data.includeMask = editInfo.undo(end).includeMask;
      data.synapseCenter = editInfo.undo(end).synapseCenter;
      data.measurePoint = editInfo.undo(end).measurePoint;
      makeMeasureMask();

      dispInfo.state = editInfo.undo(end).state;
      dispInfo.stage = editInfo.undo(end).stage;

      editInfo.undo(end) = [];

      showGUI(dispInfo.state);

    end

    if(length(editInfo.undo) > 0)
      set(handles.menuItemUndo, 'Label', ...
	  sprintf('<html><u>U</u>ndo %s</html>', editInfo.undo(end).description))
    else
      set(handles.menuItemUndo,'Label','<html><u>U</u>ndo: No history</html>')
    end

    showImage();
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function excludeRegion(source, event)

    saveUndo('add exclude region');

    stopEdit();

    excludeMask = roipoly();

    data.includeMask = data.includeMask & ~excludeMask;

    set(handles.fig,'CurrentAxes',handles.image)    
    imagesc(~data.includeMask);
    pause(2);
    detectSoma();
    dispInfo.axis = []; % Reset view
    showImage();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function clearExcludeRegions(source, event)

    saveUndo('reset exclude region');
    data.includeMask = ones(size(data.neuriteMask));

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function steerableFilter()

    % Parts of this function provided by Matthew Down
    % as well as GaussianDerivatives3D and nms.

    img = data.image(:,:,detection.morphChannel,dispInfo.curImg);

    [g0x g0y g0z g1x g1y g1z g2x g2y g2z] = ...
      GaussianDerivatives3D(detection.filterSize/data.xyRes);

    N = 4/(detection.filterSize/data.xyRes)^4;
    Rxx = -N*conv2(conv2(img, g2x, 'same'), g0y, 'same');
    Rxy = -N*conv2(conv2(img, g1x, 'same'), g1y, 'same');
    Ryy = -N*conv2(conv2(img, g0x, 'same'), g2y, 'same');

    Rg = N*conv2(conv2(img, g0x, 'same'), g0y, 'same');

    % We use inverted 2nd derivative of gaussian,
    % Theory behind Meijering 2003 uses non-inverted, hence minus
    % since eigenProperties based on Meijering

    if(detection.steerableFilterAlpha)
      alpha = detection.steerableFilterAlpha;

      [eigenValues, eigenVectors, neuriteDir, ridgeMask] = ...
        eigenProperties(Rxx+alpha*Ryy,Ryy+alpha*Rxx,Rxy*(1-alpha)); 
    else
      % Default case, no elongated filters
      [eigenValues, eigenVectors, neuriteDir, ridgeMask] = ...
        eigenProperties(Rxx,Ryy,Rxy); 
    end

    data.dirVect = squeeze(eigenVectors(:,:,1,:));
    data.rigidity = squeeze(eigenValues(:,:,2));
    data.lambdaMax = max(data.rigidity(:));

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

      [xC,yC] = meshgrid(1:data.width,1:data.height);
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

  function [g0x g0y g0z g1x g1y g1z g2x g2y g2z] = GaussianDerivatives3D(sigma)

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
    eigenProperties(Rxx,Ryy,Rxy)

    eigenValues = zeros(size(Rxx,1),size(Rxx,2),2);
    eigenVectors = zeros(size(Rxx,1),size(Rxx,2),2,2);
    neuriteDir = zeros(size(Rxx));

    if(1)
      % Homemade function to calculate eigenvalues and eigenvectors
      [eigenValues, eigenVectors] = calcEigen2x2(Rxx,Rxy,Rxy,Ryy);
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

  function cc = connectCost(pointA,pointB)

    % Taken from Meijering et al 2003
    % Note that our filter is inverted.

    % We use > 0 instead of < 0, since our filter is minus that one
    % used in Meijering
    if(data.rigidity(pointB(2),pointB(1)) > 0)
      ccLambda = 1 - data.rigidity(pointB(2),pointB(1)) / data.lambdaMax;
    else
      ccLambda = 1;
    end

    AB = pointB-pointA;
    dirAB = AB / norm(AB);

    tmpDirA = [data.dirVect(pointA(2),pointA(1),1), ...
	        data.dirVect(pointA(2),pointA(1),2)];

    wA = tmpDirA / norm(tmpDirA);
 
    tmp = abs(sum(dirAB.*wA));
    ccV = 0.5*(sqrt(1-tmp) + sqrt(1+tmp));

    cc = detection.connectCostLambda * ccLambda ...
          + (1-detection.connectCostLambda)*ccV;

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function growNeuriteFromSoma()

    % Start by seeding neurite mask with soma mask
    data.neuriteMask = data.somaMask;

    saveUndo('detect neurites')
    growNeurite(data.somaMask);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function growNeurite(neuriteSeed)

    %profile clear
    %profile on

    set(handles.calculateLabel, 'String', 'Locating neurites...');
    drawnow

    tic

    somaIdx = find(neuriteSeed);

    detection.pixelQue = somaIdx;

    loopCtr = 0;

    while(length(detection.pixelQue) > 0)      
      
      % Optimized this line by doing calculation explicit, see below
      %[yCenter,xCenter] = ind2sub(size(data.neuriteMask), ...
      %				  detection.pixelQue(1));

      xCenter = floor((detection.pixelQue(1)-1) / size(data.neuriteMask,1))+1;
      yCenter = mod(detection.pixelQue(1)-1,size(data.neuriteMask,1))+1;

      % Find valid neighours
      [neighX, neighY] = processNeighbourhood(xCenter,yCenter);

      % Make sure they are not already qued

      % Optimized this line by calculating manually...
      % neighIdx = sub2ind(size(data.neuriteMask),neighY,neighX);
      neighIdx = neighY+(neighX-1)*size(data.neuriteMask,1);

      detection.pixelQue = [detection.pixelQue(2:end); neighIdx];

      % Mark them in the mask
      data.neuriteMask(neighIdx) = 1;      

      loopCtr = loopCtr + 1;

      if(mod(loopCtr,5000) == 0)
        % figure, imagesc(data.neuriteMask), axis equal, drawnow
        fprintf('%d iterations\n', loopCtr)
        set(handles.calculateLabel, ...
	    'String', sprintf('Locating neurites... (%d)', loopCtr));
        drawnow
      end
    end
 
    % Remove spurious pixels around the neurite
    se = strel('disk',1);
    data.neuriteMask = imopen(data.neuriteMask,se);

    fprintf('%d iterations, done.\n', loopCtr)

    updateNeuriteMask();

    data.distMask = makeDistMask(data.somaMask,data.neuriteMask);

    toc

    set(handles.calculateLabel, 'String', '');

    showImage();

    %profview

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [neighX, neighY] = processNeighbourhood(xCenter,yCenter)

    nHoodX = transpose([1 1  1 0  0 -1 -1 -1] + xCenter);
    nHoodY = transpose([1 0 -1 1 -1  1  0 -1] + yCenter);

    if(xCenter == 1 | xCenter == data.width ...
       | yCenter == 1 | yCenter == data.height)

      % Only if we are at border do this more time consuming check

      validIdx = find(1 <= nHoodX & nHoodX <= data.width ...
		      & 1 <= nHoodY & nHoodY <= data.height);

      nHoodX = nHoodX(validIdx);
      nHoodY = nHoodY(validIdx);
    end

    nHoodMask = zeros(size(nHoodX));

    for iN = 1:length(nHoodX)

      % My idea was to check if connectCost was lower than a certain value
      % then if that is the case, the neighbourhood pixels are added
      % and returned. The function that calls this function then adds these
      % pixels to detection.pixelQue (if they are not already part of
      % neurite mask.

     if(data.neuriteMask(nHoodY(iN),nHoodX(iN)) == 0 ...
         & connectCost([xCenter yCenter], [nHoodX(iN) nHoodY(iN)]) ...
  	< detection.maxAddCost)
	 

        % disp('Adding pixel')
        nHoodMask(iN) = 1;    
      % else
      %  c = connectCost([xCenter yCenter], [nHoodX(iN) nHoodY(iN)])
      %  n = data.neuriteMask(nHoodY(iN),nHoodX(iN))
      end
 
    end

    neighX = nHoodX(find(nHoodMask));
    neighY = nHoodY(find(nHoodMask));

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [eigenVal, eigenVec] = calcEigen2x2(a11,a12,a21,a22)

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

  function calculateNeuriteLength()

    skelIdx = find(data.skeleton);

    tmpSkel = data.skeleton;

    distTotal = 0;

    for i = 1:length(skelIdx)
      [yC,xC] = ind2sub(size(data.skeleton),skelIdx(i));

      % Clipping so we do not get outside image
      yNeigh = min(max(yC-1:yC+1,1),size(data.skeleton,1));
      xNeigh = min(max(xC-1:xC+1,1),size(data.skeleton,2));

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

    data.neuriteLength = distTotal * data.xyRes;

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function trimSkeleton()

    branchPoints = bwmorph(data.skeleton,'branchpoints');    
    endPoints = bwmorph(data.skeleton,'endpoints');
    endPointIdx = find(endPoints);

    % We want to remove all the mini-branches smaller than the trim size.
    % To find them we remove the branch points, then calculate the locate
    % all the connected components smaller than the trim size, and if
    % they only have one branch point around them (imdilate 1), 
    % then we remove it.

    trimMask = double(data.skeleton - bwmorph(branchPoints,'dilate') > 0);

    mCC = bwconncomp(trimMask);

    for i = 1:length(mCC.PixelIdxList)
      if(length(mCC.PixelIdxList{i}) < detection.trimNeuriteSize)
        % Does it contain an endpoint?
        if(nnz(ismember(mCC.PixelIdxList{i},endPointIdx)))
	  [y,x] = ind2sub(size(data.skeleton),mCC.PixelIdxList{i});
          yAll = max(1,min([y+1;y+1;y+1;y;y;y;y-1;y-1;y-1],data.height));
          xAll = max(1,min([x+1;x;x-1;x+1;x;x-1;x+1;x;x-1],data.width));
          trimIdx = sub2ind(size(data.skeleton),yAll,xAll);

 	  %fprintf('Trimming away %d pixels.\n', length(mCC.PixelIdxList{i}))

	  % Yes it does, remove pixel from skeleton
          %data.skeleton(mCC.PixelIdxList{i}) = 0;
          data.skeleton(trimIdx) = 0;
        end    
      end
    end   

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Calculates the sholl analysis for a skeleton of the neurites
  function shollAnalysis()

    data.skeleton = bwmorph(data.neuriteMask-data.somaMask>0,'skel',inf);

    % This give the centroid of all pixels included in all somas
    somaProps = regionprops(data.somaMask,'centroid');

    % Prune away tiny pertrusions, mostly artifacts
    trimSkeleton();

    CC = bwconncomp(data.somaMask);

    if(CC.NumObjects ~= 1)
      fprintf('Detected %d somas, aborting Sholl analysis.\n', ...
	      length(somaProps))
      exportInfo.saveSholl = 0;

      % Save histogram with NaN
      maxDist = sqrt(data.width^2+data.height^2)*data.xyRes;
      data.shollEdges = 0:detection.shollBinSize:maxDist;
      data.shollDendHist = NaN*zeros(1,length(data.shollEdges)-1);

      data.shollIntMorphMean = NaN*zeros(1,length(data.shollEdges)-1);
      data.shollIntMorphStd = NaN*zeros(1,length(data.shollEdges)-1);
      data.shollIntMorphSEM = NaN*zeros(1,length(data.shollEdges)-1);
      data.shollIntSynMean = NaN*zeros(1,length(data.shollEdges)-1);
      data.shollIntSynStd = NaN*zeros(1,length(data.shollEdges)-1);
      data.shollIntSynSEM = NaN*zeros(1,length(data.shollEdges)-1);
      data.shollIntXMean = NaN*zeros(1,length(data.shollEdges)-1);
      data.shollIntXStd = NaN*zeros(1,length(data.shollEdges)-1);
      data.shollIntXSEM = NaN*zeros(1,length(data.shollEdges)-1);

      writeExportSettings();

      data.skeleton = bwmorph(data.neuriteMask-data.somaMask>0,'skel',inf);
      trimSkeleton();

      calculateNeuriteLength();

      return
    else
      exportInfo.saveSholl = 1;
      writeExportSettings();
    end

    % This must be modified to handle multiple somas
    somaCenter = somaProps.Centroid;

    [y,x] = find(data.somaMask);
    somaRadius = max(sqrt((x-somaCenter(1)).^2 + (y-somaCenter(2)).^2)) ...
		     * data.xyRes;

    [y,x] = find(data.skeleton);
    distSkeleton = sqrt((x-somaCenter(1)).^2 + (y-somaCenter(2)).^2) ...
			     * data.xyRes;

    distToSoma = zeros(size(data.neuriteMask));
    for i = 1:length(x)
      % We want to the soma radius to be defined as distance 0
      distToSoma(y(i),x(i)) = distSkeleton(i) - somaRadius;
    end

    maxDist = sqrt(data.width^2+data.height^2)*data.xyRes;
    %data.shollEdges = (round(somaRadius/5e-6)*5e-6+5e-6)...
    %			 :detection.shollBinSize:maxDist;

    data.shollEdges = 0:detection.shollBinSize:maxDist;
    data.shollDendHist = zeros(1,length(data.shollEdges)-1);
    shollHalfWidth = sqrt(2)*data.xyRes/2;

    % Soma radius defined as distance 0
    [ySyn,xSyn] = ind2sub(size(data.synapseMask),data.synapseCenter);
    synDist = sqrt((xSyn-somaCenter(1)).^2+(ySyn-somaCenter(2)).^2) ...
                * data.xyRes - somaRadius;

    data.shollIntMorphMean = zeros(1,length(data.shollEdges)-1);
    data.shollIntMorphStd = zeros(1,length(data.shollEdges)-1);
    data.shollIntMorphSEM = zeros(1,length(data.shollEdges)-1);
    data.shollIntSynMean = zeros(1,length(data.shollEdges)-1);
    data.shollIntSynStd = zeros(1,length(data.shollEdges)-1);
    data.shollIntSynSEM = zeros(1,length(data.shollEdges)-1);
    data.shollIntXMean = zeros(1,length(data.shollEdges)-1);
    data.shollIntXStd = zeros(1,length(data.shollEdges)-1);
    data.shollIntXSEM = zeros(1,length(data.shollEdges)-1);


    for i = 1:length(data.shollEdges)-1
      % Neurite Sholl analysis
      BW = (data.shollEdges(i) <= distToSoma) ...
	    & (distToSoma < data.shollEdges(i+1));
      CC = bwconncomp(BW);
      data.shollDendHist(i) = length(CC.PixelIdxList);

      % Synapse 
      synIdx = find(data.shollEdges(i) <= synDist ...
		    & synDist < data.shollEdges(i+1));

      if(~isempty(synIdx))
        data.shollIntMorphMean(i) = mean(data.synapseIntensityMorphMean(synIdx));
        data.shollIntMorphStd(i) = std(data.synapseIntensityMorphMean(synIdx));
        data.shollIntMorphSEM(i) = data.shollIntMorphStd(i) /sqrt(length(synIdx));

        data.shollIntSynMean(i) = mean(data.synapseIntensitySynMean(synIdx));
        data.shollIntSynStd(i) = std(data.synapseIntensitySynMean(synIdx));
        data.shollIntSynSEM(i) = data.shollIntSynStd(i) /sqrt(length(synIdx));

        data.shollIntXMean(i) = mean(data.synapseIntensityXMean(synIdx));
        data.shollIntXStd(i) = std(data.synapseIntensityXMean(synIdx));
        data.shollIntXSEM(i) = data.shollIntXStd(i) /sqrt(length(synIdx));
      else
        % NaN fields will be left empty when exporting to XML
        data.shollIntMorphMean(i) = NaN;
        data.shollIntMorphStd(i) = NaN;
        data.shollIntMorphSEM(i) = NaN;

        data.shollIntSynMean(i) = NaN;
        data.shollIntSynStd(i) = NaN;
        data.shollIntSynSEM(i) = NaN;

        data.shollIntXMean(i) = NaN;
        data.shollIntXStd(i) = NaN;
        data.shollIntXSEM(i) = NaN;
       end

    end

    % data.synapseDist contains the arc length
    % synDist is the shortest distance from soma to synapse...

    data.shollSynHist = histc(synDist,data.shollEdges);
    data.shollSynHist = data.shollSynHist(1:end-1);

    calculateNeuriteLength();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function calculateIntensityGradients()

    if(isempty(data.distMask))
      data.distMask = makeDistMask(data.somaMask,data.neuriteMask);
    end

    tmp = data.distMask(:);
    tmp(find(tmp == inf)) = NaN;
    maxDist = max(tmp(:));

    data.gradientEdges = ...
      0:detection.shollBinSize:(maxDist+detection.shollBinSize);

    data.morphGradientMean = NaN*zeros(size(data.gradientEdges));
    data.morphGradientStd = NaN*zeros(size(data.gradientEdges));
    data.morphGradientSEM = NaN*zeros(size(data.gradientEdges));

    data.synGradientMean = NaN*zeros(size(data.gradientEdges));
    data.synGradientStd = NaN*zeros(size(data.gradientEdges));
    data.synGradientSEM = NaN*zeros(size(data.gradientEdges));

    data.XGradientMean = NaN*zeros(size(data.gradientEdges));
    data.XGradientStd = NaN*zeros(size(data.gradientEdges));
    data.XGradientSEM = NaN*zeros(size(data.gradientEdges));

    tmpMorph = data.image(:,:,detection.morphChannel, dispInfo.curImg);
    tmpSyn = data.image(:,:,detection.synChannel, dispInfo.curImg);
    tmpX = data.image(:,:,detection.XChannel, dispInfo.curImg);

    for i = 1:length(data.gradientEdges)
      if(i == 1)
        idx = find(data.distMask == 0);
      else
        idx = find(data.gradientEdges(i-1) < data.distMask ...
		   & data.distMask <= data.gradientEdges(i));
      end

      data.morphGradientMean(i) = mean(tmpMorph(idx));
      data.morphGradientStd(i) = std(tmpMorph(idx));
      data.morphGradientSEM(i) = mean(tmpMorph(idx))/sqrt(length(idx));

      data.synGradientMean(i) = mean(tmpSyn(idx));
      data.synGradientStd(i) = std(tmpSyn(idx));
      data.synGradientSEM(i) = std(tmpSyn(idx))/sqrt(length(idx));

      data.XGradientMean(i) = mean(tmpX(idx));
      data.XGradientStd(i) = std(tmpX(idx));
      data.XGradientSEM(i) = std(tmpX(idx))/sqrt(length(idx));

    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function intensityHistograms()

    maxInt = max([data.maxRed,data.maxGreen,data.maxBlue]);

    data.intensityHistogramEdges = ...
      0:detection.intensityBinSize:(maxInt+detection.intensityBinSize-1);

    data.morphHist = histc(data.synapseIntensityMorphMean, ...
                           data.intensityHistogramEdges);

    data.synHist = histc(data.synapseIntensitySynMean, ...
                         data.intensityHistogramEdges);

    data.XHist = histc(data.synapseIntensityXMean, ...
		       data.intensityHistogramEdges);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function readExportSettings()

    exportInfo.saveMask = get(handles.saveMask,'Value');
    exportInfo.saveSholl = get(handles.saveSholl,'Value');
    exportInfo.saveIntensityHistogram = get(handles.saveIntensityHist,'Value');
    exportInfo.saveMat = get(handles.saveMat,'Value');
    exportInfo.saveSomaMeasure = get(handles.saveSomaMeasure,'Value');
    exportInfo.saveSynapseProfile = get(handles.saveSynapseProfile,'Value');
    exportInfo.saveGradient = get(handles.saveGradient,'Value');

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function writeExportSettings()

    set(handles.saveMask,'Value',exportInfo.saveMask);
    set(handles.saveSholl,'Value',exportInfo.saveSholl);
    set(handles.saveIntensityHist,'Value',exportInfo.saveIntensityHistogram);
    set(handles.saveMat,'Value',exportInfo.saveMat);
    set(handles.saveSomaMeasure,'Value',exportInfo.saveSomaMeasure);
    set(handles.saveSynapseProfile,'Value',exportInfo.saveSynapseProfile);
    set(handles.saveGradient,'Value',exportInfo.saveGradient);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function exportData(source, event)

    if(exist('source'))
      % The export function was called from a button, ask for save file
      setExportFiles();
    end

    if(isempty(data.exportXMLfile))
       % User aborted save.
      return
    end

    if(~detection.reanalyseFlag)
      progressBar = waitbar(0,sprintf('Exporting %s...', ...
				      strrep(data.fileName{dispInfo.curImg}, ...
					     '_','\_')));
    else
      % Do not show a progress bar if reanalysing, to prevent focus stealing
      progressBar = [];
    end

    barFig = gcf;
    nSteps = 6;

    readExportSettings();

    % Calculate data...
    if(exportInfo.saveSholl)
      shollAnalysis();
      exportWaitBar(1/nSteps, progressBar);
    else
      % Sholl normally calculates skeleton, but if we do not do that
      % lets do it here.
      
      data.skeleton = bwmorph(data.neuriteMask-data.somaMask>0,'skel',inf);
      trimSkeleton();

      calculateNeuriteLength();

    end

    if(exportInfo.saveIntensityHistogram)
      intensityHistograms();
      exportWaitBar(2/nSteps, progressBar);
    end

    calculateSomaMeasures();
    exportWaitBar(3/nSteps, progressBar);

    if(exportInfo.saveGradient)
      calculateIntensityGradients();
    end

    % exportMasks();

    exportToXML();
    exportWaitBar(4/nSteps, progressBar);

    % Also save the data in matlab readable format, just in case
    if(exportInfo.saveMat)
      saveData();
    end
    exportWaitBar(5/nSteps, progressBar);

    if(exportInfo.saveMask & ~isempty(progressBar))
      % Only need to return focus if we had a progressbar
      figure(handles.fig)
      saveMask();
    end

    if(~isempty(progressBar))
      % Only required if we have a progress bar
      figure(barFig)			    
    end
    exportWaitBar(6/nSteps, progressBar);

    % Increase the saved counter
    detection.numSaved = detection.numSaved + 1;

    saveConfig();

    delete(progressBar);
    % uiwait(helpdlg('All done. Tack för idag!', 'SynD'));

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setExportFiles()

    exportFile = strcat(data.fileName{dispInfo.curImg}, ...
			'-export.xml');

    curPwd = pwd;
    try
      cd(data.exportPath)
    catch
      fprintf('Unable to change to %s\n', data.exportPath)
    end

    [data.exportXMLfile, xmlExportPath] = ...
	  uiputfile('*-export.xml', ...
		    'Export synapse data', ...
		    exportFile);

    cd(curPwd);

    if(isempty(data.exportXMLfile) | data.exportXMLfile == 0)
      data.exportXMLfile = [];
      disp('Export aborted.')
      return
    end

    data.exportPath = xmlExportPath;

    if(strfind(data.exportXMLfile,'-export.xml'))
      data.exportSaveFile = ...
	strrep(data.exportXMLfile,'-export.xml','-save.mat');
      data.exportNeuriteMaskFile = ...
        strrep(data.exportXMLfile,'-export.xml','-neuritemask.tiff');
      data.exportSynapseMaskFile = ...
        strrep(data.exportXMLfile,'-export.xml','-synapsemask.tiff');
    else
      data.exportSaveFile = ...
	strcat(data.exportXMLfile,'-save.mat');
      data.exportNeuriteMaskFile = ...
        strcat(data.exportXMLfile,'-neuritemask.tiff');
      data.exportSynapseMaskFile = ...
        strcat(data.exportXMLfile,'-synapsemask.tiff');
    end


  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function saveMask()

    %%% Save neurite and soma masks

    % Get morphology image with true morphology colour
    setChannelVisible(1,0,0); 
    morphImg = getImg();

    % Show neurite mask with synapse colour
    morphImg(:,:,detection.synChannel) = data.neuriteMask;

    % Show soma with X-channel colour
    morphImg(:,:,detection.XChannel) = data.somaMask*0.5;

    imwrite(morphImg,strcat(data.exportPath,data.exportNeuriteMaskFile), ...
	    'tif', 'compression','lzw');

    %%% Save synapse picture

    setChannelVisible(0,1,0); 
    synImg = getImg();

    % Use morph channel to mark synapses
    synImg(:,:,detection.morphChannel) = 0.5*data.synapseMask;

    % Use X-channel to mark synapse centres
    tmp = zeros(data.height,data.width);
    tmp(data.synapseCenter) = 1;
    synImg(:,:,detection.XChannel) = tmp;

    imwrite(synImg,strcat(data.exportPath,data.exportSynapseMaskFile), ...
	    'tif', 'compression','lzw');


    setChannelVisible(1,1,1); 

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function exportWaitBar(progress,progressBar)
    if(~isempty(progressBar))
      waitbar(progress,progressBar);
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function saveData()

    clear old

    old.image = data.image;

    old.xyRes = data.xyRes;
    old.fileName = data.fileName;
    old.num = data.num;

    old.somaMask = data.somaMask;
    old.neuriteMask = data.neuriteMask;
    old.synapseMask = data.synapseMask;

    old.synapseCenter = data.synapseCenter;
    old.synapseArea = data.synapseArea;
    old.synapseDist = data.synapseDist;

    old.meanSynapseMorph = data.meanSynapseMorph;
    old.meanSynapseSyn = data.meanSynapseSyn;
    old.meanSynapseX = data.meanSynapseX;
    old.maxRadie = detection.maxRadie;
    old.excludeSomaSynapses = detection.excludeSomaSynapses;

    old.meanSynapseProfileDist = data.meanSynapseProfileDist;
    old.meanSynapseProfileMorph = data.meanSynapseProfileMorph;
    old.meanSynapseProfileSyn = data.meanSynapseProfileSyn;
    old.meanSynapseProfileX = data.meanSynapseProfileX;

    old.measurePoint = data.measurePoint;
    old.somaMeasureMask = data.somaMeasureMask;

    old.synapseIntensityMorphMean = data.synapseIntensityMorphMean;
    old.synapseIntensityMorphSEM = data.synapseIntensityMorphSEM;
    old.synapseIntensitySynMean = data.synapseIntensitySynMean;
    old.synapseIntensitySynSEM = data.synapseIntensitySynSEM;
    old.synapseIntensityXMean = data.synapseIntensityXMean;
    old.synapseIntensityXSEM = data.synapseIntensityXSEM;

    old.morphChannel = detection.morphChannel;
    old.synChannel = detection.synChannel;
    old.XChannel = detection.XChannel;
    old.remapChan = detection.remapChan;
    old.filterSize = detection.filterSize;

    old.totDendLength = data.neuriteLength;
    old.intensityHistogramEdges = data.intensityHistogramEdges;
    old.morphHist = data.morphHist;
    old.synHist = data.synHist;
    old.XHist = data.XHist;

    old.somaMeasureMorph = data.somaMeasureMorph;
    old.somaMeasureSyn = data.somaMeasureSyn;
    old.somaMeasureX = data.somaMeasureX;

    old.somaArea = data.somaArea;
    old.somaMajorAxisLength = data.somaMajorAxisLength;
    old.somaMinorAxisLength = data.somaMinorAxisLength;

    old.shollEdges = data.shollEdges;
    old.shollDendHist = data.shollDendHist;
    old.shollSynHist = data.shollSynHist;

    old.shollIntMorphMean = data.shollIntMorphMean;
    old.shollIntMorphStd = data.shollIntMorphStd;
    old.shollIntMorphSEM = data.shollIntMorphSEM;
    old.shollIntSynMean = data.shollIntSynMean;
    old.shollIntSynStd = data.shollIntSynStd;
    old.shollIntSynSEM = data.shollIntSynSEM;
    old.shollIntXMean = data.shollIntXMean;
    old.shollIntXStd = data.shollIntXStd;
    old.shollIntXSEM = data.shollIntXSEM;

    old.gradientEdges = data.gradientEdges;

    old.morphGradientMean = data.morphGradientMean;
    old.morphGradientStd = data.morphGradientStd; 
    old.morphGradientSEM = data.morphGradientSEM;

    old.synGradientMean = data.synGradientMean;
    old.synGradientStd = data.synGradientStd;
    old.synGradientSEM = data.synGradientSEM;

    old.XGradientMean = data.XGradientMean;
    old.XGradientStd = data.XGradientStd;
    old.XGradientSEM = data.XGradientSEM;

    old.dataImg = dispInfo.curImg;
    old.detection = detection;

    old.somaIntensityThreshold = data.somaIntensityThreshold;
    old.synapseIntensityThreshold = data.synapseIntensityThreshold;

    % Save version information also
    old.version = getVersion();

    if(data.exportSaveFile)
      fprintf('Saving data to %s\n', data.exportSaveFile);
      save(strcat(data.exportPath,data.exportSaveFile),'old')
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function importMask(source, event)

    if(isempty(data.image))
      
      disp('Unable to import mask, you must load corresponding lsm file first')
      beep

      % uiwait(warndlg('You must load the corresponding lsm file first.', ...
      %		     'Unable to import masks','modal'));

      return
    end

    curDir = pwd;

    try
      cd(data.exportPath)
    catch
      fprintf('Unable to change to %s\n', data.exportPath)
    end

    % Ask the user which file to load
    [dataFile,dataPath] = uigetfile('*-save.mat','Select old save file');

    data.exportPath = dataPath;
    cd(curDir);

    old = load(strcat(dataPath,dataFile));
    try
      old = old.old;
    catch
      fprintf('Unable to load file %s!\n', dataFile)
      return
    end


    if(~strcmp(old.fileName{1},data.fileName{1}))
      uiwait(warndlg(sprintf('Image file loaded: %s\nMask file loaded: %s', ...
			     data.fileName{1},old.fileName{1}), ...
		     'SynD : Mask import file mismatch'))

      data.somaIntensityThreshold = NaN;
    else

      % Load the old intensity threshold used
      try
        data.somaIntensityThreshold = old.somaIntensityThreshold;
      catch
        data.somaIntensityThreshold = old.detection.morphThreshold;
      end

    end

    saveUndo('loading masks')

    try
      data.neuriteMask = old.neuriteMask;
      data.somaMask = old.somaMask;
    catch
      % Old format had wrong name, revision 162...
      data.neuriteMask = old.neuriteMaskIdx;
      data.somaMask = old.somaMaskIdx;     
    end

    try
      data.somaMeasureMask = old.somaMeasureMask;
      exportInfo.saveSomaMeasure = 1;
      writeExportSettings();
    catch
      disp('No soma measure marked, ignoring.')
      data.somaMeasureMask = [];
    end

    if(data.height ~= size(data.neuriteMask,1) ...
       | data.width ~= size(data.neuriteMask,2))
      disp('Neurite and soma mask size does not match current image')

      % Clear the masks...
      data.neuriteMask = [];
      data.somaMask = [];
      activateStage('soma');
      return
    end

    if(~strcmp(data.fileName{1},old.fileName{1}))
      fprintf('Image file name: %s\nOld mask file name: %s\n', ...
	      data.fileName{1}, old.fileName{1})
      beep
    end

    updateNeuriteMask();

    data.distMask = makeDistMask(data.somaMask,data.neuriteMask);

    activateStage('synapses');

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function exportMasks()
    % !!! To be continued
    disp('exportMasks not implemented yet')
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % This exports to xml, readable by excel

  function exportToXML()

    if(isempty(data.exportXMLfile))
      return
    end

    CC = bwconncomp(data.somaMask);
    nSoma = CC.NumObjects; 

    if(nSoma > 1 & detection.singleSoma ...
       & ~detection.reanalyseFlag)
      % If we are reanalysing we are reusing an old soma mask, so the user
      % should know about this by now...

      uiwait(warndlg(sprintf('Single soma selected, but found %d somas.', ...
			     nSoma), ...
		     'Soma inconsistency error','modal'));
    end


    % Header info, to make excel feel comfortable with our data...
    docNode = makeXMLheader();

    % Create first sheet -- File info

    columnName = {'File name','Pixel size [um]', ...
		  'Morphology channel', ...
		  'Synapse channel', ...
		  'X channel', ...
                  'Number of somas', ...
		  'Soma synapses'};

    if(length(data.fileName) > 1)
      disp('Warning, export only handles one input file');
    end

    if(detection.excludeSomaSynapses)
      somaSynapsesStatus = 'excluded';
    else
      somaSynapsesStatus = 'included';
    end

    columnData = {{data.fileName{dispInfo.curImg}}, ...
		  [data.xyRes*1e6], ...
		  [detection.morphChannel], ...
		  [detection.synChannel], ...
		  [detection.XChannel], ...
                  [nSoma], ...
		  {somaSynapsesStatus}};

    makeXMLsheet(docNode,'File info', ...
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

    [y,x] = ind2sub(size(data.synapseMask),data.synapseCenter);

    % To avoid division by zero we discard those data points
    okMidx = find(data.synapseIntensityMorphMean ~= 0);
    okSidx = find(data.synapseIntensitySynMean ~= 0);

    columnData = {1:length(data.synapseCenter), ...
		  x, y, ...
		  data.synapseArea, ...
		  data.synapseIntensityMorphMean, ...
		  data.synapseIntensityMorphSEM, ...
		  data.synapseIntensitySynMean, ...
		  data.synapseIntensitySynSEM, ...
		  data.synapseIntensityXMean, ...
		  data.synapseIntensityXSEM, ...
		  data.synapseIntensitySynMean(okMidx) ...
		   ./data.synapseIntensityMorphMean(okMidx), ...
		  data.synapseIntensityXMean(okMidx) ...
		   ./data.synapseIntensityMorphMean(okMidx), ...
		  data.synapseIntensityXMean(okSidx) ...
		   ./data.synapseIntensitySynMean(okSidx) };


    makeXMLsheet(docNode,'Raw synapse data', ...
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

    nSynapses = length(data.synapseCenter);
    totLength = data.neuriteLength*1e6;


    columnData = { nSynapses, ...
		   totLength, ...
		   nSynapses/totLength, ...
		   mean(data.synapseArea), ...
		   mean(data.synapseIntensityMorphMean), ...
		   mean(data.synapseIntensitySynMean), ...  
		   mean(data.synapseIntensityXMean), ...  
		   nanmean(data.synapseIntensitySynMean(okMidx) ...
			   ./data.synapseIntensityMorphMean(okMidx)), ...
		   nanmean(data.synapseIntensityXMean(okMidx) ...
			   ./data.synapseIntensityMorphMean(okMidx)), ...
		   nanmean(data.synapseIntensityXMean(okSidx) ...
			   ./data.synapseIntensitySynMean(okSidx)) };

    makeXMLsheet(docNode,'Synapse averages', ...
		 columnName, columnData);

    if(exportInfo.saveIntensityHistogram)

      % Write histogram with channel intensities

      columnName = { 'From', 'To', ...
		     'Morph', 'Syn', 'X', ...
		     'Cum morph', 'Cum syn', 'Cum X' };
				      

      columnData = { data.intensityHistogramEdges(1:end-1), ...
		     data.intensityHistogramEdges(2:end), ...
		     data.morphHist(1:end-1), ...
		     data.synHist(1:end-1), ...
		     data.XHist(1:end-1), ...
		     cumsum(data.morphHist(1:end-1)), ...
		     cumsum(data.synHist(1:end-1)), ...
		     cumsum(data.XHist(1:end-1)) };

      makeXMLsheet(docNode,'Intensity histograms', ...
		   columnName, columnData);
    end


    if(exportInfo.saveSholl)
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

      columnData = { data.shollEdges(1:end-1)*1e6, data.shollEdges(2:end)*1e6, ...
		     data.shollDendHist, ...
		     data.shollSynHist, ...
		     data.shollIntMorphMean, ...
		     data.shollIntMorphStd, ...
		     data.shollIntMorphSEM, ...
		     data.shollIntSynMean, ...
		     data.shollIntSynStd, ...
		     data.shollIntSynSEM, ...
		     data.shollIntXMean, ...
		     data.shollIntXStd, ...
		     data.shollIntXSEM, ...
		     data.shollIntSynMean./data.shollIntMorphMean, ...
		     data.shollIntXMean./data.shollIntMorphMean, ...
		     data.shollIntXMean./data.shollIntSynMean, ...
                   };

      makeXMLsheet(docNode,'Sholl analysis', ...
		   columnName, columnData);

    end
		   

    if(exportInfo.saveSomaMeasure)
      columnName = { 'Soma morph (mean)', ...
		     'Soma morph (std)', ...
		     'Soma morph (SEM)', ...
		     'Soma syn (mean)', ...
		     'Soma syn (std)', ...
		     'Soma syn (SEM)', ...
		     'Soma X (mean)', ...
		     'Soma X (std)', ...
		     'Soma X (SEM)' };

       columnData = { mean(data.somaMeasureMorph), ...
	  	      std(data.somaMeasureMorph), ...
		      std(data.somaMeasureMorph) ...
		        /sqrt(length(data.somaMeasureMorph)), ...	 
		      mean(data.somaMeasureSyn), ...
		      std(data.somaMeasureSyn), ...		    
		      std(data.somaMeasureSyn) ...
  		        / sqrt(length(data.somaMeasureSyn)), ...       
		      mean(data.somaMeasureX), ...
		      std(data.somaMeasureX), ...		    
		      std(data.somaMeasureX) ...
		        / sqrt(length(data.somaMeasureX)), ...
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
		  mean(data.somaArea), ...
		  std(data.somaArea), ...
		  std(data.somaArea) ...
		    / sqrt(length(data.somaArea)), ...
		  mean(data.somaMajorAxisLength), ...
		  std(data.somaMajorAxisLength), ...
		  std(data.somaMajorAxisLength) ...
  		    / sqrt(length(data.somaMajorAxisLength)), ...
		  mean(data.somaMinorAxisLength), ...
		  std(data.somaMinorAxisLength), ...
		  std(data.somaMinorAxisLength) ...
  		    / sqrt(length(data.somaMinorAxisLength)), ...
                 };

    % Concatenate cell arrays
    columnName = {columnName{:},columnName2{:}};
    columnData = {columnData{:},columnData2{:}}; 

    makeXMLsheet(docNode,'Soma measures', ...
		 columnName, columnData);



    if(exportInfo.saveSynapseProfile)

      columnName = {'Distance (micrometer)', ...
		    'Morphology', ...
		    'Synapse', ...
		    'X', ...
		    'Syn/Morph', ...
		    'X/Morph', ...
		    'X/Syn' };

      columnData = { data.meanSynapseProfileDist*1e6, ...
		     data.meanSynapseProfileMorph, ...
		     data.meanSynapseProfileSyn, ...
		     data.meanSynapseProfileX, ...
		     data.meanSynapseProfileSyn ...
		     ./data.meanSynapseProfileMorph, ...
		     data.meanSynapseProfileX ...
		     ./data.meanSynapseProfileMorph, ...
		     data.meanSynapseProfileX ...
		     ./data.meanSynapseProfileSyn, ...
                   };


      makeXMLsheet(docNode,'Synapse profile', ...
		   columnName, columnData);
    end

    if(exportInfo.saveGradient)
      columnName = { 'From (mu)','To (mu)', ...
		     'Morph (mean)', 'Morph (std)', 'Morph (SEM)', ...
		     'Syn (mean)', 'Syn (std)', 'Syn (SEM)', ...
		     'X (mean)', 'X (std)', 'X (SEM)', ...	     
		     'Syn/Morph (mean)', ...
		     'X/Moph (mean)', ...
                   };

      columnData = { [0,data.gradientEdges(1:end-1)]*1e6, ...
		     data.gradientEdges*1e6,...
		     data.morphGradientMean, ...
		     data.morphGradientStd, ...
		     data.morphGradientSEM, ...
		     data.synGradientMean, ...
		     data.synGradientStd, ...
		     data.synGradientSEM, ...
		     data.XGradientMean, ...
		     data.XGradientStd, ...
		     data.XGradientSEM, ...
		     data.synGradientMean./data.morphGradientMean, ...
		     data.XGradientMean./data.morphGradientMean, ...
                   };

      makeXMLsheet(docNode,'Neurite gradients', ...
		   columnName, columnData);

    end

    % Write all to disk

    fprintf('Exporting data to %s\n',data.exportXMLfile);
    xmlwrite(strcat(data.exportPath,data.exportXMLfile),docNode);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function docNode = makeXMLheader()

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


  function makeXMLsheet(docNode, ...
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
      % docData.setAttribute('ss:Type','Number');
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

  function skipSynapseDetection(source, event)

    activateStage('analyse');
    stopEdit();

    if(max(max(data.distMask .* data.neuriteMask)) == inf)
      uiwait(errordlg(['Manually edit the neurite mask to connect the dark ' ...
                       'grey neurites to the soma or use the clean button ' ...
                       'to automatically remove unconnected components!'], ...
 	               'Unconnected components', 'modal'))
      return             
    end

    dispInfo.axis = []; % Reset view
    showGUI('analyse');

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function exportOnlyAverages(source, event)

    exportAverageFile = strcat(data.fileName{dispInfo.curImg}, ...
			       '-onlyAverage.csv');

    curPwd = pwd;
    try
      cd(data.exportPath)
    catch
      fprintf('Unable to change to %s\n', data.exportPath)
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

    data.exportPath = savePath;

    tmpMorph = data.image(:,:,detection.morphChannel, dispInfo.curImg);
    tmpSyn = data.image(:,:,detection.synChannel, dispInfo.curImg);
    tmpX = data.image(:,:,detection.XChannel, dispInfo.curImg);

    fid = fopen(strcat(data.exportPath,exportAverageFile),'w');
    fprintf(fid, ['File,Soma area,Soma M,Soma M (SEM),Soma S,Soma S (SEM),' ...
                  'Soma X,Soma X (SEM),Neurite area,Neurite M,' ...
                  'Neurite M (SEM),Neurite S,Neurite S (SEM),Neurite X,' ...
                  'Neurite X (SEM)\n']);

    fprintf(fid,'%s,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n', ...
	    data.fileName{dispInfo.curImg}, ...
	    nnz(data.somaMask)*(data.xyRes*1e6)^2, ...
	    mean(tmpMorph(find(data.somaMask))), ...
	    std(tmpMorph(find(data.somaMask))) ...
	      /sqrt(nnz(data.somaMask)), ...
	    mean(tmpSyn(find(data.somaMask))), ...
	    std(tmpSyn(find(data.somaMask))) ...
	      /sqrt(nnz(data.somaMask)), ...
	    mean(tmpX(find(data.somaMask))), ...
	    std(tmpX(find(data.somaMask))) ...
	      /sqrt(nnz(data.somaMask)), ...
	    nnz(data.neuriteMask)*(data.xyRes*1e6)^2, ...
	    mean(tmpMorph(find(data.neuriteMask))), ...
	    std(tmpMorph(find(data.neuriteMask))) ...
	      /sqrt(nnz(data.neuriteMask)), ...
	    mean(tmpSyn(find(data.neuriteMask))), ...
	    std(tmpSyn(find(data.neuriteMask))) ...
	      /sqrt(nnz(data.neuriteMask)), ...
	    mean(tmpX(find(data.neuriteMask))), ...
	    std(tmpX(find(data.neuriteMask))) ...
	      /sqrt(nnz(data.neuriteMask)));

    fclose(fid);

    saveConfig();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function nextStage(source, event)

    % Make sure the next stage is active before advancing

    % Turn of any motion handles
    set(handles.fig,'windowbuttonmotionfcn', [])

    switch(dispInfo.state)
      case 'load'
        stopEdit();

	if(dispInfo.stage >= 2)
          dispInfo.axis = []; % Reset view
	  showGUI('soma');
          if(~nnz(data.somaMask))
            % No pixels in the soma mask, try and detect it
            detectSoma();
          end
        end

      case 'soma'
        stopEdit();

	if(dispInfo.stage >= 3)
          if(nnz(data.somaMask))      
            dispInfo.axis = []; % Reset view  
	    showGUI('neurites');
            if(~nnz(data.neuriteMask))
	      detectNeurite();
            end
          else
            % There are no pixels in the soma mask...
            uiwait(errordlg('You must detect or draw a soma!', ...
			    'No soma detected', 'modal'))
            return
          end
        end

      case 'neurites'
        stopEdit();

	if(dispInfo.stage >= 4)
          % Verify that there are no unconnected neurites...
	      
          if(max(max(data.distMask .* data.neuriteMask)) == inf)
            uiwait(errordlg(['Manually edit the neurite mask to connect the dark ' ...
                             'grey neurites to the soma or use the clean button ' ...
                             'to automatically remove unconnected components!'], ...
			    'Unconnected components', 'modal'))
            return             
          end

          dispInfo.axis = []; % Reset view
	  showGUI('synapses');
          if(~nnz(data.synapseMask))
  	    detectSynapses();
          end
        end

      case 'synapses'
        stopEdit();

	if(dispInfo.stage >= 5)
          dispInfo.axis = []; % Reset view
	  showGUI('analyse');
        end

      case 'analyse'
        stopEdit();

        dispInfo.axis = []; % Reset view
        % We are already at last stage, do nothing
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function activeIcon(iconNumber)

    for i = 1:length(iconNumber)
      set(handles.allIcons(iconNumber(i)),'backgroundcolor',[1 0 0]);
    end

    restOfIcons = setdiff(1:length(handles.allIcons),iconNumber);
    for i = 1:length(restOfIcons)
      set(handles.allIcons(restOfIcons(i)),'backgroundcolor',get(gcf,'color'));
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function showGUI(GUIstate)

    makeVisible = [];
    makeInvisible = [];

    % Update icons so they are correct based on current stage
    % of analysis.
    setIcons();

    switch(GUIstate)
      case 'load'
        % The GUI shown when we are loading image and selecting
        % which frames in a stack to display.
        % Also used to select which channel are which

 	makeVisible = loadGUI;
        makeInvisible = setdiff(allGUI,loadGUI);

        if(data.num <= 1)
	  multiImgButtons = [handles.next,handles.prev,handles.num];
          makeVisible = setdiff(makeVisible,multiImgButtons);
          makeInvisible = union(makeInvisible, multiImgButtons);
        end

        dispInfo.state = 'load';  
        setChannelVisible(1,1,1);
        nextStageButtonStatus(1);

        % Used for handling the histograms
        set(handles.fig,'WindowButtonDownFcn', @mouseHandler);

        set(handles.XYres,'String', num2str(data.xyRes*1e6));

        activeIcon(1);
        uicontrol(handles.loadNeuron);

      case 'soma'
        % The GUI shown when we are locating the soma

        if(strcmp(dispInfo.state,'load'))
          % Moving from load? If so remove histogram listener
          set(handles.fig,'WindowButtonDownFcn', []);
        end

	makeVisible = somaGUI;
        makeInvisible = setdiff(allGUI,somaGUI);
        dispInfo.state = 'soma';  
        setChannelVisible(1,0,0);
        nextStageButtonStatus(2);

        dispInfo.somaColor = [1 1 1];
        dispInfo.neuriteColor = NaN;
        dispInfo.synapseColor = NaN;
        dispInfo.measurePointColor = NaN;

        activeIcon(2);
        % uicontrol(handles.detectSoma);

      case 'neurites'
        % The GUI to locate neurites and edit neurite mask

        if(strcmp(dispInfo.state,'load'))
          % Moving from load? If so remove histogram listener
          set(handles.fig,'WindowButtonDownFcn', []);
        end

	makeVisible = neuriteGUI;
        makeInvisible = setdiff(allGUI,neuriteGUI);

        % Alternate which GUI is shown depending on settings
        if(detection.useSteerableFilters)
	  cfgHandles = [handles.morphThresh, handles.morphThreshLabel];
        else
	  cfgHandles = [handles.growThresh, ...
			handles.growThreshLabel, ...
 		        handles.filterSize, ...
 		        handles.filterSizeLabel, ...
			handles.addNeurite];
        end

        makeVisible = setdiff(makeVisible,cfgHandles);
        makeInvisible = union(makeInvisible,cfgHandles);

        dispInfo.state = 'neurites';  
        setChannelVisible(1,0,0);
        nextStageButtonStatus(3);

        dispInfo.somaColor = [1 1 1];
        dispInfo.neuriteColor = [1 1 1]*0.7;
        dispInfo.synapseColor = NaN;
        dispInfo.measurePointColor = NaN;

        activeIcon(3);
        % uicontrol(handles.detectNeurite);

      case 'synapses'
        % The GUI for synapse detection

        if(strcmp(dispInfo.state,'load'))
          % Moving from load? If so remove histogram listener
          set(handles.fig,'WindowButtonDownFcn', []);
        end

	makeVisible = synapseGUI;
        makeInvisible = setdiff(allGUI,synapseGUI);

        dispInfo.state = 'synapses';  
        setChannelVisible(0,1,0);
        nextStageButtonStatus(4);

        dispInfo.somaColor = NaN;
        dispInfo.neuriteColor = NaN;
        dispInfo.synapseColor = [1 1 1];
        dispInfo.measurePointColor = NaN;

        activeIcon(4);
        % uicontrol(handles.detectSynapses);

      case 'analyse'
        % Analyse and export GUI

        % Reset viewpoint
        dispInfo.axis = [];

        if(strcmp(dispInfo.state,'load'))
          % Moving from load? If so remove histogram listener
          set(handles.fig,'WindowButtonDownFcn', []);
        end

	makeVisible = analyseGUI;
        makeInvisible = setdiff(allGUI,analyseGUI);
        dispInfo.state = 'analyse';  
        setChannelVisible(1,1,1);

        activeIcon(5);
        uicontrol(handles.exportData);

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
    set([handles.redHist, handles.greenHist, handles.blueHist],'visible','off');

    showImage();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function mouseHandler(source, event)

     if(isempty(data.image))
       % Nothing loaded yet, ignore
       return
     end

    [insideRed, xRed, yRed] = checkInside(handles.redHist);
    [insideGreen, xGreen, yGreen] = checkInside(handles.greenHist);
    [insideBlue, xBlue, yBlue] = checkInside(handles.blueHist);

    if(insideRed)
      dispInfo.scaleRed = data.maxRed / xRed;
      fprintf('Setting red scale to %d\n', dispInfo.scaleRed)
    elseif(insideGreen)
      dispInfo.scaleGreen = data.maxGreen / xGreen;
      fprintf('Setting green scale to %d\n', dispInfo.scaleGreen)
    elseif(insideBlue)
      dispInfo.scaleBlue = data.maxBlue / xBlue;
      fprintf('Setting blue scale to %d\n', dispInfo.scaleBlue)
    end
 
    showImage();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Check that we clicked inside axis of figure

  function [insideFlag, x, y] = checkInside(testHandle)

    tmpXY = get(testHandle,'CurrentPoint');
    x = tmpXY(1,1); 
    y = tmpXY(1,2);

    set(handles.fig,'CurrentAxes',testHandle)
    a = axis();

    if(a(1) <= x & x <= a(2) & a(3) <= y & y <= a(4))
      insideFlag = 1;
    else
      insideFlag = 0;
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function startZoomWheelMode(source,event)

    % disp('Zoom mode active, click to center, scroll wheel to zoom')
    set(handles.fig,'WindowScrollWheelFcn', @zoomHandler);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function zoomHandler(source, event)

    set(handles.fig,'CurrentAxes',handles.image)

    a = axis();

    xCenter = (a(1)+a(2))/2;
    yCenter = (a(3)+a(4))/2;

    aWidth = a(2)-a(1);
    aHeight = a(4)-a(3);

    % Where did the user click, within the image?
    xy = get(handles.image,'CurrentPoint');
    xCenter = xy(1,1);
    yCenter = xy(1,2);

    if(xCenter < a(1) | a(2) < xCenter ...
       | yCenter < a(3) | a(4) < yCenter)

      disp('Outside of figure, ignoring')
      % set(handles.fig,'WindowScrollWheelFcn', '');

      return
    end

    aHeight = aHeight * (1.25)^event.VerticalScrollCount;
    aWidth = aWidth * (1.25)^event.VerticalScrollCount;

    xCenter = min(max(aWidth/2+1,xCenter),data.width-aWidth/2);
    yCenter = min(max(aHeight/2+1,yCenter),data.height-aHeight/2);

    if(aWidth > data.width | aHeight > data.height)
      a = [1 data.width 1 data.height];
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

  function [nSynapse, synArea, synStd] = synapseThresholdSensitivity(threshRange)

    oldThresh = detection.synThreshold;

    nSynapse = NaN*ones(size(threshRange));
    synArea = NaN*ones(size(threshRange));
    synStd = NaN*ones(size(threshRange));

    for iT = 1:length(threshRange)
      fprintf('Thresh: %d\n', threshRange(iT))
      synStd(iT) = threshRange(iT);

      detection.synThreshold = threshRange(iT);
      detectSynapseCenters();
      calculateSynapseProperties();
      nSynapse(iT) = length(data.synapseCenter);
      synArea(iT) = nnz(data.synapseMask);
    end

    detection.synThreshold = oldThresh;

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

  function resizeFunction(source, event)

    % Resize the bitmaps for the icons
    setIcons();

    if(isempty(data.image))
      % disp('Resizing and redisplaying logo.')
      showLogo();
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setIcons()

    % We want the icons to fit the push buttons

    % Find out how big the button is
    set(handles.loadStage,'units','pixels')
    buttonPos = get(handles.loadStage,'position');
    buttonSize = buttonPos(3:4)-4;
    set(handles.loadStage,'units','normalized')

    
    % Make sure to desaturete the icons that should be

    maxStage = 5;
    stageFlag = (1:maxStage) <= dispInfo.stage;


    for i = 1:length(bilder.icon)
      bilder.iconResize{i} = imresize(bilder.icon{i},buttonSize([2 1]));
    end

    for i = 1:length(stageFlag)
      if(stageFlag(i))
	set(allIcons(i), ...
	    'Callback',allIconCallback{i}, ...
	    'Cdata', bilder.iconResize{i});
				      
      else
	set(allIcons(i), ...
	    'Callback',[], ...
	    'Cdata', desaturateIcon(bilder.iconResize{i}));
      end
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function loadStage(source, event)
    dispInfo.axis = []; % Reset view
    stopEdit();
    showGUI('load');
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function somaStage(source, event)
    dispInfo.axis = []; % Reset view
    stopEdit();
    showGUI('soma');
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function neuriteStage(source, event)
    dispInfo.axis = []; % Reset view
    stopEdit();
    showGUI('neurites');
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function synapseStage(source, event)
    dispInfo.axis = []; % Reset view
    stopEdit();
    showGUI('synapses');
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function analyseStage(source, event)
    dispInfo.axis = []; % Reset view
    stopEdit();
    showGUI('analyse');
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function activateStage(stage)

    switch(stage)
      case 'load'
        dispInfo.stage = 1;

      case 'soma'
        set(handles.loadNeuron, 'BackgroundColor', get(gcf,'color'));
        dispInfo.stage = 2;

      case 'neurites'
        set(handles.detectSoma, 'BackgroundColor', get(gcf,'color'));
        dispInfo.stage = 3;

      case 'synapses'
        set(handles.detectNeurite, 'BackgroundColor', get(gcf,'color'));
        dispInfo.stage = 4;

      case 'analyse'
        set(handles.detectSynapses, 'BackgroundColor', get(gcf,'color'));
        dispInfo.stage = 5;

    end

    showGUI(dispInfo.state);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function bwIcon = desaturateIcon(colourIcon)

    bwIcon = zeros(size(colourIcon));

    for i = 1:3
      bwIcon(:,:,i) = mean(colourIcon,3)/255.0;
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function loadConfig()

    if(exist(configFile))
      fprintf('Loading configuration from %s\n', configFile)
      cfg = load(configFile);

      try
        detection.morphChannel = cfg.old.morphChannel;
        detection.synChannel = cfg.old.synChannel;
        detection.XChannel = cfg.old.XChannel;
        detection.remapChan = cfg.old.remapChan;
        detection.enableRemapChan = cfg.old.enableRemapChan;

        detection.filterSize = cfg.old.filterSize;
        detection.singleSoma = cfg.old.singleSoma;

        data.xyRes = cfg.old.xyRes;

        detection.somaErodeRadius = cfg.old.somaErodeRadius;
        detection.morphThreshold = cfg.old.morphThreshold;
        detection.autoSomaThreshold = cfg.old.autoSomaThreshold;

        detection.measurePointRadius = cfg.old.measurePointRadius;
        detection.nMeasurePoints = cfg.old.nMeasurePoints;

        detection.synThreshold = cfg.old.synThreshold;
        detection.minSynapseSize = cfg.old.minSynapseSize;
        detection.singleSynapseRadie = cfg.old.singleSynapseRadie;
        detection.neuritePadding = cfg.old.neuritePadding;
        detection.maxRadie = cfg.old.maxRadie;
        detection.excludeSomaSynapses = cfg.old.excludeSomaSynapses;

        detection.shollBinSize = cfg.old.shollBinSize;

        exportInfo.saveMask = cfg.old.saveMask;

        % exportInfo.saveSholl = cfg.old.saveSholl;
        exportInfo.saveSholl = 1; % Always set it to save

        exportInfo.saveIntensityHistogram = cfg.old.saveIntHist;
        exportInfo.saveMat = cfg.old.saveMat;

        data.loadPath = cfg.old.loadPath;
        data.exportPath = cfg.old.exportPath;

        editInfo.fileTypeOrder = cfg.old.fileTypeOrder;

        dispInfo.figPosition = cfg.old.figPosition;

        detection.numSaved = cfg.old.numSaved;

      catch exception
        getReport(exception)
        disp('Config file belonged to old version.')
      end
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function saveSynDpath(source, event)
    disp('Saving current SynD path permanently')
    savepath
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setSynDpath()

    SyndScriptPath = SynDscriptFull(1:(max(strfind(SynDscriptFull, ...
						 SynDscript))-1));

    SynDpath = { SyndScriptPath, ...
		 strcat(SyndScriptPath,'LSM/lsm'), ...
		 strcat(SyndScriptPath,'LSM/cstruct'), ...
		 strcat(SyndScriptPath,'bilder') };

    fprintf('Set SynD work path to %s\n', SynDpath{1})

    for i = 1:length(SynDpath)
      addpath(SynDpath{i});
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function saveConfig()

    fprintf('Saving old configuration to %s\n', configFile)

    old.morphChannel = detection.morphChannel;
    old.synChannel = detection.synChannel;
    old.XChannel = detection.XChannel;
    old.remapChan = detection.remapChan;
    old.enableRemapChan = detection.enableRemapChan;

    old.filterSize = detection.filterSize;
    old.singleSoma = detection.singleSoma;

    old.xyRes = data.xyRes;

    old.somaErodeRadius = detection.somaErodeRadius;
    old.morphThreshold = detection.morphThreshold;
    old.autoSomaThreshold = detection.autoSomaThreshold;

    old.measurePointRadius = detection.measurePointRadius;
    old.nMeasurePoints = detection.nMeasurePoints;

    old.synThreshold = detection.synThreshold;
    old.minSynapseSize = detection.minSynapseSize;
    old.singleSynapseRadie = detection.singleSynapseRadie;
    old.neuritePadding = detection.neuritePadding;
    old.maxRadie = detection.maxRadie;
    old.excludeSomaSynapses = detection.excludeSomaSynapses;

    old.shollBinSize = detection.shollBinSize;

    old.saveMask = exportInfo.saveMask;
    old.saveSholl = exportInfo.saveSholl;
    old.saveIntHist = exportInfo.saveIntensityHistogram;
    old.saveMat = exportInfo.saveMat;
    old.loadPath = data.loadPath;
    old.exportPath = data.exportPath;

    old.fileTypeOrder = editInfo.fileTypeOrder;

    old.figPosition = get(handles.fig,'position');

    old.numSaved = detection.numSaved;

    try 
      save(configFile, 'old');
    catch
      fprintf('Unable to save config file %s.\n', configFile)
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function nextStageButtonStatus(dispStage)
    if(dispStage < dispInfo.stage)
      % We can advance to next stage
      set(handles.nextStage, 'BackgroundColor', [1 0.6 0]);

      % Set focus to next button
      uicontrol(handles.nextStage);
    else
      % Still not finished with this stage
      set(handles.nextStage, 'BackgroundColor', [0.8 0.8 0.8]);
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function captureKeyPress(source, event)

    try
      % Was it a keypress?
      event.Key;
    catch
      % ginput messes with scroll wheel handling, fixing it!
      startZoomWheelMode();
      zoomHandler(source, event);
      return
    end

    % Shared keys for all
    switch(event.Key)
      case 'r'
	toggleRed();
      case 'g'
        toggleGreen();
      case 'b'
        toggleBlue();
    end

    switch(event.Character)
      case '+'

	if(isempty(data.image))
          return
        end

        % Zoom in
        set(handles.fig,'CurrentAxes',handles.image)    
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

	  xCenter = min(max(aWidth/2+1,xCenter),data.width-aWidth/2);
          yCenter = min(max(aHeight/2+1,yCenter),data.height-aHeight/2);

          if(aWidth > data.width | aHeight > data.height)
	    a = [1 data.width 1 data.height];
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
	if(isempty(data.image))
          return
        end

        % Zoom out
        set(handles.fig,'CurrentAxes',handles.image)    
        a = axis;
        xCenter = mean(a(1:2));
        yCenter = mean(a(3:4));
        aWidth = min((a(2)-a(1))*2,data.width);
        aHeight = min((a(4)-a(3))*2,data.height);

	xCenter = min(max(aWidth/2+1,xCenter),data.width-aWidth/2);
        yCenter = min(max(aHeight/2+1,yCenter),data.height-aHeight/2);

        a(1) = xCenter-aWidth/2;
        a(2) = xCenter+aWidth/2;          
        a(3) = yCenter-aHeight/2;
        a(4) = yCenter+aHeight/2;
        axis(a);

    end

    % Separate keys depending on state
    switch(dispInfo.state)
      case 'load'
        switch(event.Key)
          case 'l'
	    loadNeuron();
          case 'n'
	    nextStage();
          case 'x'
            % Focus to the XY-res
	    uicontrol(handles.XYres)
        end
      case 'soma'
        switch(event.Key)
          case 'a'
            addSoma();
          case 'd'
	    detectSoma();
          case 'n'
	    nextStage();
          case 'e'
            editSoma();
          case 'm'
	    toggleMask();
          case 's'
            % Change focus
	    uicontrol(handles.somaErodeRadius);
          case 't'
	    uicontrol(handles.morphThresh);
          case 'u'
            undoLastAction();
          case 'r'
            setSomaMeasurePoints();
          case 'o'
	    addSomaMeasurePoints();
        end
     
      case 'neurites'
        switch(event.Key)
          case 'd'
	    detectNeurite();
          case 'n'
	    nextStage();
          case 'e'
            editNeurite();
          case 'c'
	    cleanMask();
          case 'a'
            if(detection.useSteerableFilters)
	      addNeurite();
            end
          case 'o'
	    killBlob();
          case 'm'
	    toggleMask();
          case 'f'
            if(detection.useSteerableFilters)
  	      uicontrol(handles.growThresh);
            else
  	      uicontrol(handles.morphThresh);
            end
          case 't'
            addThinNeurites();
          case 'u'
            undoLastAction();
          case 's'
            showSkeleton();
        end

      case 'synapses'
        switch(event.Key)
          case 'd'
	    detectSynapses();
          case 'n'
	    nextStage();
          case 'e'
            editSynapse();
          case 'm'
	    toggleMask();
          case 't'
	    uicontrol(handles.synThresh);
          case 'u'
            undoLastAction();
        end

      case 'analyse'
        switch(event.Key)
          case 'e'
	    exportData();
        end

    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
      
