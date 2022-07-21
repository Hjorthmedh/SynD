% Synapse Detector (SynD)
% 
% Johannes Hjorth
% j.j.j.hjorth@damtp.cam.ac.uk
% johannes.hjorth@cncr.vu.nl
% hjorth@kth.se
%
% Sabine Schmitz
% sabine.schmitz@cncr.vu.nl
%


function SynD_aggregate2022()

  close all

  data = [];

  conserveMemory = 1;

  groupIdx = []; % Index of group each neuron belongs to

  % Group info structure
  group = struct('Name',[], ...
		 'Members', []);

  morphLabel = 'MAP2';
  synLabel = 'VAMP2';
  XLabel = 'Munc18';

  exportFileMask = '%s-export-summary.xml';
  exportFile = [];
  exportPath = pwd;
  loadPath = pwd;

  plotDir = [];
  plotExperimentName = [];

  configFile = strrep(mfilename('fullpath'), ... 
		      mfilename(), ...
		      'SynD_aggregate_config.mat');

  versionFile = strrep(mfilename('fullpath'), ... 
		      mfilename(), ...
		      'version.txt');


  loadConfig();

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  handles.fig = figure('Name','SynD - Aggregate data (also known as BarB)', ...
		       'MenuBar','none', ...
		       'Toolbar','none', ...
		       'Position', [50 50 1150 680]);

  handles.plot = axes('Units','Pixels', 'Position', [660 250 400 400]);

  % !!! Disabled hotkeys
  %set(handles.fig,'WindowKeyPressFcn',@captureKeyPress);

  handles.selectData = uicontrol('Style','listbox', ...
				 'String', 'None', ...
				 'ToolTipString', ...
				 'Choose data to group', ...
				 'Value', [1], ...
				 'Max', 1, 'Min', 1, ...
				 'Interruptible', 'off', ...
				 'Callback', @dataListboxClicked, ...
				 'Position',[50 50 380 600]);

  handles.groupNameLabel = uicontrol('Style', 'text', ...
				     'String', 'Group name:', ...
				     'HorizontalAlignment','left', ...
				     'BackgroundColor', 0.8*[1 1 1], ...
				     'Position', [450 620 100 25]);


  handles.groupName = uicontrol('Style','edit', ...
				'String', 'WT', ...
				'HorizontalAlignment','left', ...
				'Interruptible', 'off', ...
				'Position', [450 595 150 30]);

  handles.morphNameLabel = uicontrol('Style', 'text', ...
				     'String', 'Morph marker:', ...
				     'HorizontalAlignment','left', ...
				     'BackgroundColor', 0.8*[1 1 1], ...
				     'Position', [450 560 100 25]);

  handles.morphName = uicontrol('Style','edit', ...
				'String', morphLabel, ...
				'HorizontalAlignment','left', ...
				'Interruptible', 'off', ...
				'Position', [450 535 150 30], ...
				'Callback', @updateLabels);

  handles.synNameLabel = uicontrol('Style', 'text', ...
				     'String', 'Synapse marker:', ...
				     'HorizontalAlignment','left', ...
				     'BackgroundColor', 0.8*[1 1 1], ...
				     'Position', [450 500 100 25]);

  handles.synName = uicontrol('Style','edit', ...
			      'String', synLabel, ...
			      'HorizontalAlignment','left', ...
			      'Interruptible', 'off', ...
			      'Position', [450 475 150 30], ...
			      'Callback', @updateLabels);

  handles.XNameLabel = uicontrol('Style', 'text', ...
				 'String', 'X marker:', ...
				 'HorizontalAlignment','left', ...
				 'BackgroundColor', 0.8*[1 1 1], ...
				 'Position', [450 440 100 25]);

  handles.XName = uicontrol('Style','edit', ...
			    'String', XLabel, ...
			    'HorizontalAlignment','left', ...
			    'Interruptible', 'off', ...
			    'Position', [450 415 150 30], ...
			    'Callback', @updateLabels);

  handles.makeGroup = uicontrol('Style','pushbutton', ...
				'String','<html><u>D</u>efine group</html>', ...
				'Interruptible', 'off', ...
				'Position',[450,360,100,30], ...
				'Callback', @defineGroup);

  handles.guessGroup = uicontrol('Style','pushbutton', ...
				'String','<html>Find</html>', ...
				'fontsize', 8, ...
				'Interruptible', 'off', ...
				'Position',[560,370,40,20], ...
				'Callback', @guessGroup);


  handles.makeGroup = uicontrol('Style','pushbutton', ...
				'String','<html><u>E</u>xport data</html>', ...
				'Interruptible', 'off', ...
				'Position',[450,320,100,30], ...
				'Callback', @exportData);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  handles.plotTypeLabel = uicontrol('Style','text', ...
				    'String', 'Show plot:', ...
				    'HorizontalAlignment', 'left', ...
				    'BackgroundColor', 0.8*[1 1 1], ...
				    'Position', [660 180 100 30]);

  handles.plotType = uicontrol('Style','popupmenu', ...
			       'String', {'Plot type list'}, ...
			       'fontsize', 10, ...
			       'Position', [660 160 350 30], ...
			       'Value', 1, ...
			       'Callback', @showPlot);

  handles.saveAllPlots = uicontrol('Style','pushbutton', ...
				   'String','Save all plots', ...
				   'Interruptible', 'off', ...
				   'Position',[740,120,100,30], ...
				   'Callback', @saveAllPlots);

  handles.savePlot = uicontrol('Style','pushbutton', ...
			       'String','Save plot', ...
			       'Interruptible', 'off', ...
			       'Position',[660,120,70,30], ...
			       'Callback', @savePlotCallback);

  handles.saveSummaryPlots = uicontrol('Style','pushbutton', ...
				       'String','Save summary plots', ...
				       'Interruptible', 'off', ...
				       'Position',[850,120,140,30], ...
				       'Callback', @makeAllSummaryPlots);

  saveFileType = {'eps','pdf','png'};

  handles.fileType = uicontrol('Style','popupmenu', ...
			       'String', saveFileType, ...
			       'fontsize', 10, ...
			       'Position', [1020 160 70 30], ...
			       'Value', 1);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  for iGB = 1:12
    boxX = mod(iGB-1,4);
    boxY = 2-floor((iGB-1)/4);

    handles.groupLabels(iGB) = uicontrol('Style','text', ...
					 'String',num2str(iGB), ...
					 'BackgroundColor', 0.8*[1 1 1], ...
					 'HorizontalAlignment','right', ...
					 'Fontsize', 8, ...
					 'Visible', 'off', ...
					 'Position', ...
					 [455+40*boxX 230+30*boxY 20 10]);
				       

    handles.groupBox(iGB) = uicontrol('Style','checkbox', ...
				      'BackgroundColor', 0.8*[1 1 1], ...
				      'Visible', 'off', ...
				      'Value', 0, ...
				      'Position', ...
				      [475+40*boxX 225+30*boxY 20 20], ...
				      'Callback', @selectGroups);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  handles.XlimMinLabel = uicontrol('Style', 'text', ...
				   'String', 'Xmin', ...
				   'Fontsize',7, ...
				   'HorizontalAlignment','center', ...
				   'BackgroundColor', 0.8*[1 1 1], ...
				   'Position', [450 200 40 15]);

  handles.XlimMaxLabel = uicontrol('Style', 'text', ...
				   'String', 'Xmax', ...
				   'Fontsize',7, ...
				   'HorizontalAlignment','center', ...
				   'BackgroundColor', 0.8*[1 1 1], ...
				   'Position', [490 200 40 15]);

  handles.YlimMinLabel = uicontrol('Style', 'text', ...
				   'String', 'Ymin', ...
				   'Fontsize',7, ...
				   'HorizontalAlignment','center', ...
				   'BackgroundColor', 0.8*[1 1 1], ...
				   'Position', [530 200 40 15]);

  handles.YlimMaxLabel = uicontrol('Style', 'text', ...
				   'String', 'Ymax', ...
				   'Fontsize',7, ...
				   'HorizontalAlignment','center', ...
				   'BackgroundColor', 0.8*[1 1 1], ...
				   'Position', [570 200 40 15]);


  handles.XlimMin = uicontrol('Style','edit', ...
			      'String', '', ...
			      'HorizontalAlignment','left', ...
			      'Interruptible', 'off', ...
			      'Position', [450 175 40 30]);

  handles.XlimMax = uicontrol('Style','edit', ...
			      'String', '', ...
			      'HorizontalAlignment','left', ...
			      'Interruptible', 'off', ...
			      'Position', [490 175 40 30]);

  handles.YlimMin = uicontrol('Style','edit', ...
			      'String', '', ...
			      'HorizontalAlignment','left', ...
			      'Interruptible', 'off', ...
			      'Position', [530 175 40 30]);

  handles.YlimMax = uicontrol('Style','edit', ...
			      'String', '', ...
			      'HorizontalAlignment','left', ...
			      'Interruptible', 'off', ...
			      'Position', [570 175 40 30]);

  handles.setAxis = uicontrol('Style','pushbutton', ...
			      'String','Change axis', ...
			      'Interruptible', 'off', ...
			      'Position',[450,130,100,30], ...
			      'Callback', @setAxis);

  handles.barGraphDotsText = uicontrol('Style', 'text', ...
				       'String', 'Bar graph dots:', ...
				       'Fontsize', 12, ...
				       'BackgroundColor', 0.8*[1 1 1], ...
				       'HorizontalAlignment','left', ...
				       'Position', [450 100 110 20]);

  handles.barGraphDots      = uicontrol('Style', 'checkbox', ...
					'Value', 0, 'Max', 1, ...
					'BackgroundColor', 0.8*[1 1 1], ...
					'Interruptible', 'off', ...
					'Enable', 'off', ...
					'Position', [560 100 20 20], ...
					'Callback', @setBarGraphDots);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  handles.menuFile = uimenu(handles.fig,'Label','File');
  handles.menuItemLoad = uimenu(handles.menuFile, ...
				'Label','Load data', ...
				'Interruptible','off', ...
				'Callback', @loadData);

  handles.menuItemLoad = uimenu(handles.menuFile, ...
				'Label','Export data', ...
				'Interruptible','off', ...
				'Callback', @exportData);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  handles.menuDebug = uimenu(handles.fig,'Label','Debug');
  handles.menuItemDebug = uimenu(handles.menuDebug, ...
				 'Label','Keyboard', ...
				 'Interruptible', 'off', ...
				 'Callback', @runDebug);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function runDebug(source, event)
    disp('Type return to exit debug mode')
    keyboard
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  plotInfo = struct('PlotName', [], ...
		    'PlotType', [], ...
		    'DataStr', [], ...
		    'Xlabel', [], ...
		    'Ylabel', [], ...
		    'FigName', [], ...
		    'barGraphDots', 0, ...
		    'Xlim', [], ...
		    'Ylim', []);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setPlotInfo()

    % "Basics"
    % 2) synapse number: Bar graph
    % y-axis: Mean synapse number
    % x-axis: (empty)

    plotInfo(1).Name    = 'Synapse number';
    plotInfo(1).Type    = 'bar';
    plotInfo(1).DataStr = 'nSynapses';
    plotInfo(1).Xlabel  = '';
    plotInfo(1).Ylabel  = 'Mean synapse number';
    plotInfo(1).FigName = 'synapse_number';
    plotInfo(1).barGraphDots = 1;


    % 3)Total dendrite length: Bar graph
    % y-axis: Mean total dendrite length (um)
    % x-axis: (empty)
  
    % Apologies, the totDendLength is NOT in SI units... sorry.
    % Sabine comment: it is in micrometers.
    % Update: is is in meters...
    plotInfo(2).Name    = 'Dendritic length';
    plotInfo(2).Type    = 'bar';
    plotInfo(2).DataStr = 'totDendLength*1e6';
    plotInfo(2).Xlabel  = '';
    plotInfo(2).Ylabel  = 'Mean total dendrite length (\mum)';
    plotInfo(2).FigName = 'total_dendrite_length';
    plotInfo(2).barGraphDots = 1;
  
   
    % 4) synapse number/dendrite length: Bar graph
    % y-axis: Mean synapse number per um dendrite
    % x-axis: (empty)
  
    plotInfo(3).Name    = 'Synapse number/Dendritic length';
    plotInfo(3).Type    = 'bar';
    plotInfo(3).DataStr = {'nSynapses', 'totDendLength*1e6'};
    plotInfo(3).Xlabel  = '';
    plotInfo(3).Ylabel  = 'Mean synapses number per \mum dendrite';
    plotInfo(3).FigName = 'synapse_number_per__dendrite_length';
    plotInfo(3).barGraphDots = 1;  
  
    % 5) synapse size: Bar graph
    % y-axis: Mean synapse area (um2)
    % x-axis: (empty)
  
    plotInfo(4).Name    = 'Mean synapse area';
    plotInfo(4).Type    = 'bar';
    plotInfo(4).DataStr = 'synapseArea'; % This is already in micrometers...
    plotInfo(4).Xlabel  = '';
    plotInfo(4).Ylabel  = 'Mean synapse area (\mum^2)';
    plotInfo(4).FigName = 'synapse_area';
    plotInfo(4).barGraphDots = 1;  
  

    % 6) synaptic morphology intensity: Bar graph
    % y-axis: Mean synaptic morphology intensity (a.u.)
    % x-axis: (empty)
  
    plotInfo(5).Name    = sprintf('Synaptic %s intensity', morphLabel);
    plotInfo(5).Type    = 'bar';
    plotInfo(5).DataStr = 'synapseIntensityMorphMean';
    plotInfo(5).Xlabel  = '';
    plotInfo(5).Ylabel  = sprintf('Mean synaptic %s intensity (a.u.)',morphLabel);
    plotInfo(5).FigName = sprintf('synaptic_%s_intensity', morphLabel);
    plotInfo(5).barGraphDots = 1;
  
  
    % 7) synapse intensity: Bar graph
    % y-axis: Mean synapse intensity (a.u.)
    % x-axis: (empty)
  
    plotInfo(6).Name    = sprintf('Synaptic %s intensity', synLabel);
    plotInfo(6).Type    = 'bar';
    plotInfo(6).DataStr = 'synapseIntensitySynMean';
    plotInfo(6).Xlabel  = '';
    plotInfo(6).Ylabel  = sprintf('Mean synaptic %s intensity (a.u.)', synLabel);
    plotInfo(6).FigName = sprintf('synaptic_%s_intensity', synLabel);
    plotInfo(6).barGraphDots = 1;  
  
    % 8) synaptic X intensity: Bar graph
    % y-axis: Mean synaptic X intensity (a.u.)
    % x-axis: (empty)
  
    plotInfo(7).Name    = sprintf('Synaptic %s intensity', XLabel);
    plotInfo(7).Type    = 'bar';
    plotInfo(7).DataStr = 'synapseIntensityXMean';
    plotInfo(7).Xlabel  = '';
    plotInfo(7).Ylabel  = sprintf('Mean synaptic %s intensity (a.u.)',XLabel);
    plotInfo(7).FigName = sprintf('synaptic_%s_intensity', XLabel);
    plotInfo(7).barGraphDots = 1;    
   
    % 9) synaptic synapse/morphology intensity: Bar graph
    % y-axis: Mean synaptic synapse/morphology intensity (a.u.)
    % x-axis: (empty)

    plotInfo(8).Name    = sprintf('Synaptic %s/%s intensity', synLabel,morphLabel);
    plotInfo(8).Type    = 'bar';
    plotInfo(8).DataStr = {'synapseIntensitySynMean','synapseIntensityMorphMean'};
    plotInfo(8).Xlabel  = '';
    plotInfo(8).Ylabel  = sprintf('Mean synaptic %s/%s intensity (a.u.)', ...
				  synLabel, morphLabel);
    plotInfo(8).FigName = sprintf('synaptic_ratio_of_%s_and_%s_intensity', ...
				  synLabel, morphLabel);
    plotInfo(8).barGraphDots = 1;


    % 10) synaptic X/morphology intensity: Bar graph
    % y-axis: Mean synapse intensity (a.u.)
    % x-axis: (empty)

    plotInfo(9).Name    = sprintf('Synaptic %s/%s intensity', XLabel,morphLabel);
    plotInfo(9).Type    = 'bar';
    plotInfo(9).DataStr = {'synapseIntensityXMean','synapseIntensityMorphMean'};
    plotInfo(9).Xlabel  = '';
    plotInfo(9).Ylabel  = sprintf('Mean synaptic %s/%s intensity (a.u.)', ...
				  XLabel, morphLabel);
    plotInfo(9).FigName = sprintf('synaptic_ratio_of_%s_and_%s_intensity', ...
				  XLabel, morphLabel);
    plotInfo(9).barGraphDots = 1;

 
    % 11) synaptic X/synapse intensity: Bar graph
    % y-axis: Mean synaptic X/synapse intensity (a.u.)
    % x-axis: (empty)

    plotInfo(10).Name    = sprintf('Synaptic %s/%s intensity', XLabel,synLabel);
    plotInfo(10).Type    = 'bar';
    plotInfo(10).DataStr = {'synapseIntensityXMean','synapseIntensitySynMean'};
    plotInfo(10).Xlabel  = '';
    plotInfo(10).Ylabel  = sprintf('Mean synaptic %s/%s intensity (a.u.)', ...
				   XLabel, synLabel);
    plotInfo(10).FigName = sprintf('synaptic_ratio_of_%s_and_%s_intensity', ...
				   XLabel, synLabel);
    plotInfo(10).barGraphDots = 1;

 
    % 13) somatic morphology intensity: Bar graph
    % y-axis: Mean somatic morphology intensity (a.u.)
    % x-axis: (empty)

    plotInfo(11).Name    = sprintf('Somatic %s intensity', morphLabel);
    plotInfo(11).Type    = 'bar';
    plotInfo(11).DataStr = 'somaMeasureMorph';
    plotInfo(11).Xlabel  = '';
    plotInfo(11).Ylabel  = sprintf('Mean somatic %s intensity (a.u.)', morphLabel);
    plotInfo(11).FigName = sprintf('somatic_%s_intensity',morphLabel);
    plotInfo(11).barGraphDots = 1;


    % 14) somatic synapse intensity: Bar graph
    % y-axis: Mean somatic synapse intensity (a.u.)
    % x-axis: (empty)

    plotInfo(12).Name    = sprintf('Somatic %s intensity', synLabel);
    plotInfo(12).Type    = 'bar';
    plotInfo(12).DataStr = 'somaMeasureSyn';
    plotInfo(12).Xlabel  = '';
    plotInfo(12).Ylabel  = sprintf('Mean somatic %s intensity (a.u.)', synLabel);
    plotInfo(12).FigName = sprintf('somatic_%s_intensity',synLabel);
    plotInfo(12).barGraphDots = 1;


    % 15) somatic X intensity: Bar graph
    % y-axis: Mean somatic X intensity (a.u.)
    % x-axis: (empty)

    plotInfo(13).Name    = sprintf('Somatic %s intensity', XLabel);
    plotInfo(13).Type    = 'bar';
    plotInfo(13).DataStr = 'somaMeasureX';
    plotInfo(13).Xlabel  = '';
    plotInfo(13).Ylabel  = sprintf('Mean somatic %s intensity (a.u.)',XLabel);
    plotInfo(13).FigName = sprintf('somatic_%s_intensity',XLabel);
    plotInfo(13).barGraphDots = 1;

 
    % 16) soma area: Bar graph
    % y-axis: Mean soma area (um2)
    % x-axis: (empty)

    plotInfo(14).Name    = 'Soma area';
    plotInfo(14).Type    = 'bar';
    plotInfo(14).DataStr = 'somaArea';
    plotInfo(14).Xlabel  = '';
    plotInfo(14).Ylabel  = 'Mean soma area (\mum^2)';
    plotInfo(14).FigName = 'soma_area';
    plotInfo(14).barGraphDots = 1;


    % 17) major soma axis: Bar graph
    % y-axis: Mean major soma axis (um2)
    % x-axis: (empty)

    plotInfo(15).Name    = 'Soma major axis';
    plotInfo(15).Type    = 'bar';
    plotInfo(15).DataStr = 'somaMajorAxisLength';
    plotInfo(15).Xlabel  = '';
    plotInfo(15).Ylabel  = 'Soma major axis (\mum)';
    plotInfo(15).FigName = 'soma_major_axis';
    plotInfo(15).barGraphDots = 1;

 
    % 18) minor soma axis: Bar graph
    % y-axis: Mean minor soma axis (um2)
    % x-axis: (empty)

    plotInfo(16).Name    = 'Soma minor axis';
    plotInfo(16).Type    = 'bar';
    plotInfo(16).DataStr = 'somaMinorAxisLength';
    plotInfo(16).Xlabel  = '';
    plotInfo(16).Ylabel  = 'Soma minor axis (\mum)';
    plotInfo(16).FigName = 'soma_minor_axis';
    plotInfo(16).barGraphDots = 1;


    % 19) synaptic/somatic morphology intensity: Bar graph
    % y-axis: Mean synaptic/somatic morphology intensity (a.u.)
    % x-axis: (empty)

    plotInfo(17).Name    = sprintf('Synaptic/somatic %s intensity', morphLabel);
    plotInfo(17).Type    = 'bar';
    plotInfo(17).DataStr = {'synapseIntensityMorphMean', 'somaMeasureMorph'};
    plotInfo(17).Xlabel  = '';
    plotInfo(17).Ylabel  = sprintf('Mean synaptic/somatic %s intensity', morphLabel);
    plotInfo(17).FigName = sprintf('synaptic_to_somatic_%s_intensity_ratio',morphLabel);
    plotInfo(17).barGraphDots = 1;


    % 20) synaptic/somatic synapse intensity: Bar graph
    % y-axis: Mean synaptic/somatic synapse intensity (a.u.)
    % x-axis: (empty)

    plotInfo(18).Name    = sprintf('Synaptic/somatic %s intensity', synLabel);
    plotInfo(18).Type    = 'bar';
    plotInfo(18).DataStr = {'synapseIntensitySynMean', 'somaMeasureSyn'};
    plotInfo(18).Xlabel  = '';
    plotInfo(18).Ylabel  = sprintf('Mean synaptic/somatic %s intensity', synLabel);
    plotInfo(18).FigName = sprintf('synaptic_to_somatic_%s_intensity_ratio', synLabel);
    plotInfo(18).barGraphDots = 1;


    % 21) synaptic/somatic  X intensity: Bar graph
    % y-axis: Mean synaptic/somatic X intensity (a.u.)
    % x-axis: (empty)

    plotInfo(19).Name    = sprintf('Synaptic/somatic %s intensity', XLabel);
    plotInfo(19).Type    = 'bar';
    plotInfo(19).DataStr = {'synapseIntensityXMean', 'somaMeasureX'};
    plotInfo(19).Xlabel  = '';
    plotInfo(19).Ylabel  = sprintf('Mean synaptic/somatic %s intensity', XLabel);
    plotInfo(19).FigName = sprintf('synaptic_to_somatic_%s_intensity_ratio', XLabel);
    plotInfo(19).barGraphDots = 1;


    % 22) synaptic/somatic (synapse/morphology) intensity: Bar graph
    % y-axis: Mean synaptic/somatic (synapse/morphology) intensity (a.u.)
    % x-axis: (empty)

    plotInfo(20).Name    = sprintf('Synaptic %s/somatic %s intensity', ...
				   synLabel, morphLabel);
    plotInfo(20).Type    = 'bar';
    plotInfo(20).DataStr = {'synapseIntensitySynMean', 'somaMeasureMorph'};
    plotInfo(20).Xlabel  = '';
    plotInfo(20).Ylabel  = sprintf('Mean synaptic %s/somatic %s intensity', ...
				   synLabel, morphLabel);
    plotInfo(20).FigName = sprintf('synaptic_%s_to_somatic_%s_intensity_ratio', ...
				   synLabel, morphLabel);
    plotInfo(20).barGraphDots = 1;


    % 23) synaptic/somatic (X/morphology) intensity: Bar graph
    % y-axis: Mean synaptic/somatic (X/morphology) intensity (a.u.)
    % x-axis: (empty)

    plotInfo(21).Name    = sprintf('Synaptic %s/somatic %s intensity', ...
				   XLabel, morphLabel);
    plotInfo(21).Type    = 'bar';
    plotInfo(21).DataStr = {'synapseIntensityXMean', 'somaMeasureMorph'};
    plotInfo(21).Xlabel  = '';
    plotInfo(21).Ylabel  = sprintf('Mean synaptic %s/somatic %s intensity', ...
				   XLabel, morphLabel);
    plotInfo(21).FigName = sprintf('synaptic_%s_to_somatic_%s_intensity_ratio', ...
				   XLabel, morphLabel);
    plotInfo(21).barGraphDots = 1;


    % 24) synaptic/somatic (X/synapse) intensity: Bar graph
    % y-axis: Mean synaptic/somatic (X/synapse) intensity (a.u.)
    % x-axis: (empty)

    plotInfo(22).Name    = sprintf('Synaptic %s/somatic %s intensity', ...
				   XLabel, synLabel);
    plotInfo(22).Type    = 'bar';
    plotInfo(22).DataStr = {'synapseIntensityXMean', 'somaMeasureSyn'};
    plotInfo(22).Xlabel  = '';
    plotInfo(22).Ylabel  = sprintf('Mean synaptic %s/somatic %s intensity', ...
				   XLabel, synLabel);
    plotInfo(22).FigName = sprintf('synaptic_%s_to_somatic_%s_intensity_ratio', ...
				   XLabel, synLabel);
    plotInfo(22).barGraphDots = 1;


    %"M"
    % n) cumulative histogram of synaptic morphology intensity: cumulative histogram (see supplementary Fig.4 SynD paper)
    % y-axis: cumulative probability (%)
    % x-axis: distance from soma (um)

    plotInfo(23).Name    = sprintf('Cumulative histogram of synaptic %s intensity', ...
				   morphLabel);
    plotInfo(23).Type    = 'cumhist';
    plotInfo(23).DataStr = 'morphHist';
    plotInfo(23).Xlabel  = 'Fluorescent intensity (a.u.)';
    plotInfo(23).Ylabel  = 'Cumulative probability (%)';
    plotInfo(23).FigName = sprintf('cumulative_histogram_synaptic_%s_intensity', ...
				   morphLabel);

    % "S"
    % n) cumulative histogram of synapse intensity: cumulative histogram
    % y-axis: cumulative probability (%)
    % x-axis: distance from soma (um)

    plotInfo(24).Name    = sprintf('Cumulative histogram of synaptic %s intensity', ...
				   synLabel);
    plotInfo(24).Type    = 'cumhist';
    plotInfo(24).DataStr = 'synHist';
    plotInfo(24).Xlabel  = 'Fluorescent intensity (a.u.)';
    plotInfo(24).Ylabel  = 'Cumulative probability (%)';
    plotInfo(24).FigName = sprintf('cumulative_histogram_synaptic_%s_intensity', ...
				   synLabel);

    % "X"
    % n) cumulative histogram of synaptic X intensity: cumulative histogram
    % y-axis: cumulative probability (%)
    % x-axis: distance from soma (um)

    plotInfo(25).Name    = sprintf('Cumulative histogram of synaptic %s intensity', ...
				   XLabel);
    plotInfo(25).Type    = 'cumhist';
    plotInfo(25).DataStr = 'XHist';
    plotInfo(25).Xlabel  = 'Fluorescent intensity (a.u.)';
    plotInfo(25).Ylabel  = 'Cumulative probability (%)';
    plotInfo(25).FigName = sprintf('cumulative_histogram_synaptic_%s_intensity', ...
				   XLabel);
 
    % "Sholl dend"
    % n) Sholl analysis: Histogram (connected points, see Fig3B SynD paper)
    % y-axis: number of branches
    % x-axis: distance from soma (um)

    plotInfo(26).Name    = 'Sholl analysis';
    plotInfo(26).Type    = 'hist';
    plotInfo(26).DataStr = 'DendHist';
    plotInfo(26).Xlabel  = 'Distance from soma (\mum)';
    plotInfo(26).Ylabel  = 'Mean number of branches';
    plotInfo(26).FigName = 'sholl_dendrite_analysis';


    % "Sholl synd"
    % n) Synapse localization: Histogram (as bars, see Fig3C SynD paper)
    % y-axis: number of synapses
    % x-axis: distance from soma (um)

    plotInfo(27).Name    = 'Synapse localisation';
    plotInfo(27).Type    = 'hist';
    plotInfo(27).DataStr = 'SynHist';
    plotInfo(27).Xlabel  = 'Distance from soma (\mum)';
    plotInfo(27).Ylabel  = 'Mean number of synapses';
    plotInfo(27).FigName = 'sholl_synapse_analysis';


    % "Sholl M"
    % n) Synaptic morphology intensity distribution in Sholl analysis: Histogram (as bars, see Fig3D SynD paper)
    % y-axis: synaptic morphology intensity (a.u.)
    % x-axis: distance from soma (um)

    plotInfo(28).Name    = sprintf('Synaptic %s intensity distribution in Sholl analysis', ...
				   morphLabel);
    plotInfo(28).Type    = 'hist';
    plotInfo(28).DataStr = 'IntMorphMean';
    plotInfo(28).Xlabel  = 'Distance from soma (\mum)';
    plotInfo(28).Ylabel  = sprintf('Mean synaptic %s intensity (a.u.)',morphLabel);
    plotInfo(28).FigName = sprintf('sholl_synaptic_%s_analysis',morphLabel);


    % "Sholl S"
    % n) Synapse intensity distribution in Sholl analysis: Histogram (as bars, see Fig3D SynD paper)
    % y-axis: synapse intensity (a.u.)
    % x-axis: distance from soma (um)

    plotInfo(29).Name    = sprintf('Synaptic %s intensity distribution in Sholl analysis', ...
				   synLabel);
    plotInfo(29).Type    = 'hist';
    plotInfo(29).DataStr = 'IntSynMean';
    plotInfo(29).Xlabel  = 'Distance from soma (\mum)';
    plotInfo(29).Ylabel  = sprintf('Mean synaptic %s intensity (a.u.)',synLabel);
    plotInfo(29).FigName = sprintf('sholl_synaptic_%s_analysis',synLabel);


    % "Sholl X"
    % n) Synaptic X intensity distribution in Sholl analysis: Histogram (as bars, see Fig3D SynD paper)
    % y-axis: synaptic X intensity (a.u.)
    % x-axis: distance from soma (um)

    plotInfo(30).Name    = sprintf('Synaptic %s intensity distribution in Sholl analysis', ...
				   XLabel);
    plotInfo(30).Type    = 'hist';
    plotInfo(30).DataStr = 'IntXMean';
    plotInfo(30).Xlabel  = 'Distance from soma (\mum)';
    plotInfo(30).Ylabel  = sprintf('Mean synaptic %s intensity (a.u.)',XLabel);
    plotInfo(30).FigName = sprintf('sholl_synaptic_%s_analysis',XLabel);


    % "Sholl r(S,M)"
    % n) Synaptic synapse/morphology intensity distribution in Sholl analysis: Histogram (as bars, see Fig3C SynD paper)
    % y-axis: synaptic synapse/morphology intensity (a.u.)
    % x-axis: distance from soma (um)

    plotInfo(31).Name    = sprintf('Synaptic %s/%s intensity distribution in Sholl analysis', ...
				   synLabel, morphLabel);
    plotInfo(31).Type    = 'hist';
    plotInfo(31).DataStr = {'IntSynMean', 'IntMorphMean'};
    plotInfo(31).Xlabel  = 'Distance from soma (\mum)';
    plotInfo(31).Ylabel  = sprintf('Mean synaptic %s/%s intensity (a.u.)', ...
				   synLabel, morphLabel);
    plotInfo(31).FigName = sprintf('sholl_synaptic_ratio_of_%s_%s_analysis',...
				   synLabel, morphLabel);


    % "Sholl r(X,M)"
    % n) Synaptic X/morphology intensity distribution in Sholl analysis: Histogram (as bars, see Fig3D SynD paper)
    % y-axis: synaptic X/morphology intensity (a.u.)
    % x-axis: distance from soma (um)

    plotInfo(32).Name    = sprintf('Synaptic %s/%s intensity distribution in Sholl analysis', ...
				   XLabel, morphLabel);
    plotInfo(32).Type    = 'hist';
    plotInfo(32).DataStr = {'IntXMean', 'IntMorphMean'};
    plotInfo(32).Xlabel  = 'Distance from soma (\mum)';
    plotInfo(32).Ylabel  = sprintf('Mean synaptic %s/%s intensity (a.u.)', ...
				   XLabel, morphLabel);
    plotInfo(32).FigName = sprintf('sholl_synaptic_ratio_of_%s_%s_analysis',...
				   XLabel, morphLabel);


    % "Sholl r(X,S)"
    % n) Synaptic X/synapse intensity distribution in Sholl analysis: Histogram (as bars, see Fig3D SynD paper)
    % y-axis: synaptic X/synapse intensity (a.u.)
    % x-axis: distance from soma (um)
   
    plotInfo(33).Name    = sprintf('Synaptic %s/%s intensity distribution in Sholl analysis', ...
				   XLabel, synLabel);
    plotInfo(33).Type    = 'hist';
    plotInfo(33).DataStr = {'IntXMean', 'IntSynMean'};
    plotInfo(33).Xlabel  = 'Distance from soma (\mum)';
    plotInfo(33).Ylabel  = sprintf('Mean synaptic %s/%s intensity (a.u.)', ...
				   XLabel, synLabel);
    plotInfo(33).FigName = sprintf('sholl_synaptic_ratio_of_%s_%s_analysis',...
				   XLabel, synLabel);


    % "Profile M"
    % n) (Juxta)synaptic morphology intensity profile: Histogram (connected points, basically like you do during extract)
    % y-axis: (juxta)synaptic morphology intensity (a.u.)
    % x-axis: distance from synapse centre (um)

    plotInfo(34).Name    = sprintf('(Juxta)synaptic %s intensity profile', ...
				   morphLabel);
    plotInfo(34).Type    = 'profile';
    plotInfo(34).DataStr = 'meanSynapseProfileMorph';
    plotInfo(34).Xlabel  = 'Distance from synapse centre (\mum)';
    plotInfo(34).Ylabel  = sprintf('(Juxta)synaptic %s intensity (a.u.)',morphLabel);
    plotInfo(34).FigName = sprintf('juxta_synaptic_%s_profile',morphLabel);


    % "Profile S"
    % n) (Juxta)synaptic synapse intensity profile: Histogram (connected points, basically like you do during extract)
    % y-axis: (juxta)synaptic synapse intensity (a.u.)
    % x-axis: distance from synapse centre (um)

    plotInfo(35).Name    = sprintf('(Juxta)synaptic %s intensity profile', ...
				   synLabel);
    plotInfo(35).Type    = 'profile';
    plotInfo(35).DataStr = 'meanSynapseProfileSyn';
    plotInfo(35).Xlabel  = 'Distance from synapse centre (\mum)';
    plotInfo(35).Ylabel  = sprintf('(Juxta)synaptic %s intensity (a.u.)',synLabel);
    plotInfo(35).FigName = sprintf('juxta_synaptic_%s_profile',synLabel);


    % "Profile X"
    % n) (Juxta)synaptic X intensity profile: Histogram (connected points, basically like you do during extract)
    % y-axis: (juxta)synaptic X intensity (a.u.)
    % x-axis: distance from synapse centre (um)
    
    plotInfo(36).Name    = sprintf('(Juxta)synaptic %s intensity profile', ...
				   XLabel);
    plotInfo(36).Type    = 'profile';
    plotInfo(36).DataStr = 'meanSynapseProfileX';
    plotInfo(36).Xlabel  = 'Distance from synapse centre (\mum)';
    plotInfo(36).Ylabel  = sprintf('(Juxta)synaptic %s intensity (a.u.)',XLabel);
    plotInfo(36).FigName = sprintf('juxta_synaptic_%s_profile',XLabel);


    % "Profile r(S,M)"
    % n) (Juxta)synaptic synapse/morphology intensity profile: Histogram (connected points, basically like you do during extract)
    % y-axis: (juxta)synaptic synapse/morphology intensity (a.u.)
    % x-axis: distance from synapse centre (um)

    plotInfo(37).Name    = sprintf('(Juxta)synaptic %s/%s intensity profile', ...
				   synLabel, morphLabel);
    plotInfo(37).Type    = 'profile';
    plotInfo(37).DataStr = {'meanSynapseProfileSyn', 'meanSynapseProfileMorph'};
    plotInfo(37).Xlabel  = 'Distance from synapse centre (\mum)';
    plotInfo(37).Ylabel  = sprintf('(Juxta)synaptic %s/%s intensity (a.u.)', ...
				   synLabel,morphLabel);
    plotInfo(37).FigName = sprintf('juxta_synaptic_ratio_%s_over_%s_profile', ...
				   synLabel, morphLabel);


    % "Profile r(X,M)"
    % n) (Juxta)synaptic X/morphology intensity profile: Histogram (connected points, basically like you do during extract)
    % y-axis: (juxta)synaptic X/morphology intensity (a.u.)
    % x-axis: distance from synapse centre (um)

    plotInfo(38).Name    = sprintf('(Juxta)synaptic %s/%s intensity profile', ...
				   XLabel, morphLabel);
    plotInfo(38).Type    = 'profile';
    plotInfo(38).DataStr = {'meanSynapseProfileX', 'meanSynapseProfileMorph'};
    plotInfo(38).Xlabel  = 'Distance from synapse centre (\mum)';
    plotInfo(38).Ylabel  = sprintf('(Juxta)synaptic %s/%s intensity (a.u.)', ...
				   XLabel,morphLabel);
    plotInfo(38).FigName = sprintf('juxta_synaptic_ratio_%s_over_%s_profile', ...
				   XLabel, morphLabel);



    % "Profile r(X,S)"
    % n) (Juxta)synaptic X/synapse intensity profile: Histogram (connected points, basically like you do during extract)
    % y-axis: (juxta)synaptic X/synapse intensity (a.u.)
    % x-axis: distance from synapse centre (um)
    
    plotInfo(39).Name    = sprintf('(Juxta)synaptic %s/%s intensity profile', ...
				   XLabel, synLabel);
    plotInfo(39).Type    = 'profile';
    plotInfo(39).DataStr = {'meanSynapseProfileX', 'meanSynapseProfileSyn'};
    plotInfo(39).Xlabel  = 'Distance from synapse centre (\mum)';
    plotInfo(39).Ylabel  = sprintf('(Juxta)synaptic %s/%s intensity (a.u.)', ...
				   XLabel,synLabel);
    plotInfo(39).FigName = sprintf('juxta_synaptic_ratio_%s_over_%s_profile', ...
				   XLabel, synLabel);


    plotInfo(40).Name    = sprintf('Neurite %s intensity profile', ...
				   morphLabel);
    plotInfo(40).Type    = 'gradient';
    plotInfo(40).DataStr = 'morphGradientMean';
    plotInfo(40).Xlabel  = 'Distance from soma (\mum)';
    plotInfo(40).Ylabel  = sprintf('%s intensity (a.u.)', ...
				   morphLabel);
    plotInfo(40).FigName = sprintf('neurite_%s_gradient', ...
				   morphLabel);

    plotInfo(41).Name    = sprintf('Neurite %s intensity profile', ...
				   synLabel);
    plotInfo(41).Type    = 'gradient';
    plotInfo(41).DataStr = 'synGradientMean';
    plotInfo(41).Xlabel  = 'Distance from soma (\mum)';
    plotInfo(41).Ylabel  = sprintf('%s intensity (a.u.)', ...
				   synLabel);
    plotInfo(41).FigName = sprintf('neurite_%s_gradient', ...
				   synLabel);

    plotInfo(42).Name    = sprintf('Neurite %s intensity profile', ...
				   XLabel);
    plotInfo(42).Type    = 'gradient';
    plotInfo(42).DataStr = 'XGradientMean';
    plotInfo(42).Xlabel  = 'Distance from soma (\mum)';
    plotInfo(42).Ylabel  = sprintf('%s intensity (a.u.)', ...
				   XLabel);
    plotInfo(42).FigName = sprintf('neurite_%s_gradient', ...
				   XLabel);

    plotInfo(43).Name    = sprintf('Neurite %s/%s intensity profile', ...
				   synLabel, morphLabel);
    plotInfo(43).Type    = 'gradient';
    plotInfo(43).DataStr = {'synGradientMean','morphGradientMean'};
    plotInfo(43).Xlabel  = 'Distance from soma (\mum)';
    plotInfo(43).Ylabel  = sprintf('%s/%s intensity (a.u.)', ...
				   synLabel,morphLabel);
    plotInfo(43).FigName = sprintf('neurite_ratio_%s_over_%s_gradient', ...
				   synLabel,morphLabel);

    plotInfo(44).Name    = sprintf('Neurite %s/%s intensity profile', ...
				   XLabel,morphLabel);
    plotInfo(44).Type    = 'gradient';
    plotInfo(44).DataStr = {'XGradientMean','morphGradientMean'};
    plotInfo(44).Xlabel  = 'Distance from soma (\mum)';
    plotInfo(44).Ylabel  = sprintf('%s/%s intensity (a.u.)', ...
				   XLabel,morphLabel);
    plotInfo(44).FigName = sprintf('neurite_ratio_%s_over_%s_gradient', ...
				   XLabel,morphLabel);



  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  setPlotTypeList(); 

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function loadData(source, event)
    oldDir = pwd;

    if(~isempty(loadPath))
      try
        cd(loadPath);    
      catch
        fprintf('Unable to change to %s\n', loadPath)
      end
    end

    fileType = {'*-save.mat', 'Saved neuron data'; ...
		'*.mat', 'Matlab file'; ...
		'*', 'All files'};

    [dataFile, dataPath, filterIndex] = ...
      uigetfile(fileType, ...
		'Select neuron image', ...
		'MultiSelect', 'on');

    loadPath = dataPath;

    if(~iscell(dataFile) & dataFile == 0)
      % User pressed cancel
      return;
    end

    if(~iscell(dataFile))
      dataFile = {dataFile};
    end

    for i = 1:length(dataFile)
      fprintf('Loading %s\n', dataFile{i})
      tmp = load(strcat(dataPath,dataFile{i}));

      try 
        tmp.old.neuriteMask;
      catch
        % Some old format has wrong name.
        tmp.old.neuriteMask = tmp.old.neuriteMaskIdx;
        tmp.old.somaMask = tmp.old.somaMaskIdx;
        tmp.old.neuriteMaskIdx = [];
        tmp.old.somaMaskIdx = [];
      end  
    
      try
        tmp.old.somaMeasureMorph;
      catch
        tmp.old.somaMeasureMorph = [];
      end

      try
        tmp.old.gradientEdges;
      catch
        tmp.old.gradientEdges = [];
        tmp.old.morphGradientMean = [];
        tmp.old.morphGradientStd = [];
        tmp.old.morphGradientSEM = [];
        tmp.old.synGradientMean = [];
        tmp.old.synGradientStd = [];
        tmp.old.synGradientSEM = [];
        tmp.old.XGradientMean = [];
        tmp.old.XGradientStd = [];
        tmp.old.XGradientSEM = [];
      end

      try
        tmp.old.meanSynapseProfileDist;
      catch
        tmp.old.meanSynapseProfileDist = [];
        tmp.old.meanSynapseProfileMorph = [];
        tmp.old.meanSynapseProfileSyn = [];
        tmp.old.meanSynapseProfileX = [];
      end

      % Check if it is old format, then guess info
      try
        tmp.old.detection;
      catch
        disp('No detection info, assuming neurite padding 3.')
        tmp.old.detection.neuritePadding = 3;
        tmp.old.dataImg = 1;
        tmp.old.detection.morphThreshold = NaN;
      end

      try
        % Are absolute intensity thresholds saved?
        tmp.old.somaIntensityThreshold;
      catch
        disp('Old format, some detection threshold information missing')
        tmp.old.somaIntensityThreshold = tmp.old.detection.morphThreshold;
        tmp.old.synapseIntensityThreshold = NaN;
      end

      try
        tmp.old.version = str2num(tmp.old.version);
        if(isempty(tmp.old.version))
	  disp('Analysed with an unknown version of SynD.')
	  tmp.old.version = NaN;
        else
          fprintf('Analysed with SynD revision %d.\n', tmp.old.version)
        end
      catch
 	disp('Analysed with old version, pre-revision 465.')
        tmp.old.version = 0;
      end

      CC = bwconncomp(tmp.old.somaMask);
      tmp.old.nSoma = CC.NumObjects; 

      if(conserveMemory)
        % Clear variables not used
        tmp.old.origImage = [];

        % tmp.old.somaMask = logical(tmp.old.somaMask);
        % tmp.old.neuriteMask = logical(tmp.old.neuriteMask);
        % tmp.old.synapseMask = logical(tmp.old.synapseMask);

        % Clear these ones also
        tmp.old.image = [];
        tmp.old.somaMask = [];
        tmp.old.neuriteMask = [];
        tmp.old.synapseMask = [];

      end

      % Pre-calculating, to make it easier for plots
      tmp.old.nSynapses = length(tmp.old.synapseCenter);

      if(length(tmp.old.synapseCenter) ...
	 ~= length(tmp.old.synapseIntensitySynMean))
        disp('This data file is corrupted. Please reanalyse. File not loaded.')
        beep
        continue
      end

      try
        if(isempty(data))
          data = tmp.old;
          groupIdx = 0;
        else
          data(end+1) = tmp.old;
          groupIdx(end+1) = 0;
        end
      catch e
        fprintf('Unable to import file %s, saved with different versions?\n',dataFile{i});
        % keyboard
      end  

    end

    updateListboxNames();

    cd(oldDir);

    % Check to make sure the user are not mixing

    prevFlag = data(1).detection.excludeSomaSynapses;
    allEqual = 1;

    for i = 2:length(data)
      if(prevFlag ~= data(i).detection.excludeSomaSynapses)
        allEqual = 0;
      end
    end

    if(~allEqual)
      disp('Inconsistent data, some include soma synapses, others do not.')
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function exportData(source, event)

    if(isempty(group(1).Members))
      if(isempty(data))
        disp('No data loaded')
        return
      else
        disp('No groups defined. Creating one big group.')
        oldSelect = get(handles.selectData,'Value');

        set(handles.selectData,'Value',1:length(data));
        defineGroup();
        set(handles.selectData,'Value',oldSelect);
      end
    end

    checkVersions()

    setExportFile();

    exportToXML();
    saveConfig();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Make sure the user is not using one of the versions flagged as bad
  % or mixing revisions.

  function checkVersions()

    % Any versions which are BAD should be listed here
    badVersions = [];

    % Equivalent versions should be grouped together
    versionGroups = { [1], [465, inf] }; % [startVersion, endVersion]

    allVer = cat(1,data.version);

    % Many updates does not directly affect the detection, GUI changes etc
    % Make sure they are grouped together into equivalence classes of
    % version so the user can uppgrade then getting a warning about
    % about bad data from two equivalent releases.

    for i = 1:length(versionGroups)
      idx = find(versionGroups{i}(1) <= allVer ...
		 & allVer <= versionGroups{i}(end));
      allVer(idx) = versionGroups{i}(1);
    end

    uniqueVer = unique(allVer);

    warningMsg = [];

    % Check that the user only used one version of SynD to analyse data
    if(numel(uniqueVer) > 1)
      warningMsg = 'Your data is analysed with more than one version of SynD. ';
    end


    if(isnan(getVersion()))
      warningMsg = sprintf('%s%s', warningMsg, ...
			   'Unable to determine version of SynD aggregate. ');
    end

    if(nnz(uniqueVer > getVersion()))
      warningMsg = sprintf('%s%s', warningMsg, ...
                           ['Some data files were analysed by a newer ' ...
		            'version than this version of SynD aggregate. ']);
    end

    % Make sure no data was analysed with a version flagged as bad
    if(nnz(ismember(uniqueVer,badVersions)))
      warningMsg = sprintf('%s%s', warningMsg, ...
			  ['Some of your data is analysed with a version ' ...
			   'of SynD flagged as bad. ']);
    end

    % Recommend that the user reanalyses data that was done in a version
    % prior to version control, or at the least be very careful
    if(numel(find(uniqueVer == 0)))
      warningMsg = sprintf('%s%s', warningMsg, ...
			  ['Some of your data is analysed with an old ' ...
			   'version without version information stored. ' ...
                           'It is important that all files used the ' ...
			   'same version. ']);
    end

    % Found some unknown version, version.txt missing...
    if(nnz(isnan(uniqueVer)))
      warningMsg = sprintf('%s%s', warningMsg, ...
			  ['Some of your data is analysed with an unknown ' ...
		           'version. Is version.txt missing? ']);
    end

    % If there are some warnings, show them to the user!
    if(~isempty(warningMsg))
      warningMsg = sprintf('%s%s', warningMsg, ...
			  ['Please reanalyse with the latest version of ' ...
			   'SynD to ensure data integrity. You can use ' ...
                           'the reanalyse function in SynD extract.']);

      uiwait(warndlg(warningMsg, 'No mixing!','modal'));
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function dataListboxClicked(source, event)

    %disp('Data listbox clicked.')

    if(isempty(data))
      loadData();
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setExportFile()
    
    oldPath = pwd;
    try
      cd(exportPath);
    catch
      fprintf('Unable to change to %s\n', char(exportPath))
    end

    groupNameList = group(1).Name;
    for i=2:length(group)
      groupNameList = sprintf('%s-%s', groupNameList,group(i).Name);
    end

    exportFile = sprintf(exportFileMask,groupNameList);

    [exportFile, exportPath] = ...
	  uiputfile('*-export-summary.xml', ...
		    'Export aggregated data', ...
		    exportFile);

    cd(oldPath);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function exportToXML()

    if(isempty(exportFile) | exportFile == 0)
      disp('No export file selected')
      return
    end

    if(isempty(group(1).Members))
      uiwait(warndlg('No groups specified','Export error','modal'));
      return
    end

    % Header info, to make excel feel comfortable with our data...
    docNode = makeXMLheader();

    [columnName, columnData] = makeInfoSheet();

    disp('Creating synapse average sheet')
    makeXMLsheet(docNode, 'Info', columnName, columnData);

    for iG = 1:length(group)
      [columnName, columnData] = aggregateMeanSynapses(group(iG).Members);
      makeXMLsheet(docNode,sprintf('Basics_%s',group(iG).Name), ...
		   columnName, columnData);
    end

    disp('Creating intensity sheets')

    % Write histogram with channel intensities
    for iG = 1:length(group)
      [columnName, columnData] = ...
      aggregateIntensityHist(group(iG).Members,'morph');
      makeXMLsheet(docNode, ....
                   sprintf('M_%s', ...
			   group(iG).Name), ...
		   columnName, columnData);
    end

    for iG = 1:length(group)
      [columnName, columnData] = ...
	       aggregateIntensityHist(group(iG).Members,'syn');
      makeXMLsheet(docNode, ...
    	           sprintf('S_%s', ...
			   group(iG).Name), ...
		 columnName, columnData);
    end

    for iG = 1:length(group)
      [columnName, columnData] = aggregateIntensityHist(group(iG).Members,'X');
      makeXMLsheet(docNode, ...
		   sprintf('X_%s', ...
			   group(iG).Name), ...
		   columnName, columnData);
    end

    % Sholl 
    disp('Creating sholl analysis sheet')

    for iG = 1:length(group)
      [columnName, columnData] = aggregateSholl(group(iG).Members,'DendHist');
      makeXMLsheet(docNode, ...
		   sprintf('Sholl dend_%s', group(iG).Name), ...
		 columnName, columnData);
    end

    for iG = 1:length(group)
      [columnName, columnData] = aggregateSholl(group(iG).Members,'SynHist');
      makeXMLsheet(docNode, ...
		   sprintf('Sholl syn_%s', group(iG).Name), ...
		   columnName, columnData);
    end

    for iG = 1:length(group)
      [columnName, columnData] = ...
	       aggregateSholl(group(iG).Members,'IntMorphMean');
      makeXMLsheet(docNode, ...
		   sprintf('Sholl M_%s', group(iG).Name), ...
		 columnName, columnData);
    end

    for iG = 1:length(group)
      [columnName, columnData] = aggregateSholl(group(iG).Members,'IntSynMean');
      makeXMLsheet(docNode, ...
		   sprintf('Sholl S_%s', group(iG).Name), ...
		 columnName, columnData);
    end

    for iG = 1:length(group)
      [columnName, columnData] = aggregateSholl(group(iG).Members,'IntXMean');
      makeXMLsheet(docNode, ...
		   sprintf('Sholl X_%s', group(iG).Name), ...
		 columnName, columnData);
    end

    % Fractions of intensities (aggregateSholl then gets extra argument)

    for iG = 1:length(group)
      [columnName, columnData] = aggregateSholl(group(iG).Members, ...
						'IntSynMean', ...
						'IntMorphMean');
      makeXMLsheet(docNode, ...
		   sprintf('Sholl r(S,M)_%s', group(iG).Name), ...
		 columnName, columnData);
    end

    for iG = 1:length(group)
      [columnName, columnData] = aggregateSholl(group(iG).Members, ...
						'IntXMean', ...
						'IntMorphMean');
      makeXMLsheet(docNode, ...
		   sprintf('Sholl r(X,M)_%s', group(iG).Name), ...
		 columnName, columnData);
    end

    for iG = 1:length(group)
      [columnName, columnData] = aggregateSholl(group(iG).Members, ...
						'IntXMean', ...
						'IntSynMean');
      makeXMLsheet(docNode, ...
		   sprintf('Sholl r(X,S)_%s', group(iG).Name), ...
		 columnName, columnData);
    end

    disp('Creating synapse profile sheets')

    for iG = 1:length(group)
      [columnName, columnData] = aggregateSynapseProfile(group(iG).Members, ...
							 'Morph');
      makeXMLsheet(docNode, ...
		   sprintf('Profile M_%s', group(iG).Name), ...
		   columnName, columnData);
    end

    for iG = 1:length(group)
      [columnName, columnData] = aggregateSynapseProfile(group(iG).Members, ...
							 'Syn');
      makeXMLsheet(docNode, ...
		   sprintf('Profile S_%s', group(iG).Name), ...
		   columnName, columnData);
    end

    for iG = 1:length(group)
      [columnName, columnData] = aggregateSynapseProfile(group(iG).Members, ...
							 'X');
      makeXMLsheet(docNode, ...
		   sprintf('Profile X_%s', group(iG).Name), ...
		   columnName, columnData);
    end

    for iG = 1:length(group)
      [columnName, columnData] = aggregateSynapseProfile(group(iG).Members, ...
							 'Syn','Morph');
      makeXMLsheet(docNode, ...
		   sprintf('Profile r(S,M)_%s', group(iG).Name), ...
		   columnName, columnData);
    end

    for iG = 1:length(group)
      [columnName, columnData] = aggregateSynapseProfile(group(iG).Members, ...
							 'X','Morph');
      makeXMLsheet(docNode, ...
		   sprintf('Profile r(X,M)_%s', group(iG).Name), ...
		   columnName, columnData);
    end

    for iG = 1:length(group)
      [columnName, columnData] = aggregateSynapseProfile(group(iG).Members, ...
							 'X','Syn');
      makeXMLsheet(docNode, ...
		   sprintf('Profile r(X,S)_%s', group(iG).Name), ...
		   columnName, columnData);
    end

    for iG = 1:length(group)
      [columnName, columnData] = aggregateNeuriteGradient(group(iG).Members, ...
							 'morph');

      makeXMLsheet(docNode, ...
		   sprintf('Gradient M_%s', group(iG).Name), ...
		   columnName, columnData);
    end

    for iG = 1:length(group)
      [columnName, columnData] = aggregateNeuriteGradient(group(iG).Members, ...
							 'syn');
      makeXMLsheet(docNode, ...
		   sprintf('Gradient S_%s', group(iG).Name), ...
		   columnName, columnData);

    end

    for iG = 1:length(group)
      [columnName, columnData] = aggregateNeuriteGradient(group(iG).Members, ...
							 'X');
      makeXMLsheet(docNode, ...
		   sprintf('Gradient X_%s', group(iG).Name), ...
		   columnName, columnData);

    end

    for iG = 1:length(group)
      [columnName, columnData] = aggregateNeuriteGradient(group(iG).Members, ...
							  'syn','morph');
      makeXMLsheet(docNode, ...
		   sprintf('Gradient r(S,M)_%s', group(iG).Name), ...
		   columnName, columnData);

    end

    for iG = 1:length(group)
      [columnName, columnData] = aggregateNeuriteGradient(group(iG).Members, ...
							  'X','morph');
      makeXMLsheet(docNode, ...
		   sprintf('Gradient r(X,M)_%s', group(iG).Name), ...
		   columnName, columnData);

    end

    [columnName, columnData] = aggregateParameters();
    makeXMLsheet(docNode, 'Parameters', columnName, columnData);


    % Write all to disk

    fprintf('Exporting data to %s\n',exportFile);
    xmlwrite(strcat(exportPath,exportFile),docNode);

    uiwait(helpdlg('Export done. Time for table tennis!', 'SynD'));

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
        disp('Warning, XML write called with matrix, this is bad! Report it.')
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

  function [columnName, columnData] = makeInfoSheet()

    columnName = { 'Abbreviation', 'Description' };
    shortName = {};
    fullName = {};

    for i = 1:length(group)

      shortName{i,1} = sprintf('Basics_%s',group(i).Name);
      fullName{i,1} = ...
	sprintf('Basic information for all synapses in group %s', group(i).Name);

      shortName{i,2} = sprintf('M_%s',group(i).Name);
      fullName{i,2} = ...
	sprintf('Average intensity of morphology channel for synapses in group %s', group(i).Name);

      shortName{i,3} = sprintf('S_%s',group(i).Name);
      fullName{i,3} = ...
	sprintf('Average intensity of synapse channel for synapses in group %s', group(i).Name);

      shortName{i,4} = sprintf('X_%s',group(i).Name);
      fullName{i,4} = ...
	sprintf('Average intensity of synapse channel for synapses in group %s', group(i).Name);

      shortName{i,5} = sprintf('Sholl dend_%s',group(i).Name);
      fullName{i,5} = ...
	sprintf('Sholl analysis of dendrites in group %s', group(i).Name);

      shortName{i,6} = sprintf('Sholl syn_%s',group(i).Name);
      fullName{i,6} = ...
	sprintf('Sholl analysis of synapses in group %s', group(i).Name);

      shortName{i,7} = sprintf('Sholl M_%s',group(i).Name);
      fullName{i,7} = ...
	sprintf('Average intensity of morphology channel for synapses at different distances from soma in group %s', group(i).Name);

      shortName{i,8} = sprintf('Sholl S_%s',group(i).Name);
      fullName{i,8} = ...
	sprintf('Average intensity of synapse channel for synapses at different distances from soma in group %s', group(i).Name);

      shortName{i,9} = sprintf('Sholl X_%s',group(i).Name);
      fullName{i,9} = ...
	sprintf('Average intensity of X channel for synapses at different distances from soma in group %s', group(i).Name);

      shortName{i,10} = sprintf('Sholl r(S,M)_%s',group(i).Name);
      fullName{i,10} = ...
	sprintf('Average ratio of synapse channel to morphology channel for synapses at different distances from soma in group %s', group(i).Name);

      shortName{i,11} = sprintf('Sholl r(X,M)_%s',group(i).Name);
      fullName{i,11} = ...
	sprintf('Average ratio of X channel to morphology channel for synapses at different distances from soma in group %s', group(i).Name);

      shortName{i,12} = sprintf('Sholl r(X,S)_%s',group(i).Name);
      fullName{i,12} = ...
	sprintf('Average ratio of X channel to synapse channel for synapses at different distances from soma in group %s', group(i).Name);

      shortName{i,13} = sprintf('Profile M_%s',group(i).Name);
      fullName{i,13} = ...
	sprintf('Radial profile of morphology channel around single synapses in group %s', group(i).Name);

      shortName{i,14} = sprintf('Profile S_%s',group(i).Name);
      fullName{i,14} = ...
	sprintf('Radial profile of synapse channel around single synapses in group %s', group(i).Name);

      shortName{i,15} = sprintf('Profile X_%s',group(i).Name);
      fullName{i,15} = ...
	sprintf('Radial profile of X channel around single synapses in group %s', group(i).Name);

      shortName{i,16} = sprintf('Profile r(S,M)_%s',group(i).Name);
      fullName{i,16} = ...
	sprintf('Radial profile of ratio of synapse channel to morphology channel around single synapses in group %s', group(i).Name);

      shortName{i,17} = sprintf('Profile r(X,M)_%s',group(i).Name);
      fullName{i,17} = ...
	sprintf('Radial profile of ratio of X channel to morphology channel around single synapses in group %s', group(i).Name);

      shortName{i,18} = sprintf('Profile r(X,M)_%s',group(i).Name);
      fullName{i,18} = ...
	sprintf('Radial profile of ratio of X channel to synapse channel around single synapses in group %s', group(i).Name);

      shortName{i,19} = sprintf('Gradient M_%s',group(i).Name);
      fullName{i,19} = ...
	sprintf('Morphology channel intensity gradient along neurite, distance measured from soma, in group %s (0 to 0 column is soma)', group(i).Name);

      shortName{i,20} = sprintf('Gradient S_%s',group(i).Name);
      fullName{i,20} = ...
	sprintf('Synapse channel intensity gradient along neurite, distance measured from soma, in group %s (0 to 0 column is soma)', group(i).Name);

      shortName{i,21} = sprintf('Gradient X_%s',group(i).Name);
      fullName{i,21} = ...
	sprintf('X channel intensity gradient along neurite, distance measured from soma, in group %s (0 to 0 column is soma)', group(i).Name);

    end

    columnData = { reshape(shortName,numel(shortName),1), ...
		   reshape(fullName,numel(fullName),1) };

    columnData{1}{end+1} = 'Parameters';
    columnData{2}{end+1} = 'Parameters used for neurite and synapse detection';

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [columnName, columnData] = aggregateParameters()

    columnName = {'File name', ...
		  'Version', ...
		  'Pixel size (micrometers)', ...
		  'Morph channel', ...
		  'Syn channel', ...
		  'X channel', ...
		  'Wiener filter size (pixels)', ...
		  'Soma threshold', ...
		  'Soma erode radius (pixels)', ...
		  'Exclude soma synapses', ...
		  'Measure point radius (pixels)', ...
		  'Default #soma measures', ...
		  'Steerable filters', ...
		  'Filter size (micrometers)', ...
		  'Max add cost', ...
		  'Connect cost lambda', ...
		  'Steerable filter alpha', ...
		  'Min protrusion length (micrometers)', ...
		  'Background filtering', ...
		  'Smoothing radius (pixels)', ...
		  'Synapse threshold (std)', ...
		  'Synapse threshold (intensity, applied to smoothed and possibly BG removed image)', ...
		  'Min synapse size (micrometers2)', ...
		  'Neurite padding (micrometers)', ...
		  'Single synapse kernel radie (micrometers)', ...
		  'Max synapse profile radie (pixels)', ...
		  'Trim neurite size (pixels)', ...
		  'Sholl bin size', ...
		  'Intensity bin size', ...
                  };

    try
      tmpDet = cat(1,data.detection);
    catch
      disp(['If you see this then Johannes updated the data.detection ' ...
            'structure in SynD_extract.m without properly checking that ' ...
  	    'SynD_aggregate.m can handle it. Email him at hjorth@kth.se!'])
      keyboard
    end

    allFileNames = {};
    for i = 1:length(data)
      allFileNames{i} = data(i).fileName{1};
    end

    columnData = { allFileNames, ...
		   cat(1,data.version), ...
		   cat(1,data.xyRes)*1e6, ...
		   cat(1,data.morphChannel), ...
		   cat(1,data.synChannel), ...
		   cat(1,data.XChannel), ...
		   cat(1,tmpDet.wienerSize), ...
		   cat(1,data.somaIntensityThreshold), ...
		   cat(1,tmpDet.somaErodeRadius), ...
		   yesNo(cat(1,data.excludeSomaSynapses)), ...
		   cat(1,tmpDet.measurePointRadius), ...
		   cat(1,tmpDet.nMeasurePoints), ...
		   yesNo(cat(1,tmpDet.useSteerableFilters)), ...
		   cat(1,tmpDet.filterSize)*1e6, ...
		   cat(1,tmpDet.maxAddCost), ...
		   cat(1,tmpDet.connectCostLambda), ...
		   cat(1,tmpDet.steerableFilterAlpha), ...
		   cat(1,tmpDet.minProtrusionLength)*1e6, ...
		   yesNo(cat(1,tmpDet.backgroundFiltering)), ...
		   cat(1,tmpDet.smoothingRadius), ...
		   numOrAuto(cat(1,tmpDet.synThreshold)), ...
		   cat(1,data.synapseIntensityThreshold), ...
		   cat(1,tmpDet.minSynapseSize)*1e12, ...
		   cat(1,tmpDet.neuritePadding)*1e6, ...
		   numOrAuto(cat(1,tmpDet.singleSynapseRadie)*1e6), ...
		   cat(1,tmpDet.maxRadie), ...
		   cat(1,tmpDet.trimNeuriteSize), ...
		   cat(1,tmpDet.shollBinSize), ...
		   cat(1,tmpDet.intensityBinSize), ...
                  };

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function status = yesNo(flag)
 
    status = {};

    for iS = 1:length(flag)
      if(flag(iS))
	status{iS} = 'yes';
      else
        status{iS} = 'no';
      end
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function status = numOrAuto(flag)

    status = {};

    for iS = 1:length(flag)
      if(isnan(flag(iS)))
	status{iS} = 'Auto';
      else
        status{iS} = num2str(flag(iS));
      end
    end


  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [columnName, columnData] = aggregateMeanSynapses(dataIdx)

    columnName = {'File name', ...
		  'Number of synapses', ...
		  'Total dendritic length (micrometer)', ...
		  'Synapses/micrometer', ...
		  'Mean synapse area (micrometer^2)', ...
		  'SEM synapse area (micrometer^2)', ...
		  'Mean morph channel', ...
		  'SEM morph channel', ...
		  'Mean synapse channel', ...
		  'SEM synapse channel', ...
		  'Mean X channel', ...
		  'SEM X channel', ...
		  'Mean syn/morph', ...
		  'Mean X/morph', ...
		  'Mean X/syn', ...
 	          'Number of somas' };

    % clear fileNames nSynapses totLength meanSynArea 
    % clear meanSynMorphInt meanSynSynInt meanSynXInt
    % clear fracXM fracXM fracXS

    for i = 1:length(dataIdx)
      idx = dataIdx(i);
      fileNames{i} = data(idx).fileName{1};
      nSynapses(i) = length(data(idx).synapseCenter);
      totLength(i) = data(idx).totDendLength*1e6;

      meanSynArea(i) = mean(data(idx).synapseArea);
      semSynArea(i) = std(data(idx).synapseArea) ...
 	               / sqrt(length(data(idx).synapseArea));

      meanSynMorphInt(i) = mean(data(idx).synapseIntensityMorphMean);
      semSynMorphInt(i) = std(data(idx).synapseIntensityMorphMean) ...
	                   / sqrt(length(data(idx).synapseIntensityMorphMean));

      meanSynSynInt(i) = mean(data(idx).synapseIntensitySynMean);
      semSynSynInt(i) = std(data(idx).synapseIntensitySynMean) ...
	                  / sqrt(length(data(idx).synapseIntensitySynMean));

      meanSynXInt(i) = mean(data(idx).synapseIntensityXMean);
      semSynXInt(i) = std(data(idx).synapseIntensityXMean) ...
	                  / sqrt(length(data(idx).synapseIntensityXMean));

      % To avoid division by zero we discard data points that are zero
      okMidx = find(data(idx).synapseIntensityMorphMean ~= 0);
      okSidx = find(data(idx).synapseIntensitySynMean ~= 0);

      fracSM(i) = mean(data(idx).synapseIntensitySynMean(okMidx) ...
		       ./data(idx).synapseIntensityMorphMean(okMidx), 'omitnan');
      fracXM(i) = mean(data(idx).synapseIntensityXMean(okMidx) ...
		       ./data(idx).synapseIntensityMorphMean(okMidx), 'omitnan');
      fracXS(i) = mean(data(idx).synapseIntensityXMean(okSidx) ...
		       ./data(idx).synapseIntensitySynMean(okSidx), 'omitnan');

      nSoma(i) = data(idx).nSoma;

    end

    nSynPerLen = nSynapses ./ totLength;

    % At the end Sabine wanted mean,std,sem,N
    fileNames{end+1} = 'Mean';
    fileNames{end+1} = 'Std';
    fileNames{end+1} = 'SEM';
    fileNames{end+1} = 'N';

    % makeSummary adds mean, std and SEM of the column to the end
    nSynapses       = [nSynapses,       makeSummary(nSynapses)];
    totLength       = [totLength,       makeSummary(totLength)];
    nSynPerLen      = [nSynPerLen,      makeSummary(nSynPerLen)];
    meanSynArea     = [meanSynArea,     makeSummary(meanSynArea)];
    semSynArea      = [semSynArea,      makeSummary(semSynArea)];
    meanSynMorphInt = [meanSynMorphInt, makeSummary(meanSynMorphInt)];
    semSynMorphInt  = [semSynMorphInt,  makeSummary(semSynMorphInt)];
    meanSynSynInt   = [meanSynSynInt,   makeSummary(meanSynSynInt)];
    semSynSynInt    = [semSynSynInt,    makeSummary(semSynSynInt)];
    meanSynXInt     = [meanSynXInt,     makeSummary(meanSynXInt)];
    semSynXInt      = [semSynXInt,      makeSummary(semSynXInt)];
    fracSM          = [fracSM,          makeSummary(fracSM)];
    fracXM          = [fracXM,          makeSummary(fracXM)];
    fracXS          = [fracXS,          makeSummary(fracXS)];
    nSoma           = [nSoma,           makeSummary(nSoma)];

    columnData = { fileNames, ...
		   nSynapses, ...
		   totLength, ...
		   nSynPerLen, ...
		   meanSynArea, ...
		   semSynArea, ...
		   meanSynMorphInt, ...
		   semSynMorphInt, ...
		   meanSynSynInt, ...  
		   semSynSynInt, ...  
		   meanSynXInt, ...  
		   semSynXInt, ...  
		   fracSM, ...
		   fracXM, ...
		   fracXS, ...
                   nSoma };

    columnName{end+1} = 'Soma area (mean)';
    columnName{end+1} = 'Soma Major Axis (mean)';
    columnName{end+1} = 'Soma Minor Axis (mean)';

    allArea = [];
    allMajor = [];
    allMinor = [];

    % These measures are always present
    for i = 1:length(dataIdx)
      allArea(i)  = mean(data(dataIdx(i)).somaArea, 'omitnan');
      allMajor(i) = mean(data(dataIdx(i)).somaMajorAxisLength, 'omitnan');
      allMinor(i) = mean(data(dataIdx(i)).somaMinorAxisLength, 'omitnan');

    end

    columnData{end+1} = [allArea,  makeSummary(allArea)];
    columnData{end+1} = [allMajor, makeSummary(allMajor)];
    columnData{end+1} = [allMinor, makeSummary(allMinor)];

    % These measures are optional
    if(~isempty(data(dataIdx(1)).somaMeasureMorph))
      columnName{end+1} = 'Soma M (mean)';
      columnName{end+1} = 'Soma S (mean)';
      columnName{end+1} = 'Soma X (mean)';
      columnName{end+1} = 'Synapse M/Soma M (mean)';
      columnName{end+1} = 'Synapse S/Soma S (mean)';
      columnName{end+1} = 'Synapse X/Soma X (mean)';
      columnName{end+1} = 'Synapse (S/M) / Soma (S/M) (mean)';
      columnName{end+1} = 'Synapse (X/M / Soma (X/M) (mean)';
      columnName{end+1} = 'Synapse (X/S) / Soma (X/S) (mean)';

      allSomaMorph = [];
      allSomaSyn = [];
      allSomaX = [];

      allFracMorph = [];
      allFracSyn = [];
      allFracX = [];
      allFracSynMorph = [];
      allFracXMorph = [];
      allFracXSyn = [];

      for i = 1:length(dataIdx)
	allSomaMorph(i) = mean(data(dataIdx(i)).somaMeasureMorph, 'omitnan');
	allSomaSyn(i)   = mean(data(dataIdx(i)).somaMeasureSyn, 'omitnan');
	allSomaX(i)     = mean(data(dataIdx(i)).somaMeasureX, 'omitnan');

        allFracMorph(i) = mean(data(dataIdx(i)).synapseIntensityMorphMean) ...
	                  / allSomaMorph(i);
        allFracSyn(i)   = mean(data(dataIdx(i)).synapseIntensitySynMean) ...
	                  / allSomaSyn(i);
        allFracX(i) = mean(data(dataIdx(i)).synapseIntensityXMean)/allSomaX(i);

      end

      allFracSynMorph = allFracSyn ./ allFracMorph;
      allFracXMorph = allFracX ./ allFracMorph;
      allFracXSyn = allFracX ./ allFracSyn;

      columnData{end+1} = [allSomaMorph,    makeSummary(allSomaMorph)];
      columnData{end+1} = [allSomaSyn,      makeSummary(allSomaSyn)];
      columnData{end+1} = [allSomaX,        makeSummary(allSomaX)];
      columnData{end+1} = [allFracMorph,    makeSummary(allFracMorph)];
      columnData{end+1} = [allFracSyn,      makeSummary(allFracSyn)];
      columnData{end+1} = [allFracX,        makeSummary(allFracX)];
      columnData{end+1} = [allFracSynMorph, makeSummary(allFracSynMorph)];
      columnData{end+1} = [allFracXMorph,   makeSummary(allFracXMorph)];
      columnData{end+1} = [allFracXSyn,     makeSummary(allFracXSyn)];

    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [columnName, columnData] = aggregateIntensityHist(dataIdx,histType)

    columnName = { 'From', 'To' };

    e1 = data(dataIdx(1)).intensityHistogramEdges(2) ...
      - data(dataIdx(1)).intensityHistogramEdges(1);

    for i = 2:length(dataIdx)
      e2 = data(dataIdx(i)).intensityHistogramEdges(2) ...
	   - data(dataIdx(i)).intensityHistogramEdges(1);

      if(e1 ~= e2)
	uiwait(errordlg(sprintf('Different intensity bin sizes in %s and %s', ...
				data(dataIdx(1)).fileName{1}, ...
				data(dataIdx(i)).fileName{1}), ...
			'SynD: Histogram edges inconsistent'))
	columnName = {'Inconsistent histogram edges'};
        columnData = {[-1]};
        return
      end

    end

    [maxRows,maxRowsIdx] = maxVectLen('intensityHistogramEdges');
    maxRows = maxRows - 1; % There are n+1 edges when there is n intervalls

    columnData = { data(dataIdx(1)).intensityHistogramEdges(1:end-1), ...
		   data(dataIdx(1)).intensityHistogramEdges(2:end) };

    columnSum = zeros(maxRows,1);
    
    for i = 1:length(dataIdx)
      idx = dataIdx(i);
      columnName{end+1} = data(idx).fileName{1};
      columnData{end+1} = eval(sprintf('data(idx).%sHist(1:end-1)',histType));

      nRows = length(columnData{end});

      % Add them together
      for j = 1:nRows
        columnSum(j) = columnSum(j) + columnData{end}(j);
      end

      if(nnz(data(maxRowsIdx).intensityHistogramEdges(1:nRows) ...
	     - data(idx).intensityHistogramEdges(1:nRows) > 1e-9))
        disp('Histogram edges are different - talk to Johannes')
        keyboard
      end

    end

    columnName{end+1} = 'Sum';
    columnData{end+1} = columnSum;

    csum = cumsum(columnSum);

    columnName{end+1} = 'Cum sum';
    columnData{end+1} = csum;

    columnName{end+1} = '%';
    columnData{end+1} = csum / csum(end) * 100;

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % countTypeDenominator is optional
  %

  function [columnName, columnData] = aggregateSholl(dataIdx, ...
						     countType, ...
						     countTypeDenominator)
    columnName = { 'From (mu)', 'To (mu)' };

    e1 = data(dataIdx(1)).shollEdges(2) - data(dataIdx(1)).shollEdges(1);
    for i = 2:length(dataIdx)
      e2 = data(dataIdx(i)).shollEdges(2) - data(dataIdx(i)).shollEdges(1);
 
      if(e1 ~= e2)
	uiwait(errordlg(sprintf('Different sholl bin sizes in %s and %s', ...
				data(dataIdx(1)).fileName{1}, ...
				data(dataIdx(i)).fileName{1}), ...
			'SynD: Histogram edges inconsistent'))

	columnName = {'Inconsistent sholl bin size'};
        columnData = {[-1]};
        return
      end
    end

    maxRows = 0;
    maxRowsIdx = 0;

    for i = 1:length(dataIdx)
      curRows = length(data(dataIdx(i)).shollEdges)-1;
      if(curRows > maxRows)
	maxRows = curRows;
        maxRowsIdx = i;
      end
    end

    columnData = { data(dataIdx(maxRowsIdx)).shollEdges(1:end-1)*1e6, ...
		   data(dataIdx(maxRowsIdx)).shollEdges(2:end)*1e6 };

    allData = NaN*zeros(maxRows,length(dataIdx));

    for i = 1:length(dataIdx)
      idx = dataIdx(i);
      columnName{end+1} = data(idx).fileName{1};

      if(~exist('countTypeDenominator'))
        columnData{end+1} = eval(sprintf('data(idx).sholl%s', ...
	 	 	 	         countType));
      else
        % This allows us to export ratios
        columnData{end+1} = eval(sprintf('data(idx).sholl%s./data(idx).sholl%s', ...
	 	 	 	         countType, countTypeDenominator));

      end

      allData(1:numel(columnData{end}),i) = columnData{end}(:);

    end

    columnName{end+1} = 'Mean';
    columnData{end+1} = mean(allData,2,'omitnan');

    columnName{end+1} = 'Std';
    columnData{end+1} = std(allData,[],2, 'omitnan');

    columnName{end+1} = 'SEM';
    nOk = sum(~isnan(allData),2);
    columnData{end+1} = std(allData,[],2,'omitnan') ./ sqrt(nOk);

    columnName{end+1} = 'N';
    columnData{end+1} = nOk;

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [columnName, columnData] = aggregateSynapseProfile(dataIdx, ...
							      profileType, ...
							      profileTypeDenom)
    columnName = { 'Distance (micrometer)' };

    [maxLen, maxIdx] = maxVectLen('meanSynapseProfileDist', dataIdx);

    if(maxLen == 0)
      columnName = { 'No data' };
      columnData = { [-1] };
      return
    end

    origEdges = data(maxIdx).meanSynapseProfileDist;
    columnData = { data(maxIdx).meanSynapseProfileDist *1e6};

    allData = NaN*zeros(maxLen,length(dataIdx));

    for i = 1:length(dataIdx)
      if(length(origEdges) ...
	 ~= length(data(dataIdx(i)).meanSynapseProfileDist) ...
         | nnz(abs(origEdges-data(dataIdx(i)).meanSynapseProfileDist)>1e-9))
	uiwait(errordlg(sprintf(['Different distance edges in %s and %s. ' ...
				 'Is the pixel size different ' ...
                                 '(in micrometers) ' ...
				 'for the two files?'], ...
				data(maxIdx).fileName{1}, ...
				data(dataIdx(i)).fileName{1}), ...
			'SynD: Inconsistent synapse profile data.'))
	columnName = {'Inconsistent distance points'};
        columnData = {[-1]};
        return
      end

      columnName{end+1} = data(dataIdx(i)).fileName{1};

      if(~exist('profileTypeDenom'))
        tmp = eval(sprintf('data(dataIdx(i)).meanSynapseProfile%s', ...
			   profileType));
        columnData{end+1} = tmp;

        allData(1:length(tmp),i) = tmp;

      else
	tmp = eval(sprintf(['data(dataIdx(i)).meanSynapseProfile%s' ...
			    './data(dataIdx(i)).meanSynapseProfile%s'], ...
			   profileType, profileTypeDenom));

        columnData{end+1} = tmp;
        allData(1:length(tmp),i) = tmp;

      end
    end

    columnName{end+1} = 'Mean';
    columnData{end+1} = mean(allData,2,'omitnan');

    columnName{end+1} = 'Std';
    columnData{end+1} = std(allData,[],2,'omitnan');

    columnName{end+1} = 'SEM';
    nOk = sum(~isnan(allData),2);
    columnData{end+1} = std(allData,[],2,'omitnan') ./ sqrt(nOk);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [columnName, columnData] = aggregateNeuriteGradient(dataIdx, ...
							       gradientType, ...
							       gradientTypeDenom)

    columnName = { 'From (micrometer)', 'To (micrometer)' };

    [maxLen, maxIdx] = maxVectLen('gradientEdges', dataIdx);

    if(maxLen == 0)
      columnName = { 'No data' };
      columnData = { [-1] };
      return
    end

    origEdges = data(maxIdx).gradientEdges;
    columnData = { [0, data(maxIdx).gradientEdges(1:end-1)] *1e6 };
    columnData{2} = data(maxIdx).gradientEdges*1e6;

    allData = NaN*zeros(maxLen,length(dataIdx));

    for i = 1:length(dataIdx)
      tmp = data(dataIdx(i)).gradientEdges;
      % We can get rounding errors, hence why I check with 1e-9 rather than 0
      if(nnz(abs(origEdges(1:length(tmp))-tmp)>1e-9))
	uiwait(errordlg(sprintf('Different gradient edges in %s and %s', ...
				data(maxIdx).fileName{1}, ...
				data(dataIdx(i)).fileName{1}), ...
			'SynD: Inconsistent neurite gradient data'))
	columnName = {'Inconsistent distance points'};
        columnData = {[-1]};
        keyboard
        return
      end

      columnName{end+1} = data(dataIdx(i)).fileName{1};

      if(exist('gradientTypeDenom'))
        tmp = eval(sprintf(['data(dataIdx(i)).%sGradientMean' ...
                            './data(dataIdx(i)).%sGradientMean'] , ...
			   gradientType,gradientTypeDenom));

      else
        tmp = eval(sprintf('data(dataIdx(i)).%sGradientMean', ...
			   gradientType));
      end

      columnData{end+1} = tmp;

      allData(1:length(tmp),i) = tmp;

    end

    columnName{end+1} = 'Mean';
    columnData{end+1} = mean(allData,2, 'omitnan');

    columnName{end+1} = 'Std';
    columnData{end+1} = std(allData,[],2, 'omitnan');

    columnName{end+1} = 'SEM';
    nOk = sum(~isnan(allData),2);
    columnData{end+1} = std(allData,[],2,'omitnan') ./ sqrt(nOk);

    columnName{end+1} = 'N';
    columnData{end+1} = nOk;


  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function enableGroupGUI()

    nMax = length(handles.groupBox);

    activeGroups = 1:length(group);
    inactiveGroups = setdiff(1:nMax,activeGroups);

    for i = 1:nMax
      if(ismember(i,activeGroups))
	set(handles.groupLabels(i),'String', group(i).Name)
	set(handles.groupLabels(i),'visible','on')
	  set(handles.groupBox(i),'visible','on', 'Value', 1)
      else
	set(handles.groupLabels(i),'visible','off')
        set(handles.groupBox(i),'visible','off', 'Value', 0)
      end
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function groupID = getGroupsForPlot()

    groupVisible = zeros(1,length(group));

    for i = 1:length(groupVisible)
      if(get(handles.groupBox(i),'Value'))
        groupVisible(i) = 1;
      end
    end

    groupID = find(groupVisible);

    if(isempty(groupID))
      % Our plot commands does not handle no groups selected
      groupID = 1:length(group);

      % Turn on everyone
      enableGroupGUI();
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  function defineGroup(source, event)

    if(isempty(data))
      disp('No data loaded, can not define group.')
      return
    end

    selectIdx = get(handles.selectData,'Value');
    newGroupName = get(handles.groupName,'String');
    morphLabel = get(handles.morphName,'String');
    synLabel = get(handles.synName,'String');
    XLabel = get(handles.XName,'String');


    if(length(newGroupName) > 16)
      uiwait(warndlg(['You are using a long group name, excel allows ' ...
                      'a very limited length on work sheet names. ' ...
		      ' This might cause problems with reading XML file.'], ...
		     'Potential Excel problem'));
    end

    if(isempty(newGroupName))
      uiwait(warndlg('You must give the group a name','Input error'));
      return
    end

    if(isempty(selectIdx))
      uiwait(warndlg('No neurons selected','Input error'));
      return
    end

    % Remove selected neurons from other group if there
    for i = 1:length(group)
      group(i).Members = setdiff(group(i).Members,selectIdx);
    end

    % Clear empty groups, obs iterate from the end, since list shortens
    for i=length(group):-1:1
      if(isempty(group(i).Members))
        group(i) = [];
      end
    end

    oldGroupId = 0;

    % Check if group already exists
    for i=1:length(group)
      if(strcmp(group(i).Name,newGroupName))
        % Put the neurons in the existing group
        oldGroupId = i;
      end
    end

    %Add neurons to group
    if(oldGroupId)
      group(i).Members = union(group(i).Members,selectIdx);
      groupIdx(selectIdx) = oldGroupId;
    else
      tmpGroup.Name = newGroupName;
      tmpGroup.Members = selectIdx;
      group(end+1) = tmpGroup;

      groupIdx(selectIdx) = length(group);
    end

    updateListboxNames();

    setPlotTypeList();

    enableGroupGUI();

    showPlot();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function updateLabels(source, event)

    morphLabel = get(handles.morphName,'String');
    synLabel = get(handles.synName,'String');
    XLabel = get(handles.XName,'String');

    setPlotTypeList();
    showPlot();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function guessGroup(source, event)

    selectIdx = [];

    subStrIdx = strfind(get(handles.selectData,'String'), ...
			get(handles.groupName,'String'));

    for i = 1:length(subStrIdx)
      if(~isempty(subStrIdx{i}))
	selectIdx = [selectIdx;i];
      end
    end

    set(handles.selectData,'Value',selectIdx);
     
    defineGroup();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function updateListboxNames()

    listNames = {};

    for i = 1:length(data)
      if(groupIdx(i))
	listNames{i} = sprintf('%s (%s)',  ...
			       data(i).fileName{1}, ...
			       group(groupIdx(i)).Name);
      else
	listNames{i} = data(i).fileName{1};
      end

    end

    set(handles.selectData,'String',listNames);
    set(handles.selectData,'min',1,'max',length(data));

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setAxis(source, event)

    XlimMin = str2num(get(handles.XlimMin,'String'));
    XlimMax = str2num(get(handles.XlimMax,'String'));
    YlimMin = str2num(get(handles.YlimMin,'String'));
    YlimMax = str2num(get(handles.YlimMax,'String'));

    plotId = get(handles.plotType,'Value');

    plotInfo(plotId).Xlim = [XlimMin XlimMax];
    plotInfo(plotId).Ylim = [YlimMin YlimMax];

    showPlot();
 
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setAxisGUI(axisRange)

    set(handles.XlimMin,'String',num2str(axisRange(1)));
    set(handles.XlimMax,'String',num2str(axisRange(2)));
    set(handles.YlimMin,'String',num2str(axisRange(3)));
    set(handles.YlimMax,'String',num2str(axisRange(4)));

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function loadConfig()
    if(exist(configFile))
      try
        cfg = load(configFile);

        exportPath = cfg.old.exportPath;
        loadPath   = cfg.old.loadPath;

        morphLabel = cfg.old.morphLabel;
        synLabel   = cfg.old.synLabel;
        XLabel     = cfg.old.XLabel;

      catch e
        % getReport(e)
        disp('Failed to read config file.')
      end
    else
      disp('Could not find config file.')
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function saveConfig()
    clear old
    old.exportPath = exportPath;
    old.loadPath = loadPath;
    old.morphLabel = morphLabel;
    old.synLabel = synLabel;
    old.XLabel = XLabel;

    save(configFile,'old'); 
    disp('Saving config.')
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function summary = makeSummary(vect)

    if(min(size(vect)) ~= 1)
      disp('Make summary called with non-vector. Talk to Johannes.')
      keyboard
    end

    summary(1) = mean(vect,'omitnan');
    summary(2) = std(vect,'omitnan');
    summary(3) = std(vect,'omitnan')/sqrt(nnz(~isnan(vect)));
    summary(4) = nnz(~isnan(vect)); % Count the number of non-NaN:s

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function captureKeyPress(source, event)

    switch(event.Key)
      case 'd'
	defineGroup();
      case 'e'
        exportData();
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % If colIdx is specified, all elements will be set to that colour
  % if it is not set, then the markers will get different colours
  % depending on which set they belong to

  function setMarkerColor(plotHandle, colIdx)

    col = colormap('lines');

    for ip = 1:length(plotHandle)

      if(exist('colIdx'))
        if(length(colIdx) == length(plotHandle))
	  faceCol = col(colIdx(ip),:);
        else
  	  faceCol = col(colIdx,:);
        end
      else
	faceCol = col(ip,:);
      end

      set(plotHandle(ip), ...
          'markeredgecolor', [0 0 0], ...
	  'markerfacecolor', faceCol);

    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setFaceColor(plotHandle, colIdx)

    col = colormap('lines');

    for ip = 1:length(plotHandle)

      if(exist('colIdx'))
        if(length(colIdx) == length(plotHandle))
  	  faceCol = col(colIdx(ip),:);
        else
  	  faceCol = col(colIdx,:);
        end
      else
	faceCol = col(ip,:);
      end

      set(plotHandle(ip), ...
          'edgecolor', [0 0 0], ...
	  'facecolor', faceCol);

    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setLineColor(plotHandle, colIdx)

    col = colormap('lines');

    for ip = 1:length(plotHandle)

      if(exist('colIdx'))
        if(length(colIdx) == length(plotHandle))
  	  lineCol = col(colIdx(ip),:);
        else
  	  lineCol = col(colIdx,:);
        end
      else
	lineCol = col(ip,:);
      end

      set(plotHandle(ip), ...
	  'color', lineCol);

    end
  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function col = getGroupColor(groupID)

    cMap = colormap('lines');
    col = cMap(groupID,:);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  function setErrorbarWidth(plotHandle)

    h = get(plotHandle,'children');
    xp = get(h(2),'xdata');

    pAxes = get(plotHandle,'parent');
    pXlim = get(pAxes,'xlim');
    set(pAxes,'units','centimeters')
    pPos = get(pAxes,'position');

    if(strcmpi(get(gca,'Units'),'centimeters'))
      % disp('Units: cm')
      errorbarWidth = 0.1;
    else
      % This should have given us 7.5 mm error bar width
      errorbarWidth = 0.75 * (pXlim(2)-pXlim(1)) / (pPos(1)-pPos(2));
    end 

    idxNaN = find(isnan(xp));

    xp(4:9:end) = xp(1:9:end) - errorbarWidth/2;
    xp(5:9:end) = xp(1:9:end) + errorbarWidth/2;
    xp(7:9:end) = xp(1:9:end) - errorbarWidth/2;
    xp(8:9:end) = xp(1:9:end) + errorbarWidth/2;

    % Those that did not have error bars before should not get it now either
    xp(idxNaN) = NaN;

    set(h(2),'xdata',xp);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % If dataStr is a cell array, first cell is nominator 2nd denominator

  function makeBarGraph(plotId, saveName, subplotFlag)

    allData = {};
    barLabels = {};

    % User selected groups to plot
    plotGroupID = getGroupsForPlot();

    dataPoint    = NaN*zeros(length(plotGroupID),1);
    dataPointSEM = NaN*zeros(length(plotGroupID),1);
    nElem        = NaN*zeros(length(plotGroupID),1);

    for barNum = 1:length(plotGroupID)

      groupId = plotGroupID(barNum);
      
      barLabels{barNum} = group(groupId).Name;
      allData{barNum} = NaN*zeros(length(group(groupId).Members),1);

      for i = 1:length(group(groupId).Members)
        memberIdx = group(groupId).Members(i);

        if(iscell(plotInfo(plotId).DataStr))
          % We got a nominator and denominator to take care of
          tmpA = eval(sprintf('data(memberIdx).%s', ...
			      plotInfo(plotId).DataStr{1}));
          tmpB = eval(sprintf('data(memberIdx).%s', ...
			      plotInfo(plotId).DataStr{2}));

          if(numel(tmpA) == numel(tmpB) & numel(tmpA) > 0)
             % Calculate all ratios then average them together
             % But make sure to avoid division by zero
	     tmpB(find(tmpB == 0)) = NaN;
             tmp = mean(tmpA ./ tmpB,'omitnan');
	  else
	    % For soma and synapse measure ratio, there are not
	    % the same number of data points, thus average first
  	    tmp = mean(tmpA,'omitnan') / mean(tmpB,'omitnan');
          end
        else
          tmp = eval(sprintf('data(memberIdx).%s', ...
			     plotInfo(plotId).DataStr));
          if(numel(tmp) > 1)
	    tmp = mean(tmp, 'omitnan');
          end
        end

	if(numel(tmp) == 0)
          tmp = NaN;
          fprintf('Missing %s for %s.\n', ...
		  plotInfo(plotId).Name, group(groupId).Name)
        end

        try
          allData{barNum}(i) = tmp;
        catch e
  	  getReport(e)
  	  disp('This should not happen, talk to Johannes.')
          keyboard
        end
      end

      dataPoint(barNum) = mean(allData{barNum},'omitnan');
      nElem(barNum) = nnz(~isnan(allData{barNum}));
      dataPointSEM(barNum) = std(allData{barNum},'omitnan') ...
  	                      / sqrt(nElem(barNum));

    end

    if(exist('saveName'))
      if(~exist('subplotFlag'))
        altFig = figure('visible','off');
        set(altFig,'paperunits','centimeters')
        set(altFig,'papertype','A4')
        ax = axes();
        set(ax,'Units','centimeter');
        set(ax,'Position', [2 2 3.3 3.1])
      end
    else
      saveName = [];
      set(handles.fig,'CurrentAxes',handles.plot);
    end

    nIdx = 1:length(plotGroupID);

    p = errorbar(nIdx,dataPoint,dataPointSEM, ...
		 'linestyle','none', ...
		 'linewidth',1, ...
		 'marker','none', ...
		 'color', [0 0 0]);

    setErrorbarWidth(p);

    hold on

    for i = 1:length(plotGroupID)
      b = bar(nIdx(i),dataPoint(i),'facecolor',[0 0 0]);
      setFaceColor(b,plotGroupID(i));

      if(plotInfo(plotId).barGraphDots)
        pm = plot(i+0.3,allData{i},'o','markersize',3,'linewidth',0.75);
        setMarkerColor(pm,plotGroupID(i));
      end
    end

    ylabel(plotInfo(plotId).Ylabel,'fontname','Arial','fontsize', 8)

    if(~isempty(saveName))
      title(sprintf('%s_%s', plotInfo(plotId).Name, saveName), ...
	    'fontname','Arial','fontsize', 8, 'interpreter', 'none')
    else
      title(sprintf('%s', plotInfo(plotId).Name), ...
	    'fontname','Arial','fontsize', 8, 'interpreter', 'none')
    end

    set(gca, 'xtick',nIdx)
    set(gca, 'xticklabel', barLabels)
    set(gca, 'fontsize',7, 'fontname','Arial')
    set(gca, 'linewidth', 1, 'tickdir', 'out')
    set(gca, 'ticklength', get(gca,'ticklength')*3);

    a = axis;
    a(1) = 0.5;
    a(2) = length(plotGroupID) + 0.5;
    a(3) = 0;

    if(~isempty(plotInfo(plotId).Ylim))
      a(3:4) = plotInfo(plotId).Ylim;
    end

    axis(a);
    setAxisGUI(a);

    textY = (a(4)-a(3))*0.05+a(3);
    for i = 1:length(nElem)
      text(i,textY,sprintf('n=%d',nElem(i)), ...
	   'fontname','Arial', ...
	   'fontsize', 7, ...
	   'horizontalalignment', 'center', ...
	   'color', [1 1 1]);
    end

    box off
    hold off

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % maxIdx is the real index to the data, not a index into subIdx

  function [maxLen, maxIdx] = maxVectLen(vectName,subIdx)

    if(~exist('subIdx'))
      subIdx = 1:length(data);
    end

    maxLen = 0;
    maxIdx = 0;

    for iV = 1:length(subIdx)
       vectLen = eval(sprintf('length(data(subIdx(iV)).%s)', vectName));
       if(vectLen > maxLen)
	 maxLen = vectLen;
         maxIdx = subIdx(iV);     
       end
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  function makeCumHist(plotId, saveName, subplotFlag)

    allData = {};
    allEdges = {};

    [maxLenALL, maxIdxALL] = maxVectLen('intensityHistogramEdges');

    graphLegend = {};

    % User selected groups to plot
    plotGroupID = getGroupsForPlot();

    cumHistEdges = {};
    cumHistData = {};
    cumHistSEM = {};

    for histNum = 1:length(plotGroupID)

      groupId = plotGroupID(histNum);
      memberIdx = group(groupId).Members;

      [maxLen,maxIdx] = maxVectLen('intensityHistogramEdges',memberIdx);

      % Initiate memory
      allData{histNum} = NaN*zeros(length(memberIdx), maxLen);
      allEdges{histNum} = allData{histNum};

      for i = 1:length(memberIdx)
        tmpData = eval(sprintf('data(memberIdx(i)).%s', ...
			       plotInfo(plotId).DataStr));
        tmpEdges = data(memberIdx(i)).intensityHistogramEdges;

        nData = length(tmpData);
        nEdges = length(tmpEdges);

        allData{histNum}(i,1:nData) = tmpData;	  
        allEdges{histNum}(i,1:nEdges) = tmpEdges;

        % This one will keep getting overwritten, but since the edges 
        % for a group should be the same (we check that later)
        cumHistEdges{histNum}(1:nEdges) = tmpEdges;

      end

      % !!! Use diff to check all edges are similar
      if(nnz(diff(allEdges{histNum},1,1) > 1e-9))
        fprintf('Error, different edges for group %s. Plot incorrect.\n', ...
		group(groupId).Name)
        beep
        uiwait(errordlg(sprintf(['Error, different edges for group %s. ' ...
                                 'Plot incorrect.\n'], ...
				group(groupId).Name), ...
			'SynD: Cumulative histogram edges inconsistent'))


      end

      graphLegend{histNum} = group(groupId).Name;



      cumHistData{histNum}  = mean(cumsum(allData{histNum},2)...
				   ./ repmat(sum(allData{histNum},2), ...
					     1,size(allData{histNum},2)), ...
				   1, 'omitnan');
      cumHistSEM{histNum}   = std(cumsum(allData{histNum},2),[],1) ...
	                          ./ sqrt(sum(~isnan(allData{histNum}),1), 'omitnan');


    end

    if(exist('saveName'))
      if(~exist('subplotFlag'))
        altFig = figure('visible','off');
        set(altFig,'paperunits','centimeters')
        set(altFig,'papertype','A4')
        ax = axes();
        set(ax,'Units','centimeter');
        set(ax,'Position', [2 2 3.3 3.1])
      end
    else
      saveName = [];
      set(handles.fig,'CurrentAxes',handles.plot);
    end
    
    p = [];
    for i = 1:length(cumHistData)
      nData = length(cumHistData{i});
      p(i) = stairs(cumHistEdges{i}(1:nData),cumHistData{i}*100);
      hold on
    end
    hold off

    legend(p,graphLegend,'location','best')
    xlabel(plotInfo(plotId).Xlabel,'fontsize',8,'fontname','arial')
    ylabel(plotInfo(plotId).Ylabel,'fontsize',8,'fontname','arial')

    if(~isempty(saveName))
      title(sprintf('%s_%s', plotInfo(plotId).Name, saveName), ...
	    'fontname','Arial','fontsize', 8, 'interpreter', 'none')
    else
      title(sprintf('%s', plotInfo(plotId).Name), ...
	    'fontname','Arial','fontsize', 8, 'interpreter', 'none')
    end

    for i = 1:length(plotGroupID)
      set(p(i),'color',getGroupColor(plotGroupID(i)));
    end

    set(gca,'fontsize',7,'fontname','Arial')
    set(gca,'linewidth', 1, 'tickdir', 'out')
    set(gca,'ticklength',get(gca,'ticklength')*3);

    box off
    hold off

    % !!! We can add error bars here if we want

    a = axis; 
    a(3) = 0;

    if(~isempty(plotInfo(plotId).Xlim))
      a(1:2) = plotInfo(plotId).Xlim;
    end

    if(~isempty(plotInfo(plotId).Ylim))
      a(3:4) = plotInfo(plotId).Ylim;
    end

    axis(a);
    setAxisGUI(a);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function makeHist(plotId, saveName, subplotFlag)

    allData = {};
    allCenters = {};
    graphLegend = {};

    [maxLenALL,maxIdxALL] = maxVectLen('shollEdges');
    maxLenALL = maxLenALL - 1; % N intervalls have N+1 edges

    % User selected groups to plot
    plotGroupID = getGroupsForPlot();

    histData = NaN*ones(maxLenALL,length(plotGroupID));
    histSTD = NaN*ones(maxLenALL,length(plotGroupID));

    for histNum = 1:length(plotGroupID)

      groupId = plotGroupID(histNum);
      memberIdx = group(groupId).Members;

      [maxLen,maxIdx] = maxVectLen('shollEdges',memberIdx);
      maxLen = maxLen - 1;

      allData{histNum} = NaN*zeros(maxLen, length(memberIdx));
      allCenters{histNum} = NaN*zeros(maxLen, length(memberIdx));

      for i = 1:length(memberIdx)	
	if(iscell(plotInfo(plotId).DataStr))
	  tmp = eval(sprintf(['data(memberIdx(i)).sholl%s' ...
			      './data(memberIdx(i)).sholl%s'], ...
			     plotInfo(plotId).DataStr{1}, ...
			     plotInfo(plotId).DataStr{2}));


        else 
          tmp = eval(sprintf('data(memberIdx(i)).sholl%s', ...
			     plotInfo(plotId).DataStr));
        end

	if(isempty(tmp))
	  allData{histNum}(:,i) = NaN;
        else
          allData{histNum}(1:length(tmp),i) = tmp;
        end

        tmp = (data(memberIdx(i)).shollEdges(1:end-1) ...
	       + data(memberIdx(i)).shollEdges(2:end))/2;

        if(isempty(tmp))
 	  allCenters{histNum}(:,i) = NaN;
        else
          allCenters{histNum}(1:length(tmp),i) = tmp*1e6; % To get micrometers
        end
      end

      if(nnz(diff(allCenters{histNum},1,2) > 1e-9))
        fprintf('Error, different edges for group %s. Plot incorrect.\n', ...
		group(groupId).Name)
        beep
      end

      graphLegend{histNum} = group(groupId).Name;

      tmpData = mean(allData{histNum},2, 'omitnan');
      histData(1:length(tmpData),histNum) = tmpData;
      nElem = sum(~isnan(allData{histNum}),2);

      % !!! Just zero or one data point, then do not show errorbars for them
      nElem(find(nElem < 2)) = NaN;

      histDataSEM(1:length(tmpData),histNum) ...
        = std(allData{histNum},[],2, 'omitnan') ./ sqrt(nElem);

      histCenters(1:length(tmpData),histNum) = allCenters{histNum}(:,1);

    end

    if(exist('saveName'))
      if(~exist('subplotFlag'))
        altFig = figure('visible','off');
        set(altFig,'paperunits','centimeters')
        set(altFig,'papertype','A4')
        ax = axes();
        set(ax,'Units','centimeter');
        set(ax,'Position', [2 2 7.9 3.1])
      end
    else
      saveName = [];
      set(handles.fig,'CurrentAxes',handles.plot);
    end

    goodIdx = find(sum(~isnan(histData),2));

    if(length(plotInfo(plotId).Xlim) > 1)
      rangeIdx = find(plotInfo(plotId).Xlim(1) <= histCenters ...
		      & histCenters <= plotInfo(plotId).Xlim(2));

      useIdx = intersect(goodIdx, rangeIdx);
    else
      useIdx = goodIdx;
    end

    pe = errorbar(histCenters(useIdx,:), ...
		  histData(useIdx,:), ...
		  histDataSEM(useIdx,:), ...
		  'linestyle', 'none', ...
		  'marker', 'none'); 

    setLineColor(pe, plotGroupID);

    hold on
    % Draw markers seperately so we get them on top
    p = plot(histCenters(useIdx,:), ...
	     histData(useIdx,:), ...
	     'o', ...
	     'linewidth',0.75, ...
	     'markersize',3);     

    setMarkerColor(p, plotGroupID);


    legend(p,graphLegend,'location','best')
    xlabel(plotInfo(plotId).Xlabel,'fontname','Arial','fontsize', 8)
    ylabel(plotInfo(plotId).Ylabel,'fontname','Arial','fontsize', 8)

    if(~isempty(saveName))
      title(sprintf('%s_%s', plotInfo(plotId).Name, saveName), ...
	    'fontname','Arial','fontsize', 8, 'interpreter', 'none')
    else
      title(sprintf('%s', plotInfo(plotId).Name), ...
	    'fontname','Arial','fontsize', 8, 'interpreter', 'none')
    end

    set(gca, 'fontsize', 7, 'fontname', 'Arial')
    set(gca, 'linewidth', 1, 'tickdir', 'out')
    set(gca,'ticklength',get(gca,'ticklength')*3);

    box off
    hold off

    halfBinWidth = (histCenters(2,1)-histCenters(1,1))/2;
    a = axis; 
    if(~isnan(halfBinWidth))
      a(1) = 0 - halfBinWidth; 
    end
    a(3) = 0;

    if(~isempty(plotInfo(plotId).Xlim))
      a(1:2) = plotInfo(plotId).Xlim;
    end

    if(~isempty(plotInfo(plotId).Ylim))
      a(3:4) = plotInfo(plotId).Ylim;
    end

    axis(a);
    setAxisGUI(a);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function makeProfileGraph(plotId,saveName,subplotFlag)

    allData = {};
    allCenters = {};
    graphLegend = {};

    [maxLenALL,maxIdxALL] = maxVectLen('meanSynapseProfileDist');

    % User selected groups to plot
    plotGroupID = getGroupsForPlot();

    profileCenter = NaN*zeros(maxLenALL, length(plotGroupID));
    profileData = NaN*zeros(maxLenALL, length(plotGroupID));
    profileDataSEM = NaN*zeros(maxLenALL, length(plotGroupID));

    for profNum = 1:length(plotGroupID)

      groupId = plotGroupID(profNum);
      memberIdx = group(groupId).Members;

      [maxLen,maxIdx] = maxVectLen('meanSynapseProfileDist',memberIdx);

      allData{profNum} = NaN*zeros(maxLen, length(memberIdx));
      allCenters{profNum} = NaN*zeros(maxLen, length(memberIdx));

      for i = 1:length(memberIdx)
	if(iscell(plotInfo(plotId).DataStr))
 	  tmp = eval(sprintf('data(memberIdx(i)).%s./data(memberIdx(i)).%s', ...
			     plotInfo(plotId).DataStr{1}, ...
			     plotInfo(plotId).DataStr{2}));
        else
 	  tmp = eval(sprintf('data(memberIdx(i)).%s', ...
			     plotInfo(plotId).DataStr));
        end

	if(~isempty(tmp))
          allData{profNum}(1:length(tmp),i) = tmp;
          allCenters{profNum}(1:length(tmp),i) ...
	    = data(memberIdx(i)).meanSynapseProfileDist;
        end

      end

      if(nnz(diff(allCenters{profNum},1,2) > 1e-9))
        fprintf('Error, different edges for group %s. Plot incorrect.\n', ...
		group(groupId).Name)
        beep
      end

      graphLegend{profNum} = group(groupId).Name;

      if(~isempty(allData{profNum}))
        profileCenter(:,profNum) ...
	  = allCenters{profNum}(:,find(memberIdx == maxIdx,1));
        profileData(:,profNum) = mean(allData{profNum},2, 'omitnan');
        nElem = sum(~isnan(allData{profNum}),2);

        profileDataSEM(:,profNum) = std(allData{profNum},[],2, 'omitnan') ...
                                    ./ sqrt(nElem);
      end

    end

    if(exist('saveName'))
      if(~exist('subplotFlag'))
        altFig = figure('visible','off');
        set(altFig,'paperunits','centimeters')
        set(altFig,'papertype','A4')
        ax = axes();
        set(ax,'Units','centimeter');
        set(ax,'Position', [2 2 3.3 3.1])
      end
    else
      saveName = [];
      set(handles.fig,'CurrentAxes',handles.plot);
    end

    if(~isempty(profileCenter))
      goodIdx = find(sum(~isnan(profileData),2));;
      halfBinWidth = 1e6*(profileCenter(2,1)-profileCenter(1,1))/2;


      if(length(plotInfo(plotId).Xlim) > 1)
        rangeIdx = find(plotInfo(plotId).Xlim(1) <= profileCenter*1e6 ...
			& profileCenter*1e6 <= plotInfo(plotId).Xlim(2));

        useIdx = intersect(goodIdx, rangeIdx);
      else
        useIdx = goodIdx;
      end

      pE = errorbar(profileCenter(useIdx,:)*1e6, ...
		    profileData(useIdx,:), ...
		    profileDataSEM(useIdx,:), ...
		    'linestyle','none', ...
		    'linewidth',1, ...
		    'marker','none', ...
		    'color', [0 0 0]);

      setLineColor(pE, plotGroupID);

      hold on
      p = plot(profileCenter(useIdx,:)*1e6, ...
	       profileData(useIdx,:), ...
	       'marker', 'o', ...
	       'linewidth',0.75, ...
	       'markersize',3);

      setMarkerColor(p, plotGroupID);
      setLineColor(p, plotGroupID);

      legend(p,graphLegend,'location','best')
    else
      p = [];
      pE = [];
      halfBinWidth = 0;
    end

    xlabel(plotInfo(plotId).Xlabel,'fontname','Arial','fontsize', 8)
    ylabel(plotInfo(plotId).Ylabel,'fontname','Arial','fontsize', 8)

    if(~isempty(saveName))
      title(sprintf('%s_%s', plotInfo(plotId).Name, saveName), ...
	    'fontname','Arial','fontsize', 8, 'interpreter', 'none')
    else
      title(sprintf('%s', plotInfo(plotId).Name), ...
	    'fontname','Arial','fontsize', 8, 'interpreter', 'none')
    end

    set(gca,'fontsize',7,'fontname','Arial')
    set(gca, 'linewidth', 1, 'tickdir', 'out')
    set(gca,'ticklength',get(gca,'ticklength')*3);

    a = axis; 
    a(1) = 0-halfBinWidth;
    a(3) = 0;

    if(~isempty(plotInfo(plotId).Xlim))
      a(1:2) = plotInfo(plotId).Xlim;
    end

    if(~isempty(plotInfo(plotId).Ylim))
      a(3:4) = plotInfo(plotId).Ylim;
    end

    axis(a);
    setAxisGUI(a);

    for i = 1:length(pE)
      setErrorbarWidth(pE(i));
    end
    setMarkerColor(p);

    box off
    hold off

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function makeGradientPlot(plotId, saveName, subplotFlag)

    allData = {};
    allEdges = {};

    [maxLenALL, maxIdxALL] = maxVectLen('gradientEdges');

    if(maxLenALL == 0)
      disp('No data for gradients.')
      return
    end

    % User selected groups to plot
    plotGroupID = getGroupsForPlot();

    gradientData = NaN*zeros(length(plotGroupID), maxLenALL);
    gradientSEM = NaN*zeros(length(plotGroupID), maxLenALL);
    gradientEdges = NaN*zeros(length(plotGroupID), maxLenALL);

    for gradNum = 1:length(plotGroupID)

      groupId = plotGroupID(gradNum);
      memberIdx = group(groupId).Members;

      [maxLen,maxIdx] = maxVectLen('gradientEdges',memberIdx);

      % Initiate memory
      allData{gradNum} = NaN*zeros(length(memberIdx), maxLen);
      allEdges{gradNum} = NaN*zeros(length(memberIdx), maxLen);

      for i = 1:length(memberIdx)

        if(iscell(plotInfo(plotId).DataStr))
          % We got a nominator and denominator to take care of
          tmpA = eval(sprintf('data(memberIdx(i)).%s', ...
			      plotInfo(plotId).DataStr{1}));
          tmpB = eval(sprintf('data(memberIdx(i)).%s', ...
			      plotInfo(plotId).DataStr{2}));

          tmpB(find(tmpB == 0)) = NaN; % Avoid division by zero
          tmpData = tmpA ./ tmpB;
        else
          tmpData = eval(sprintf('data(memberIdx(i)).%s', ...
  	 		         plotInfo(plotId).DataStr));
        end

        tmpEdges = data(memberIdx(i)).gradientEdges*1e6;

        nData = length(tmpData);
        nEdges = length(tmpEdges);

        allData{gradNum}(i,1:nData) = tmpData;	  
        allEdges{gradNum}(i,1:nEdges) = tmpEdges;

      end

      if(nnz(diff(allEdges{gradNum},1,1) > 1e-9))
        fprintf('Error, different edges for group %s. Plot incorrect.\n', ...
		group(groupId).Name)
        beep
      end

      graphLegend{gradNum} = group(groupId).Name;

      gradientEdges(gradNum,1:maxLen) = data(maxIdx).gradientEdges*1e6;
      gradientData(gradNum,1:maxLen) = mean(allData{gradNum},1, 'omitnan');
      gradientSEM(gradNum,1:maxLen) = ...
	std(allData{gradNum},[],1, 'omitnan') ./ sqrt(sum(~isnan(allData{gradNum})));

    end

    % Initiate the plot
    if(exist('saveName'))
      if(~exist('subplotFlag'))
        altFig = figure('visible','off');
        set(altFig,'paperunits','centimeters')
        set(altFig,'papertype','A4')
        ax = axes();
        set(ax,'Units','centimeter');
        set(ax,'Position', [2 2 3.3 3.1])
      end

      SEMalpha = 1; % Saving with transparency takes a looong time...
    else
      saveName = [];
      set(handles.fig,'CurrentAxes',handles.plot);
      SEMalpha = 0.3;
    end

    hold off   
    p = stairs(transpose(gradientEdges(:,1:end-1)), ...
	       transpose(gradientData(:,2:end)),'linewidth',1);

    setLineColor(p,plotGroupID);

    hold on

    % Points for the soma data
    pe = errorbar(gradientEdges(:,1),gradientData(:,1),gradientSEM(:,1),'k.');

    for i = 1:size(gradientEdges,1)
      p2(i,1) = plot(gradientEdges(i,1),gradientData(i,1),'.', 'markersize',20);
      set(p2(i,1),'color',get(p(i),'color'));
      graphLegend{end+1} = sprintf('%s soma', graphLegend{i});
    end

    % Mark standard error of the mean
    for i = 1:size(gradientEdges,1)
      [x1,y1] = stairs(gradientEdges(i,1:end-1), ...
		       gradientData(i,2:end)+gradientSEM(i,2:end));
      [x2,y2] = stairs(transpose(gradientEdges(i,1:end-1)), ...
		       transpose(gradientData(i,2:end) ...
				 -gradientSEM(i,2:end)));

      i1 = find(~isnan(x1) & ~isnan(y1));
      i2 = find(~isnan(x2) & ~isnan(y2));

      x1 = x1(i1); y1 = y1(i1);
      x2 = x2(i2); y2 = y2(i2);

      patchCol = get(p(i),'color');
      if(SEMalpha == 1)
        % In case of no transparency, just change shade slightly
        patchCol = patchCol*0.3 +0.7*[1 1 1];
      end
      p3(i,1) = patch([x1;x2(end:-1:1)],[y1;y2(end:-1:1)],patchCol);
      set(p3(i,1), 'facealpha', SEMalpha, 'edgealpha', SEMalpha);
      set(p3(i,1),'edgecolor',patchCol);
    end

    P = get(p(1),'parent');
    c = get(P,'children');

    % Move the SEM patches to the back (ie first in list of children)
    set(P,'children', [c(find(~ismember(c,p3))); p3]);

    hold off

    % We could add errorbars if we want.

    legend([p; p2],graphLegend,'location','best')

    xlabel(plotInfo(plotId).Xlabel,'fontsize',8,'fontname','arial')
    ylabel(plotInfo(plotId).Ylabel,'fontsize',8,'fontname','arial')

    
    if(~isempty(saveName))
      title(sprintf('%s_%s', plotInfo(plotId).Name, saveName), ...
	    'fontname','Arial','fontsize', 8, 'interpreter', 'none')
    else
      title(sprintf('%s', plotInfo(plotId).Name), ...
	    'fontname','Arial','fontsize', 8, 'interpreter', 'none')
    end

    set(gca,'fontsize',7,'fontname','Arial')
    set(gca, 'linewidth', 1, 'tickdir', 'out')
    set(gca,'ticklength',get(gca,'ticklength')*3);

    box off
    hold off

    a = axis;
    a(1) = 0;
    a(2) = data(maxIdxALL).gradientEdges(end)*1e6;
    a(3) = 0;

    if(~isempty(plotInfo(plotId).Xlim))
      a(1:2) = plotInfo(plotId).Xlim;
    end

    if(~isempty(plotInfo(plotId).Ylim))
      a(3:4) = plotInfo(plotId).Ylim;
    end

    axis(a);
    setAxisGUI(a);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setPlotTypeList()
    setPlotInfo();

    plotTypeList = {};

    for i = 1:length(plotInfo)
      plotTypeList{i} = plotInfo(i).Name;
    end

    set(handles.plotType,'String',plotTypeList)

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function setBarGraphDots(source, event)
    plotId = get(handles.plotType,'Value');

    plotInfo(plotId).barGraphDots = get(handles.barGraphDots,'Value'); 

    showPlot();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function selectGroups(source, event)

    showPlot();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function showPlot(source, event)

    if(isempty(data) | isempty(group(1).Members))
      disp('You need to load data and define groups to show plots.')
      return
    end
   
    plotId = get(handles.plotType,'Value');
    fprintf('Showing %s\n', plotInfo(plotId).Name)

    switch(plotInfo(plotId).Type)
      case 'bar'
	makeBarGraph(plotId);
      case 'cumhist'
        makeCumHist(plotId);
      case 'hist'
        makeHist(plotId);
      case 'sholl'
        makeHist(plotId);
      case 'profile'
        makeProfileGraph(plotId);
      case 'gradient'
        makeGradientPlot(plotId);
      otherwise
        fprintf('Unknown plot type %s\n', plotInfo(plotId).Type)
    end

    if(~isempty(plotInfo(plotId).barGraphDots))
      set(handles.barGraphDots, ...
	  'enable','on', ...
	  'value', plotInfo(plotId).barGraphDots)
    else
      set(handles.barGraphDots, ...
	  'enable','off', ...
	  'value', 0)
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function savePlotCallback(source, event)

    plotId = get(handles.plotType,'Value');

    [plotDir, plotExperimentName] = savePlot(plotId,plotDir,plotExperimentName);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function saveAllPlots(source, event)

    progressBar = waitbar(0,'Saving figures...');

    [plotDir, plotExperimentName] = savePlot(1,plotDir,plotExperimentName);
    waitbar(1/length(plotInfo),progressBar);

    for plotId = 2:length(plotInfo)
      [plotDir, plotExperimentName] = ...
	 savePlot(plotId, plotDir, plotExperimentName);
      waitbar(plotId/length(plotInfo),progressBar);
    end

    delete(progressBar);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [plotDir, plotExperimentName] = savePlot(plotId, ...
						    plotDir, ...
						    plotExperimentName)

    fprintf('Saving %s\n', plotInfo(plotId).Name)

    groupNameList = group(1).Name;
    for i=2:length(group)
      groupNameList = sprintf('%s-%s', groupNameList,group(i).Name);
    end

    % Get file type.
    fileType = saveFileType{get(handles.fileType,'Value')};

    if(~exist('plotDir') | isempty(plotDir))
      answer = inputdlg('Experiment name', ...
			'Name:', 1, ...
			{ groupNameList });

      plotExperimentName = answer{1};

      plotDir = uigetdir(exportPath);

    end

    switch(plotInfo(plotId).Type)
      case 'bar'
	makeBarGraph(plotId, plotExperimentName);
      case 'cumhist'
        makeCumHist(plotId, plotExperimentName);
      case 'hist'
        makeHist(plotId, plotExperimentName);
      case 'sholl'
        makeHist(plotId, plotExperimentName);
      case 'profile'
        makeProfileGraph(plotId, plotExperimentName);
      case 'gradient'
        makeGradientPlot(plotId, plotExperimentName);
      otherwise
        fprintf('Unknown plot type %s\n', plotInfo(plotId).Type)
        return
    end


    fullFigName = sprintf('%s/%s_%s_SynD.%s', ...
			  plotDir, ...
			  plotInfo(plotId).FigName, ...
			  plotExperimentName, fileType);

    switch(fileType)
      case 'pdf'
	saveas(gcf,fullFigName,'pdf');
      case 'eps'
	saveas(gcf,fullFigName,'psc2');
      case 'png'
	saveas(gcf,fullFigName,'png');
      otherwise
        if(~isempty(fileType))
          fprintf('Unknown file type: %s\n', fileType)
        end
    end

    set(handles.fig,'CurrentAxes',handles.plot);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function makeSummaryPlots(plotIDs, summaryFileName, width, height)

    groupNameList = group(1).Name;
    for i=2:length(group)
      groupNameList = sprintf('%s-%s', groupNameList,group(i).Name);
    end

    if(~exist('plotDir') | isempty(plotDir))
      answer = inputdlg('Experiment name', ...
			'Name:', 1, ...
			{ groupNameList });

      plotExperimentName = answer{1};

      plotDir = uigetdir(exportPath);

    end

    altFig = figure('visible','off');
    set(altFig,'paperunits','centimeters')
    set(altFig,'papertype','A4')

    for yPlot = 1:size(plotIDs,1)
      for xPlot = 1:size(plotIDs,2)

 	plotId = plotIDs(yPlot,xPlot);

        if(isnan(plotId))
          % Skip this plot
          continue
        end


        if(isempty(plotInfo(plotIDs(1)).Xlabel))
          wPad = 1.4;
        else
          wPad = 1.9;
        end

        hPad = 1.6;


        ax = axes();
        set(ax,'Units','centimeter');
        set(ax,'Position', ...
            [2+(xPlot-1)*(width+wPad), ...
             2+(size(plotIDs,1)-yPlot)*(height+hPad), ...
	     width, height])

        switch(plotInfo(plotId).Type)
          case 'bar'
     	    makeBarGraph(plotId, plotExperimentName,'subplot');
          case 'cumhist'
            makeCumHist(plotId, plotExperimentName,'subplot');
          case 'hist'
            makeHist(plotId, plotExperimentName,'subplot');
          case 'sholl'
            makeHist(plotId, plotExperimentName,'subplot');
          case 'profile'
            makeProfileGraph(plotId, plotExperimentName,'subplot');
          case 'gradient'
            makeGradientPlot(plotId, plotExperimentName,'subplot');
          otherwise
            fprintf('Unknown plot type %s\n', plotInfo(plotId).Type)
          end
      end
    end

    fileType = saveFileType{get(handles.fileType,'Value')};

    fullFigName = sprintf('%s/%s_%s_SynD.%s', ...
			  plotDir, ...
			  summaryFileName, ...
			  plotExperimentName, fileType);
    switch(fileType)
      case 'pdf'
	saveas(gcf,fullFigName,'pdf');
      case 'eps'
	saveas(gcf,fullFigName,'psc2');
      case 'png'
	saveas(gcf,fullFigName,'png');
      otherwise
        if(~isempty(fileType))
          fprintf('Unknown file type: %s\n', fileType)
        end
    end

    set(handles.fig,'CurrentAxes',handles.plot);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function makeAllSummaryPlots(source, event)

    progressBar = waitbar(0,'Saving summary figures...');

    makeSummaryPlots([2 1 3 4; 14 15 16 NaN],'summary-basics',3.3,3.1);
    waitbar(1/7,progressBar);

    makeSummaryPlots([34 35 36; 37 38 39],'summary-juxta',3.3,3.1);
    waitbar(2/7,progressBar);

    makeSummaryPlots([26 27; 28 29; ...
		      30 NaN],'summary-Sholl',7.9,3.1);
    waitbar(3/7,progressBar);

    makeSummaryPlots([31 32; 33 NaN],'summary-Sholl-ratios',7.9,3.1);
    waitbar(4/7,progressBar);

    makeSummaryPlots([5 6 7; 8 9 10; 23 24 25],'summary-intensity1',3.3,3.1);
    waitbar(5/7,progressBar);

    makeSummaryPlots([11 12 13; 17 18 19; 20 21 22],'summary-intensity2',3.3,3.1);
    waitbar(6/7,progressBar);

    makeSummaryPlots([40 41 42; 43 44 NaN],'summary-gradients',3.3,3.1);
    waitbar(7/7,progressBar);

    delete(progressBar);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function SynDversion = getVersion()

    try
     fid = fopen(versionFile,'r');
      SynDversion = fgets(fid);
      if(SynDversion(end) == char(10))
        SynDversion = str2num(SynDversion(1:end-1));
      else
	SynDversion = str2num(SynDversion);
      end
      fclose(fid);

      if(isempty(SynDversion))
	SynDversion = NaN;
      end
    catch
      SynDversion = NaN;
      disp('No version information available.')
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
