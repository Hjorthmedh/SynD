%
% SynD - Imports tiff stacks where each series is divided into one image
%        per channel and Z-depth. The program reads in the channel info
%        by looking at the colour map of the tiff image.
%
% Johannes Hjorth, 2011
% Marielle Deurloo 
%

function SynD_importLeicaTiffStack(source, event)

  close all

  data.fileName = {};
  data.origImage = {};
  data.height = 0;
  data.width = 0;
  data.collapsedImage = [];
  data.channel = [];

  data.loadPath = pwd;

  configFile = 'SynD_stack_config.mat';
  loadConfig();

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  handles.fig = figure('Name','SynD collapse stack', ...
		       'MenuBar','none', ...
		       'Toolbar','figure', ...
		       'Position', [50 50 1150 680]);

  handles.image = axes('Units','Pixels', ...
		       'Position', [50 50 600 600]);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  handles.menuFile = uimenu(handles.fig,'Label','File');

  handles.menuItemLoad = uimenu(handles.menuFile, ...
				'Label','Load Tiff image(s)', ...
				'Interruptible','off', ...
				'Callback', @loadFile);

  handles.menuItemExport = uimenu(handles.menuFile, ...
				  'Label', 'Export collapsed image', ...
				  'Interruptible', 'off', ...
				  'Callback', @exportImage);

  handles.menuItemRestart = uimenu(handles.menuFile, ...
				   'Label','Restart', ...
				   'Interruptible', 'off', ...
				   'Callback', @SynD_importLeicaTiffStack);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  handles.menuDebug = uimenu(handles.fig,'Label','Debug');
  handles.menuItemDebug =  uimenu(handles.menuDebug, ...
                                      'Label','Keyboard', ...
                                      'Interruptible', 'off', ...
                                      'Callback', @keyboard);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function runDebug(source, event)
    disp('Type return to exit debug mode')
    keyboard
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function loadFile(source, event)

    data.tiffName = selectFile({'*.tif*','Tiff image'}, ...
			       'Select tiff image');

    if(~iscell(data.tiffName))
      data.tiffName = {data.tiffName};
    end

    for j = 1:length(data.tiffName)
      loadTiffFile(strcat(data.loadPath, data.tiffName{j}));
    end

    collapseImage();
    showImage();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function fileName = selectFile(fileMask, selectTitle)

    curPwd = pwd;
    try
      cd(data.loadPath);
    catch
      fprintf('Unable to change to %s\n', data.loadPath)
    end

    [fileName, filePath] = uigetfile(fileMask, ...
				     selectTitle, ...
				     'Multiselect','on');
    cd(curPwd);

    if(~iscell(fileName) & fileName == 0)
      % User pressed cancel
      fileName = [];
      return;
    end

    data.loadPath = filePath;

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function loadTiffFile(fileName)

    if(isempty(fileName))
      return
    end

    tiffInfo = imfinfo(fileName);

    if(isempty(data.origImage))
      data.height = tiffInfo.Height;
      data.width = tiffInfo.Width;
    elseif(data.height ~= tiffInfo.Height ...
 	   | data.width ~= tiffInfo.Width)
      fprintf('Dimension mismatch: H1: %d W1: %d, H2: %d W2: %d\n', ...
	      data.height, data.width, ...
	      tiffInfo.Height, tiffInfo.Width)
       return
    end

    data.fileName{end+1} = fileName;

    for i = 1:length(tiffInfo)
      tmp = imread(fileName,i);

      try
        fprintf('%s - channel %d\n', fileName, find(sum(tiffInfo(i).Colormap)))

        data.channel(end+1) = find(sum(tiffInfo(i).Colormap));
      catch e
        getReport(e)
        disp('This type of Tiff image is not supported, talk to Johannes.')
        disp('Type dbquit to exit debug mode')
        keyboard
      end

      data.origImage{end+1} = tmp;

    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function newImage = expandImage(oldImage)

    newImage = zeros(size(oldImage,1),size(oldImage,2),3,'uint16');

    for i = 1:size(oldImage,3)
      newImage(:,:,i) = oldImage(:,:,i);
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function showImage()

    set(handles.fig,'CurrentAxes',handles.image);

    if(size(data.collapsedImage,3) == 1)
      imagesc(data.collapsedImage)
    else
      image(normaliseImage(expandImage(data.collapsedImage)));
    end 

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function normImg = normaliseImage(img)
    normImg = zeros(size(img),'double');
    for ci = 1:size(img,3)
      normImg(:,:,ci) = double(img(:,:,ci)) / double(max(max(img(:,:,ci))));
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function collapseImage()

    uChan = unique(data.channel);

    data.collapsedImage = zeros(data.height,data.width, ...
				length(uChan), ...
				'uint16');

    for i = 1:length(uChan)
      idx = find(data.channel == uChan(i));

      for j = 1:length(idx)
        data.collapsedImage(:,:,i) = max(data.collapsedImage(:,:,i), ...
					 uint16(data.origImage{idx(j)}));
      end
    end

  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function exportImage(source, event)

    if(isempty(data.collapsedImage))
      disp('No data loaded')
      return
    end

    curPwd = pwd;
    cd(data.loadPath)

    exportFile = strrep(data.fileName{1},'.tif','-stacked.tif');
    if(length(exportFile) == length(data.fileName{1}))
      disp('Oops, making sure we do not overwrite original data.')
	exportFile = strcat(data.fileName{1},'-stacked.tif');
    end

    [exportFile, exportPath] = ...
	  uiputfile('*.tif', ...
		    'Export stacked image', ...
		    exportFile);

    cd(curPwd);

    if(isempty(exportFile) | exportFile == 0)
      exportFile = [];
      disp('Export aborted.')
      return
    end

    fprintf('Writing %s.\n',exportFile)
    imwrite(expandImage(data.collapsedImage),strcat(exportPath,exportFile), ...
	    'tif','compression','lzw');

    saveConfig();

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function loadConfig()
    if(exist(configFile))
      fprintf('Loading configuration from %s\n', configFile)
      cfg = load(configFile);

      try
        data.loadPath = cfg.old.loadPath;
      catch e
        getReport(e);
        disp('Unable to read config file.')

      end

    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function saveConfig()

    fprintf('Saving old configuration to %s\n', configFile)
    old.loadPath = data.loadPath;

    try 
      save(configFile, 'old');
    catch
      fprintf('Unable to save config file %s.\n', configFile)
    end

  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
