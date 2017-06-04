function SynD_alignTiff(source, event)

  close all, format compact
  addpath(strcat(pwd,'/LSM'));

  data.height = 0;
  data.width = 0;

  data.fileName = {};
  data.origImage = [];
  data.alignedImage = [];
  data.numChan = 0;

  data.loadPath = pwd;

  detection.alignOffset = zeros(3,2);
  detection.alignState = 'Auto';

  dispInfo.showRed = 1;
  dispInfo.showGreen = 1;
  dispInfo.showBlue = 1;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  handles.fig = figure('Name','SynD Tiff alignment', ...
		       'MenuBar','none', ...
		       'Toolbar','figure', ...
		       'Position', [50 50 1150 680]);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  handles.image = axes('Units','Pixels', ...
		       'Position', [50 50 600 600]);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  handles.alignButtons = uibuttongroup('Units','Pixels', ...
				       'Position', [700 500 200 110], ...
				       'SelectionChangeFcn', @selectAlign, ...
				       'Interruptible','off', ...
				       'SelectedObject', []);
			       
  handles.autoAlign = uicontrol('Style','radio', ...
				'Parent', handles.alignButtons, ...
				'String', 'Auto align', ...
				'Position', [20 70 150 20]);

  handles.align2 = uicontrol('Style','radio', ...
				 'Parent', handles.alignButtons, ...
				 'String', 'Align channel 2', ...
				 'Position', [20 40 150 20]);

  handles.align3 = uicontrol('Style','radio', ...
			       'Parent', handles.alignButtons, ...
			       'String', 'Align channel 3', ...
			       'Position', [20 10 150 20]);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  handles.red = uicontrol('Style', 'checkbox', ...
			  'Position', [716 450 20 20], ... 
			  'BackgroundColor', 0.8*[1 1 1], ...
			  'TooltipString', ...
			  'Show red channel', ...
			  'Value', 1, ...
			  'Callback', @toggleRed);

  handles.green = uicontrol('Style', 'checkbox', ...
			    'Position', [770 450 20 20], ...
			    'BackgroundColor', 0.8*[1 1 1], ...
			    'TooltipString', ...
			    'Show green channel', ...
			    'Value', 1, ...
			    'Callback', @toggleGreen);

  handles.blue = uicontrol('Style', 'checkbox', ...
			   'Position', [820 450 20 20], ... 
			   'BackgroundColor', 0.8*[1 1 1], ...
			   'TooltipString', ...
			   'Show blue channel', ...
			   'Value', 1, ...
			   'Callback', @toggleBlue);

  handles.redLabel = uicontrol('Style','text', ...
			       'String', 'R', ...
			       'Position', [700 450 20 20], ... 
			       'BackgroundColor', 0.8*[1 1 1], ...
			       'FontSize', 12, ...
			       'HorizontalAlignment', 'left');

  handles.greenLabel = uicontrol('Style','text', ...
				 'String', 'G', ...
				 'Position', [750 450 20 20], ... 
				 'BackgroundColor', 0.8*[1 1 1], ...
				 'FontSize', 12, ...
				 'HorizontalAlignment', 'left');

  handles.blueLabel = uicontrol('Style','text', ...
				'String', 'B', ...
				'Position', [800 450 20 20], ... 
				'BackgroundColor', 0.8*[1 1 1], ...
				'FontSize', 12, ...
				'HorizontalAlignment', 'left');


  handles.infoText = uicontrol('Style','text', ...
			       'String', 'When in manual editing mode, use arrows to align the channel. Shift+arrow moves 5 pixels instead of 1 pixel.', ...
			       'Position', [700 370 200 100], ... 
			       'BackgroundColor', 0.8*[1 1 1], ...
			       'FontSize', 10, ...
			       'HorizontalAlignment', 'left');


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  set([handles.fig, handles.image, ...
       handles.alignButtons, ...
       handles.autoAlign, ...
       handles.align2, ...
       handles.align3], 'units','normalized')

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  handles.menuFile = uimenu(handles.fig,'Label','File');

  handles.menuItemLoad = uimenu(handles.menuFile, ...
				'Label','Load Tiff image(s)', ...
				'Interruptible','off', ...
				'Callback', @loadFile);

  handles.menuItemAlign = uimenu(handles.menuFile, ...
				 'Label','Align image', ...
				 'Interruptible','off', ...
				 'Callback', @alignImages);


  handles.menuItemExportAligned = uimenu(handles.menuFile, ...
					 'Label', 'Export aligned image', ...
					 'Interruptible', 'off', ...
					 'Callback', @exportAlignedImage);

  handles.menuItemRestart = uimenu(handles.menuFile, ...
				   'Label','Restart (to do another stack)', ...
				   'Interruptible','off', ...
				   'Callback', @SynD_alignTiff);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  handles.menuDebug = uimenu(handles.fig,'Label','Debug');
  handles.menuItemDebug =  uimenu(handles.menuDebug, ...
                                      'Label','Keyboard', ...
                                      'Interruptible', 'off', ...
                                      'Callback', @runDebug);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function runDebug(source, event)
    disp('Type return to exit debug mode')
    keyboard
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

  function toggleGreen(source, event)
    if(~exist('source') | source ~= handles.green)
      set(handles.green,'Value',~get(handles.green,'Value'));
    end

    dispInfo.showGreen = get(handles.green,'Value');
    showImage();
  end

  function toggleBlue(source, event)
    if(~exist('source') | source ~= handles.blue)
      set(handles.blue,'Value',~get(handles.blue,'Value'));
    end

    dispInfo.showBlue = get(handles.blue,'Value');
    showImage();
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function loadTiffFile(fileName)

    if(isempty(fileName))
      tiffData = [];
      return
    end

    tiffInfo = imfinfo(fileName);
    data.height = tiffInfo.Height;
    data.width = tiffInfo.Width;

    data.fileName{end+1} = fileName;

    for i = 1:length(tiffInfo)
      tmp = imread(fileName,i);

      if(isempty(data.origImage))
        data.origImage = tmp;
      else
        data.origImage(:,:,end+1:end+size(tmp,3)) = tmp;
      end
    end

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
				     'Multiselect','off');
    cd(curPwd);

    if(~iscell(fileName) & fileName == 0)
      % User pressed cancel
      fileName = [];
      return;
    end

    data.loadPath = filePath;

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function loadFile(source, event)

    data.tiffName = selectFile({'*.tif*','Tiff image'}, ...
			       'Select tiff image');

    loadTiffFile(strcat(data.loadPath, data.tiffName));

    alignImages();
    set(handles.alignButtons,'SelectedObject',handles.autoAlign)

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Using FFT to align images
  %
  % http://stackoverflow.com/questions/2050510/matlabs-fminsearch-function
  %

  function alignImages(source, event)

    if(size(data.origImage,3) <= 1)
      disp('You need to load at least two images.')
      return
    end

    tic

    data.alignedImage = zeros(size(data.origImage),'uint16');
    data.alignedImage(:,:,1) = data.origImage(:,:,1);

    for i = 2:size(data.origImage,3)

      [output Greg] = dftregistration(fft2(data.origImage(:,:,1)), ...
				      fft2(data.origImage(:,:,i)), ...
				      1);

      data.alignedImage(:,:,i) = uint16(abs(ifft2(Greg)));
      data.alignOffset(i,:) = output(3:4);

    end

    toc

    showImage();

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

    if(isempty(data.alignedImage))
      return
    end

    tmp = zeros(data.height,data.width,3);

    if(dispInfo.showRed & size(data.alignedImage,3) >= 1)
      tmp(:,:,1) = data.alignedImage(:,:,1);
    end

    if(dispInfo.showGreen & size(data.alignedImage,3) >= 2)
      tmp(:,:,2) = data.alignedImage(:,:,2);
    end

    if(dispInfo.showBlue & size(data.alignedImage,3) >= 3)
      tmp(:,:,3) = data.alignedImage(:,:,3);
    end

    set(handles.fig,'CurrentAxes',handles.image);
    image(normaliseImage(tmp))

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function normImg = normaliseImage(img)
    normImg = zeros(size(img),'double');
    for ci = 1:size(img,3)
      normImg(:,:,ci) = double(img(:,:,ci)) / double(max(max(img(:,:,ci))));
    end

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function selectAlign(source, event)

    switch(get(handles.alignButtons,'SelectedObject'))
      case handles.autoAlign
        detection.alignState = 'Auto';
	set(handles.fig,'WindowKeyPressFcn',[]);
	alignImages();
      case handles.align2
        detection.alignState = 'Chan2'; % Morphology marker
	set(handles.fig,'WindowKeyPressFcn',@alignKeyHandler);
      case handles.align3
        detection.alignState = 'Chan3'; % Synapse marker
	set(handles.fig,'WindowKeyPressFcn',@alignKeyHandler);
      otherwise
        disp('Unknown align option.')
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function alignKeyHandler(source, event)

    if(size(data.origImage,3) <= 1)
      disp('Need at least 2 channels to align')
      return
    end

    dx = 0;
    dy = 0;

    switch(event.Key)
      case 'leftarrow'
	dx = -1;
      case 'rightarrow'
        dx = 1;
      case 'uparrow'
        dy = -1;
      case 'downarrow'
        dy = 1;
    end

    % Hold down shift to make it go 5x as fast
    if(strcmp(event.Modifier,'shift'))
      dx = dx * 5;
      dy = dy * 5;
    end

    switch(detection.alignState)
      case 'Chan2'
	detection.alignOffset(2,:) = detection.alignOffset(2,:) + [dx dy];
      case 'Chan3'
        detection.alignOffset(3,:) = detection.alignOffset(3,:) + [dx dy];
    end

    data.alignedImage = zeros(data.height,data.width,3,'uint16');
       
    data.alignedImage(:,:,1) = data.origImage(:,:,1);

    for i = 2:min(size(data.origImage,3),3)
      data.alignedImage(:,:,i) = shiftImage(data.origImage(:,:,i), ...
					    detection.alignOffset(i,:));
    end

    set(handles.fig,'CurrentAxes',handles.image);
    a = axis;
    showImage();
    axis(a);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function newImg = shiftImage(oldImg,shift)

    wShift = shift(1);
    hShift = shift(2);

    if(hShift < 0)
      hIdxOld = (abs(hShift)+1):size(oldImg,1);
      hIdxNew = 1:(size(oldImg,1)-abs(hShift));
    else
      hIdxOld = 1:(size(oldImg,1)-abs(hShift));
      hIdxNew = (abs(hShift)+1):size(oldImg,1);
    end

    if(wShift < 0)
      wIdxOld = (abs(wShift)+1):size(oldImg,2);
      wIdxNew = 1:(size(oldImg,2)-abs(wShift));
    else
      wIdxOld = 1:(size(oldImg,2)-abs(wShift));
      wIdxNew = (abs(wShift)+1):size(oldImg,2);
    end

    newImg = oldImg; % Make sure we retain same data type
    newImg(:) = 0;
    newImg(hIdxNew,wIdxNew) = oldImg(hIdxOld,wIdxOld);

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function exportAlignedImage(source, event)

    if(isempty(data.alignedImage))
      disp('All images are not loaded')
      return
    end

    curPwd = pwd;
    cd(data.loadPath)

    exportFile = strrep(data.fileName{1},'.tif','-aligned.tif');
    if(length(exportFile) == length(data.fileName{1}))
      disp('Oops, making sure we do not overwrite original data.')
	exportFile = strcat(data.fileName{1},'-aligned.tif');
    end

    [exportFile, exportPath] = ...
	  uiputfile('*.tif', ...
		    'Export aligned images', ...
		    exportFile);

    cd(curPwd);

    if(isempty(exportFile) | exportFile == 0)
      exportFile = [];
      disp('Export aborted.')
      return
    end

    fprintf('Writing %s.\n',exportFile)
    imwrite(expandImage(data.alignedImage),strcat(exportPath,exportFile), ...
	    'tif','compression','lzw');

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
