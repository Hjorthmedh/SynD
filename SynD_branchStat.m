% Script written by Johannes Hjorth for Kate Marler
%

function SynD_branchStat()

  close all

  % Post processing tool for SynD that investigates how close
  % synapses are to branch points.
  
  data = struct([]);
  
  detection.trimNeuriteSize = 12;  
  
  lazyLoad = 0; % Should we load a default test file, for faster testing...
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function loadData(source, event)
   
    if(lazyLoad & ~exist('source'))
      % Lazy loading while testing code
      %dataPath = '/Users/hjorth/matlab/DATA/Sabine/';
      %dataFile = '100316_nicecellzoomc_16xaverage.lsm-save.mat';
      
      dataPath = '/Users/hjorth/DATA/Kate/Cambridge/';
      dataFile = 'combine arl8 red arl 8 green phal 4.tif-save.mat';
    else
      [dataFile,dataPath] = uigetfile('*-save.mat');
    end
      
    if(dataFile == 0)
      disp('User pressed cancel.')
      return
    end
    
    tmp = load(strcat(dataPath,dataFile));
    
    data = tmp.old;
    data.height = size(data.image,1);
    data.width = size(data.image,2);
    
    % Calculate skeleton
    data.skeleton = bwmorph(data.neuriteMask-data.somaMask>0,'skel',inf);
    
    % Trim away tiny branches that do not belong
    trimSkeleton();
    synapseBranchDist();
    % synapseBranchDistOLD();
    
    intensityCorrelation();
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [dClosest,neighIdx] = closestNeighDist(pointA,pointB)
  
    if(~exist('pointB'))
      % We do self comparison
      pointB = pointA;
      clearDiag = 1;
    else
      clearDiag = 0;
    end
  
    xA = repmat(pointA(:,1),1,size(pointB,1));
    yA = repmat(pointA(:,2),1,size(pointB,1));
    xB = repmat(transpose(pointB(:,1)),size(pointA,1),1);
    yB = repmat(transpose(pointB(:,2)),size(pointA,1),1);

    d = sqrt((xA-xB).^2+(yA-yB).^2);
    
    if(clearDiag)
      d = d + eye(size(d))*max(d(:));
    end
    
    [dClosest,neighIdx] = min(d,[],2);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function branchPoints = findBranchPoints()

    % Locate all branch points
    branchPoints = bwmorph(data.skeleton,'branchpoints');    

    % Sometimes a branch is detected by multiple adjacent pixels. 
    % and other times small inhomogenities lead to multiple nearby
    % false branches.
    
    % Remove branch points that are too close together
    [yBran,xBran] = find(branchPoints);

    xB = repmat(xBran,1,length(yBran));
    yB = repmat(yBran,1,length(xBran));
    
    d = sqrt((xB-transpose(xB)).^2 + (yB-transpose(yB)).^2);

    cutOffDist = 4;
    [nA,nB] = find(0 < d & d < cutOffDist);
    
    try
    for i = 1:length(nA)
      if(data.image(yBran(nA(i)),xBran(nA(i)),data.morphChannel) ...
         > data.image(yBran(nB(i)),xBran(nB(i)),data.morphChannel))
        branchPoints(yBran(nB(i)),xBran(nB(i))) = 0;
      else
        branchPoints(yBran(nA(i)),xBran(nA(i))) = 0;
      end
    end
    catch e
      getReport(e)
      keyboard
    end

    endpoints = bwmorph(data.skeleton,'endpoints');
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Allow the user to specify the synapse centers, that way we can
  % test the null hypothesis with randomly placed synapses.
  
  function [synapseDistToBranch,isoSkelDistToBranch] = ...
      findClosestBranchPoint(synapseCenterIdx)
  
    % Find skeleton pixels
    skelIdx = find(data.skeleton);
  
    % Calculate arc length for every pixel in the skeleton
    distMask = makeDistMask(data.somaMask,data.neuriteMask);
    skelDist = distMask(skelIdx);
    
    % Remove branch points (and one pixel around them) from skeleton
    branchPoints = findBranchPoints();
    isolatedBranchMask = double(data.skeleton - bwmorph(branchPoints,'dilate') > 0);
    
    % Label the now isolated branches
    isolatedBranchMaskID = bwlabel(isolatedBranchMask);
    
    % Find the neurite pixel closest to every synapse (euclidean distance)
    [synY,synX] = ind2sub(size(data.synapseMask),synapseCenterIdx);
    
    % Find all the pixels of the isolated pieces, and their arc
    % length distance to the soma, and their ID.
    isolatedIdx = find(isolatedBranchMaskID);
    [isoY,isoX] = ind2sub(size(isolatedBranchMaskID),isolatedIdx);
    isolatedID = isolatedBranchMaskID(isolatedIdx);
    
    [dClosest, neighIdx] = closestNeighDist([synX,synY],[isoX,isoY]);
    
    if(max(dClosest) > 20)
      disp('WARNING, max(dClosest) is very large. Talk to Johannes!')
    end
    
    %%% Verify neighIDx
    
    if(1)
      figure
      imagesc(isolatedBranchMaskID)
      hold on
   
      cmap = colormap;
      [ySyn,xSyn] = ind2sub(size(data.synapseMask), ...
                            synapseCenterIdx);
    
      for i = 1:length(synapseCenterIdx)
        
        for j = 1:3
          c(j) = interp1(linspace(min(isolatedID(:)), ...
                                  max(isolatedID(:)),size(cmap,1)), ...
                         cmap(:,j), ...
                         isolatedID(neighIdx(i)));
        end
        
        plot(xSyn(i),ySyn(i),'*','color',c)
        
      end
    end
    
    
    % The isolated branch ID of every synapse, and their distance
    % to the soma.
    synapseID = isolatedBranchMaskID(isolatedIdx(neighIdx));
    synapseDistToSoma = distMask(isolatedIdx(neighIdx));
    
   
    % Find the branch pixel fatherest from the synapse (arc length)
    % this branch can be either closer, or further away from the
    % soma.
    synapseDistToBranch = NaN*zeros(size(synapseCenterIdx));
    endPoints = bwmorph(data.skeleton,'endpoints');
        
    for i = 1:length(synapseCenterIdx)
      idx = find(synapseID(i) == isolatedID);
      pixelIdx = isolatedIdx(idx);
      
      % Is the most distal point an end point?
      maxDist = max(distMask(pixelIdx));
      minDist = min(distMask(pixelIdx));

      switch(nnz(endPoints(pixelIdx)))
       case 0
        synapseDistToBranch(i) = min(maxDist - synapseDistToSoma(i), ...
                                     synapseDistToSoma(i) - minDist);
       case 1
        if(minDist <= data.xyRes) % This branch borders the soma
          synapseDistToBranch(i) = maxDist - synapseDistToSoma(i);          
        else
          synapseDistToBranch(i) = synapseDistToSoma(i) - minDist;
        end
       case 2
        % Two end points, no branch point. Ignoring.
        synapseDistToBranch(i) = NaN;
       otherwise
        synapseDistToBranch(i) = NaN;
        fprintf('Multiple nearby branch points, confused. %d pixels ignored\n', ...
                length(pixelIdx))
      end

      if(nnz(endPoints(pixelIdx)))
        % Segment contains an end point, only calculate distance to
        % proximal end
        synapseDistToBranch(i) = synapseDistToSoma(i) - minDist;
      else
        % The distal point is also a branch point
        synapseDistToBranch(i) = min(maxDist - synapseDistToSoma(i), ...
                                     synapseDistToSoma(i) - minDist);
      end
           
    end
        
    % Also for every pixel in the isolated branch, calculate the distance to
    % the closest branch point (this can be more proximal or more distal).
    allID = unique(isolatedBranchMaskID);
   
    isoSkelDistToBranch = NaN*zeros(size(isolatedIdx));
    
    for i = 1:length(allID)
      idx = find(allID(i) == isolatedID);
      pixelIdx = isolatedIdx(idx);
      
      if(nnz(~isnan(isoSkelDistToBranch(idx))))
        disp('This should never happen, internal error! Talk to Johannes.')
        keyboard
      end

      % Is the most distance point an end point?
      maxDist = max(distMask(pixelIdx));
      minDist = min(distMask(pixelIdx));

      switch(nnz(endPoints(pixelIdx)))
       case 0
        % No end points, do it normally
        isoSkelDistToBranch(idx) = min(maxDist - distMask(pixelIdx), ...
                                       distMask(pixelIdx) - minDist);                                
       case 1
        % One end point, is the soma falsely detected as end point
        if(minDist <= data.xyRes)
          % It was the soma
          isoSkelDistToBranch(idx) = maxDist - distMask(pixelIdx);
        else
          % No soma just one normal end point
          isoSkelDistToBranch(idx) = distMask(pixelIdx) - minDist;        
        end
       case 2
         % Two end points, one is the soma... so no branch
         isoSkelDistToBranch(idx) = NaN;
       otherwise
         fprintf('(Iter 2) Multiple nearby branch points, confused. %d pixels ignored\n', ...
                 length(pixelIdx))
         isoSkelDistToBranch(idx) = NaN;      
      end
      
    end

    %% !!! Plot membership ID for the pixels in the isolated
    %neurites....
    
    figure
    tmp = zeros(size(data.skeleton));
    tmp(isolatedIdx) = isoSkelDistToBranch;
    imagesc(tmp)
    axis equal
    colorbar
    hold on

    [bY,bX] = find(branchPoints);
    plot(bX,bY,'wo')

    [eY,eX] = find(bwmorph(data.skeleton,'endpoints'));
    plot(eX,eY,'wv')   
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function synapseBranchDist()

    [synapseDistToBranch,isoSkelDistToBranch] = ...
        findClosestBranchPoint(data.synapseCenter); 
    
    
    figure
    branchPoints = findBranchPoints();
    tmp = zeros(data.height,data.width,3);
    for i = 1:3
      tmp(:,:,i) = 0.2*data.neuriteMask + data.skeleton*0.5;
    end
      
    image(tmp);
    axis equal

    % Plot branch points as circles
    [bY,bX] = find(branchPoints);
    hold on
    plot(bX,bY,'bo')

    [eY,eX] = find(bwmorph(data.skeleton,'endpoints'));
    plot(eX,eY,'bv')
        
    % Plot synapses as coloured stars, coded by branch distance
    cmap = colormap('autumn');
    
    hold on
    [ySyn,xSyn] = ind2sub(size(data.synapseMask), ...
                          data.synapseCenter);

    for i = 1:length(data.synapseCenter)
      
      for j = 1:3
        c(j) = interp1(linspace(min(synapseDistToBranch), ...
                            max(synapseDistToBranch),size(cmap,1)), ...
                       cmap(:,j), ...
                       synapseDistToBranch(i));
      end
      
      plot(xSyn(i),ySyn(i),'*','color',c)

    end
    % colorbar % this didnt use the colours of the stars, only the
    % background image
    hold off
    
    figure
    maxEdge = ceil((max(union(synapseDistToBranch,isoSkelDistToBranch))*1e6)/10)*10;
    edges = 0:2:maxEdge;
    nSyn = histc(synapseDistToBranch*1e6,edges);
    nSkel = histc(isoSkelDistToBranch*1e6,edges);
    nSynSkel = nSyn./nSkel;
    
    subplot(3,1,1)
    bar(edges,nSyn)
    set(gca,'fontsize',12)
    ylabel('#Synapses','fontsize',14)
    title('Synapses','fontsize',14)
    subplot(3,1,2)
    bar(edges,nSkel)
    set(gca,'fontsize',12)
    ylabel('#Pixels','fontsize',14)
    title('Skeleton','fontsize',14)
    subplot(3,1,3)
    bar(edges,nSynSkel)
    set(gca,'fontsize',12)    
    ylabel('#Synapses/#Skeleton','fontsize',14)
    xlabel('Distance to closest branch (\mum)','fontsize',14)
    title('Normalised data','fontsize',14)
       
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % This old function calculated the euclidean distance to the
  % branch point, rather than the arc length distance.
  
  function synapseBranchDistOLD()
    % Find all branch points
    branchPoints = findBranchPoints();
    [yBran,xBran] = find(branchPoints);

    [ySkel,xSkel] = find(data.skeleton);
    
    fprintf('%d branch points, %d synapses\n', ...
            nnz(branchPoints), ...
            length(data.synapseCenter))
    
    % Find all synapse centers
    [ySyn,xSyn] = ind2sub(size(data.synapseMask),data.synapseCenter);
    
    % Calculate distance to closest branch point for every synapse
    dSyn = closestNeighDist([xSyn,ySyn],[xBran,yBran])*data.xyRes;
 
    % Calculate distance to closest branch point for every point in
    % skeleton
    dSkel = closestNeighDist([xSkel,ySkel],[xBran,yBran])*data.xyRes;
    
    figure
    tmp = data.neuriteMask + data.skeleton*5 + 10*branchPoints;
    tmp(data.synapseCenter) = 20;
    imagesc(tmp);
    figure
    maxEdge = ceil((max(union(dSyn,dSkel))*1e6)/10)*10;
    edges = 0:2:maxEdge;
    nSyn = histc(dSyn*1e6,edges);
    nSkel = histc(dSkel*1e6,edges);
    nSynSkel = nSyn./nSkel;
    
    subplot(3,1,1)
    bar(edges,nSyn)
    set(gca,'fontsize',12)
    ylabel('#Synapses','fontsize',14)
    title('Synapses','fontsize',14)
    subplot(3,1,2)
    bar(edges,nSkel)
    set(gca,'fontsize',12)
    ylabel('#Pixels','fontsize',14)
    title('Skeleton','fontsize',14)
    subplot(3,1,3)
    bar(edges,nSynSkel)
    set(gca,'fontsize',12)    
    ylabel('#Synapses/#Skeleton','fontsize',14)
    xlabel('Distance to closest branch (\mum)','fontsize',14)
    title('Normalised data','fontsize',14)
    
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % !! Lazy load...
  loadData();


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function [longestNeurite,nEndPoints] = longestNeurite()

    distMask = makeDistMask(data.somaMask,data.neuriteMask);
    endPoints = bwmorph(data.skeleton,'endpoints');
    
    longestNeurite = max(distMask(endPoints));
    nEndPoints = nnz(endPoints);
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Helper functions taken from SynD_extract.m
  
  function trimSkeleton(nIter)

    oldSkel = data.skeleton;
  
    if(~exist('nIter'))
      nIter = 3;
    else
      if(nIter < 1)
        return
      end
    end
  
    fprintf('Trimming skeleton: %d\n',nIter)
    
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
    
    % UPDATE --- Finally we recalulate the end points
    % and see if there are any end points neighbouring branch points
    % if so remove those end points.
    
    ep = bwmorph(data.skeleton,'endpoints');
    
    [epY,epX] = find(ep);
    [bpY,bpX] = find(branchPoints);
    
    d = sqrt((kron(epY,ones(1,length(bpY))) ...
              - kron(ones(length(epY),1),transpose(bpY))).^2 ...
             + (kron(epX,ones(1,length(bpX))) ...
              - kron(ones(length(epX),1),transpose(bpX))).^2);
          
    dMin = min(d,[],2);
    
    epNeigh = find(dMin <= sqrt(2));
    
    if(~isempty(epNeigh))
      fprintf('Found %d neighbours.\n', length(epNeigh))
      data.skeleton(sub2ind(size(data.skeleton),...
                    epY(epNeigh),epX(epNeigh))) = 0; 
    end
          
    % Remove non-endpoint neighbours
    data.skeleton = bwmorph(data.skeleton,'skel',inf);
    
    fprintf('Modifying %d pixels\n', nnz(oldSkel - data.skeleton))
      
    % Previously tried with multiple iterations
    % trimSkeleton(nIter-1);
    
  end  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Helper function taken from SynD_extract.m

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
  
  function intensityCorrelation()
    
    disp('Only analysing first image')
    synChan = data.image(:,:,data.synChannel);
    XChan = data.image(:,:,data.XChannel);

    maskIdx = find(data.neuriteMask);
    synChan = synChan(maskIdx);
    XChan = XChan(maskIdx);
    
    figure
    plot(synChan(:),XChan(:),'.k');
    xlabel('Synapse channel intensity')
    ylabel('XChannel intensity')
    
  end
  
  % keyboard  
  
end