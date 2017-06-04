function Melina_branchStats()

  detection.trimNeuriteSize = 12;  
  
  try
    curPath = pwd();
    old = load('SynD_config.mat');
    filePath = old.old.exportPath;
    cd(filePath)
  catch
    disp('Unable to load SynD configuration file')
  end
  
  [fileName,filePath] = uigetfile('Select a mat file to load','*.mat');
  cd(curPath);
  
  tmp = load(strcat(filePath,fileName));
  data = tmp.old;
  data.height = size(data.image,1);
  data.width = size(data.image,2);
  clear tmp;

  data.skeleton = bwmorph(data.neuriteMask-data.somaMask>0,'skel',inf);
  trimSkeleton();
  
  distMask = makeDistMask(data.somaMask,data.neuriteMask);
  endPoints = bwmorph(data.skeleton,'endpoints');

  % Remove the end points that are linking to the soma
  pad = strel('disk',10);
  endPoints = endPoints & ~imdilate(data.somaMask,pad);

  % Remove end points that are too close together
  minDist = 10;
  
  [ey,ex] = find(endPoints);
  for i = 1:numel(ex)
    for j = i+1:numel(ex)
      d = sqrt((ex(i)-ex(j))^2 + (ey(i)-ey(j))^2);
      if(d < minDist)
        if(distMask(ey(i),ex(i)) < distMask(ey(j),ex(j)))
          endPoints(ey(i),ex(i)) = 0;
        else
          endPoints(ey(j),ex(j)) = 0;        
        end
      end
    end
  end
  
  [ey,ex] = find(endPoints);

  longestNeuriteLength = max(distMask(endPoints));
  nEndPoints = nnz(endPoints);
    
  fprintf('Longest neurite is %f micrometer\n', ...
          longestNeuriteLength*1e6)
  fprintf('%d end points found.\n', nEndPoints)
  
  figure
  imagesc(data.neuriteMask + data.somaMask + data.skeleton + endPoints*2)
  hold on
  plot(ex,ey,'wo')
  
  title(sprintf('%d end points, longest neurite %.3f micrometer', ...
                nEndPoints, longestNeuriteLength*1e6))
  
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
 
  
end