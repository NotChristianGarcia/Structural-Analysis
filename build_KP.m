function [K,P] = build_KP(nnode,node_coor,nelem,elemdata,nforce,forcedata)
    % Outline % D.O.F 1 = x and 2 = y
    % Create global & local K's
    % Create EA/L Matrix
    % Create connectivity Matrix
    % Throw it all in global K
    K = zeros(2*nnode); % making global K matrix with size of nnode*2 by nnode*2 (for x and y of each node)
	P = zeros(2*nnode,1);
    localK = zeros(4,4,nelem); % local, by element k matrices
    EAdL = zeros(nelem,1); % EA divided by L matrix
    locConnect = zeros(nelem,4); % connectivity matrix for each elem.
    gloConnect = zeros(nnode,2); % global connectivity, each node with 2 layers for x/y (i/j) directions
	
	for node = 1:nnode
		gloConnect(node,:) = (2*node-1):2*node;
	end
	
	for iterateForce = 1:nforce
		P((2*forcedata(iterateForce,1)-2+forcedata(iterateForce,2)),1) = forcedata(iterateForce,3);
	end
	
	for elem = 1:nelem
		x1 = node_coor(elemdata(elem,1),1);
		y1 = node_coor(elemdata(elem,1),2);
		x2 = node_coor(elemdata(elem,2),1);
		y2 = node_coor(elemdata(elem,2),2);
		EA = elemdata(elem,4)*elemdata(elem,3);
		L =  (12*sqrt((x2-x1)^2 + (y2-y1)^2));
		EAdL = EA/L; % in inch
		theta = atan((y2-y1)/(x2-x1)); % arctan of opposite length over adjacent length
		
		localK(:,:,elem) = [cos(theta)^2, cos(theta)*sin(theta), -cos(theta)^2, -cos(theta)*sin(theta); ...
							cos(theta)*sin(theta), sin(theta)^2, -cos(theta)*sin(theta), -sin(theta)^2; ...
							-cos(theta)^2, -cos(theta)*sin(theta), cos(theta)^2, cos(theta)*sin(theta); ...
							-cos(theta)*sin(theta), -sin(theta)^2, cos(theta)*sin(theta), sin(theta)^2;];
		localK(:,:,elem) = localK(:,:,elem) * EAdL;
		
		locConnect(elem,:) = [gloConnect(elemdata(elem,1),:),gloConnect(elemdata(elem,2),:)];
		
		for localki = 1:4
			for localkj = 1:4
				K(locConnect(elem,localki),locConnect(elem,localkj)) = ...
					K(locConnect(elem,localki),locConnect(elem,localkj)) + localK(localki,localkj,elem);
			end
		end
	end
end