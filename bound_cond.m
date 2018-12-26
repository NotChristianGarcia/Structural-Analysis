function [K,P] = bound_cond(K,P,ndisp,dispdata,M)
	nnode = size(K,1);
	gloConnect = zeros(nnode,2); % global connectivity, each node with 2 layers for x/y (i/j) directions
	for node = 1:nnode
		gloConnect(node,:) = (2*node-1):2*node;
	end
	
	for iterateDisp = 1:ndisp
		K(gloConnect(dispdata(iterateDisp,1),dispdata(iterateDisp,2)),gloConnect(dispdata(iterateDisp,1),dispdata(iterateDisp,2))) ...
			= K(gloConnect(dispdata(iterateDisp,1),dispdata(iterateDisp,2)),gloConnect(dispdata(iterateDisp,1),dispdata(iterateDisp,2))) + M;
		P((2*dispdata(iterateDisp,1)-2+dispdata(iterateDisp,2)),1) = dispdata(iterateDisp,3)*M;
	end
end