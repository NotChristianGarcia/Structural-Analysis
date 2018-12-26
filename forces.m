function [Ps,F,Sigma] = forces(nnode,node_coor,nelem,elemdata,ndisp,dispdata,U,M,K,P)
% P = KU
% M = Penalty boy
% Sigma = Stress = Force/Area
% Force = EA/L * displacement
	Pdir = K*U;
	Pdir = Pdir - P;
	F = zeros(nelem,1);
	Sigma = zeros(nelem,1);
	unitV = zeros(nelem,2); % 1 = i, 2 = j
	Uij = zeros(nelem,2);
	Ps = zeros(nnode,2);
	for node = 1:nnode
		Ps(node,:) = [Pdir(2*node-1,1),Pdir(2*node,1)];
	end
	
	for elem = 1:nelem
		x1 = node_coor(elemdata(elem,1),1);
		y1 = node_coor(elemdata(elem,1),2);
		x2 = node_coor(elemdata(elem,2),1);
		y2 = node_coor(elemdata(elem,2),2);
		
		ui = U(elemdata(elem,1)*2 - 1,1);
		vi = U(elemdata(elem,1)*2,1);
		uj = U(elemdata(elem,2)*2 - 1,1);
		vj = U(elemdata(elem,2)*2,1);
		
		EA = elemdata(elem,4)*elemdata(elem,3);
		L =  (12*sqrt((x2-x1)^2 + (y2-y1)^2));
		EAdL = EA/L; % in inch
		
		hypo = sqrt((x2-x1)^2 +(y2-y1)^2);
		unitV(elem,:) = [(x2-x1)/hypo,(y2-y1)/hypo];
		
		Uij(elem,:) = [uj-ui,vj-vi];
		
		D = dot(Uij(elem,:),unitV(elem,:));
		
		F(elem) = D * EAdL;
		Sigma(elem) = F(elem)/elemdata(elem,3);
	end
end