function [ FNo, FNAo] = AR_calcFN_wrl(V,F)
% Calculate facet normals and facet normal areas for any set of 
% vertices and facets (useful for shape scaling)
% Necessary because .wrl shape files index from 0 instead of 1
% 
% Input:
%  V  -- number of vertices x 3 array of vertex positions
%  F  -- number of faces x 3 array of face indices

% Output:
% FNA: Facet normals
% FNAo: Facet normal areas


sizeF=size(F);
for ind=1:sizeF(1)
    v1=V(F(ind,1)+1,:);
    v2=V(F(ind,2)+1,:);
    v3=V(F(ind,3)+1,:);
    a=v1-v3;
    b=v2-v3;
    FN(ind,:)=cross(a,b);
    FNA(ind)=norm(FN(ind,:))/2;
    FN(ind,:)=FN(ind,:)/(2*FNA(ind));
    
end

%result.FN=FN;
%result.FNA=FNA';

FNo = FN;
FNAo = FNA';
end
