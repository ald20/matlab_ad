function [results] = ShadowingGeometry(V,F,FN)
%Generate 'above horizon and facing away' map for 

%% Initiate tables and constants


sizeFN = size(FN);
sizeV = size(V);

vismap = zeros(sizeFN(1));
logic = true(sizeFN(1));

Dtab=zeros(sizeFN(1));
Zangle=zeros(sizeFN(1));

Vcentr=(V(F(:,1),:)+V(F(:,2),:)+V(F(:,3),:))/3;

%% Calculate the 'above horizon and facing awawy' map 

for point = 1:sizeFN(1)
    
vismap(point,:) = zeros(sizeFN(1),1);
vismap(point,point) = 4;


Vfac=repmat(V(F(point,1),:),[sizeFN(1) 1]);

diffN = sqrt(sum((Vcentr-Vfac).^2,2));

Dtab(point,:) = (FN(point,:)*(Vcentr-Vfac)')';
Zangle(point,:) =  Dtab(point,:)'./diffN;


logic(point,:)=(Dtab(point,:)>0);

Vc= Vcentr(point,:);
FNp=FN(point,:);

P0(:,:,point) = repmat(Vcentr(point,:),sizeFN(1),1)-V(F(:,1),:); 


for i = 1:sizeFN
    

    
    if (~logic(point,i))
        continue
    end
        Vci=Vcentr(i,:);
        FNi = FN(i,:);
        viewfac = dot(FNp,Vci-Vc)*dot(FNi,Vc-Vci);
        if (viewfac>=0)
            logic(point,i)=false;
        end

end

end

%% Calculate v and u's for all facets 

v = V(F(:,2),1:3) - V(F(:,1),1:3);
u = V(F(:,3),1:3) - V(F(:,1),1:3);
uv = dot(u,v,2);
uu = dot(u,u,2);
vv = dot(v,v,2);
denom = uv.^2-uu.*vv;


%% Output shape parameters

results.Vcentr=Vcentr;

results.vismap=vismap;
results.logic=logic;

results.Dtab=Dtab;
results.Zangle=Zangle;
results.v=v;
results.u=u;
results.uv=uv./denom;
results.uu=uu./denom;
results.vv=vv./denom;
results.denom=1.0./denom;
results.P0=P0;

end