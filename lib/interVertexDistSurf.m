function [Mean SD Dist] = interVertexDist(surf1)

% Mean, SD: Global mean and stadard deviation of intervertex distance
% intervertex distance on each vertex
Dist=0;
for i=1:length(surf1.coord(1,:))
    Tri_1=find(surf1.tri(:,1)==i);
    Tri_2=find(surf1.tri(:,2)==i);
    Tri_3=find(surf1.tri(:,3)==i);
    
    Dist11=dist3(surf1.coord(:,i)', surf1.coord(:,surf1.tri(Tri_1,2))');
    Dist12=dist3(surf1.coord(:,i)', surf1.coord(:,surf1.tri(Tri_1,3))');
    Dist21=dist3(surf1.coord(:,i)', surf1.coord(:,surf1.tri(Tri_2,1))');
    Dist22=dist3(surf1.coord(:,i)', surf1.coord(:,surf1.tri(Tri_2,3))');
    Dist31=dist3(surf1.coord(:,i)', surf1.coord(:,surf1.tri(Tri_3,1))');
    Dist32=dist3(surf1.coord(:,i)', surf1.coord(:,surf1.tri(Tri_3,2))');
    Dist(i)=(sum(Dist11)+sum(Dist12)+sum(Dist21)+sum(Dist22)+sum(Dist31)+sum(Dist32)) ...
        / (length(Dist11)+length(Dist12)+length(Dist21)+length(Dist22)+length(Dist31)+length(Dist32));
    if(mod(i,4000)==0)
        i/400
    end
end

Mean=mean(Dist)
SD=std(Dist)
% Dist=SurfStatSmooth(Dist, surf1, 20);
% SurfStatView(Dist, surf1);
% SurfStatColLim([0,3]);

% h=gcf
% exportfigbo(h, fileimage, 'png', 10);
