function  surf_outskel = MakeVertexCorrespond_surf_skel(surf1, SegLaplVol, DispXVol, DispYVol, DispZVol)
surf_orig=surf1;
SegLapl=SurfStatReadVol1(SegLaplVol);
DispX=SurfStatReadVol1(DispXVol);
DispY=SurfStatReadVol1(DispYVol);
DispZ=SurfStatReadVol1(DispZVol);

% Initialization
vert_ROI=1:length(surf1.coord(1,:));
one = ones(length(surf1.coord(1,:)),1);
Vol = DispX;
surf1.edg=SurfStatEdg(surf1);
length(vert_ROI)
iter=1;
deform_index=ones(1,length(vert_ROI));
Edg = SurfStatEdg(surf1);
% Deformation using laplacian vector fields
while sum(deform_index) > length(deform_index)/5 && iter < 75
    vert_ROI2 = vert_ROI(deform_index>0.5);
    length(vert_ROI2);
    s2volcoord = ((double(surf1.coord')-one*Vol.origin)./(one*Vol.vox(1:3))+1)';
%     s2volcoord_tmp = s2volcoord;
    DefX = interpn(DispX.data, s2volcoord(1,:), s2volcoord(2,:), s2volcoord(3,:), 'linear',0);
    DefY = interpn(DispY.data, s2volcoord(1,:), s2volcoord(2,:), s2volcoord(3,:), 'linear',0);
    DefZ = interpn(DispZ.data, s2volcoord(1,:), s2volcoord(2,:), s2volcoord(3,:), 'linear',0);
    Def = [ DefX' DefY' DefZ'];
    
    if iter == 1
        Def_old =Def;
    end
    if iter > 20
        if mod(iter,3) == 0 %iter == 6 || iter == 9 || iter == 12  || iter == 15 || iter == 18 || iter == 21 || iter == 24
            Curv = SurfStatComputeCurv_abs(surf1, 1:length(surf1.coord(1,:)));
            Curv1 = SurfStatSmooth(Curv, surf1, 3);
            vertEdge = find(Curv1<0.07);    
        end
        surf_tmp=surf1;
        surf_tmp.coord(:,vertEdge)=surf_tmp.coord(:,vertEdge) + (Def(vertEdge,:)' * 0.2);
    else
        surf_tmp=surf1;
        surf_tmp.coord=surf_tmp.coord + (Def' * 0.2);
        
    end    
    s2volcoord = ((double(surf_tmp.coord')-one*Vol.origin)./(one*Vol.vox(1:3))+1)';
    
    SegValSurf = interpn(SegLapl.data, s2volcoord(1,:), s2volcoord(2,:), s2volcoord(3,:), 'linear',0);
    vertInside = find(SegValSurf<1.1);
    vv1 = normalize3d(surf_orig.coord' - surf_tmp.coord');    
     [Meantmp SD_tmp Dist_tmp] = interVertexDistSurf(surf_tmp);
    if iter>2
        testCoord=s2volcoord+(vv1' * 0.4);
        test_Class=interpn(SegLapl.data, testCoord(1,:), testCoord(2,:), testCoord(3,:), 'linear');
        for Vi=vertInside
            if test_Class(Vi)<1.1
                surf_tmp.coord(:,Vi)=surf_tmp.coord(:,Vi)+(vv1(Vi,:)' * 0.4);
%                 Vi;
            end
        end
    end
     [Mean SD Dist] = interVertexDistSurf(surf1);
     Cand = find(Dist - Dist_tmp>0.1);
     %surf_tmp.coord(1,:)=SurfStatSmoothRegion(surf_tmp.coord(1,:), surf1, 2, Cand, Edg);
     %surf_tmp.coord(2,:)=SurfStatSmoothRegion(surf_tmp.coord(2,:), surf1, 2, Cand, Edg);
     %surf_tmp.coord(3,:)=SurfStatSmoothRegion(surf_tmp.coord(3,:), surf1, 2, Cand, Edg);

%     s2volcoord_tmp = ((double(surf_tmp.coord')-Vol.origin)./Vol.vox(1:3)+1)';        
    
    surf1=surf_tmp;
    
    Def_old = Def;
    
%     deform_index(find(feature(:,2) == 1)) = 0;
    if mod(iter,8) == 0 && iter < 32
        surf1.coord(1,:)=SurfStatSmooth(surf1.coord(1,:), surf1, 1);
        surf1.coord(2,:)=SurfStatSmooth(surf1.coord(2,:), surf1, 1);
        surf1.coord(3,:)=SurfStatSmooth(surf1.coord(3,:), surf1, 1);
    end
    
    if mod(iter,4) == 0 && iter > 32 
        surf1.coord(1,:)=SurfStatSmooth(surf1.coord(1,:), surf1, 1);
        surf1.coord(2,:)=SurfStatSmooth(surf1.coord(2,:), surf1, 1);
        surf1.coord(3,:)=SurfStatSmooth(surf1.coord(3,:), surf1, 1);
    end
    
    iter=iter+1
%     if iter==2
%         Curv=SurfStatComputeCurv_abs(surf1, 1:length(surf1.coord(1,:)));
%         Curv1=SurfStatSmooth(Curv, surf1, 5);
%         Curv1=zscore(Curv1);
%     end
    
end

surf_outskel=surf1;
