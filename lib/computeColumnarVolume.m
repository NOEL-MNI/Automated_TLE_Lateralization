function vertex_col_vol = computeColumnarVolume(outerSurfFile, bladeSurfFile)

    bladeSurf=SurfStatReadSurf1(bladeSurfFile);
    outerSurf=SurfStatReadSurf1(outerSurfFile);

    tlink = dist3(outerSurf.coord',bladeSurf.coord')';

    columnarVol = zeros(3, length(bladeSurf.tri));
    for t = 1:length(bladeSurf.tri)

        % Compute area on skeleton
        L1_Blade = dist3(bladeSurf.coord(:,bladeSurf.tri(t,1))',bladeSurf.coord(:,bladeSurf.tri(t,2))');
        L2_Blade = dist3(bladeSurf.coord(:,bladeSurf.tri(t,2))',bladeSurf.coord(:,bladeSurf.tri(t,3))');
        L3_Blade = dist3(bladeSurf.coord(:,bladeSurf.tri(t,1))',bladeSurf.coord(:,bladeSurf.tri(t,3))');
        Peri_Blade = (L1_Blade + L2_Blade + L3_Blade)/2;
        Area_Blade = sqrt(Peri_Blade * (Peri_Blade - L1_Blade) * (Peri_Blade - L2_Blade) * (Peri_Blade - L3_Blade));

        % Compute area on outer shell
        L1_Outer = dist3(outerSurf.coord(:,outerSurf.tri(t,1))',outerSurf.coord(:,outerSurf.tri(t,2))');
        L2_Outer = dist3(outerSurf.coord(:,outerSurf.tri(t,2))',outerSurf.coord(:,outerSurf.tri(t,3))');
        L3_Outer = dist3(outerSurf.coord(:,outerSurf.tri(t,1))',outerSurf.coord(:,outerSurf.tri(t,3))');
        Peri_Outer = (L1_Outer + L2_Outer + L3_Outer)/2;
        Area_Outer = sqrt(Peri_Outer * (Peri_Outer - L1_Outer) * (Peri_Outer - L2_Outer) * (Peri_Outer - L3_Outer));

        % Get Columnar Volume
        Area = (Area_Blade + Area_Outer)/2;
        columnarVol(:,t) = Area*tlink(bladeSurf.tri(t,:));

    end

    % Get Vertex-Wise Columnar Volume
    columnarVol_Vector = [columnarVol(1,:),columnarVol(2,:),columnarVol(3,:)];
    triangle_Indices   = [bladeSurf.tri(:,1);bladeSurf.tri(:,2);bladeSurf.tri(:,3)];
    vertex_count = hist(triangle_Indices,unique(triangle_Indices));
    vertex_col_vol = accumarray(triangle_Indices, columnarVol_Vector');
    vertex_col_vol = vertex_col_vol ./ vertex_count';
    
end
