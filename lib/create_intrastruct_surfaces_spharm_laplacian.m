function create_intrastruct_surfaces_spharm_laplacian(surf_obj, out_intra_surface_file, label_file)  
% label_file: subfields file
% out_intra_surface_file: file name without .obj: ex) output_surf
% Then, this function will produce output_surf_interUpper.obj,
% output_surf_skel.obj, output_surf_interLow.obj, output_surf_boundary.obj
[ll Rand1]= system(['echo $$'])


inter_surf_obj=strcat(out_intra_surface_file, '_inter.obj');
bound_surf_obj=strcat(out_intra_surface_file, '_surf.obj');
skel_surf_obj=strcat(out_intra_surface_file, '_skelFinal.obj');
[rr ll] = system(['ls ' inter_surf_obj]);

if length(strrep(ll , inter_surf_obj,'')) > 5
    disp('%% convert obj to vtk %%');
    disp('');
    bound_surf=surf_obj;
    bound_surf_vtk=strcat('/tmp/tmp_MCsurf_', Rand1, 'TS.vtk')

    disp('%% Read MRI volume %%');
    disp(''); 
    surf1=SurfStatReadSurf1(bound_surf);
    % convert obj to vtk
    save_surface_vtk(surf1, strcat('/tmp/tmp_MCsurf_', Rand1, 'TS.vtk'), 'ASCII')

    % define output file names
    skel_surf_vtk=strcat('/tmp/tmp_MCsurf_', Rand1, '_skel.vtk');
    medialSurfTemp=strcat('/tmp/tmp_MCsurf_', Rand1, '_skel.obj');
    medialSurfVol=strcat('/tmp/tmp_MCsurf_', Rand1, '_skel.mnc');
    OuterSurfVol=strcat('/tmp/tmp_MCsurf_', Rand1, '_outer.mnc');
    OuterSurfVol_pad=strcat('/tmp/tmp_MCsurf_', Rand1, '_outer_pad.mnc');
    SegforLaplVol=strcat('/tmp/tmp_MCsurf_', Rand1, '_segforlapl.mnc');
    SegforLaplVol2=strcat('/tmp/tmp_MCsurf_', Rand1, '_segforlapl_smooth.mnc');
    SegforLaplVol3=strcat('/tmp/tmp_MCsurf_', Rand1, '_segforlapl_smooth_close.mnc');
    SegforLaplVol4=strcat('/tmp/tmp_MCsurf_', Rand1, '_segforlapl_smooth_close2.mnc');
    LaplaceVol=strcat('/tmp/tmp_MCsurf_', Rand1, '_laplace');
    DispXVol=strcat('/tmp/tmp_MCsurf_', Rand1, '_laplace_GradX.mnc');
    DispYVol=strcat('/tmp/tmp_MCsurf_', Rand1, '_laplace_GradY.mnc');
    DispZVol=strcat('/tmp/tmp_MCsurf_', Rand1, '_laplace_GradZ.mnc');
    
    % Skeletonization using voronoi (external x64 binary)
    % Output is vtk which has a skeleton surface with irregular
    % triangulation and with many secondary branches 
    disp('%% Skeletonization %%');
    disp('');
    [ll rr] = system(['cmrep_vskel -Q qvoronoi -p 2.0 -e 3 -q 50 ' ... 
              bound_surf_vtk ' ' skel_surf_vtk])

    disp('%% Load the skeleton from VTK mesh %%');
    disp('');
    m=vtk_polydata_read(skel_surf_vtk);
    rad=vtk_get_point_data(m,'Radius');

    disp('%% Convert vtk to surfstat format %%');
    disp('');
    surf_med1.coord=m.points'
    surf_med1.tri = cell2mat(m.cells.polygons)'
    surf_med1.normal=m.point_data(1).data';
    Prune_ratio=m.point_data(3).data;
    Rad=m.point_data(1).data;
    GeoDist=m.point_data(2).data;
    % SurfStatWriteSurf1('skel_inter1.obj',surf_med1);
    
    % Prune secondary braches using a combination of parameters including radius, geodesic
    % distance and prune ratio.
    disp('%% Prune secondary branches %%');
    disp(''); 
    vertex_to_remove=Prune_ratio<4 | Rad==0 | GeoDist<1;
    Vrm=find(vertex_to_remove==1);
    surf_med2=surf_med1;
    Rad2=Rad;
    GeoDist2=GeoDist;
    Prune_ratio2=Prune_ratio;
    for i=1:sum(vertex_to_remove)
       surf_med2.coord(:, Vrm(i))=[];
       Rad2(Vrm(i))=[];
       Prune_ratio2(Vrm(i))=[];
       GeoDist2(Vrm(i))=[];
       surf_med2.normal(:,Vrm(i))=[];
       tt1=findn(surf_med2.tri == Vrm(i));
       if ~isempty(tt1)
           tt1=unique(tt1(:,1));
           ss1=size(tt1);
           for j=1:ss1(1)
               surf_med2.tri(tt1(j,1),:)=[];
               tt1(j:ss1(1))=tt1(j:ss1(1))-1;
           end
       end
               surf_med2.tri(surf_med2.tri(:,1)>Vrm(i),1)=surf_med2.tri(surf_med2.tri(:,1)>Vrm(i),1)-1;
               surf_med2.tri(surf_med2.tri(:,2)>Vrm(i),2)=surf_med2.tri(surf_med2.tri(:,2)>Vrm(i),2)-1;
               surf_med2.tri(surf_med2.tri(:,3)>Vrm(i),3)=surf_med2.tri(surf_med2.tri(:,3)>Vrm(i),3)-1;

       Vrm(i+1:length(Vrm))=Vrm(i+1:length(Vrm))-1;
    end
    % SurfStatWriteSurf1('skel_inter2.obj',surf_med2);
    % figure; SurfStatView1(Prune_ratio2, surf_med2); cameramenu
    
    % More smoothing
    disp('%% Smooth skeleton shape %%');
    disp('');
    TTT_X=SurfStatSmooth(surf_med2.coord(1,:), surf_med2, 2);
    TTT_Y=SurfStatSmooth(surf_med2.coord(2,:), surf_med2, 2);
    TTT_Z=SurfStatSmooth(surf_med2.coord(3,:), surf_med2, 2);
    surf_med3=surf_med2;
    surf_med3.coord = [TTT_X; TTT_Y; TTT_Z];
    surf_med3.normal = computeVertexWiseNormal(surf_med3);
    % figure; SurfStatView1(Prune_ratio2, surf_med3); cameramenu
    
    disp('%% Generate Laplacian field %%');
    disp('');
    SurfStatWriteSurf1(medialSurfTemp, surf_med3);
    system(['surface_mask2 -binary_mask ' label_file ' ' surf_obj ' ' OuterSurfVol]);
    d.dim 	= fliplr(miinquire(OuterSurfVol,'imagesize')');
    system(['mincreshape -clobber -dimrange  xspace=-3,' num2str(d.dim(1)+6) ' -dimrange yspace=-3,' num2str(d.dim(2)+6) ' -dimrange zspace=-3,' num2str(d.dim(3)+6) ' ' OuterSurfVol ' ' OuterSurfVol_pad]);
    system(['scan_object_to_volume ' OuterSurfVol_pad ' ' medialSurfTemp ' ' medialSurfVol]);
    system(['minccalc -clobber -expr "if(A[0]==1 && A[1]==0) out=2 else if(A[1]>0) out=1 else out=3" ' OuterSurfVol_pad ' ' medialSurfVol ' '  SegforLaplVol]);
    
    % Extend the edges of skeleton along its curvature to the boundary
    % remove the voxel touching to the SPHARM surface so that we ensure all
    % voxels making up of the skeleton to be inside the structure and to be
    % safe against PVE
    skel_edge_extend(SegforLaplVol, SegforLaplVol2);
    system(['mincmorph -clobber -close ' SegforLaplVol2 ' ' SegforLaplVol3]);
    system(['minccalc -expr "if(A[0]==1) out=1 else out=A[1]" ' SegforLaplVol3 ' ' SegforLaplVol2 ' ' SegforLaplVol4])
    
    % Generating Laplacian field within the volume define by the SPHARM
    % surface (outer boundary) and skeleton surface (inner boundary)
    system(['my_MincLaplaceDist -i ' SegforLaplVol4 ' -o ' LaplaceVol ' -like ' SegforLaplVol4]);
    system(['mincblur -clobber -fwhm 1 ' DispXVol ' ' strrep(DispXVol, '.mnc', '')]);
    system(['mincblur -clobber -fwhm 1 ' DispYVol ' ' strrep(DispYVol, '.mnc', '')]);
    system(['mincblur -clobber -fwhm 1 ' DispZVol ' ' strrep(DispZVol, '.mnc', '')]);
    system(['mincblur -clobber -fwhm 1 ' SegforLaplVol4 ' ' strrep(SegforLaplVol4, '.mnc', '')]);
    
    system(['minccalc -clobber -expr "if(A[0]<0.04 && A[0]>-0.04) out=A[1] else out=A[0]" ' DispXVol ' ' strrep(DispXVol, '.mnc', '_blur.mnc') ' ' strrep(DispXVol, '.mnc', '_new.mnc')]);
    system(['minccalc -clobber -expr "if(A[0]<0.04 && A[0]>-0.04) out=A[1] else out=A[0]" ' DispYVol ' ' strrep(DispYVol, '.mnc', '_blur.mnc') ' ' strrep(DispYVol, '.mnc', '_new.mnc')]);
    system(['minccalc -clobber -expr "if(A[0]<0.04 && A[0]>-0.04) out=A[1] else out=A[0]" ' DispZVol ' ' strrep(DispZVol, '.mnc', '_blur.mnc') ' ' strrep(DispZVol, '.mnc', '_new.mnc')]);

    % Deforming the SPHARM surface on to the skeleton, allowing the
    % transposition of the spharm point correspondence to the skeleton
    disp('%% Make correspondences between outer-surf and skeleton %%');
    disp('');
    surf_skel_2side = MakeVertexCorrespond_surf_skel(surf1, strrep(SegforLaplVol4, '.mnc', '_blur.mnc'), strrep(DispXVol, '.mnc', '_new.mnc'), strrep(DispYVol, '.mnc', '_new.mnc'), strrep(DispZVol, '.mnc', '_new.mnc'));

    
    disp('%% Generate t-links %%');
    disp('');
    T_link = (surf_skel_2side.coord - surf1.coord);


    disp('%% Generate intermediate surface %%');
    disp('');
    surf_inter=surf1;
    surf_inter.coord = (surf1.coord + surf_skel_2side.coord) /2;

    SurfStatWriteSurf1(bound_surf_obj, surf1);
    SurfStatWriteSurf1(inter_surf_obj, surf_inter);
    SurfStatWriteSurf1(skel_surf_obj, surf_skel_2side);
else
    surf1 = SurfStatReadSurf1(bound_surf_obj);
    surf_inter = SurfStatReadSurf1(inter_surf_obj);
    surf_skel_2side = SurfStatReadSurf1(skel_surf_obj);
%     DistVertex=dist3(surf_skel_2side.coord', surf1.coord')';

end    


