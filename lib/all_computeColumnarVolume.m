struct = 'CA';
outDir = '/host/scarus/local_raid/benoit/03_Experiments/PatchSurfHybrid_CorrectedSurfaces_IncludePatients/07_Blades/BladeSampling/ColumnVolume';
inDir  = strcat('/host/scarus/local_raid/benoit/03_Experiments/PatchSurfHybrid_CorrectedSurfaces_IncludePatients/07_Blades/',struct,'/Blades');

[~,tmp] = system(['ls ',inDir]);
listCases = textscan(tmp,'%s');
listCases = listCases{1};

for i = 1:length(listCases)

    disp(listCases{i});
    
    outerSurfFile = strcat(inDir,'/',listCases{i},'/TLE_',listCases{i},'.obj');
    bladeSurfFile = strcat(inDir,'/',listCases{i},'/TLE_',listCases{i},'_skelFinal.obj');
    columnVolFile = strcat(outDir,'/TLE_',listCases{i},'_ColVol.txt');
    
    if exist(outerSurfFile,'file') && exist(bladeSurfFile,'file') && ~exist(columnVolFile,'file')
    
        vertex_col_vol = computeColumnarVolume(outerSurfFile,bladeSurfFile);
        SurfStatWriteData(columnVolFile,vertex_col_vol);
        
    end

end