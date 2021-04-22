function skel_edge_extend(inVolFile, outVolFile)
skelVol=SurfStatReadVol1(inVolFile);
padSiz=30;
skelVol1=skelVol;
skelVol1.data=padarray(round(skelVol.data), [padSiz padSiz padSiz]);
% siz=size(skelVol1.data);
k_rad=17;
TTT=findn(skelVol1.data==2);
skelVol2=skelVol;
iter=1;
old_L=1000000;
while iter<80 && length(TTT(:,1))<old_L;
    new_L=length(TTT(:,1))
    for nn=1:length(TTT(:,1))
        i=TTT(nn,1);
        j=TTT(nn,2);
        k=TTT(nn,3);

        parcel=skelVol1.data(i-k_rad:i+k_rad, j-k_rad:j+k_rad, k-k_rad:k+k_rad);
 
        % Check if at least one neighbor voxel of the parcel's center is part of the medial
        % surface
        if sum(parcel(k_rad:k_rad+2, k_rad+1, k_rad+1)==1) || sum(parcel(k_rad+1, k_rad:k_rad+2, k_rad+1)==1) || sum(parcel(k_rad+1, k_rad+1, k_rad:k_rad+2)==1)
           surfV=findn(parcel==1);
           Mean1 = mean(surfV);
           %STD1= std(surfV);
           Disp1 = [k_rad+1 k_rad+1 k_rad+1] - Mean1;
           PCA1=pca(surfV);
           Angle=acos(dot(PCA1(:,3),Disp1) / (norm(PCA1(:,3))*norm(Disp1))) * 180 / pi;
           if abs(Angle) > 83 && abs(Angle) < 97 && norm(Disp1)>2
               if sum(parcel(k_rad:k_rad+2, k_rad+2, k_rad+1)==1)<3 && sum(parcel(k_rad:k_rad+2, k_rad, k_rad+1)==1)<3 && sum(parcel(k_rad:k_rad+2, k_rad+1, k_rad+2)==1)<3 && sum(parcel(k_rad:k_rad+2, k_rad+1, k_rad)==1)<3
                   if sum(parcel(k_rad+2, k_rad:k_rad+2, k_rad+1)==1)<3 && sum(parcel(k_rad, k_rad:k_rad+2, k_rad+1)==1)<3 && sum(parcel(k_rad+1, k_rad:k_rad+2, k_rad+2)==1)<3 && sum(parcel(k_rad+1, k_rad:k_rad+2, k_rad)==1)<3
                       if sum(parcel(k_rad+1, k_rad+2, k_rad:k_rad+2)==1)<3 && sum(parcel(k_rad+1, k_rad, k_rad:k_rad+2)==1)<3 && sum(parcel(k_rad+2, k_rad+1, k_rad:k_rad+2)==1)<3 && sum(parcel(k_rad, k_rad+1, k_rad:k_rad+2)==1)<3
                          if sum(sum(sum(parcel(k_rad:k_rad+2, k_rad:k_rad+2, k_rad:k_rad+2)>2)))==0
                             skelVol2.data(i-padSiz, j-padSiz, k-padSiz)=1;
                          end
                       end
                   end
               end
           end

        end
    end
    skelVol1.data=padarray(round(skelVol2.data), [padSiz padSiz padSiz]);
    skelVol1.data(skelVol1.data==4)=1;
    old_L=new_L;
    TTT=findn(skelVol1.data==2);
    iter=iter+1;
end

skelVol2.file_name=outVolFile;
SurfStatWriteVol1(skelVol2);
