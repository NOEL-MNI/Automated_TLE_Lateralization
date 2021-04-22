function K_mean=SurfStatComputeCurv(surf, vertex)
K_mean=zeros(1, length(vertex));
vv2=1;
for vv=vertex
    tri=findn(surf.tri==vv);
    normals=computeMeshTriangleNormals3(surf, tri);
    siz = size(normals);
    Norm=zeros(siz(2), siz(2));
    for i=1:siz(2)
        for j=1:siz(2)
            if i~=j
                Norm(i,j)=normals(:,i)'*normals(:,j);
            end
        end
    end
    N_max=max(max(Norm));
    k1=find(Norm>0);
    N_min=min(min(Norm(k1)));

    K_mean(vv2)=1-(N_max+N_min)/2;
    vv2=vv2+1;
end


