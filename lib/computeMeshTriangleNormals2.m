function normals = computeMeshTriangleNormals(surf)
% function normals = computeMeshTriangleNormals(surf)

normals = zeros(3,length(surf.tri));
for t = 1 : length(surf.tri)
    tp0        = surf.coord(:,surf.tri(t,1));
    tp1        = surf.coord(:,surf.tri(t,2));
    tp2        = surf.coord(:,surf.tri(t,3));

    normal = cross((tp0 - tp1),(tp2 - tp1));
    normals(:,t) = normalize3d(normal')';
end