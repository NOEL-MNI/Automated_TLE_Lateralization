function normals = computeMeshTriangleNormals3(surf, tri)
% function normals = computeMeshTriangleNormals(surf)

normals = zeros(3,length(tri));
for t = 1:length(tri)
    tp0        = surf.coord(:,surf.tri(tri(t),1));
    tp1        = surf.coord(:,surf.tri(tri(t),2));
    tp2        = surf.coord(:,surf.tri(tri(t),3));

    normal = cross((tp0 - tp1),(tp2 - tp1));
    normals(:,t) = normalize3d(normal')';
end