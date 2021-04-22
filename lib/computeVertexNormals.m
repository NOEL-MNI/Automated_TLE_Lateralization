function normal = computeVertexNormals(surf, vertex)
% Norm=computeMeshTriangleNormals2(surf);
normal=surf.normal;
l=length(surf.tri(:,1));
for i=vertex
tri=find(surf.tri==i);
tri2=tri-floor((tri-1)/l)*l;
normalV=computeMeshTriangleNormals3(surf, tri2);
normal(:,i)=mean(normalV');
normal(:,i)=normalize3d(normal(:,i)')';
end
