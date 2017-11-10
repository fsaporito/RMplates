function [] = meshbuilder(h)

n = length(h);

for i=1:n
    [xv,yv,vertices,edges,endpoints,boundary,boundedges] = meshgen(h(i));
    meshname = ['./meshes/mesh' num2str(h(i)) '.mat'];
    meshplot(xv,yv,endpoints,meshname,i);
    save(meshname,'xv','yv','vertices','edges','endpoints','boundary','boundedges');
end
