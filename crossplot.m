function crossplot(filename)

load(filename);

i = IG(2:end/2, :);
g = IG(end/2 +2:end, :);

figure;
scatter(i(:),g(:));

C = opCurvelet2d(size(i,1), size(i, 2));


Ci = C*i(:);
Cg = C*g(:);

figure;
scatter(Ci, Cg)




end