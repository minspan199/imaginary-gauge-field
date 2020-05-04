

N = 20;


IntensityNear(N, N) = 0;
IntensityNear(1:2:N - 1, 1:2:N - 1) = 1;
save('imageNear_meshgridX.mat', 'IntensityNear');

figure
imagesc(IntensityNear);

IntensityNear(N, N) = 0;
IntensityNear(1:2:N - 1, :) = 1;
IntensityNear(:, 1:2:N - 1) = 1;
save('imageNear_meshgrid.mat', 'IntensityNear');

IntensityNear = ones(N, N);
save('imageNear_uniform.mat', 'IntensityNear');
