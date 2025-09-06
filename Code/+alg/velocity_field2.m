function [velocities] = velocity_field2(velocityField,points)
%VELOCITY_FIELD2 - Interpolates 2D velocity field to given particle positions
%
%   velocities = velocity_field2(velocityField,points)
%
%   INPUT:
%     - velocityField : [Nx × Nz × 2] velocity grid (x- and z- components)
%     - points        : [M×3] scatterer positions (x,z used; y ignored)
%
%   OUTPUT:
%     - velocities    : [M×3] velocity vectors interpolated at query points
%
%   Note:
%     - Velocity field is assumed 2D (x-z plane), y-components set to 0
%     - Out-of-bound points are assigned zero velocity

[m,n,~]=size(velocityField);
x_grid = 1:m; % x-axis grid indices
z_grid = 1:n; % z-axis grid indices

% Define interpolants for vx and vz
Fvx = griddedInterpolant({x_grid, z_grid}, velocityField(:,:,1), 'linear', 'none');
Fvz = griddedInterpolant({x_grid, z_grid}, velocityField(:,:,2), 'linear', 'none');

% Extract query coordinates
points_x = points(:,1);
points_z = points(:,3);

% Interpolate velocities
vx = Fvx(points_x, points_z);
vz = Fvz(points_x, points_z);

% Assign zeros to out-of-bound interpolations
vx(isnan(vx)) = 0;
vz(isnan(vz)) = 0;

% y-component = 0
vy = zeros(size(vx));

% Combine components
velocities = [vx, vy, vz];
end