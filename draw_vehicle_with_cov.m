% Name        : draw_vehicle_with_cov (X, P, robotSize)
% Description : Draws a triangle at the pose X and the 95% uncertainty
%               ellipse of N(X,P) (only position).
% Input       : X - Pose (x,y,o)'
%               P - Pose covariance
%               robotSize - The size of the triangle.
function draw_vehicle_with_cov (X, P, robotSize)
    draw_vehicle(X, robotSize); hold on;
    draw_ellipse(X, P, 'r');
end