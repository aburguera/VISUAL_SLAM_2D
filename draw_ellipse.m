% Name        : draw_ellipse(X,P,theColor)
% Description : Draws a 95% uncertainty ellipse of a Normal with mean X and
%               covariance P. The distribution must be bivariate.
% Input       : X - The mean. One column and two or more rows, though only
%                   X(1:2,1) will be used.
%               P - The covariance. If larger than 2x2, only P(1:2,1:2)
%                   will be considered.
%               theColor - Valid Matlab color spec.
function draw_ellipse(X,P,theColor)
    tita=linspace(0,2*pi,20);
    theCircle=[cos(tita);sin(tita)];
    [V,D]=eig(full(P(1:2,1:2)));
    ejes=sqrt(9.2103*diag(D));
    tmp=(V*diag(ejes))*theCircle;
    hp=line(tmp(1,:)+X(1),tmp(2,:)+X(2));
    set(hp,'Color',theColor);
    set(hp,'LineWidth',1.5);
return;