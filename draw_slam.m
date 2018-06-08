% Name        : draw_slam(X,P,theLoops,drawRobot,drawLoops)
% Description : Plots SLAM output
% Input       : X         - Mean of the state vector
%               P         - Covariance of the state vector
%               drawRobot - If <>0 plot robot and associated uncertainty
%                           ellipse at current pose and some intermediate
%                           poses.
%               drawLoops - If <>0 plots the detected loops.
function draw_slam(X,P,theLoops,drawRobot,drawLoops)
    Xcur=zeros(3,1);
    Pcur=zeros(3);
    numNodes=size(X,1)/3;
    Xh=zeros(3,numNodes);
    for i=1:numNodes
        i0=i*3-2;
        i1=i*3;
        [Xcur,Pcur]=compose_references(Xcur,X(i0:i1,1),Pcur,P(i0:i1,i0:i1));
        Xh(:,i)=Xcur;
        if drawRobot && mod(i,10)==0
            draw_vehicle_with_cov(Xcur,Pcur,100);
            hold on;
        end;
    end;
    if drawRobot
        draw_vehicle_with_cov(Xcur,Pcur,100);
        hold on;
    end;
    if drawLoops
        if size(theLoops,2)>0
            theX=[Xh(1,theLoops(1,:));Xh(1,theLoops(2,:))];
            theY=[Xh(2,theLoops(1,:));Xh(2,theLoops(2,:))];
            plot(theX,theY,'b');
            hold on;
        end;
    end;
    plot(Xh(1,:),Xh(2,:),'k','LineWidth',2);
return;
