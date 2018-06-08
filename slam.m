% Name        : [Xslam,Pslam,theLoops]=slam(odoData,drawRobot,drawLoops,onlinePlot,finalPlot)
% Description : Performs 2D SLAM.
% Input       : odoData    - Output of compute_odometry
%               drawRobot  - If <>0, the robot and the uncertainty ellipse
%                            are plotted at the last estimated pose and at
%                            some intermediate poses.
%               drawLoops  - If <>0, the loop closures are plotted.
%               onlinePlot - SLAM output is plotted during execution.
%               finalPlot  - SLAM output is plotted at the end.
% Output      : Xslam      - State vector (mean)
%               Pslam      - State vector (covariance)
%               theLoops   - Detected loops, in the format required by
%                            draw_slam.
% Note        : The SLAM state vector is a set of relative motions from
%               one frame to the next so that X(i*3-2:i*3) denotes the
%               motion from frame i-1 to frame i.
% Note        : If you use this software, please cite the paper stated in
%               README.TXT
% Author      : Antoni Burguera
%               antoni.burguera@uib.es
function [Xslam,Pslam,theLoops]=slam(odoData,drawRobot,drawLoops,onlinePlot,finalPlot)
    % Initialize plotting/logging resources
    theLoops=[];
    if onlinePlot
        figure;
    end;

    % Parameters determined experimentally.
    TwoSigmaPos=2;
    TwoSigmaAngle=3*pi/180;
    Podo=[(TwoSigmaPos/2)^2,0,0;0,(TwoSigmaPos/2)^2,0;0,0,(TwoSigmaAngle/2)^2];
    loopMaxDistance=500;
    
    % Initialize mean (Xslam) and covariance (Pslam) of the state vector.
    Xslam=[];
    Pslam=[];

    % For each odometric estimate
    for curFrame=1:size(odoData,2)

        % State augmentation: add the current odometric estimate into the
        % state vector.
        Xslam=[Xslam;odoData(curFrame).X];
        Pslam(end+1:end+3,end+1:end+3)=Podo;
        
        % Search loop closings
        [theIndexes,theMeasurements]=search_loops(Xslam,odoData,loopMaxDistance);

        % Store the loop information (just for plotting/logging purposes)
        theLoops=[theLoops [theIndexes;zeros(1,size(theIndexes,2))+curFrame]];
        
        % State vector update (i.e. graph optimization) according to detected loops.
        [Xslam,Pslam]=slam_update(Xslam,Pslam,Podo/10,theMeasurements,theIndexes);

        % Online plot
        if onlinePlot
            cla;
            draw_slam(Xslam,Pslam,theLoops,drawRobot,drawLoops);
            axis equal;
            drawnow;
        end;
    end;
    
    % Final plot
    if finalPlot && ~onlinePlot
        figure;
        draw_slam(Xslam,Pslam,theLoops,drawRobot,drawLoops);
        axis equal;
    end;
return;

% Name        : [Xslam,Pslam]=slam_update(Xslam,Pslam,Ploop,actualMeasurements,frameIndexes)
% Description : Updates the state vector (i.e. optimizes the graph) when a
%               loop is detected.
% Input       : Xslam      - State vector (mean)
%               Pslam      - State vector (covariance)
%               Ploop      - Error of a loop closure modeled as a
%                            covariance.
%               actualMeasurements - 3xN matrix of loop closures. Each
%                            column is (x,y,o)' transformation.
%               frameIndexes - The indexes of the frames that close the
%                            loop with the current one.
% Output      : Xslam      - State vector (mean)
%               Pslam      - State vector (covariance)
function [Xslam,Pslam]=slam_update(Xslam,Pslam,Ploop,actualMeasurements,frameIndexes)
    predictedMeasurements=zeros(3,size(frameIndexes,2));
    for i=1:size(frameIndexes,2)
        X=zeros(3,1);
        for j=frameIndexes(:,i):(size(Xslam,1)/3)-1,
            [X,~]=compose_references(X,Xslam((j+1)*3-2:(j+1)*3,1),[],[]);
        end;
        predictedMeasurements(:,i)=X;
    end;
    for j=1:10,  
        H=[];  
        theInnovation=zeros(size(actualMeasurements,2)*3,1);
        Rsm=zeros(size(actualMeasurements,2)*3);
        Xtmp=Xslam;
        Ptmp=Pslam;
        for i=1:size(actualMeasurements,2),
            [Htmp]=compute_observation_jacobian(Xtmp, predictedMeasurements(:,i), frameIndexes(:,i)+1);
            H=[H;Htmp];
            theMeasurement=actualMeasurements(:,i);
            thePrediction=predictedMeasurements(:,i);
            theDifference=theMeasurement-thePrediction;
            theDifference(3,1)=normalize(theDifference(3,1));
            i1=i*3;
            i0=i1-2;
            theInnovation(i0:i1,1)=theDifference;
            Rsm(i0:i1,i0:i1)=Ploop;
        end;       
        if (size(actualMeasurements,2)>0),
            Ptmp=Pslam-Pslam*H'*(inv(H*Pslam*H'+Rsm))*H*Pslam;
            Xtmp=Xtmp+Ptmp*H'*(inv(Rsm))*(theInnovation)-Ptmp/Pslam*(Xtmp-Xslam);              
        end;
    end;
    Xslam=Xtmp;
    Pslam=Ptmp;          
return;

% Name        : H=compute_observation_jacobian(Xslam,hk,startIndex)
% Description : Computes the Jacobian matrix corresponding to an
%               observation (i.e. loop closure)
% Input       : Xslam      - State vector (mean)
%               hk         - Predicted loop closures (obtained from the
%                            state vector).
%               startIndex - Frame index where the loop starts.
% Output      : H          - The Jacobian matrix.
function H=compute_observation_jacobian(Xslam,hk,startIndex)
  % Build the first part of the Jacobian, which is composed of zeros
  H=zeros(3,3*(startIndex-1));
  % The rest of items
  glk=zeros(3,1);
  % Initialize c=cos(0), s=sin(0)
  c=1;
  s=0;
  for i=startIndex:size(Xslam,1)/3,
    % Compute glk
    Xs=Xslam(i*3-2:i*3,1);
    glk=[glk(1)+Xs(1)*c-Xs(2)*s;glk(2)+Xs(1)*s+Xs(2)*c;glk(3)+Xs(3)];
    % Precompute some parameters
    c=cos(glk(3));
    s=sin(glk(3));
    p1=[-glk(1)*c-glk(2)*s+hk(1)*c+hk(2)*s;
        glk(1)*s-glk(2)*c-hk(1)*s+hk(2)*c];
    p2=glk(3)-Xslam(i*3,1);
    sp=sin(p2);
    cp=cos(p2);
    % Store the partial Jacobian into the output Jacobian matric
    H=[H [1,0,-p1(1)*s-p1(2)*c;0,1,p1(1)*c-p1(2)*s;0,0,1]*[cp,-sp,0;sp,cp,0;0,0,1]];    
  end;
return;

% Name        : [theIndexes,theMeasurements]=search_loops(Xslam,odoData,maxDistance)
% Description : Searches loop closures.
% Input       : Xslam       - State vector (mean)
%               odoData     - Output of compute_odometry
%               maxDistance - Frames which are (accordint to the state
%                             vector) farther than this distances are not
%                             considered as candidate loop closures.
% Output      : theIndexes  - Indexes within odoData of frames that close a
%                             loop with the last one.
%               theMeasurements - 3xN matrix. Each column is (x,y,o)' and
%                             denotes the estimated loop roto-translation
%                             (from previous to current frame).
function [theIndexes,theMeasurements]=search_loops(Xslam,odoData,maxDistance)
    numPoses=size(Xslam,1)/3;
    X=Xslam(end-2:end,1);
    theIndexes=[];
    theMeasurements=[];
    for i=numPoses-2:-1:1
        curSLAM=Xslam((i+1)*3-2:(i+1)*3,1);
        [X,~]=compose_references(curSLAM,X,[],[]);
        if sqrt(X(1:2)'*X(1:2))<maxDistance
            [matches,~]=vl_ubcmatch(odoData(i).d,odoData(numPoses).d);
            if size(matches,2)>10
                [Z,fail]=ransac_estimate_motion(odoData(i).f(1:2,matches(1,:)),odoData(numPoses).f(1:2,matches(2,:)),1000,5,10,.75);
                if ~fail
                    theIndexes=[theIndexes i];
                    theMeasurements=[theMeasurements Z];
                end;
            end;
        end;
    end;
return;