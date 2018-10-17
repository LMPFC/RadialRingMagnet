clear
clc
close all

rPipe = 0.0127;         % Pipe radius [m]
lPipe = 0.2;            % Pipe test length [m]
radialRes = 16;    % radial resolution
lengthRes = 16;    % length resolution
thetaRes = 8;       % theta resolution (ideally multiple of 2)

pipeRMeshCyl = linspace(0,rPipe,radialRes);                    % Create radial points for pipe mesh
pipeZMeshCyl = linspace(-lPipe/2,lPipe/2,lengthRes);                    % Create lengthwise points for pipe mesh
pipeThetaMeshCyl = linspace(0,2*pi-2*pi/thetaRes,thetaRes);    % Create theta points for pipe mesh

% Pipe is set up such that in cylindrical coordinates the z-axis is the
% same as the x-axis in cartesian. Pipe flow travels in the positive
% x-direction, while the vertical is given as the y-direction, and
% correspondingly the horizontal-transverse is given as the z-direction.
% The mesh is laid out corersponding to (r,theta,z).

pipeXMesh = ones(radialRes,thetaRes,lengthRes);
pipeYMesh = pipeXMesh;
pipeZMesh = pipeXMesh;

% Create x-coordinates in cartesian
for zIter = 1:lengthRes
    pipeXMesh(:,:,zIter) = pipeZMeshCyl(zIter)*pipeXMesh(:,:,zIter);
end

% Create y-coordinates in cartesian
for thetaIter = 1:thetaRes
    for rIter = 1:radialRes
        pipeYMesh(rIter,thetaIter,:) = pipeRMeshCyl(rIter)*sin(pipeThetaMeshCyl(thetaIter));
    end
end

% Create z-coordinates in cartesian
for thetaIter = 1:thetaRes
    for rIter = 1:radialRes
        pipeZMesh(rIter,thetaIter,:) = pipeRMeshCyl(rIter)*cos(pipeThetaMeshCyl(thetaIter));
    end
end

% Un-comment the following code to plot out the mesh

% figure()
% hold on
% for rIter = 1:radialRes
%     for thetaIter = 1:thetaRes
%         for zIter = 1:lengthRes
%             plot3(pipeXMesh(rIter,thetaIter,zIter),pipeYMesh(rIter,thetaIter,zIter),pipeZMesh(rIter,thetaIter,zIter),'bo')
%         end
%     end
% end
% hold off


% Create r/z values for ring magnet--distance from center of magnet to the
% pipe center is a necessary parameter, called rFlowEdge

rFlowEdge = 0.3;    % Distance from pipe center to magnet center

magRMesh = ones(radialRes,thetaRes,zRes);
magZMesh = ones(radialRes,thetaRes,zRes);

% Create r-coordinates for ring magnet
for rIter = 1:radialRes
    for thetaIter = 1:thetaRes
        for zIter = 1:lengthRes
            magRMesh(rIter,thetaIter,zIter) = sqrt((rFlowEdge+pipeZMesh(rIter,thetaIter,zIter))^2+pipeXMesh(rIter,thetaIter,zIter)^2);
        end
    end
end

% Create z-coordinates for ring magnet
magZMesh = pipeYMesh;








