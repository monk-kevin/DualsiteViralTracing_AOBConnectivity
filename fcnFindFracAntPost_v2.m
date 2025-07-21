function [FracAntPost,RotMatStruct] = fcnFindFracAntPost_v2(T,RotMatStruct)
%input T = CellCounterResults Table, plot_flag (plots orig, rot, rot+shift data)
%       > Note! T should only contain one type of cell
%output = FracAntPost where:
%         1st column 0 is more POSTERIOR and 1 is more ANTERIOR.
%         2nd column 0 is more VENTRAL and 1 is more DORSAL
%               edited 250312

%function begins
%default NO plotting
plot_flag = 0;
NT_flag = 0;
if nargin==1
    NT_flag = 1;
    RotMatStruct = struct;
end
xData = T.X(T.Type<5)';yData = -T.Y(T.Type<5)';
if NT_flag == 1
    %setting up for rotation if neurotrace
    Ant = [T.X(T.Type==5),-T.Y(T.Type==5)]; %taken from CellCounterResults file
    Post = [T.X(T.Type==6),-T.Y(T.Type==6)];
    xLandmark = [Ant(1) Post(1)];
    yLandmark = [Ant(2) Post(2)];
    RotMatStruct.xLandmark = xLandmark;
    RotMatStruct.yLandmark = yLandmark;

    vLandmark = [xLandmark;yLandmark];
    %sets center of rotation around POSTERIOR landmark
    x_center = xLandmark(2);
    y_center = yLandmark(2);
    %output to structure
    RotMatStruct.x_center = x_center;
    RotMatStruct.y_center = y_center;
    
    M = (yLandmark(2) - yLandmark(1))./(xLandmark(2) - xLandmark(1)); %finding slope of line
    ANG = atand(M); %calc angle of line to be rotated
    DEG = 360 - ANG; %calc angle of rotation
    theta = DEG*pi/180; %convert to rad
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)]; %rotation matrix
    %output to structure
    RotMatStruct.R = R;
    
    %rotating Ant Post landmarks to normalize axis
    centerLandmark = repmat([x_center;y_center],1,length(xLandmark));    
    voLandmark = R*(vLandmark-centerLandmark)+centerLandmark; %rotating line: first shifts points in plane so that center of R is at origin; applying rotation at origin; shift points back
    xRotLandmark = voLandmark(1,:);
    yRotLandmark = voLandmark(2,:);
    if round(diff(yRotLandmark),10)~=0
        error('Rotation did not work! deltaY ~= 0')
    end
    %output to structure
    RotMatStruct.voLandmark = voLandmark;

    
    % applying rotation to actual labeled data
    vData = [xData;yData];
    centerData = repmat([x_center;y_center],1,length(xData));
    voData = R*(vData-centerData)+centerData; %rotating labeled data based on line above
    xRotData = voData(1,:);
    yRotData = voData(2,:);

    %put POSTERIOR landmark X val as 0
    % OffsetVal = 0 - xRotLandmark(2);
    OffsetVal = 0 - min(xRotData);
    xRotShiftLandmark = abs(xRotLandmark + OffsetVal);
    xRotShiftData = abs(xRotData + OffsetVal);
    %output offset to struct
    RotMatStruct.OffsetVal = OffsetVal;
    RotMatStruct.xRotShiftLandmark = xRotShiftLandmark;
    RotMatStruct.yRotLandmark = yRotLandmark;
    % adjusted 250312
    %After all of that, normalize to mins and maxes (note that y values
    %must be brought above 0 to make dorsal 1 and ventral 0
    % AP normalization:
    xRotShiftDataMin = min(xRotShiftData);xRotShiftDataMax = max(xRotShiftData);
    xRotShiftDataNorm = (xRotShiftData - xRotShiftDataMin) ./ (xRotShiftDataMax - xRotShiftDataMin);
    %output to structure
    RotMatStruct.xRotShiftDataMin = xRotShiftDataMin;
    RotMatStruct.xRotShiftDataMax = xRotShiftDataMax;



    % normalizing DV based on neurotrace min/max DV
    yRotDataMin = min(yRotData);
    yRotDataPos = yRotData - yRotDataMin;
    yRotDataNorm = yRotDataPos ./ max(yRotDataPos);
    %output to structure
    RotMatStruct.MinMaxScale = [yRotDataMin max(yRotDataPos)];

else % rotate other data based on rotation of NT

    %apply rotation to labeled data based on rotation matrix from NT
    vData = [xData;yData];
    centerData = repmat([RotMatStruct.x_center;RotMatStruct.y_center],1,length(xData));
    voData = RotMatStruct.R*(vData-centerData)+centerData; %rotating labeled data based on line above
    xRotData = voData(1,:);
    yRotData = voData(2,:);

    %shift X coordinates based on NT offset value
    xRotShiftData = abs(xRotData + RotMatStruct.OffsetVal);
    
    %normalize AP based on NT
    xRotShiftDataNorm = (xRotShiftData - RotMatStruct.xRotShiftDataMin) ./ (RotMatStruct.xRotShiftDataMax - RotMatStruct.xRotShiftDataMin);

    %normalizing DV based on NT min/max DV
    yRotDataPos = yRotData - RotMatStruct.MinMaxScale(1);
    yRotDataNorm = yRotDataPos ./ RotMatStruct.MinMaxScale(2);

end



FracAntPost(:,1) = xRotShiftDataNorm';
FracAntPost(:,2) = yRotDataNorm';

%plotting original, rotated, and rotated+shifted lines and data pints
if plot_flag==1
    clf
    hold on
    plot(RotMatStruct.xLandmark(1),RotMatStruct.yLandmark(1),'k*')
    plot(RotMatStruct.xLandmark(2),RotMatStruct.yLandmark(2),'ks')
    plot(RotMatStruct.xLandmark,RotMatStruct.yLandmark,'k-')
    plot(xData,yData,'k.')
    %     plot(xRotLandmark,yRotLandmark,'ro-')
    %     plot(xRotData,yRotData,'r.')
    plot(RotMatStruct.xRotShiftLandmark(1),RotMatStruct.yRotLandmark(1),'b*')
    plot(RotMatStruct.xRotShiftLandmark(2),RotMatStruct.yRotLandmark(2),'bs')
    plot(RotMatStruct.xRotShiftLandmark,RotMatStruct.yRotLandmark,'b-')
    plot(xRotShiftData,yRotData,'b.')

    %inset for normalized values
    LgFigPsn = gca().Position;
    axes('Position',[LgFigPsn(1) + .025*LgFigPsn(3) LgFigPsn(2) + .025*LgFigPsn(4) LgFigPsn(3)*.25 LgFigPsn(4)*.25])
    plot(FracAntPost(:,1),FracAntPost(:,2),'.','color',[0 0 0.5])
    set(gca,'xtick',[0 1],'ytick',[0 1])

end


% CellCounterResults_
%adding data as new row to table - FORTHCOMING
