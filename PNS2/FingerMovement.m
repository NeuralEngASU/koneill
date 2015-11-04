%% Finger Movement Logic Script



%% Finger Movements
% Maximum finger flex (degrees)
flexMax = [0, 20, 60, 90;...    % Thumb
           0, 90, 90, 60;...    % Index
           0, 90, 90, 60;...    % Middle
           0, 90, 90, 60;...    % Ring
           0, 90, 90, 60];      % Little
          
% Minimum finger flex (degrees)
flexMin = [0, -10, 0, 0;...    % Thumb
           0, -10, 0, 0;...    % Index
           0, -10, 0, 0;...    % Middle
           0, -10, 0, 0;...    % Ring
           0, -10, 0, 0];      % Little      
          
% finger rest (degrees)
% There will be some discussion on these angles
restPos = [0, 10, 10, 10;...    % Thumb
           0, 20, 30, 10;...    % Index
           0, 20, 30, 10;...    % Middle
           0, 20, 30, 10;...    % Ring
           0, 20, 30, 10];      % Little      
          
% Max ab/duction
abdMax = [10, 10, 10, 10;...    % Thumb
           10, 20, 30, 10;...    % Index
           10, 20, 30, 10;...    % Middle
           10, 20, 30, 10;...    % Ring
           10, 20, 30, 10];      % Little   
       
% Min ab/duction
abdMin = [-10, 10, 10, 10;...    % Thumb
           -10, 20, 30, 10;...    % Index
           -10, 20, 30, 10;...    % Middle
           -10, 20, 30, 10;...    % Ring
           -10, 20, 30, 10];      % Little   

%% Linear Transform

% Trial length: 3 seconds with 75 frames per second
time = linspace(0,3, 3*75);

trans = [linspace(0,1, 103), ones(1,19), linspace(1,0, 103)];

% Calculate the difference matrix.
diffMat = flexMax-restPos;

% Initialize transMat. the container for the joing angles for all times for
% a trial
transMat = zeros(5,4);

% loop over all Time
for i = 1:length(time)
    transMat(:,:,i) = trans(i).* diffMat + restPos;
end % END FOR

derp = transMat(2,2,:);
derp = derp(:);
derp1 = [ones(75,1)*20; derp; ones(75, 1)*20];;
plot([ones(50,1)*20; derp; ones(50, 1)*20])


%% Sinusoidal Transform(s)

% Trial length: 3 seconds with 75 frames per second
time = linspace(0,3, 3*75);

% Ease on / Ease off

% Create sinusoid that starts from zero, rises to 1, then glides to 0
trans = -cosd(linspace(0,360,3*75)) + 1;
trans = trans*0.5;

% Calculate a ratio matrix. How many restPos does it take to make a flexMax
divMat = flexMax./restPos;
divMat(isnan(divMat)) = 0;

% Calculate the difference matrix.
diffMat = flexMax-restPos;

% Initialize transMat. the container for the joing angles for all times for
% a trial
transMat = zeros(5,4);

% loop over all Time
for i = 1:length(time)
    transMat(:,:,i) = trans(i).* diffMat + restPos;
end % END FOR

derp = transMat(2,2,:);
derp = derp(:);
derp2 = [ones(75,1)*20; derp; ones(75, 1)*20];;

plot([ones(50,1)*20; derp; ones(50, 1)*20])


%% Linear on / Ease off

% Create sinusoid that starts from zero, rises to 1, then glides to 0
trans = -cosd(linspace(90,270,3*75));

% Calculate a ratio matrix. How many restPos does it take to make a flexMax
divMat = flexMax./restPos;
divMat(isnan(divMat)) = 0;

% Calculate the difference matrix.
diffMat = flexMax-restPos;

% Initialize transMat. the container for the joing angles for all times for
% a trial
transMat = zeros(5,4);

% loop over all Time
for i = 1:length(time)
    transMat(:,:,i) = trans(i).* diffMat + restPos;
end % END FOR

derp = transMat(2,2,:);
derp = derp(:);
derp3 = [ones(75,1)*20; derp; ones(75, 1)*20];

plot([ones(7,1)*20; derp; ones(50, 1)*20])

%%

time = linspace(-1,4, 5*75);

hold on
plot(time, derp1, 'r', 'LineWidth', 2.75)
plot(time, derp2, 'g', 'LineWidth', 2.75)
plot(time, derp3, 'b', 'LineWidth', 2.75)

legend({'Linear', 'Full Sinus', 'Half Sinus'})

title('Finger Movement Options')
xlabel('Time, seconds')
ylabel('Joint Angle, degrees')

xlim([time(1), time(end)]);
ylim([0,100])
hold off

% EOF