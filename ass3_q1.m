% =========
% ass3_q1.m
% =========
%
% This assignment will introduce you to the idea of first building an
% occupancy grid then using that grid to estimate a robot's motion using a
% particle filter.
% 
% There are two questions to complete (5 marks each):
%
%    Question 1: code occupancy mapping algorithm 
%    Question 2: see ass3_q2.m
%
% Fill in the required sections of this script with your code, run it to
% generate the requested plot/movie, then paste the plots into a short report
% that includes a few comments about what you've observed.  Append your
% version of this script to the report.  Hand in the report as a PDF file
% and the two resulting AVI files from Questions 1 and 2.
%
% requires: basic Matlab, 'gazebo.mat'
%
% T D Barfoot, January 2016
%
clear all;

% set random seed for repeatability
rng(1);

% ==========================
% load the dataset from file
% ==========================
%
%    ground truth poses: t_true x_true y_true theta_true
% odometry measurements: t_odom v_odom omega_odom
%           laser scans: t_laser y_laser
%    laser range limits: r_min_laser r_max_laser
%    laser angle limits: phi_min_laser phi_max_laser
%
load gazebo.mat;

% =======================================
% Question 1: build an occupancy grid map
% =======================================
%
% Write an occupancy grid mapping algorithm that builds the map from the
% perfect ground-truth localization.  Some of the setup is done for you
% below.  The resulting map should look like "ass2_q1_soln.png".  You can
% watch the movie "ass2_q1_soln.mp4" to see what the entire mapping process
% should look like.  At the end you will save your occupancy grid map to
% the file "occmap.mat" for use in Question 2 of this assignment.

% allocate a big 2D array for the occupancy grid
ogres = 0.05;                   % resolution of occ grid (cell size side, divide for x and y pixels)
ogxmin = -7;                    % minimum x value
ogxmax = 8;                     % maximum x value
ogymin = -3;                    % minimum y value
ogymax = 6;                     % maximum y value
ognx = (ogxmax-ogxmin)/ogres;   % number of cells in x direction
ogny = (ogymax-ogymin)/ogres;   % number of cells in y direction
oglo = zeros(ogny,ognx);        % occupancy grid in log-odds format
ogp = zeros(ogny,ognx);         % occupancy grid in probability format

% precalculate some quantities
numodom = size(t_odom,1);
npoints = size(y_laser,2);
angles = linspace(phi_min_laser, phi_max_laser,npoints);
dx = ogres*cos(angles);
dy = ogres*sin(angles);

% interpolate the noise-free ground-truth at the laser timestamps
t_interp = linspace(t_true(1),t_true(numodom),numodom);
x_interp = interp1(t_interp,x_true,t_laser);
y_interp = interp1(t_interp,y_true,t_laser);
theta_interp = interp1(t_interp,theta_true,t_laser);
omega_interp = interp1(t_interp,omega_odom,t_laser);
  
% set up the plotting/movie recording
vid = VideoWriter('ass2_q1.avi');
open(vid);
figure(1);
clf;
pcolor(ogp);
colormap(1-gray);
shading('flat');
axis equal;
axis off;
M = getframe;
writeVideo(vid,M);

ALPHA = 1; %occupied cell
BETA = 1; %unoccupied cell
omega_max = 0.1; %only plot if angular velocity is less than 0.1
d_scan = -0.1; %10cm distance behind origin of robot, x distance??
x_i = NaN(npoints,1); %initialize map points as NaNs for NaN scans 
y_i = NaN(npoints,1);
nscans = size(1:5:size(t_laser,1));
x_start = NaN(nscans); 
y_start = NaN(nscans);

NUM_PTS_OBSTACLE = 2; %+1 num points to consider as obstacles

% loop over laser scans (every fifth)
for i=1:5:size(t_laser,1)
    
    % ------insert your occupancy grid mapping algorithm here------

    %calculate starting coord for each laser scan, (inertial frame -> pixel coord)
    x_start = round((x_interp(i) + d_scan*cos(theta_interp(i)) - ogxmin)/ogres);
    y_start = round((y_interp(i) + d_scan*sin(theta_interp(i)) - ogymin)/ogres);
    
    % loop over laser scan angles (columns, associated to angles)
    for j=1:npoints
        %Check if scans are NaN and time stamp interpolation errors and r_min_laser
        if omega_max > abs(omega_interp(i)) && ~isnan(y_laser(i,j)) && y_laser(i,j) > r_min_laser
            
            %set new log odds for each scan
            new_oglo = zeros(ogny,ognx);
            
            %boolean for max laser distance
            max_laser_dist = 0;
            if y_laser(i,j) > r_max_laser
                max_laser_dist = 1; %true
                r_laser = 10; %set laser distance as max laser distance
            else
                r_laser = y_laser(i,j);
            end
            
            %r_laser = y_laser(i,j);
            
            %current frame (laser on robot)
            x_c = r_laser*cos(angles(j));
            y_c = r_laser*sin(angles(j))';

            %robot frame
            x_r = x_c + d_scan;
            y_r = y_c;

            %inertial frame, need to rotate using 2D rotation matrix
            x_i(j) = x_interp(i) + cos(theta_interp(i))*x_r - sin(theta_interp(i))*y_r;
            y_i(j) = y_interp(i) + sin(theta_interp(i))*x_r + cos(theta_interp(i))*y_r;
            
            %end coords of laser scan, (inertial frame -> pixel coord)
            x_end = round((x_i(j) - ogxmin)/ogres);
            y_end = round((y_i(j) - ogymin)/ogres);
            
            %raytrace - use built in bresenham line optimized algorithm
            %[x y]=bresenham(x1,y1,x2,y2)
            coord = bresenham(x_start, y_start, x_end, y_end);
            x_coords = coord(:,1); %start to end
            y_coords = coord(:,2); %start to end
         
            if max_laser_dist == 1
                ind = sub2ind(size(oglo), y_coords, x_coords);
                new_oglo(ind) = -BETA;
            else
                %unoccupied cells, find linear matrix indices on rows and columns to assign -BETA for logodds
                ind = sub2ind(size(oglo), y_coords(1:end-NUM_PTS_OBSTACLE), x_coords(1:end-NUM_PTS_OBSTACLE));
                new_oglo(ind) = -BETA;

                %occupied cells, find linear matrix indices on rows and columns to assign ALPHA for logodds
                ind = sub2ind(size(oglo), y_coords(end-NUM_PTS_OBSTACLE+1:end), x_coords(end-NUM_PTS_OBSTACLE+1:end));
                new_oglo(ind) = ALPHA;
            end
            
            %update log odds
            oglo = oglo + new_oglo;
            
            %clipping
            oglo(find(oglo < -BETA*5)) = -BETA*5;
            
        end
    end
    
    %update probabilities by converting to probabilities
    %check for too large log values >100
    ogp(find(ogp > 100)) = 100;
    ogp = exp(oglo)./(1 + exp(oglo));
    
    % ------end of your occupancy grid mapping algorithm-------

    % draw the map
    clf;
    pcolor(ogp);
    colormap(1-gray);
    shading('flat');
    axis equal;
    axis off;

    % draw the robot
    hold on;
    x = (x_interp(i)-ogxmin)/ogres;
    y = (y_interp(i)-ogymin)/ogres;
    th = theta_interp(i);
    r = 0.15/ogres;
    set(rectangle( 'Position', [x-r y-r 2*r 2*r], 'Curvature', [1 1]),'LineWidth',2,'FaceColor',[0.35 0.35 0.75]);
    set(plot([x x+r*cos(th)]', [y y+r*sin(th)]', 'k-'),'LineWidth',2);
    
    %debugging
    %set(plot( x_end, y_end, 'b.' ),'MarkerSize',5);
    %set(plot( x_start, y_start, 'r.' ),'MarkerSize',20);
    
    % save the video frame
    M = getframe;
    writeVideo(vid,M);
    
    pause(0.1);
    
end

close(vid);
print -dpng ass2_q1.png

%saved populated occupancy map
save occmap.mat ogres ogxmin ogxmax ogymin ogymax ognx ogny oglo ogp;

