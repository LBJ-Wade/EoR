% Script to rotate tile beams to the station tile orientation.

beams = dlmread('tile_beams_unrotated.txt');
angles = dlmread('angles.txt');
angles = angles * 180.0/pi;

%[l0,m0,~] = sph2cart(beams(:,1)*(pi/180), beams(:,2)*(pi/180), 1);
%scatter(l0, m0);
stations = [dir('CS*'); dir('RS*')];

assert(length(stations) == length(angles));

nstations = length(angles);
for i = 1:nstations
    beams_rotated = beams;
    beams_rotated(:,1) = beams_rotated(:,1) + angles(i);
    
    stationDir = stations(i).name;
    beamsPath = [stationDir '/tile/permitted_beams.txt'];
    %delete(beamsPath)
    dlmwrite(beamsPath, beams_rotated);
    
    %[l,m,~] = sph2cart(beams_rotated(:,1)*(pi/180), beams_rotated(:,2)*(pi/180), 1);
    %figure(1);
    %scatter(l,m, 'filled');
    %drawnow;
    %pause(0.2);
    
end
