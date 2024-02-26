clear all
close all

% Relevant Constants
G = 6.6743e-11; 
M = 2e30; % star mass
m1 = 6e27; % planet 1 mass
m2 = 5e24; % planet 2 mass
earth_year = 365.25*24*60*60; 
R = 0.5*(sqrt(2*G*M)*earth_year/pi)^(2/3); % earth radius


%Time
tmax = 4 * earth_year; 
clockmax = 80000; 
dt = tmax/clockmax; 

%Planet 1 IC
theta = pi * 5/4; % <----- FREE, ANGLE -------------------
w_mag = sqrt(G*M/R);
z = R * [cos(theta), sin(theta)]; %Position Vector
w = w_mag * [-sin(theta), cos(theta)]; %Velocity Vector

%Planet 2 IC
n = 3; % <----- FREE, FREQUENCY --------------------------
ecc = 0.9; % <----- FREE, ECCENTRICITY -------------------
a = R * n^(-2 / 3); %semimajor axis
c = ecc * a; %focal distance
b = sqrt(a^2 - c^2); % semiminor axis  
x = [c + a, 0]; %Position Vector
v = [0, -sqrt((a-c) / (a+c) * (G * M / a))]; %Velocity Vector




% Arrays to store trajectory:
tsave = zeros(1,clockmax);
x1save = zeros(1,clockmax);
y1save = zeros(1,clockmax);
x2save = zeros(1,clockmax);
y2save = zeros(1,clockmax);

energy_save = zeros(1,clockmax); %energy


% Initialization for animation:
plot(0,0,'k*','linewidth',6) % Sun at center
hold on

% Planet Handles
hp1 = plot(z(1),z(2),'bo','linewidth',3); %blue 
hp2 = plot(x(1),x(2),'mo','linewidth',3); %magenta

% Trajectory Handles
ht1 = plot(z(1),z(2),'b','linewidth',1);
ht2 = plot(x(1),x(2),'m','linewidth',1); 


% Axes
a = 2 * R;
axis equal 
axis([-a,a,-a,a])
axis manual 


% Main Loop
for clock = 1:clockmax
    t = clock * dt;

    % Intermediate Values
    D1 = z / norm(z)^3;
    D2 = x / norm(x)^3;
    D3 = (x - z) / norm(x - z)^3;
    
    % Update velocities
    w = w + dt * G*(m2 * D3 - M * D1);
    v = v + dt * G*(m1 * D3 - M * D2);

    % Update Positions (using new velocities)
    z = z + dt * w;
    x = x + dt * v;
    
    % Save Data
    tsave(clock) = t;
    x1save(clock) = z(1);
    y1save(clock) = z(2);
    x2save(clock) = x(1);
    y2save(clock) = x(2);

    % Energy
    kin_1 = m1 * norm(w)^2 / 2;
    kin_2 = m1 * norm(v)^2 / 2;
    pot_1 = -2 * G * M * m1 / norm(z);
    pot_2 = -2 * G * M * m2 / norm(x);
    pot_3 = -2 * G * m1 * m2 / norm(x-z);
    energy = kin_1 + kin_2 + pot_1 + pot_2 + pot_3;
    energy_save(clock) = energy;

    % Update Handles
    hp1.XData = z(1);
    hp1.YData = z(2);
    ht1.XData = x1save(1:clock);
    ht1.YData = y1save(1:clock);

    hp2.XData = x(1);
    hp2.YData = x(2);
    ht2.XData = x2save(1:clock);
    ht2.YData = y2save(1:clock);

    drawnow limitrate

end

error = 0;
for i = 2: clockmax
    error = error + abs(energy_save(i) - energy_save(i-1));
end

%Error Calculation
avg_error = error / clockmax;
disp('Average Error:');
disp(avg_error);
net_energy = energy_save(clockmax) - energy_save(1);
disp('Net energy chage: ');
disp(net_energy);





