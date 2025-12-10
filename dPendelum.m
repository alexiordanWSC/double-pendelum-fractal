function flip_time = dPendelum(z1_0, z3_0)
% Returns the first time the lower pendulum angle z3(t) crosses pi.
% If no crossing occurs: flip_time = NaN.

% parameters of rods 
m1  = 1;                       
m2  = 1;                        
L1  = 0.1;                        
L2  = 0.1;                       
g   = 9.798; % 0.798 is weird but is measured gravity in Denver, CO

u2 = m2 / (m2 + m1);

% time parameters
tmin = 0;
tmax = 20;
h = 0.01;
t = tmin:h:tmax;
M = length(t);

% allocate arrays
z1 = zeros(1, M);
z2 = zeros(1, M);
z3 = zeros(1, M);
z4 = zeros(1, M);

% initial conditions
z1(1) = z1_0;
z2(1) = 0;
z3(1) = z3_0;
z4(1) = 0;

% double pendelum equations
A11 = @(z1,z3) L1;
A12 = @(z1,z3) u2*L2*cos(z1-z3);
A21 = @(z1,z3) L1*cos(z1-z3);
A22 = @(z1,z3) L2;

B1 = @(z1,z3,z4) -g*sin(z1) - u2*L2*sin(z1-z3).*z4.^2;
B2 = @(z1,z2,z3) -g*sin(z3) + L1*sin(z1-z3).*z2.^2;

D = @(z1,z3) A11(z1,z3).*A22(z1,z3) - A12(z1,z3).*A21(z1,z3);

f1 = @(z2) z2;
f3 = @(z4) z4;

f2 = @(z1,z2,z3,z4) ...
    ( A22(z1,z3)*B1(z1,z3,z4) - A12(z1,z3)*B2(z1,z2,z3) ) ./ D(z1,z3);

f4 = @(z1,z2,z3,z4) ...
    ( -A21(z1,z3)*B1(z1,z3,z4) + A11(z1,z3)*B2(z1,z2,z3) ) ./ D(z1,z3);

% rungekutta 4 integration
for k = 1:M-1

    % k1
    k1_1 = f1(z2(k));
    k1_2 = f2(z1(k),z2(k),z3(k),z4(k));
    k1_3 = f3(z4(k));
    k1_4 = f4(z1(k),z2(k),z3(k),z4(k));

    % k2
    z1m = z1(k) + 0.5*h*k1_1;
    z2m = z2(k) + 0.5*h*k1_2;
    z3m = z3(k) + 0.5*h*k1_3;
    z4m = z4(k) + 0.5*h*k1_4;

    k2_1 = f1(z2m);
    k2_2 = f2(z1m,z2m,z3m,z4m);
    k2_3 = f3(z4m);
    k2_4 = f4(z1m,z2m,z3m,z4m);

    % k3
    z1m = z1(k) + 0.5*h*k2_1;
    z2m = z2(k) + 0.5*h*k2_2;
    z3m = z3(k) + 0.5*h*k2_3;
    z4m = z4(k) + 0.5*h*k2_4;

    k3_1 = f1(z2m);
    k3_2 = f2(z1m,z2m,z3m,z4m);
    k3_3 = f3(z4m);
    k3_4 = f4(z1m,z2m,z3m,z4m);

    % k4
    z1m = z1(k) + h*k3_1;
    z2m = z2(k) + h*k3_2;
    z3m = z3(k) + h*k3_3;
    z4m = z4(k) + h*k3_4;

    k4_1 = f1(z2m);
    k4_2 = f2(z1m,z2m,z3m,z4m);
    k4_3 = f3(z4m);
    k4_4 = f4(z1m,z2m,z3m,z4m);

    % update
    z1(k+1) = z1(k) + (h/6)*(k1_1 + 2*k2_1 + 2*k3_1 + k4_1);
    z2(k+1) = z2(k) + (h/6)*(k1_2 + 2*k2_2 + 2*k3_2 + k4_2);
    z3(k+1) = z3(k) + (h/6)*(k1_3 + 2*k2_3 + 2*k3_3 + k4_3);
    z4(k+1) = z4(k) + (h/6)*(k1_4 + 2*k2_4 + 2*k3_4 + k4_4);
end

% flip detection
t_min_check = 0.05;      
target = pi;
escape_tol = 0.05;    % must move 0.05 rad away from pi first

flip_time = NaN;

start_idx = find(t >= t_min_check, 1, 'first');

% to prevent initial conditions starting near pi
% kinda doesn't work
escaped = false;

for k = start_idx:(M-1)

    % check if escaped from the "stuck" zone
    if ~escaped
        if abs(z3(k) - target) > escape_tol
            escaped = true;
        else
            continue; 
        end
    end

    % find crossing of pi
    v1 = z3(k)   - target;
    v2 = z3(k+1) - target;

    % exact hit
    if v1 == 0
        flip_time = t(k);
        return;
    end

    % crossing, necessary because we almost never actually hit pi
    if v1 * v2 < 0
        alpha = abs(v1) / (abs(v1) + abs(v2));
        flip_time = t(k) + alpha*(t(k+1) - t(k));
        return;
    end
end
