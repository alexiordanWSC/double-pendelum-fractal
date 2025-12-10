
th1_0 = pi/2;
th2_0 = pi/2;
N = 100;  % dimension in pixels of grid
th1_vals = linspace(th1_0 - 0.001, th1_0 + 0.001, N);  
th2_vals = linspace(th2_0 - 0.001, th2_0 + 0.001, N);  

flip_times = NaN(N, N);  

% loop
for i = 1:N
    for j = 1:N
        flip_times(i,j) = dPendelum(th1_vals(i), th2_vals(j));
    end
end

% plot fractal
figure;
imagesc(th1_vals, th2_vals, flip_times);
axis xy;
xlabel('\theta_1 initial (rad)');
ylabel('\theta_2 initial (rad)');
title('Double Pendulum Flip-Time Fractal');
colorbar;
colormap(turbo);
caxis([0 20]);    
