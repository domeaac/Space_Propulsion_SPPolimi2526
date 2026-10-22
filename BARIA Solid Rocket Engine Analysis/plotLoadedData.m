%% PLOTTING THE LOADED DATA

% PLOT 3 pressure curves for each of the 9 beatches firings tests
% figure;
for i = 1:9

    p_mat_i = data.(batch{i});              % pressures of i-th batch

    p_L_i = p_mat_i(:, 1);                  % LP [bar]
    p_M_i = p_mat_i(:, 2);                  % MP [bar]
    p_H_i = p_mat_i(:, 3);                  % HP [bar]

    % time vector: sampling frequency f = 1 kHz, so sampling time dt = 0.001s = 1ms
    N_t_i = length(p_L_i);
    time_vec_i = 0 : 0.001 : (N_t_i - 1)*0.001;

    % %subplot 3x3
    % subplot(3,3,i)
    % plot(time_vec_i, p_L_i); hold on;
    % plot(time_vec_i, p_M_i);
    % plot(time_vec_i, p_H_i);
    % 
    % title(['Batch ', num2str(i)])
    % xlabel('Time [s]')
    % ylabel('Pressure [bar]')

end

% PLOT a single batch
i = 1;
p_mat = data.(batch{i});
N_t = length(p_mat(:, 1));
time_vec = 0 : 0.001 : (N_t - 1)*0.001;

figure;
plot(time_vec, p_mat(:, 1)); hold on;
plot(time_vec, p_mat(:, 2));
plot(time_vec, p_mat(:, 3));

title(['Batch', num2str(i)])
xlabel('Time [s]')
ylabel('Pressure [bar]')

% Zoom on ignition transient (before steady-state transient)
i = 1;
p_mat = data.(batch{i});
N_t = length(p_mat(:, 1));
time_vec = 0 : 0.001 : (N_t - 1)*0.001;

figure;
plot(time_vec(1:40), p_mat(1:40, 1)); hold on;
plot(time_vec(1:40), p_mat(1:40, 2));
plot(time_vec(1:40), p_mat(1:40, 3));

title(['Ignition transient - Batch', num2str(i)])
xlabel('Time [s]')
ylabel('Pressure [bar]')


% % Zoom on batch 7 ignition transient (before steady-state transient)
% i = 7;
% p_mat = data.(batch{i});
% N_t = length(p_mat(:, 1));
% time_vec = 0 : 0.001 : (N_t - 1)*0.001;
% 
% figure;
% plot(time_vec(1:150), p_mat(1:150, 1)); hold on;
% plot(time_vec(1:150), p_mat(1:150, 2));
% plot(time_vec(1:150), p_mat(1:150, 3));
% 
% title(['Ignition transient - Batch', num2str(i)])
% xlabel('Time [s]')
% ylabel('Pressure [bar]')
