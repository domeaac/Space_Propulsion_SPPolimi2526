
% Identify p_max,
% NOTE: that shall not be in the transient, so if it is found in the
% first 200 points

p_eff_vec = zeros(9,3);
rb_vec    = zeros(9,3);
c_vec     = zeros(27,1);
toll = 0.05;

figure;
c1 = []; c2 = []; c3 = [];
for i = 1 : 9 % batches
    for j = 1 : 3 % 3 motors for each batch

        p_mat_i = data.(batch{i}); % pressures of i-th batch

        % switch for pressure trace 7
        if i == 7
            if j == 2
                p = p_mat_i(:,3);
            elseif j == 3
                p = p_mat_i(:,2);
            elseif j == 1
                p = p_mat_i(:,1);
            end
        else
            p = p_mat_i(:, j); % [bar]
        end

        % time vector:
        N_t = length(p);
        time_vec = 0 : 0.001 : (N_t - 1)*0.001;
        [p_max, idx_max] = max(p(200:end)); % find max pressure
        idx_max = idx_max+199;

        % find action interval times tA, tD:

        % tA:
        p_i = inf;
        k = 0;
        while p_i > toll*p_max
            k = k+1;
            p_i = p(idx_max - k);
        end
        idxA = idx_max - k;
        tA = time_vec(idxA);

        % tD:
        p_i = inf;
        k = 0;
        while p_i > toll*p_max
            k = k+1;
            p_i = p(idx_max + k);
        end
        idxD = idx_max + k;
        tD = time_vec(idxD);


        % Integration p(t) in [tA, tD]
        intg = trapz(time_vec(idxA:idxD), p(idxA:idxD));
        p_mean = intg/(tD-tA);
        p_ref = p_mean/2;

        % 2) BURNING TIME INTERVAL
        % similar to above
        % find burning interval times tB, tC:

        % tB:
        p_i = inf;
        k = 0;
        while p_i > p_ref
            k = k+1;
            p_i = p(idx_max - k);
        end
        idxB = idx_max - k;
        tB = time_vec(idxB);

        % tC:
        p_i = inf;
        k = 0;
        while p_i > p_ref
            k = k+1;
            p_i = p(idx_max + k);
        end
        idxC = idx_max + k;
        tC = time_vec(idxC);

        % Integral in [tB, tC]
        intg = trapz(time_vec(idxB:idxC), p(idxB:idxC));
        p_mean = intg/(tC-tB);

        % p_eff is the mean of p(t) in the burning time interval
        p_eff = p_mean;

        % burning rate:
        rb = wt/(tC-tB);

        % save data
        p_eff_vec (i,j) = p_eff;
        rb_vec (i,j) = rb;

        % Computation of c*
        % Integral in [tB, tC]
        c = (intg*1e5)*(A_t(j)*1e-6)/M_tot; % m/s ;
        c_vec(j+(i-1)*3) = c;

        if j == 1 
            c1 = [c1; c];
        else
            if j == 2
            c2 = [c2; c];
            else 
            c3 = [c3; c];
            end
        end

        % plots
        % subplot 3x9
        %if j == 2  % only for the MP, per non sovrapporre troppe linee. Se li vuoi vedere tutti 27, togli questo if
            subplot(3,3,i)
            plot(time_vec, p); hold on;
            plot([tA,tD], ones(1,2)* p_ref,'--');
            plot([tB,tC], ones(1,2)* p_eff,'--');
            plot(tA*ones(2,1), [0,p_ref],'--');
            plot(tD*ones(2,1), [0,p_ref],'--');
            plot(tB*ones(2,1), [0,p_eff],'--');
            plot(tC*ones(2,1), [0,p_eff],'--');

            title(['Batch ', num2str(i)])
            xlabel('Time [s]')
            ylabel('Pressure [bar]')
        %end

    end

    hold off;

end


rb_vec;
rb_LP = mean(rb_vec(:,1));
rb_MP = mean(rb_vec(:,2)); % 7.6888 mm/s
rb_HP = mean(rb_vec(:,3));

p_eff_vec;
p_eff_LP = mean(p_eff_vec(:,1));
p_eff_MP = mean(p_eff_vec(:,2));
p_eff_HP = mean(p_eff_vec(:,3));

% rb = 7.591 mm/s ok for MP, ok for solid propellants
% expected solutions:
% LP : 6 - 6.5 mm/s , 30 - 40 bar
% MP : 7.5 mm/s     , 40 - 50 bar
% HP : 8.5 - 9.5 m/s, 70 - 85 bar

p_val = [p_eff_vec(:,1); p_eff_vec(:,2); p_eff_vec(:,3)];
r_val = [rb_vec(:,1); rb_vec(:,2); rb_vec(:,3)];