static_years = [2015, 2016, 2021, 2023];
mag_years = [2021, 2022, 2023, 2024];

t_mag_all = [];
Bx_all = []; By_all = []; Bz_all = [];
X_all = []; Y_all = []; Z_all = [];

for y = mag_years
    var_name = ['MVN_MAG_' num2str(y)];
    if exist(var_name, 'var')
        m = eval(var_name);
        md = m.MAG_Data;
        t_mag_all = [t_mag_all; md.Time];
        Bx_all = [Bx_all; double(md.B_MSO(:,1))];
        By_all = [By_all; double(md.B_MSO(:,2))];
        Bz_all = [Bz_all; double(md.B_MSO(:,3))];
        X_all = [X_all; double(md.Pos_MSO(:,1)) / 3389.5];
        Y_all = [Y_all; double(md.Pos_MSO(:,2)) / 3389.5];
        Z_all = [Z_all; double(md.Pos_MSO(:,3)) / 3389.5];
    end
end

dt_seconds = 8;
t_start = min(t_mag_all);
t_end = max(t_mag_all);
t_binned_edges = (t_start : (dt_seconds/86400) : t_end)';

[~, ~, mag_bin_idx] = histcounts(t_mag_all, t_binned_edges);
num_bins = numel(t_binned_edges) - 1;
t_grid = 0.5 * (t_binned_edges(1:end-1) + t_binned_edges(2:end));

Bx_grid = nan(num_bins, 1); By_grid = nan(num_bins, 1); Bz_grid = nan(num_bins, 1);
X_grid  = nan(num_bins, 1); Y_grid  = nan(num_bins, 1); Z_grid  = nan(num_bins, 1);

for b = 1:num_bins
    m_mask = (mag_bin_idx == b);
    if any(m_mask)
        Bx_grid(b) = mean(Bx_all(m_mask), 'omitnan');
        By_grid(b) = mean(By_all(m_mask), 'omitnan');
        Bz_grid(b) = mean(Bz_all(m_mask), 'omitnan');
        X_grid(b)  = mean(X_all(m_mask), 'omitnan');
        Y_grid(b)  = mean(Y_all(m_mask), 'omitnan');
        Z_grid(b)  = mean(Z_all(m_mask), 'omitnan');
    end
end

t_ion_raw = []; Vx_raw = []; Vy_raw = []; Vz_raw = []; n_raw = []; Ti_raw = [];
DENSITY_FIELD = '';

for y = static_years
    var_name = ['STATIC_MOM_' num2str(y)];
    if exist(var_name, 'var')
        s = eval(var_name);
        sd = s.STATIC_Data;
    else
        matpath = ['STATIC_MOM_' num2str(y) '.mat'];
        if ~isfile(matpath)
            continue;
        end
        temp_load = load(matpath);
        sd = temp_load.STATIC_Data;
    end

    if isempty(DENSITY_FIELD)
        candidates = {'den_protons_high', 'density', 'density_protons', 'n_p', 'n_tot', 'density_tot', 'abundance_protons'};
        for c = 1:numel(candidates)
            if isfield(sd, candidates{c})
                DENSITY_FIELD = candidates{c};
                break;
            end
        end
    end

    if ~isempty(DENSITY_FIELD) && isfield(sd, DENSITY_FIELD)
        t_ion_raw = [t_ion_raw; sd.Time];
        Vx_raw = [Vx_raw; double(sd.vel_xyz_protons_high(:,1))];
        Vy_raw = [Vy_raw; double(sd.vel_xyz_protons_high(:,2))];
        Vz_raw = [Vz_raw; double(sd.vel_xyz_protons_high(:,3))];
        n_this = double(sd.(DENSITY_FIELD));
        if size(n_this, 2) > 1, n_this = n_this(:,1); end
        n_raw = [n_raw; n_this];

        if isfield(sd, 'temperature_protons')
            t_this = double(sd.temperature_protons);
            if size(t_this, 2) > 1, t_this = t_this(:,1); end
            Ti_raw = [Ti_raw; t_this];
        else
            Ti_raw = [Ti_raw; nan(size(n_this))];
        end
    end
end

[t_ion_raw, idx_sort] = sort(t_ion_raw);
Vx_raw = Vx_raw(idx_sort); Vy_raw = Vy_raw(idx_sort); Vz_raw = Vz_raw(idx_sort);
n_raw = n_raw(idx_sort); Ti_raw = Ti_raw(idx_sort);

[t_ion_u, idx_u] = unique(t_ion_raw);
Vx_u = Vx_raw(idx_u); Vy_u = Vy_raw(idx_u); Vz_u = Vz_raw(idx_u);
n_u = n_raw(idx_u); Ti_u = Ti_raw(idx_u);

Vx_grid = interp1(t_ion_u, Vx_u, t_grid, 'linear', NaN);
Vy_grid = interp1(t_ion_u, Vy_u, t_grid, 'linear', NaN);
Vz_grid = interp1(t_ion_u, Vz_u, t_grid, 'linear', NaN);
n_grid  = interp1(t_ion_u, n_u,  t_grid, 'linear', NaN);
Ti_grid = interp1(t_ion_u, Ti_u, t_grid, 'linear', NaN);

Ex_grid = -1e-3 * (Vy_grid .* Bz_grid - Vz_grid .* By_grid);
Ey_grid = -1e-3 * (Vz_grid .* Bx_grid - Vx_grid .* Bz_grid);
Ez_grid = -1e-3 * (Vx_grid .* By_grid - Vy_grid .* Bx_grid);
E_mag_grid = sqrt(Ex_grid.^2 + Ey_grid.^2 + Ez_grid.^2);

B_mag_grid = sqrt(Bx_grid.^2 + By_grid.^2 + Bz_grid.^2);
V_mag_grid = sqrt(Vx_grid.^2 + Vy_grid.^2 + Vz_grid.^2);
Flux_grid  = n_grid .* V_mag_grid;
P_dyn_grid = (1.6726e-6) * n_grid .* (V_mag_grid.^2);

Va_grid = 21.8 * (B_mag_grid ./ sqrt(max(n_grid, 1e-4)));
Ma_grid = V_mag_grid ./ max(Va_grid, 1e-4);

valid_mask = ~isnan(Ex_grid) & ~isnan(Ey_grid) & ~isnan(Ez_grid) & ...
             (E_mag_grid < 50) & (abs(Ex_grid) < 20) & (abs(Ey_grid) < 20) & (abs(Ez_grid) < 20);

X = X_grid(valid_mask); Y = Y_grid(valid_mask); Z = Z_grid(valid_mask);
Bx = Bx_grid(valid_mask); By = By_grid(valid_mask); Bz = Bz_grid(valid_mask); Bm = B_mag_grid(valid_mask);
Vx = Vx_grid(valid_mask); Vy = Vy_grid(valid_mask); Vz = Vz_grid(valid_mask); Vm = V_mag_grid(valid_mask);
Ex = Ex_grid(valid_mask); Ey = Ey_grid(valid_mask); Ez = Ez_grid(valid_mask); Em = E_mag_grid(valid_mask);
n = n_grid(valid_mask); Flux = Flux_grid(valid_mask); P_dyn = P_dyn_grid(valid_mask);
Ti = Ti_grid(valid_mask); Ma = Ma_grid(valid_mask);

x_edges = -2:0.2:3; y_edges = -2:0.2:2; z_edges = -2:0.2:2;

plot_quad_map(X, Y, Z, Bx, By, Bz, Bm, x_edges, y_edges, z_edges, 'XY', 'Figure 1: B Map: XY', 'B_x (nT)', 'B_y (nT)', 'B_z (nT)', '|B| (nT)', [-8 8], [0 35]);
plot_quad_map(X, Y, Z, Bx, By, Bz, Bm, x_edges, y_edges, z_edges, 'XZ', 'Figure 2: B Map: XZ', 'B_x (nT)', 'B_y (nT)', 'B_z (nT)', '|B| (nT)', [-6 4], [0 22]);
plot_quad_map(X, Y, Z, Bx, By, Bz, Bm, x_edges, y_edges, z_edges, 'YZ', 'Figure 3: B Map: YZ', 'B_x (nT)', 'B_y (nT)', 'B_z (nT)', '|B| (nT)', [-2 3], [0 22]);

plot_single_map(X, Y, Z, n, x_edges, y_edges, z_edges, 'XY', 'Figure 4: n Map: XY', 'n (cm^{-3})', [0 10], 'mean');
plot_single_map(X, Y, Z, n, x_edges, y_edges, z_edges, 'XZ', 'Figure 5: n Map: XZ', 'n (cm^{-3})', [0 10], 'mean');
plot_single_map(X, Y, Z, n, x_edges, y_edges, z_edges, 'YZ', 'Figure 6: n Map: YZ', 'n (cm^{-3})', [0 10], 'mean');

plot_quad_map(X, Y, Z, Vx, Vy, Vz, Vm, x_edges, y_edges, z_edges, 'XY', 'Figure 7: V Map: XY', 'V_x (km/s)', 'V_y (km/s)', 'V_z (km/s)', '|V| (km/s)', [-400 100], [0 450]);
plot_quad_map(X, Y, Z, Vx, Vy, Vz, Vm, x_edges, y_edges, z_edges, 'XZ', 'Figure 8: V Map: XZ', 'V_x (km/s)', 'V_y (km/s)', 'V_z (km/s)', '|V| (km/s)', [-400 0], [0 450]);
plot_quad_map(X, Y, Z, Vx, Vy, Vz, Vm, x_edges, y_edges, z_edges, 'YZ', 'Figure 9: V Map: YZ', 'V_x (km/s)', 'V_y (km/s)', 'V_z (km/s)', '|V| (km/s)', [-450 0], [100 450]);

plot_single_map(X, Y, Z, Flux, x_edges, y_edges, z_edges, 'XY', 'Figure 10: Flux Map: XY', 'Flux  (cm^{-3} \cdot km/s)', [0 2800], 'mean');
plot_single_map(X, Y, Z, Flux, x_edges, y_edges, z_edges, 'XZ', 'Figure 11: Flux Map: XZ', 'Flux  (cm^{-3} \cdot km/s)', [0 2500], 'mean');
plot_single_map(X, Y, Z, Flux, x_edges, y_edges, z_edges, 'YZ', 'Figure 12: Flux Map: YZ', 'Flux  (cm^{-3} \cdot km/s)', [0 2500], 'mean');

plot_quad_map(X, Y, Z, Ex, Ey, Ez, Em, x_edges, y_edges, z_edges, 'XY', 'Figure 13: E Map: XY', 'E_x (mV/m)', 'E_y (mV/m)', 'E_z (mV/m)', '|E| (mV/m)', [-8 8], [0 12]);
plot_quad_map(X, Y, Z, Ex, Ey, Ez, Em, x_edges, y_edges, z_edges, 'XZ', 'Figure 14: E Map: XZ', 'E_x (mV/m)', 'E_y (mV/m)', 'E_z (mV/m)', '|E| (mV/m)', [-8 8], [0 12]);
plot_quad_map(X, Y, Z, Ex, Ey, Ez, Em, x_edges, y_edges, z_edges, 'YZ', 'Figure 15: E Map: YZ', 'E_x (mV/m)', 'E_y (mV/m)', 'E_z (mV/m)', '|E| (mV/m)', [-8 8], [0 12]);

upstream_mask = (X > 1.2) & ~isnan(P_dyn) & (P_dyn > 0) & (P_dyn < 10);
P_up = P_dyn(upstream_mask);
figure('Name', 'Figure 16: Solar Wind Dynamic Pressure', 'Color', 'w', 'Position', [150 150 800 600]);
histogram(P_up, 'BinWidth', 0.15, 'Normalization', 'probability', 'FaceColor', [0.45 0.62 0.82], 'EdgeColor', 'k');
hold on;
med_p = median(P_up, 'omitnan'); mean_p = mean(P_up, 'omitnan'); yl = ylim;
line([med_p med_p], yl, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
line([mean_p mean_p], yl, 'Color', 'g', 'LineStyle', '--', 'LineWidth', 2);
text(med_p - 0.05, yl(2)*0.75, ['Median = ' num2str(round(med_p*1000)/1000) ' nPa'], 'Color', 'r', 'Rotation', 90, 'FontWeight', 'bold', 'FontSize', 11);
text(mean_p + 0.05, yl(2)*0.75, ['Mean = ' num2str(round(mean_p*1000)/1000) ' nPa'], 'Color', 'g', 'Rotation', 90, 'FontWeight', 'bold', 'FontSize', 11);
xlabel('P_{sw,dy} (nPa)', 'FontSize', 12); ylabel('Occurrence Rate', 'FontSize', 12);
title('Solar Wind Dynamic Pressure Distribution (Upstream)', 'FontSize', 13, 'FontWeight', 'bold');
grid on; set(gca, 'FontSize', 11, 'LineWidth', 1);

plot_single_map(X, Y, Z, ones(size(X)), x_edges, y_edges, z_edges, 'XY', 'Figure 17: Occupancy Map: XY', 'Sample Counts', [0 500], 'count');
plot_single_map(X, Y, Z, ones(size(X)), x_edges, y_edges, z_edges, 'XZ', 'Figure 18: Occupancy Map: XZ', 'Sample Counts', [0 500], 'count');
plot_single_map(X, Y, Z, ones(size(X)), x_edges, y_edges, z_edges, 'YZ', 'Figure 19: Occupancy Map: YZ', 'Sample Counts', [0 500], 'count');

plot_single_map(X, Y, Z, Ti, x_edges, y_edges, z_edges, 'XY', 'Figure 20: Ion Temperature: XY', 'T_i (eV)', [0 200], 'mean');
plot_single_map(X, Y, Z, Ti, x_edges, y_edges, z_edges, 'XZ', 'Figure 21: Ion Temperature: XZ', 'T_i (eV)', [0 200], 'mean');
plot_single_map(X, Y, Z, Ti, x_edges, y_edges, z_edges, 'YZ', 'Figure 22: Ion Temperature: YZ', 'T_i (eV)', [0 200], 'mean');

plot_single_map(X, Y, Z, Ma, x_edges, y_edges, z_edges, 'XY', 'Figure 23: Mach Number: XY', 'M_A', [0 10], 'mean');
plot_single_map(X, Y, Z, Ma, x_edges, y_edges, z_edges, 'XZ', 'Figure 24: Mach Number: XZ', 'M_A', [0 10], 'mean');
plot_single_map(X, Y, Z, Ma, x_edges, y_edges, z_edges, 'YZ', 'Figure 25: Mach Number: YZ', 'M_A', [0 10], 'mean');

figure('Name', 'Figure 26: 3D Orbit Trajectories', 'Color', 'w', 'Position', [100 100 800 700]);
plot3(X, Y, Z, '.', 'Color', [0.2 0.4 0.8], 'MarkerSize', 2); hold on;
[sx, sy, sz] = sphere(50);
surf(sx, sy, sz, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
xlabel('X_{MSO} (R_m)'); ylabel('Y_{MSO} (R_m)'); zlabel('Z_{MSO} (R_m)');
title('MAVEN Spacecraft Orbit Coverage (MSO)', 'FontSize', 12, 'FontWeight', 'bold');
axis equal; grid on; view(3);

figure('Name', 'Figure 27: MAG Time Series', 'Color', 'w', 'Position', [100 100 900 500]);
plot(t_grid, Bx_grid, 'r', t_grid, By_grid, 'g', t_grid, Bz_grid, 'b', t_grid, B_mag_grid, 'k');
xlabel('Time'); ylabel('B (nT)'); title('Magnetic Field Vector Time Series');
legend('B_x', 'B_y', 'B_z', '|B|'); grid on; datetick('x', 'yyyy-mm', 'keeplimits');

figure('Name', 'Figure 28: STATIC Time Series', 'Color', 'w', 'Position', [100 100 900 500]);
yyaxis left; plot(t_grid, n_grid, 'b'); ylabel('Density n (cm^{-3})');
yyaxis right; plot(t_grid, V_mag_grid, 'r'); ylabel('Bulk Speed |V| (km/s)');
xlabel('Time'); title('STATIC Plasma Moments Time Series'); grid on; datetick('x', 'yyyy-mm', 'keeplimits');

figure('Name', 'Figure 29: Electric Field Time Series', 'Color', 'w', 'Position', [100 100 900 500]);
plot(t_grid, Ex_grid, 'r', t_grid, Ey_grid, 'g', t_grid, Ez_grid, 'b', t_grid, E_mag_grid, 'k');
xlabel('Time'); ylabel('E (mV/m)'); title('Convective Electric Field Time Series');
legend('E_x', 'E_y', 'E_z', '|E|'); grid on; datetick('x', 'yyyy-mm', 'keeplimits');

function plot_quad_map(X_data, Y_data, Z_data, V1, V2, V3, Vm, x_e, y_e, z_e, plane, fig_t, t1, t2, t3, tm, c_lim, m_lim)
    figure('Name', fig_t, 'Color', 'w', 'Position', [100 100 1200 800]);

    if strcmp(plane, 'XY')
        P1 = X_data; P2 = Y_data; e1 = x_e; e2 = y_e; x_lbl = 'X_{MSO} (R_m)'; y_lbl = 'Y_{MSO} (R_m)';
    elseif strcmp(plane, 'XZ')
        P1 = X_data; P2 = Z_data; e1 = x_e; e2 = z_e; x_lbl = 'X_{MSO} (R_m)'; y_lbl = 'Z_{MSO} (R_m)';
    else
        P1 = Y_data; P2 = Z_data; e1 = y_e; e2 = z_e; x_lbl = 'Y_{MSO} (R_m)'; y_lbl = 'Z_{MSO} (R_m)';
    end

    [g1, g2] = meshgrid(e1, e2);

    c1 = bin_2d_proj(P1, P2, V1, e1, e2, 'mean');
    c2 = bin_2d_proj(P1, P2, V2, e1, e2, 'mean');
    c3 = bin_2d_proj(P1, P2, V3, e1, e2, 'mean');
    cm = bin_2d_proj(P1, P2, Vm, e1, e2, 'mean');

    subplot(2,2,1); render_panel_proj(g1, g2, c1, t1, c_lim, x_lbl, y_lbl, plane);
    subplot(2,2,2); render_panel_proj(g1, g2, c2, t2, c_lim, x_lbl, y_lbl, plane);
    subplot(2,2,3); render_panel_proj(g1, g2, c3, t3, c_lim, x_lbl, y_lbl, plane);
    subplot(2,2,4); render_panel_proj(g1, g2, cm, tm, m_lim, x_lbl, y_lbl, plane);
end

function plot_single_map(X_data, Y_data, Z_data, Val_data, x_e, y_e, z_e, plane, fig_t, map_t, c_lim, mode_type)
    figure('Name', fig_t, 'Color', 'w', 'Position', [200 200 700 600]);

    if strcmp(plane, 'XY')
        P1 = X_data; P2 = Y_data; e1 = x_e; e2 = y_e; x_lbl = 'X_{MSO} (R_m)'; y_lbl = 'Y_{MSO} (R_m)';
    elseif strcmp(plane, 'XZ')
        P1 = X_data; P2 = Z_data; e1 = x_e; e2 = z_e; x_lbl = 'X_{MSO} (R_m)'; y_lbl = 'Z_{MSO} (R_m)';
    else
        P1 = Y_data; P2 = Z_data; e1 = y_e; e2 = z_e; x_lbl = 'Y_{MSO} (R_m)'; y_lbl = 'Z_{MSO} (R_m)';
    end

    [g1, g2] = meshgrid(e1, e2);
    grid_c = bin_2d_proj(P1, P2, Val_data, e1, e2, mode_type);

    render_panel_proj(g1, g2, grid_c, map_t, c_lim, x_lbl, y_lbl, plane);
end

function grid_val = bin_2d_proj(p1_data, p2_data, val_data, e1, e2, mode_type)
    grid_val = nan(numel(e2)-1, numel(e1)-1);
    [~, ~, bin_1] = histcounts(p1_data, e1);
    [~, ~, bin_2] = histcounts(p2_data, e2);

    for i = 1:(numel(e1)-1)
        for j = 1:(numel(e2)-1)
            idx = (bin_1 == i) & (bin_2 == j);
            if any(idx)
                if strcmp(mode_type, 'count')
                    grid_val(j, i) = sum(idx);
                else
                    grid_val(j, i) = mean(val_data(idx), 'omitnan');
                end
            end
        end
    end
end

function render_panel_proj(G1, G2, C_data, title_str, c_lims, x_lbl, y_lbl, plane)
    pcolor(G1(1:end-1, 1:end-1), G2(1:end-1, 1:end-1), C_data);
    shading flat;
    colormap(jet);
    caxis(c_lims);
    colorbar;
    hold on;

    th = linspace(0, 2*pi, 200);
    plot(cos(th), sin(th), 'k-', 'LineWidth', 2);

    if ~strcmp(plane, 'YZ')
        L_bs = 2.04; e_bs = 0.9; x0_bs = 0.78;
        r_bs = L_bs ./ (1 + e_bs * cos(th));
        plot(x0_bs + r_bs .* cos(th), r_bs .* sin(th), 'r--', 'LineWidth', 1.5);

        L_mpb = 0.96; e_mpb = 0.9; x0_mpb = 0.78;
        r_mpb = L_mpb ./ (1 + e_mpb * cos(th));
        plot(x0_mpb + r_mpb .* cos(th), r_mpb .* sin(th), 'm--', 'LineWidth', 1.5);
    end

    title(title_str, 'FontSize', 11, 'FontWeight', 'bold');
    xlabel(x_lbl);
    ylabel(y_lbl);
    axis equal;
    xlim([min(G1(:)) max(G1(:))]);
    ylim([min(G2(:)) max(G2(:))]);
    set(gca, 'FontSize', 10, 'LineWidth', 1);
end
