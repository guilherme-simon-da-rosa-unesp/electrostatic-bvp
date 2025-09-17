function finite_difference_solution(a_cm, b_cm, V_0, N_x, N_y, num_max_iter, max_rel_diff, plot_iterations)
% FINITE_DIFFERENCE_SOLUTION  Solves the electric potential in a
% rectangular covered trough using a 2D iterative finite difference method
% and plots the results.
%
% DESCRIPTION:
%   This function computes the electric scalar potential Phi_e(x,y) in a
%   rectangular cross-section metallic tube ('covered trough') with
%   Dirichlet boundary conditions. Three sides are held at zero potential
%   and the top face has a sinusoidal voltage profile. The solution is
%   computed iteratively using the finite difference method.
%
% SYNTAX:
%   finite_difference_solution(a_cm, b_cm, V_0, N_x, N_y, num_max_iter,
%   max_rel_diff, plot_iterations)
%
% INPUTS:
%   a_cm            - Width of the trough [cm]
%   b_cm            - Height of the trough [cm]
%   V_0             - Maximum voltage at the top boundary [V]
%   N_x             - Number of grid points along the y-axis
%   N_y             - Number of grid points along the x-axis
%   num_max_iter    - Maximum number of iterations
%   max_rel_diff    - Maximum relative difference to stop iterations
%   plot_iterations - Logical: true to plot potential at each iteration
%
% EXAMPLE:
%   finite_difference_solution(16, 10, 8, 33, 21, 500, 1e-3, true)
%
% AUTHOR:
%   Guilherme S. Rosa - Last modified: 2025-09-15

%% Default parameters if not provided
if nargin < 8
    a_cm = 16; 
    b_cm = 10;
    V_0 = 8;
    N_x = 33;
    N_y = 21;
    num_max_iter = 500;
    max_rel_diff = 1e-3;
    plot_iterations = true;
end

% Convert cm to meters
cm = 1e-2;
a = a_cm * cm;
b = b_cm * cm;

% Grid generation
x = linspace(0, a, N_x);
y = linspace(0, b, N_y);
[X, Y] = meshgrid(x, y);

% Discrete points for interior nodes
x_points = 2:N_x-1;
y_points = 2:N_y-1;

delta_x = (a)/(N_x-1);
delta_y = (b)/(N_y-1);

% Initialize potential matrix
Phi_e = zeros(N_x, N_y);

% Dirichlet boundaries
Phi_e(1,:) = 0;          % bottom
Phi_e(N_x,:) = 0;        % top in x direction
Phi_e(:,1) = 0;          % left

% Top boundary sinusoidal voltage
Phi_e(:,N_y) = V_0 .* sin(linspace(0, pi, N_x));

%% Plot setup if plotting iterations
if plot_iterations
    figure;
    hold on;
    fig = imagesc(x/cm, y/cm, Phi_e');  % Initial plot
    set(gca,'YDir','normal');     % so Y increases upwards
    colormap(parula);
    clim([0 V_0]);                % set color limits from 0 to V_0    
    axis([0 a/cm 0 b/cm]);
    xlabel('Position $x$ (cm)', 'Interpreter', 'LaTeX');
    ylabel('Position $y$ (cm)', 'Interpreter', 'LaTeX');
    h_title = title('Electric scalar potential $\Phi_e$ (V)', 'Interpreter', 'LaTeX');
    axis equal;
    hcb = colorbar;
    hcb.TickLabelInterpreter = 'latex';
    xticks(linspace(0, a, 5)/cm);
    yticks(linspace(0, b, 5)/cm);
    xlim([0 a]./cm);
    ylim([0 b]./cm);
    grid on; box on;
    format_fig();
end

Phi_e_new = Phi_e;

%% Iterative solution
for ind_iter = 1:num_max_iter
    if delta_x == delta_y
        % Simple averaging for square grid
        Phi_e_new(x_points, y_points) = 0.25 * ( ...
            Phi_e(x_points+1, y_points) + Phi_e(x_points-1, y_points) + ...
            Phi_e(x_points, y_points+1) + Phi_e(x_points, y_points-1) );
    else
        % Weighted for non-square grid
        Phi_e_new(x_points, y_points) = 0.5 / (1/delta_x^2 + 1/delta_y^2) * ( ...
            (Phi_e(x_points+1, y_points) + Phi_e(x_points-1, y_points))/delta_x^2 + ...
            (Phi_e(x_points, y_points+1) + Phi_e(x_points, y_points-1))/delta_y^2 );
    end

    % Enforce boundary conditions
    Phi_e_new(1,:) = 0;
    Phi_e_new(N_x,:) = 0;
    Phi_e_new(:,1) = 0;
    Phi_e_new(:,N_y) = V_0 .* sin(linspace(0, pi, N_x));

    % Plot iteration if required
    if plot_iterations
        set(h_title, ...
        'String', sprintf('Electric scalar potential $\\Phi_e$ (V), iteration %d', ind_iter), ...
        'Interpreter', 'LaTeX');
        set(fig, 'CData', Phi_e');

        drawnow;
        pause(0.1);
    end

    % Check convergence
    rel_diff = max(abs(Phi_e_new - Phi_e) ./ Phi_e_new, [], 'all');
    if rel_diff < max_rel_diff
        fprintf('Converged after %d iterations, relative difference = %.3e\n', ind_iter, rel_diff);
        break;
    end

    Phi_e = Phi_e_new;
end

fprintf('Total iterations: %d, final relative difference: %.3e\n', ind_iter, rel_diff);

% Plotting
figure;
hold on;

imagesc(x/cm, y/cm, Phi_e');  % map your physical coordinates
set(gca,'YDir','normal');     % so Y increases upwards
colormap(parula);
clim([0 V_0]);                % set color limits from 0 to V_0

% Overlay dashed contour lines with labels
[C, h] = contour(X/cm, Y/cm, Phi_e', 0:1:V_0, 'k--', 'LineWidth', 1);
clabel(C, h, 'FontSize', 16, 'LabelSpacing', 300, 'Interpreter','latex');

xlabel('Position $x$ (cm)', 'Interpreter', 'LaTeX');
ylabel('Position $y$ (cm)', 'Interpreter', 'LaTeX');
title(sprintf('Electric scalar potential $\\Phi_e$ (V), iteration %d', ind_iter), 'Interpreter', 'LaTeX');
axis equal;
colorbar('TickLabelInterpreter', 'latex');
grid on;
box on;
xticks(linspace(0, a, 5) ./ cm);
yticks(linspace(0, b, 5) ./ cm);

% Apply formatting
format_fig();

end

%% Helper function for figure formatting
function format_fig()
    % Apply consistent figure formatting for IEEE
    scale = 2;

    % Set figure size to IEEE single-column dimensions (8.85 x 4 cm)
    set(gcf, 'Units', 'centimeters', 'Position', 2.*[1 1 8.85 4]);

    % Update all font sizes and interpreters
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 8*scale);
    set(findall(gcf, '-property', 'Interpreter'), 'Interpreter', 'LaTeX');
    set(findall(gca, '-property', 'FontSize'), 'FontSize', 8*scale);
    set(findall(gca, '-property', 'Interpreter'), 'Interpreter', 'LaTeX');
    set(findall(gcf, 'type', 'text'), 'FontSize', 8*scale);

    % Use LaTeX for tick labels
    set(gca, 'TickLabelInterpreter', 'latex');

    set(gca, 'LooseInset', [0.05 0.05 0.05 0.05]);
end
