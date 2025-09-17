function analytical_solution(a_cm, b_cm, V_0, N_x, N_y)
% ANALYTICAL_SOLUTION  Plot the analytical solution of the scalar electric
% potential in a rectangular covered trough with Dirichlet boundary
% conditions.
%
% DESCRIPTION:
%   This function computes and plots the electric scalar potential 
%   Phi_e(x, y) governed by Laplace's equation in a rectangular a x b
%   cross-section metallic tube ('covered trough'). Three sides are held at
%   zero potential, and the top face is excited with a sinusoidal voltage
%   profile.
%
%   The solution is derived using the method of separation of variables.
%
% SYNTAX:
%   analytical_solution(a_cm, b_cm, V_0, N_x, N_y)
%
% INPUTS:
%   a_cm  - Width of the trough in centimeters [cm]
%   b_cm  - Height of the trough in centimeters [cm]
%   V_0   - Maximum voltage at the top boundary [V]
%   N_x   - Number of grid points along the y-axis used for plotting
%   N_y   - Number of grid points along the x-axis used for plotting
%
% NOTES:
%   If no inputs are provided, the function uses default values:
%     a_cm = 16 cm, b_cm = 10 cm, V_0 = 8 V, N_x = 100, N_y = 100
%
% EXAMPLE:
%   analytical_solution(16, 10, 8, 100, 100)
%   analytical_solution()  % Uses default values
%
% AUTHOR:
%   Guilherme S. Rosa - Last modified: 2025-09-15

%% Default parameters if no inputs are given
if nargin < 5
    a_cm = 16;   % Width [cm]
    b_cm = 10;   % Height [cm]
    V_0 = 8;     % Max voltage [V]
    N_x = 100;    % Points for plotting along y
    N_y = 100;    % Points for plotting along x
end

% Convert dimensions to meters
cm = 1e-2;
a = a_cm * cm;
b = b_cm * cm;

% Grid generation
x = linspace(0, a, N_y);
y = linspace(0, b, N_x);
[X, Y] = meshgrid(x, y);

% Analytical potential (method of separation of variables)
Phi_e = @(x, y) V_0 ./ sinh(pi/a * b) .* sin(pi/a .* x) .* sinh(pi/a .* y);
Phi_e_mesh = Phi_e(X, Y);

% Plotting
figure;
hold on;

% 2D color map
imagesc(x/cm, y/cm, Phi_e_mesh);  % transpose to match axes
set(gca,'YDir','normal');         % correct Y-axis direction
colormap(parula);
clim([0 V_0]);                    % set color limits from 0 to V_0

% Overlay dashed contour lines with labels
[C, h] = contour(X/cm, Y/cm, Phi_e_mesh, 0:1:V_0, 'k--', 'LineWidth', 1);
clabel(C, h, 'FontSize', 16, 'LabelSpacing', 300, 'Interpreter', 'latex');

xlabel('Position $x$ (cm)', 'Interpreter', 'LaTeX');
ylabel('Position $y$ (cm)', 'Interpreter', 'LaTeX');
title('Electric scalar potential $\Phi_e$ (V)', 'Interpreter', 'LaTeX');
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
