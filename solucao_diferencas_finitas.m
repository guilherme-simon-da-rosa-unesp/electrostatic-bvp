function solucao_diferencas_finitas(a_cm, b_cm, V_0, N_x, N_y, num_max_iter, max_rel_diff, plot_iterations)
% SOLUCAO_DIFERENCAS_FINITAS  Resolve o potencial escalr elétrico em uma
% calha coberta usando método iterativo 2D de diferenças finitas e plota os
% resultados.
%
% DESCRIÇÃO:
%   Esta função calcula o potencial escalar elétrico Phi_e(x,y) em um tubo
%   metálico com seção transversal retangular ('canal coberto') com
%   condições de contorno de Dirichlet. Três lados são mantidos em
%   potencial zero e a face superior tem um perfil senoidal de tensão. A
%   solução é obtida iterativamente usando o método das diferenças finitas.
%
% SINTAXE:
%   solucao_diferencas_finitas(a_cm, b_cm, V_0, N_x, N_y, num_max_iter,
%   max_rel_diff, plot_iterations)
%
% ENTRADAS:
%   a_cm            - Largura do canal [cm] 
%   b_cm            - Altura do
%   canal [cm] V_0  - Tensão máxima na face superior [V] 
%   N_x             - Número de pontos na direção y
%   N_y             - Número de pontos na direção x
%   num_max_iter    - Número máximo de iterações
%   max_rel_diff    - Diferença relativa máxima para parar iterações
%   plot_iterations - Lógico: true para plotar o potencial a cada iteração
%
% EXEMPLO:
%   solucao_diferencas_finitas(16, 10, 8, 33, 21, 500, 1e-3, true)
%
% AUTOR:
%   Guilherme S. Rosa - Última modificação: 15/09/2025

%% Parâmetros padrão, caso não sejam fornecidos
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

% Conversão de cm para metros
cm = 1e-2;
a = a_cm * cm;
b = b_cm * cm;

% Geração da malha
x = linspace(0, a, N_x);
y = linspace(0, b, N_y);
[X, Y] = meshgrid(x, y);

% Pontos internos discretos
x_points = 2:N_x-1;
y_points = 2:N_y-1;

delta_x = a/(N_x-1);
delta_y = b/(N_y-1);

% Inicializa matriz do potencial
Phi_e = zeros(N_x, N_y);

% Condições de contorno de Dirichlet
Phi_e(1,:) = 0;          % inferior
Phi_e(N_x,:) = 0;        % superior na direção x
Phi_e(:,1) = 0;          % esquerda

% Tensão senoidal na face superior
Phi_e(:,N_y) = V_0 .* sin(linspace(0, pi, N_x));

%% Configuração do plot caso plotar iterações
if plot_iterations
    figure;
    hold on;
    fig = imagesc(x/cm, y/cm, Phi_e');  % Plot inicial
    set(gca,'YDir','normal');           % Y cresce para cima
    colormap(parula);
    clim([0 V_0]);                      % escala de cor de 0 a V_0    
    axis([0 a/cm 0 b/cm]);
    xlabel("Posi\c{c}\~{a}o $x$ (cm)", 'Interpreter', 'LaTeX');
    ylabel("Posi\c{c}\~{a}o $y$ (cm)", 'Interpreter', 'LaTeX');
    h_title = title("Potencial el\'{e}trico escalar $\Phi_e$ (V)", 'Interpreter', 'LaTeX');
    axis equal;
    hcb = colorbar;
    hcb.TickLabelInterpreter = 'latex';
    xticks(linspace(0, a, 5)/cm);
    yticks(linspace(0, b, 5)/cm);
    xlim([0 a]./cm);
    ylim([0 b]./cm);
    grid on; box on;
    formatar_figura();
end

Phi_e_new = Phi_e;

%% Solução iterativa
for ind_iter = 1:num_max_iter
    if delta_x == delta_y
        % Média simples para malha quadrada
        Phi_e_new(x_points, y_points) = 0.25 * ( ...
            Phi_e(x_points+1, y_points) + Phi_e(x_points-1, y_points) + ...
            Phi_e(x_points, y_points+1) + Phi_e(x_points, y_points-1) );
    else
        % Média ponderada para malha não quadrada
        Phi_e_new(x_points, y_points) = 0.5 / (1/delta_x^2 + 1/delta_y^2) * ( ...
            (Phi_e(x_points+1, y_points) + Phi_e(x_points-1, y_points))/delta_x^2 + ...
            (Phi_e(x_points, y_points+1) + Phi_e(x_points, y_points-1))/delta_y^2 );
    end

    % Reforçar condições de contorno
    Phi_e_new(1,:) = 0;
    Phi_e_new(N_x,:) = 0;
    Phi_e_new(:,1) = 0;
    Phi_e_new(:,N_y) = V_0 .* sin(linspace(0, pi, N_x));

    % Plotar iteração se necessário
    if plot_iterations
        set(h_title, ...
            'String', [ strcat( "Potencial escalar el\'{e}trico $\Phi_e$ (V), itera\c{c}\~{a}o ", num2str(ind_iter) ) ], ...
            'Interpreter', 'LaTeX');
        set(fig, 'CData', Phi_e');

        drawnow;
        pause(0.1);
    end

    % Verificar convergência
    rel_diff = max(abs(Phi_e_new - Phi_e) ./ Phi_e_new, [], 'all');
    if rel_diff < max_rel_diff
        fprintf('Convergência após %d iterações, diferença relativa = %.3e\n', ind_iter, rel_diff);
        break;
    end

    Phi_e = Phi_e_new;
end

fprintf('Total de iterações: %d, diferença relativa final: %.3e\n', ind_iter, rel_diff);

%% Plot final
figure;
hold on;

imagesc(x/cm, y/cm, Phi_e');  % mapa de coordenadas físicas
set(gca,'YDir','normal');     % Y cresce para cima
colormap(parula);
clim([0 V_0]);                % escala de cor de 0 a V_0

% Linhas de contorno com rótulos
[C, h] = contour(X/cm, Y/cm, Phi_e', 0:1:V_0, 'k--', 'LineWidth', 1);
clabel(C, h, 'FontSize', 16, 'LabelSpacing', 300, 'Interpreter','latex');

xlabel("Posi\c{c}\~{a}o $x$ (cm)", 'Interpreter', 'LaTeX');
ylabel("Posi\c{c}\~{a}o $y$ (cm)", 'Interpreter', 'LaTeX');
title([ strcat( "Potencial escalar el\'{e}trico $\Phi_e$ (V), itera\c{c}\~{a}o ", num2str(ind_iter) ) ], 'Interpreter', 'LaTeX');
axis equal;
colorbar('TickLabelInterpreter', 'latex');
grid on;
box on;
xticks(linspace(0, a, 5) ./ cm);
yticks(linspace(0, b, 5) ./ cm);

% Aplicar formatação
formatar_figura();

end

%% Função auxiliar para formatação de figuras
function formatar_figura()
    % Aplica formatação consistente para padrões IEEE
    escala = 2;

    % Define tamanho da figura (coluna simples: 8.85 x 4 cm)
    set(gcf, 'Units', 'centimeters', 'Position', 2.*[1 1 8.85 4]);

    % Atualiza tamanho de fontes e interpretadores
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 8*escala);
    set(findall(gcf, '-property', 'Interpreter'), 'Interpreter', 'LaTeX');
    set(findall(gca, '-property', 'FontSize'), 'FontSize', 8*escala);
    set(findall(gca, '-property', 'Interpreter'), 'Interpreter', 'LaTeX');
    set(findall(gcf, 'type', 'text'), 'FontSize', 8*escala);

    % Usa LaTeX para rótulos dos eixos
    set(gca, 'TickLabelInterpreter', 'latex');

    set(gca, 'LooseInset', [0.05 0.05 0.05 0.05]);
end
