function solucao_analitica(a_cm, b_cm, V_0, N_x, N_y)
% SOLUCAO_ANALITICA  Plota a solução analítica do potencial escalar
% elétrico em um canal retangular metálico com condições de contorno de
% Dirichlet.
%
% DESCRIÇÃO:
%   Esta função calcula e plota o potencial escalar elétrico
%   Phi_e(x, y), governado pela equação de Laplace, em um tubo metálico
%   com seção transversal retangular de dimensões a x b.
%   Três lados estão mantidos em potencial nulo, e a face superior
%   é excitada com um perfil senoidal de tensão.
%
%   A solução é obtida usando o método da separação de variáveis.
%
% SINTAXE:
%   solucao_analitica(a_cm, b_cm, V_0, N_x, N_y)
%
% ENTRADAS:
%   a_cm  - Largura do canal em centímetros [cm]
%   b_cm  - Altura do canal em centímetros [cm]
%   V_0   - Tensão máxima aplicada na face superior [V]
%   N_x   - Número de pontos na direção y usados para o gráfico
%   N_y   - Número de pontos na direção x usados para o gráfico
%
% OBSERVAÇÕES:
%   Se nenhuma entrada for fornecida, a função utiliza valores padrão:
%     a_cm = 16 cm, b_cm = 10 cm, V_0 = 8 V, N_x = 100, N_y = 100
%
% EXEMPLO:
%   solucao_analitica(16, 10, 8, 100, 100)
%   solucao_analitica()  % Usa os valores padrão
%
% AUTOR:
%   Guilherme S. Rosa - Última modificação: 15/09/2025

%% Parâmetros padrão, caso nenhuma entrada seja fornecida
if nargin < 5
    a_cm = 16;   % Largura [cm]
    b_cm = 10;   % Altura [cm]
    V_0 = 8;     % Tensão máxima [V]
    N_x = 100;   % Pontos para gráfico ao longo de y
    N_y = 100;   % Pontos para gráfico ao longo de x
end

% Conversão de unidades para metros
cm = 1e-2;
a = a_cm * cm;
b = b_cm * cm;

% Geração da malha
x = linspace(0, a, N_y);
y = linspace(0, b, N_x);
[X, Y] = meshgrid(x, y);

% Solução analítica (método da separação de variáveis)
Phi_e = @(x, y) V_0 ./ sinh(pi/a * b) .* sin(pi/a .* x) .* sinh(pi/a .* y);
Phi_e_malha = Phi_e(X, Y);

% Plotagem
figure;
hold on;

% Mapa de cores 2D
imagesc(x/cm, y/cm, Phi_e_malha);  % transpor para ajustar aos eixos
set(gca,'YDir','normal');          % corrigir direção do eixo Y
colormap(parula);
clim([0 V_0]);                     % limitar escala de cor de 0 a V_0

% Linhas de contorno com rótulos
[C, h] = contour(X/cm, Y/cm, Phi_e_malha, 0:1:V_0, 'k--', 'LineWidth', 1);
clabel(C, h, 'FontSize', 16, 'LabelSpacing', 300, 'Interpreter', 'latex');

xlabel("Posi\c{c}\~{a}o $x$ (cm)", 'Interpreter', 'LaTeX');
ylabel("Posi\c{c}\~{a}o $y$ (cm)", 'Interpreter', 'LaTeX');
title("Potencial el\'{e}trico escalar $\Phi_e$ (V)", 'Interpreter', 'LaTeX');
axis equal;
colorbar('TickLabelInterpreter', 'latex');
grid on;
box on;
xticks(linspace(0, a, 5) ./ cm);
yticks(linspace(0, b, 5) ./ cm);

% Aplicar formatação
formatar_figura();

end

%% Função auxiliar para formatação do gráfico
function formatar_figura()
    % Aplica formatação consistente ao gráfico para padrões IEEE
    escala = 2;

    % Define o tamanho da figura (coluna simples: 8.85 x 4 cm)
    set(gcf, 'Units', 'centimeters', 'Position', 2.*[1 1 8.85 4]);

    % Atualiza tamanho das fontes e interpretadores
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 8*escala);
    set(findall(gcf, '-property', 'Interpreter'), 'Interpreter', 'LaTeX');
    set(findall(gca, '-property', 'FontSize'), 'FontSize', 8*escala);
    set(findall(gca, '-property', 'Interpreter'), 'Interpreter', 'LaTeX');
    set(findall(gcf, 'type', 'text'), 'FontSize', 8*escala);

    % Usa LaTeX para os rótulos dos eixos
    set(gca, 'TickLabelInterpreter', 'latex');
    set(gca, 'LooseInset', [0.05 0.05 0.05 0.05]);
end
