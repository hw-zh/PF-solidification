function plot_scheil_path(pathNames, xRange, yRange)
% plot_scheil_path
% -------------------------------------------------------------------------
% 功能：
%   从当前工作文件夹中自动读取第一个 Excel 文件，
%   根据 Thermo-Calc / Scheil 导出的列数据，绘制凝固路径图。
%
% 输入：
%   pathNames : 你想画的 Scheil 路径名称，写成 cell 或 string，例如：
%               {'LIQUID','LIQUID + BCC_A2','LIQUID + BCC_A2 + FCC_A1'}
%
%   xRange    : 横坐标范围，例如 [0 1]
%               不输入则自动
%
%   yRange    : 纵坐标范围，例如 [1400 1440]
%               不输入则自动
%
% 用法示例：
%   plot_scheil_path({'LIQUID','LIQUID + BCC_A2','LIQUID + BCC_A2 + FCC_A1'})
%
%   plot_scheil_path({'LIQUID','LIQUID + BCC_A2','LIQUID + BCC_A2 + FCC_A1'}, [0 1], [1400 1440])
%
% 说明：
%   1. 当前代码按照你截图中的 Excel 结构写：
%         A:B, C:D, E:F, G:H, I:J   -> 平衡分段
%         K:L                       -> LIQUID
%         M:N                       -> LIQUID + BCC_A2
%         O:P                       -> LIQUID + BCC_A2 + FCC_A1
%
%   2. 如果你后续导出的 Excel 列顺序变了，只需修改下面的 eqPairs 和
%      scheilDefs 即可。
% -------------------------------------------------------------------------

    % =========================
    % 1. 输入处理
    % =========================
    if nargin < 1 || isempty(pathNames)
        error(['请输入要绘制的路径名称，例如：', newline, ...
               'plot_scheil_path({''LIQUID'',''LIQUID + BCC_A2'',''LIQUID + BCC_A2 + FCC_A1''})']);
    end

    if nargin < 2
        xRange = [];
    end

    if nargin < 3
        yRange = [];
    end

    if isstring(pathNames)
        pathNames = cellstr(pathNames);
    end

    if ischar(pathNames)
        pathNames = {pathNames};
    end

    % =========================
    % 2. 自动读取当前文件夹中的第一个 Excel
    % =========================
    files = [dir('*.xlsx'); dir('*.xls')];

    if isempty(files)
        error('当前文件夹中没有找到 Excel 文件。');
    end

    excelFile = files(1).name;
    fprintf('当前读取的 Excel 文件：%s\n', excelFile);

    % 读数值矩阵
    A = readmatrix(excelFile);

    % =========================
    % 3. 根据你当前截图定义列组
    % =========================
    % 平衡路径：前5组（A:B, C:D, E:F, G:H, I:J）
    eqPairs = [
        1 2;
        3 4;
        5 6;
        7 8;
        9 10
    ];

    % Scheil 路径定义：
    % 每一行格式：{名称, x列, y列, 颜色}
    scheilDefs = {
        'LIQUID',                    11, 12, [0.00 0.45 0.74];   % 蓝
        'LIQUID + BCC_A2',          13, 14, [1.00 0.00 0.00];   % 红
        'LIQUID + BCC_A2 + FCC_A1',15, 16, [0.00 0.50 0.00]    % 绿
    };

    % =========================
    % 4. 绘图
    % =========================
    figure('Color', 'w');
    hold on;

    % ---- 先画平衡线（灰色虚线）----
    for k = 1:size(eqPairs,1)
        xcol = eqPairs(k,1);
        ycol = eqPairs(k,2);

        if xcol <= size(A,2) && ycol <= size(A,2)
            x = A(:, xcol);
            y = A(:, ycol);

            idx = ~isnan(x) & ~isnan(y);

            if any(idx)
                plot(x(idx), y(idx), ':', ...
                    'Color', [0.55 0.55 0.55], ...
                    'LineWidth', 1.5, ...
                    'HandleVisibility', 'off');
            end
        end
    end

    % 为图例补一个“平衡”
    plot(nan, nan, ':', ...
        'Color', [0.55 0.55 0.55], ...
        'LineWidth', 1.5, ...
        'DisplayName', '平衡');

    % ---- 再画你指定的 Scheil 路径 ----
    foundPath = false(size(pathNames));

    for i = 1:numel(pathNames)
        name_i = pathNames{i};

        for j = 1:size(scheilDefs,1)
            defName = scheilDefs{j,1};

            if strcmp(name_i, defName)
                xcol = scheilDefs{j,2};
                ycol = scheilDefs{j,3};
                c    = scheilDefs{j,4};

                if xcol <= size(A,2) && ycol <= size(A,2)
                    x = A(:, xcol);
                    y = A(:, ycol);

                    idx = ~isnan(x) & ~isnan(y);

                    if any(idx)
                        plot(x(idx), y(idx), '-', ...
                            'Color', c, ...
                            'LineWidth', 1.5, ...
                            'DisplayName', defName);

                        foundPath(i) = true;
                    end
                end
            end
        end
    end

    % =========================
    % 5. 坐标轴与图例格式
    % =========================
   xlabel('Mole fraction of solid', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 16, ...
    'FontWeight', 'bold');

ylabel('Temperature [^{\circ}C]', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 16, ...
    'FontWeight', 'bold');

    set(gca, ...
        'FontName', 'Times New Roman', ...
        'FontSize', 14, ...
        'LineWidth', 1, ...
        'Box', 'on');

    legend('Location', 'northeast', ...
        'FontName', 'SimSun', ...
        'FontSize', 11, ...
        'Box', 'on');

    grid off;

    % =========================
    % 6. 坐标范围
    % =========================
    if ~isempty(xRange)
        xlim(xRange);
    end

    if ~isempty(yRange)
        ylim(yRange);
    end

    hold off;

    % =========================
    % 7. 提示哪些路径没找到
    % =========================
    for i = 1:numel(pathNames)
        if ~foundPath(i)
            warning('未找到路径：%s', pathNames{i});
        end
    end

end