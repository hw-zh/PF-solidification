function plot_phase_fraction_from_excel(phaseList, xRange, xStep)
% plot_phase_fraction_from_excel
% -------------------------------------------------------------------------
% 功能：
%   自动读取当前工作文件夹中的第一个 Excel 文件，
%   按“每两列一组（温度-相分数）”的格式绘制相分数曲线。
%
% 输入：
%   phaseList : 想画的相名，例如：
%               {'LIQUID','BCC_A2','FCC_A1'}
%
%   xRange    : 横坐标范围，例如：
%               [1675 1715]
%               不输入则自动
%
%   xStep     : 横坐标刻度间隔，例如：
%               5
%               不输入则使用默认刻度
%
% 用法：
%   plot_phase_fraction_from_excel({'LIQUID','BCC_A2','FCC_A1'})
%   plot_phase_fraction_from_excel({'LIQUID','BCC_A2','FCC_A1'}, [1675 1715])
%   plot_phase_fraction_from_excel({'LIQUID','BCC_A2','FCC_A1'}, [1675 1715], 5)
% -------------------------------------------------------------------------

    % =========================
    % 1. 输入检查
    % =========================
    if nargin < 1 || isempty(phaseList)
        error('请输入相名列表，例如：plot_phase_fraction_from_excel({''LIQUID'',''BCC_A2'',''FCC_A1''})');
    end

    if nargin < 2
        xRange = [];
    end

    if nargin < 3
        xStep = [];
    end

    if isstring(phaseList)
        phaseList = cellstr(phaseList);
    end

    if ischar(phaseList)
        phaseList = {phaseList};
    end

    if ~iscell(phaseList)
        error('phaseList 必须是 cell、string 或字符。');
    end

    phaseList = cellfun(@char, phaseList, 'UniformOutput', false);

    % =========================
    % 2. 自动读取当前文件夹里的第一个 Excel
    % =========================
    files = [dir('*.xlsx'); dir('*.xls')];

    if isempty(files)
        error('当前文件夹中没有找到 Excel 文件。');
    end

    excelFile = files(1).name;
    fprintf('当前读取的 Excel 文件：%s\n', excelFile);

    T = readtable(excelFile);
    ncol = width(T);

    % =========================
    % 3. 配色
    % =========================
    colors = [ ...
        0.00 0.45 0.74;   % 蓝
        1.00 0.00 0.00;   % 红
        0.00 0.00 0.00;   % 黑
        0.35 0.00 0.45;   % 深紫
        0.00 0.50 0.00;   % 深绿
        0.00 0.45 0.45;   % 深青
        0.63 0.08 0.18];  % 深红

    % =========================
    % 4. 绘图
    % =========================
    figure('Color', 'w');
    hold on;

    nplot = 0;
    foundPhase = false(size(phaseList));

    % 每两列一组：奇数列 = 温度，偶数列 = 相分数
    for k = 1:2:ncol-1
        x = T{:, k};
        y = T{:, k+1};
        phaseName = T.Properties.VariableNames{k+1};

        idxPhase = find(strcmp(phaseList, phaseName), 1);

        if ~isempty(idxPhase)
            idx = ~isnan(x) & ~isnan(y);

            if any(idx)
                nplot = nplot + 1;
                c = colors(mod(nplot-1, size(colors,1)) + 1, :);

                plot(x(idx), y(idx), ...
                    'LineWidth', 1.5, ...
                    'Color', c, ...
                    'DisplayName', phaseName);

                foundPhase(idxPhase) = true;
            end
        end
    end

    % =========================
    % 5. 坐标轴格式
    % =========================
    xlabel('Temperature [K]', ...
        'FontName', 'Times New Roman', ...
        'FontSize', 16, ...
        'FontWeight', 'bold');

    ylabel('Phase fraction', ...
        'FontName', 'Times New Roman', ...
        'FontSize', 16, ...
        'FontWeight', 'bold');

    set(gca, ...
        'FontName', 'Times New Roman', ...
        'FontSize', 14, ...
        'LineWidth', 1, ...
        'Box', 'on');

    legend('Location', 'northeast', ...
        'FontName', 'Times New Roman', ...
        'FontSize', 11, ...
        'Box', 'on');

    grid off;

    % =========================
    % 6. 横坐标范围与刻度
    % =========================
    if ~isempty(xRange)
        xlim(xRange);
    end

    if ~isempty(xStep)
        if isempty(xRange)
            xl = xlim;
            xticks(xl(1):xStep:xl(2));
        else
            xticks(xRange(1):xStep:xRange(2));
        end
    end

    hold off;

    % =========================
    % 7. 提示哪些相没找到
    % =========================
    for i = 1:numel(phaseList)
        if ~foundPhase(i)
            warning('未在 Excel 中找到相：%s', phaseList{i});
        end
    end

end