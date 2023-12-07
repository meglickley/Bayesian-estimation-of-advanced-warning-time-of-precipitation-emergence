function savefig(FigHandle, str, figure_width,figure_height)

export_ppi = 600; % % number of pixels per inch *in exported file*; only relevant if exporting as pixel-based figure (such as png, tiff)

% WE ALSO NEED
screen_ppi = 72; % inherent property of Matlab. CHANGING THIS WILL MAKE FIGURE LOOK INCONSISTENT BETWEEN SCREEN AND EXPORT!

% DERIVED PROPERTIES (cannot be changed; for info only)
screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

print(gcf, '-dpdf', str);