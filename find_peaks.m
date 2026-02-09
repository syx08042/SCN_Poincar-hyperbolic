function [meanAmp, varAmp] = find_peaks(v, t, n_cells)

amplitude = nan(n_cells, 1);
for i = 1:n_cells
    x = v(:,i);
    [pks, locs_pks] = findpeaks(x, t, 'MinPeakProminence', 0.05);
    [trs, locs_trs] = findpeaks(-x, t, 'MinPeakProminence', 0.05);
    trs = -trs;
    if ~isempty(pks) && ~isempty(trs)
        amplitude(i) = (mean(pks) - mean(trs)) / 2;
    end
end
meanAmp = mean(amplitude, 'omitnan');
varAmp = var(amplitude, 'omitnan');
end
