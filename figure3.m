
clear all
close all
AA=readtable('16k-Event_Age-All-Speleothems.xlsx')
LW=1.5
FS=12


try
        if exist('cols','var') && size(cols,2) >= 4
        cols(:,2) = cols(:,2) / 2; 
        cols(:,3) = cols(:,3) / 2; 
        AA{:,3} = cols(:,2);
        AA{:,4} = cols(:,3);
    else
            if width(AA) >= 4
            AA{:,3} = AA{:,3} ./ 2;
            AA{:,4} = AA{:,4} ./ 2;
        end
    end
catch
    end



cols = AA{:,2:5}; 
event_age = cols(:,4);
lower_unc = cols(:,2)/1000;
upper_unc = cols(:,3)/1000;


if ~isnumeric(event_age) || ~isnumeric(lower_unc) || ~isnumeric(upper_unc)
    event_age = str2double(string(AA{:,2}));
    lower_unc = str2double(string(AA{:,3}));
    upper_unc = str2double(string(AA{:,4}));
end


valid = ~(isnan(event_age) | isnan(lower_unc) | isnan(upper_unc));
event_age = event_age(valid);
lower_unc = lower_unc(valid);
upper_unc = upper_unc(valid);


figure;
errorbar(1:numel(event_age), event_age, lower_unc, upper_unc, 'o');
xlabel('Sample index');
ylabel('Event age');
title('Event age with lower and upper uncertainties');
grid on;


legend('Event Age with Uncertainties');
set(gca, 'XTick', 1:numel(event_age), 'XTickLabel', 1:numel(event_age));
xlim([0 numel(event_age) + 1]);
saveas(gcf, 'event_age_plot.png');


nSamples = numel(event_age);
nGrid = 201;



sn_x = cell(nSamples,1);
sn_pdf = cell(nSamples,1);
for i = 1:nSamples
    mu = event_age(i);
    sL = lower_unc(i);
    sR = upper_unc(i);
    sL = max(sL, eps);
    sR = max(sR, eps);
    mu = max(mu, eps);
    r = sR / sL;
    alpha = sign(log(r)) * min(10, abs(3*log(r))); 
    omega = 0.5*(sL + sR);
    omega = max(omega, eps);
    delta = alpha / sqrt(1 + alpha^2);
    xi = mu - omega * delta * sqrt(2/pi);
    left = xi - 8*omega;
    right = xi + 8*omega;
    x = linspace(left, right, nGrid);
    t = (x - xi) ./ omega;
    phi = exp(-0.5 * t.^2) / sqrt(2*pi);
    Phi = 0.5 * (1 + erf( (alpha .* t) / sqrt(2) ));
    pdf = (2./omega) .* phi .* Phi;
    pdf = max(pdf, 0);
    dx = x(2)-x(1);
    s = sum(pdf)*dx;
    if s > 0
        pdf = pdf / s;
    end
    sn_x{i} = x;
    sn_pdf{i} = pdf;
end



figure;
set(gcf, 'Position',  [0, 0, 500, 800])
subplot(2,1,1)
hold on;
colors = lines(nSamples);
offsets = (0:nSamples-1) * max(cellfun(@max, sn_pdf))*0.3;

for i = nSamples:-1:1
    plot(sn_x{i}, sn_pdf{i} + offsets(i), 'Color', colors(mod(i-1,size(colors,1))+1,:),'LineWidth',LW);
end
xlabel('Event age');
ylabel('Skew-normal PDF (offset)');
hold off;



hold on;
for i = 1:nSamples
    xmean = event_age(i);
    xl = xmean - lower_unc(i);
    xr = xmean + upper_unc(i);
    y = offsets(i) + max(sn_pdf{i}); 
        plot([xl, xr], [y, y], 'k-', 'LineWidth', LW,'Color',colors(mod(i-1,size(colors,1))+1,:));
    plot(xmean, y, 'ko', 'MarkerFaceColor','k', 'MarkerSize',6);
end
hold off;
grid
legend(flip({'Klang','Furong','Haozhu','Qingtian','Hulu','Zhangjia','La Vallina'}))


set(gca,'FontSize',FS,'FontName','Arial','Box','off')
xlim([15.9 16.4])
xlabel('Event age (ka) estimated from RAMPFIT')


subplot(2,1,2)
all_left = min(cellfun(@(x) x(1), sn_x));
all_right = max(cellfun(@(x) x(end), sn_x));
ngrid_pool = 1001;
xpool = linspace(all_left, all_right, ngrid_pool);
pdfs_interp = zeros(nSamples, ngrid_pool);
for i = 1:nSamples
    pdfs_interp(i,:) = interp1(sn_x{i}, sn_pdf{i}, xpool, 'linear', 0);
end
w = ones(nSamples,1) / nSamples;
pooled_pdf = w' * pdfs_interp;
dxp = xpool(2)-xpool(1);
pooled_pdf = pooled_pdf / (sum(pooled_pdf)*dxp + eps);



figure(gcf); hold on;
max_offset = offsets(end) + max(cellfun(@max, sn_pdf));
max_offset=0
plot(xpool, pooled_pdf + max_offset + 0.05*max_offset, 'k-', 'LineWidth', LW);
hold off;


log_pdfs = zeros(nSamples, ngrid_pool);
for i = 1:nSamples
    p = pdfs_interp(i,:);
    p(p<=0) = realmin;
    log_pdfs(i,:) = log(p);
end
log_geom = (w' * log_pdfs); 
geom_unnorm = exp(log_geom);
geom_pdf = geom_unnorm / (sum(geom_unnorm)*dxp + eps);


figure(gcf); hold on;
plot(xpool, geom_pdf + max_offset + 0.05*max_offset, 'r--', 'LineWidth', LW);



cdfs = cumsum(pdfs_interp, 2) * dxp;
pgrid = linspace(0,1,ngrid_pool);
quantiles = zeros(nSamples, ngrid_pool);
for i = 1:nSamples
    cdf_i = cdfs(i,:);
    [cdf_u, ia] = unique(cdf_i);
    x_u = xpool(ia);
    quantiles(i,:) = interp1(cdf_u, x_u, pgrid, 'linear', 'extrap');
end
avg_quantile = mean(quantiles, 1);

[~, ia] = unique(avg_quantile);
xp_v = pgrid(ia);
xq_v = avg_quantile(ia);
inv_q = interp1(xq_v, xp_v, xpool, 'linear', 0); 
vincent_pdf = gradient(inv_q, xpool);
vincent_pdf = max(vincent_pdf, 0);
vincent_pdf = vincent_pdf / (sum(vincent_pdf)*dxp + eps);


figure(gcf); hold on;
plot(xpool, vincent_pdf + max_offset + 0.05*max_offset, 'b-.', 'LineWidth', LW);
legend('Linear pooling','Geometric pooling', 'Quantile averaging')
grid
hold off;

xlim([15.9 16.4])
set(gca,'FontSize',FS,'FontName','Arial','Box','off')
xlabel('Event age (ka) estimated from RAMPFIT')
ylabel('Consolidated PDF from 7 speleothems')
fname = 'Figure3_pooled_pdfs.png';
print(gcf, fname, '-dpng', '-r300');


mean_linear = sum(xpool .* pooled_pdf) * dxp;
mean_geometric = sum(xpool .* geom_pdf) * dxp;
mean_vincent = sum(xpool .* vincent_pdf) * dxp;
var_linear = sum(((xpool - mean_linear).^2) .* pooled_pdf) * dxp;
var_geometric = sum(((xpool - mean_geometric).^2) .* geom_pdf) * dxp;
var_vincent = sum(((xpool - mean_vincent).^2) .* vincent_pdf) * dxp;
std_linear=var_linear.^0.5;
std_geometric=var_geometric.^0.5;
std_vincent=var_vincent.^0.5;
fprintf('Linear pooled: mean = %.4f, std = %.6f\n', mean_linear, std_linear);
fprintf('Geometric pooled: mean = %.4f, std = %.6f\n', mean_geometric, std_geometric);
fprintf('Vincentized: mean = %.4f, std= %.6f\n', mean_vincent, std_vincent);