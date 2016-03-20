vote = load('voting-yi-shan.txt');
load Metric;
simMetric = Metric(vote == 1);
notSimMetric = Metric(vote == 2);

close all;
ax = 5.4 : 0.05 : max(Metric);
simHist = hist(simMetric, ax);
notSimHist = hist(notSimMetric, ax);
plot(ax, simHist, 'b-');
hold on 
plot(ax, notSimHist, 'r-');