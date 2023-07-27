function [] = betterBoxplot(plot_location,data,color,dotsize,marker,linewidth,outliers)

% BETTERBOXPLOT will plot a single box plot with scatter at a coordinate
% 
% the final input is a boolean indicating whether or not to show outliers.
%	The default is to show outliers. 

% edited: EHS20201202

% default inputs
if nargin<2
    fprintf('\nRequired input arguments:\n    1) plot location\n     2) data')
elseif nargin==2
    color = [0 0 0];
    linewidth = 3;
    dotsize = 2;
	marker ='o';
	outliers = true; 
elseif nargin==3
    dotsize = 2;
    linewidth = 1;
	marker = 'o';
outliers = true; 
elseif nargin==4
    linewidth = 1;
	marker = 'o';
	outliers = true; 
elseif nargin==5
    linewidth = 1;
	outliers = true; 
elseif nargin==6
	outliers=true;
end

% just turning data into a column vector
data = data(:);
if size(data,2)>size(data,1), data = data'; end

% plotting
ca = get(gca);
hold on

% plot scatter of the data\
if outliers
	boxplot(data,'positions',plot_location,'BoxStyle','filled','Colors',color)
else
	boxplot(data,'positions',plot_location,'BoxStyle','filled','Colors',color,'Symbol','w.')
end
scatter(randsBtwInterval(length(data),plot_location-0.1,plot_location+0.1),data,dotsize,color,'filled',marker)

hold off




% [20200304] the genesis: The plots from violinPlot are too cartoony, and
% appear to overestimate the extent of the discrtibution (yikes!), so I've
% decided to make my own, while also scatterplotting the data. 



function r = randsBtwInterval(n,a,b)

% RANDSBTWINTERVAL outputs n random numbers between a and b

r = (b-a).*rand(n,1) + a;
