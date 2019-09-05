function plotfitdist(data,distname,par,dtype,varargin)
%PLOTFITDIST plots data and fitted distribution
%	PLOTFITDIST(DATA,NAME,PAR,DTYPE) plots the probability density 
%	function (specified with name and parameters), and the histogram 
%	of values in vector DATA scaled by total area using 'trapz'. 
%
%	NAME (char. string) is the name of the distribution
%	(for examples of accepted names see 'HELP PDF'). PAR is a 
%	vector parameters of the distribution, the size of PAR (1 to 3) 
%	must match the requirements of function PDF for the specified 
%	distribution. DTYPE is a character string, either 'cont' for continuous 
%	distributions or 'disc' for discrete ones. Continuous 
%	distributions are plotted as a line, while discrete ones are 
%	plotted as an histogram superimposed on the data histogram.
%
%	If NAME is 'binomial' an additional parameter is required:
%		PLOTFITDIST(...,'ntrials',NTRIALS) NTRIALS is a scalar
%		specifying the number of trials. If NTRIALS is a vector (i.e.
%		different number of trials for each value in DATA) the plot cannot be
%		made due to mismatch in the number of simulated data and number of
%		trials.
%
%	No input control is made, so the user is responsible for providing
%	correct parameters.


% Arguments
for j= 1:2:length(varargin)
	string= lower(varargin{j});
	switch string(1:3)
		case 'ntr'
			ntrials= varargin{j+1};
		case 'pre'
			if ~isempty(varargin{j+1}) distname= varargin{j+1}; end
		otherwise
			error('Unknown argument name');
	end
end


% Histogram for Discrete/Continuous
if strcmp(dtype,'cont')
	x = min(data):range(data)/100:max(data);
	[bincount,binpos] = hist(data,min(100,numel(data)/5));
else
	x = unique(data);
	[bincount,binpos] = hist(data,x);
end


% Calculate predicted
switch numel(par)
	case 1
		if strcmp('binomial',distname)
			y = pdf('bino',x,ntrials(1),par(1));
		else
			y = pdf(distname,x,par(1));
		end
	case 2
		y = pdf(distname,x,par(1),par(2));
	case 3
		y = pdf(distname,x,par(1),par(2),par(3));
end


% Plot
figure; hold on
% binwidth = binpos(2) - binpos(1);
% histarea = binwidth*sum(bincount);
bincount= bincount/trapz(binpos,bincount); % scaled frequencies
data= bar(binpos,bincount,'FaceColor',[.8 .8 .8],'EdgeColor',[1 1 1],'BarWidth',1);
if strcmp('cont',dtype)
	model= plot(x,y,'r','LineWidth',2);
else
	model= bar(x,y,'FaceColor',[1 0 0],'EdgeColor','none','BarWidth',0.5);
end
xlabel('Data'); ylabel('PDF')
legend([data,model],'Data',distname); legend('boxoff');









