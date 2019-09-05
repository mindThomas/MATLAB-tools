function F= fitmethis(data,varargin)
% FITMETHIS finds best-fitting distribution
%	F= fitmethis(X) fits all distributions available in MATLAB's function 
%	MLE to data in vector X and returns them ordered by some criterion: 
%	LogLikelihood or Akaike. Either 20 continuous or 5 discrete distributions 
%	are used based on user input or the type of data supplied (see below)
%
%	The function returns a structure array with fields: 
%	name: name of the distribution (see HELP MLE for a list)
%	par:  vector of parameter estimates (1, 2 or 3 values)
%	ci:   matrix of confidence limits, one column per parameter
%	LL:   Log-Likelihood of the data
%	aic:  Akaike Information Criterion
%
%	Optional arguments can be specified as name/value pairs. Argument 
%	names are case insensitive and partial matches are allowed.
%
%	Name			Value
%	'dtype'		character string that specifies if the data are continuous 
% 				('cont') or discrete ('disc'). If missing, the function 
% 				decides the data are discrete if all values of X are 
% 				natural numbers.
%
%	'ntrials'	Specifies number of trials for the binomial distribution. 
% 				NTRIALS must be either a scalar or a vector of the same 
% 				size as X. If missing the binomial is not fitted.
%
%	'figure'		Either 'on' (default), or 'off'. If 'on' a plot of the 
% 				data and the best fitting is produced (scaled to match the 
%				data). Requires aditional function 'plotfitdist'.
%
%	'alpha'		A value between 0 and 1 specifying a confidence level
%				for CI of 100*(1-alpha) (default is 0.05).
%
%	'criterion'	Criterion to use to order the fits. It can be: 'LL' for
%				Log-Likelihood (default), or 'AIC' for Akaike.
%	'output'		If set to 'off' supresses output to the command window.
%				Default 'on'
%	'pref'		Preferred distribution to plot
%
%	If X contains negative values, only the Normal distribution is fitted.
%	Also, if X contains values > 1 the Beta distribution is not fitted. If X
%	contains 0 some distributions are not fitted.
%
%	Requires Statistics Toolbox
%
%	Example:
%	x= gamrnd(5,0.5,1000,1);
%	F= fitmethis(x);
%
%   Name      Par1       Par2        Par3         LogL
%  gamma   4.947e+00   5.034e-01               -1.461e+03
%    gev  -1.126e-02   8.965e-01   1.982e+00   -1.463e+03
%	  ...


warning off


% Defaults & Storage
dtype= {};
ntrials= [];
fig= 'on';
alpha= 0.05;
criterion= 'LL';
output= 'on';
F= struct('name',{},'par',[],'ci',[],'LL',[],'aic',[]);
prefdist= [];


% Arguments
for j= 1:2:length(varargin)
	string= lower(varargin{j});
	switch string(1:min(3,length(string)))
		case 'dty'
			dtype= varargin{j+1};
		case 'ntr'
			ntrials= varargin{j+1};
		case 'fig'
			fig= varargin{j+1};
		case 'alp'
			alpha= varargin{j+1};
		case 'cri'
			criterion= varargin{j+1};
		case 'out'
			output= varargin{j+1};
		case 'pre'
			prefdist= varargin{j+1}
		otherwise
			error('Unknown argument name');
	end
end



% Distributions
Cdist= {'normal'; 'exponential'; 'gamma'; 'logistic'; ...
		  'tlocationscale';...
		  'uniform'; 'ev'; 'rayleigh'; 'gev'; 'beta'; ...
		  'nakagami'; 'rician'; 'inversegaussian'; 'birnbaumsaunders'; ...
		  'gp'; 'loglogistic'; 'lognormal'; 'weibull'};
mustbepos= 11;

Ddist= {'binomial'; 'nbin'; 'unid';'geometric';'poisson'};



% Try determine data type: Discrete or Continuous (avoid 0)
if isempty(dtype)
	if isempty(find(1- (data+1)./(fix(data)+1), 1)) 
		dtype= 'disc';
	else
		dtype= 'cont';
	end
end



% Fit all... 
switch dtype(1:4)

	% Continuous
	case 'cont'
	for j= 1:numel(Cdist)
		% If negative values, only fit normal
		if min(data) < 0
			[phat,pci]= mle(data,'distribution','normal','alpha',alpha);
			F(j).name= Cdist{j};
			F(j).par= phat;
			F(j).ci=  pci;
			pdfv= pdf('normal',data,F(j).par(1),F(j).par(2));
			F(j).LL=  sum(log(pdfv(pdfv>0 & ~isinf(pdfv))));
			F(j).aic= 2*2- 2*F(j).LL;
			break

		% Check: if values > 1 for Beta, do nothing
		elseif strcmp('beta',Cdist{j}) && max(data) > 1
			F(j).name= 'beta';
			F(j).LL= -Inf;
			F(j).aic= Inf;

		% Check: if values > 0 for some distr. (they are sorted), do nothing
		elseif  j >= mustbepos && min(data) == 0
			F(j).name= Cdist{j};
			F(j).LL= -Inf;
			F(j).aic= Inf;

		% Any other case do the fit ...
		else
			[phat,pci]= mle(data,'distribution',Cdist{j},'alpha',alpha);
			F(j).name= Cdist{j};
			F(j).par= phat;
			F(j).ci=  pci;
			if numel(F(j).par) == 1
				pdfv= pdf(F(j).name,data,F(j).par(1));
			elseif numel(F(j).par) == 2
				pdfv= pdf(F(j).name,data,F(j).par(1),F(j).par(2));
			else
				pdfv= pdf(F(j).name,data,F(j).par(1),F(j).par(2),F(j).par(3));
			end
			F(j).LL=  sum(log(pdfv(pdfv>0 & ~isinf(pdfv))));
			F(j).aic= 2*numel(F(j).par)- 2*F(j).LL;
		end

	end


	% Discrete
	case 'disc'
	for j= 1:numel(Ddist)

		% Binomial needs number of trials
		if strcmp('binomial',Ddist{j}) 
			F(j).name= 'binomial';
			if isempty(ntrials) || (numel(ntrials) > 1 && numel(data) ~= numel(ntrials))
				F(j).LL= -Inf;
				F(j).aic= Inf;
			else
				[phat,pci]= mle(data,'ntrials',ntrials,'distribution','binomial','alpha',alpha);
				F(j).par= phat;
				F(j).ci=  pci;
				pdfv= pdf('bino',data,ntrials,F(j).par(1));
				F(j).LL=  sum(log(pdfv(pdfv>0 & ~isinf(pdfv))));
			end
		else 
			[phat,pci]= mle(data,'distribution',Ddist{j},'alpha',alpha);
			F(j).name= Ddist{j};
			F(j).par= phat;
			F(j).ci=  pci;
			if numel(F(j).par) == 1
				pdfv= pdf(F(j).name,data,F(j).par(1));
			elseif numel(F(j).par) == 2
				pdfv= pdf(F(j).name,data,F(j).par(1),F(j).par(2));
			else
				pdfv= pdf(F(j).name,data,F(j).par(1),F(j).par(2),F(j).par(3));
			end
			F(j).LL=  sum(log(pdfv(pdfv>0 & ~isinf(pdfv))));
			F(j).aic= 2*numel(F(j).par)- 2*F(j).LL;
		end
	end

end


% Order by criterion
switch criterion
	case 'LL'
		index= sortrows([(1:size(F,2))',[F.LL]'],-2);
	case 'AIC'
		index= sortrows([(1:size(F,2))',[F.aic]'],2);
end
F= F(index(:,1));


% Nice screen output
if strcmp('on',output)
	fprintf('\n\t\t\t\tName\t\tPar1\t\tPar2\t\tPar3\t\tLogL\t\tAIC\n')
	for j= 1:size(F,2)
		switch numel(F(j).par)
			case 1
				fprintf('%20s \t%10.3e \t\t\t\t\t\t\t%10.3e \t%10.3e\n',F(j).name,F(j).par,F(j).LL,F(j).aic)
			case 2
				fprintf('%20s \t%10.3e \t%10.3e \t\t\t\t%10.3e \t%10.3e\n',F(j).name,F(j).par,F(j).LL,F(j).aic)
			case 3
				fprintf('%20s \t%10.3e \t%10.3e \t%10.3e \t%10.3e \t%10.3e\n',F(j).name,F(j).par,F(j).LL,F(j).aic)
		end
	end
end


% Choose the preferred distr. for plotting (if specified)
best= 1;
if ~isempty(prefdist)
	best= find(strcmp(cellstr(prefdist),cellstr({F.name})));
	if isempty(best)
		best= 1;
		fprintf('\n *** Warning: the preferred distr. is not in the list of fitted \n');
	else
		fprintf('\n *** Warning: Plotting preferred distr. [%s] NOT BEST FIT \n',F(best).name);
	end
end


% Plot data histogram & best fit
if strcmp('on',fig)
	if strcmp('binomial',F(best).name)
		plotfitdist(data,F(best).name,F(best).par,dtype,'ntrials',ntrials(best));
	else
		plotfitdist(data,F(best).name,F(best).par,dtype);
	end
end






