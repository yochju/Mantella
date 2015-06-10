function newstate = cmaes_iteration(state)

% MODICATIONS:
% 1. Exclude stopflag={'fitness','stagnation',...}, because F(x) of surrogate model changes over time.
% lin 950: if ~isempty(stopflag) || time.saving < 0.05 * time.nonoutput || countiter == 100F
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD STATE
[fitfun,xstart,insigma,inopts,varargin,cmaVersion,definput,defopts,flg_future_setting,nargin,input,opts,counteval,countevalNaN,irun ...
    flgresume, xmean,N,numberofvariables,lambda0,popsize,lambda,lambda_last,stopFitness,stopMaxFunEvals,stopMaxIter,stopFunEvals,stopIter, ... 
    stopTolFun,stopTolHistFun,stopOnStagnation,stopOnWarnings,flgreadsignals,flgWarnOnEqualFunctionValues,flgEvalParallel,stopOnEqualFunctionValues, ... 
    arrEqualFunvals, flgDiagonalOnly, flgActiveCMA,noiseHandling,noiseMinMaxEvals,noiseAlphaEvals,noiseCallback,flgdisplay,flgplotting,verbosemodulo, ... 
    flgscience,flgsaving,strsaving,flgsavingfinal, savemodulo,savetime,time,maxdx,mindx,lbounds,ubounds,stopTolX,stopTolUpX,sigma,pc,diagD,diagC,B,BD, ... 
    C,fitness,bnd,out,startseed,chiN,countiter,outiter,filenameprefix,filenames, lambda_hist,mu,weights,mueff,cc,cs,ccov1,ccovmu,ccov1_sep, ccovmu_sep, ... 
    damps,noiseReevals,noiseAlpha,noiseEpsilon,noiseTheta,noisecum,noiseCutOff,arx, arxvalid,tries,noiseS,noiseSS,noiseN,xold,zmean,fmean,ps,neg, ... 
    stopflag,noiseX,iterplotted,arfitness,Xnew_sorted,invsqrtC,stop] = cmaes_loadState(state);
 stopflag = {}; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD STATE

% set internal parameters
  if countiter == 0 || lambda ~= lambda_last
    if countiter > 0 && floor(log10(lambda)) ~= floor(log10(lambda_last)) ...
          && flgdisplay
      disp(['  lambda = ' num2str(lambda)]);
      lambda_hist(:,end+1) = [countiter+1; lambda];
    else
      lambda_hist = [countiter+1; lambda]; 
    end
    lambda_last = lambda;
    % Strategy internal parameter setting: Selection  
    mu = myeval(opts.ParentNumber); % number of parents/points for recombination
    if strncmp(lower(opts.RecombinationWeights), 'equal', 3)
      weights = ones(mu,1); % (mu_I,lambda)-CMA-ES
    elseif strncmp(lower(opts.RecombinationWeights), 'linear', 3)
      weights = mu+0.5-(1:mu)'; 
    elseif strncmp(lower(opts.RecombinationWeights), 'superlinear', 3)
      weights = log(mu+0.5)-log(1:mu)'; % muXone array for weighted recombination
                                        % qqq mu can be non-integer and
                                        % should become ceil(mu-0.5) (minor correction)
    else
      error(['Recombination weights to be "' opts.RecombinationWeights ...
             '" is not implemented']);
    end
    mueff=sum(weights)^2/sum(weights.^2); % variance-effective size of mu
    weights = weights/sum(weights);     % normalize recombination weights array
    if mueff == lambda
      error(['Combination of values for PopSize, ParentNumber and ' ...
             ' and RecombinationWeights is not reasonable']);
    end
    
    % Strategy internal parameter setting: Adaptation
    cc = myeval(opts.CMA.ccum); % time constant for cumulation for covariance matrix
    cs = myeval(opts.CMA.cs); 

    % old way TODO: remove this at some point
    % mucov = mueff;   % size of mu used for calculating learning rate ccov
    % ccov = (1/mucov) * 2/(N+1.41)^2 ... % learning rate for covariance matrix
    %        + (1-1/mucov) * min(1,(2*mucov-1)/((N+2)^2+mucov)); 

    % new way
    if myevalbool(opts.CMA.on) 
      ccov1 = myeval(opts.CMA.ccov1); 
      ccovmu = min(1-ccov1, myeval(opts.CMA.ccovmu));
    else
      ccov1 = 0;
      ccovmu = 0;
    end
    
    % flgDiagonalOnly = -lambda*4*1/ccov; % for ccov==1 it is not needed
    % 0 : C will never be diagonal anymore
    % 1 : C will always be diagonal
    % >1: C is diagonal for first iterations, set to 0 afterwards
    if flgDiagonalOnly < 1
      flgDiagonalOnly = 0; 
    end
    if flgDiagonalOnly
      ccov1_sep = min(1, ccov1 * (N+1.5) / 3); 
      ccovmu_sep = min(1-ccov1, ccovmu * (N+1.5) / 3); 
    elseif N > 98 && flgdisplay && countiter == 0
      disp('consider option DiagonalOnly for high-dimensional problems');
    end

    % ||ps|| is close to sqrt(mueff/N) for mueff large on linear fitness
    %damps = ... % damping for step size control, usually close to one 
    %    (1 + 2*max(0,sqrt((mueff-1)/(N+1))-1)) ... % limit sigma increase
    %    * max(0.3, ... % reduce damps, if max. iteration number is small
    %          1 - N/min(stopMaxIter,stopMaxFunEvals/lambda)) + cs; 
    damps = myeval(opts.CMA.damps); 

    %qqq hacking of a different parameter setting, e.g. for ccov or damps,
    %  can be done here, but is not necessary anymore, see opts.CMA. 
    % ccov1 = 0.0*ccov1; disp(['CAVE: ccov1=' num2str(ccov1)]);
    % ccovmu = 0.0*ccovmu; disp(['CAVE: ccovmu=' num2str(ccovmu)]);
    % damps = inf*damps; disp(['CAVE: damps=' num2str(damps)]);
    % cc = 1; disp(['CAVE: cc=' num2str(cc)]);

  end

  countiter = countiter + 1; 

  % Generate and evaluate lambda offspring
 
  fitness.raw = repmat(NaN, 1, lambda);

  % parallel evaluation
  if flgEvalParallel
      arz = randn(N,lambda);

      if ~flgDiagonalOnly
        arx = repmat(xmean, 1, lambda) + sigma * (BD * arz); % Eq. (1)
      else
        arx = repmat(xmean, 1, lambda) + repmat(sigma * diagD, 1, lambda) .* arz; 
      end

      % You may handle constraints here. You may either resample
      % arz(:,k) and/or multiply it with a factor between -1 and 1
      % (the latter will decrease the overall step size) and
      % recalculate arx accordingly. Do not change arx or arz in any
      % other way.
 
      if ~bnd.isactive
        arxvalid = arx;
      else
        arxvalid = xintobounds(arx, lbounds, ubounds);
      end
      % You may handle constraints here.  You may copy and alter
      % (columns of) arxvalid(:,k) only for the evaluation of the
      % fitness function. arx and arxvalid should not be changed.
      fitness.raw = feval(fitfun, arxvalid, varargin{:}); %mantella: no idea what varargin{:} is supposed to do... probably nothing
      countevalNaN = countevalNaN + sum(isnan(fitness.raw));
      counteval = counteval + sum(~isnan(fitness.raw)); 
  end

  % non-parallel evaluation and remaining NaN-values
  % set also the reevaluated solution to NaN
  fitness.raw(lambda + find(isnan(fitness.raw(1:0)))) = NaN;  %mantella: 0 was noiseReevals before.
  for k=find(isnan(fitness.raw)), 
    % fitness.raw(k) = NaN; 
    tries = 0;
    % Resample, until fitness is not NaN
    while isnan(fitness.raw(k))
      if k <= lambda  % regular samples (not the re-evaluation-samples)
        arz(:,k) = randn(N,1); % (re)sample

        if flgDiagonalOnly  
          arx(:,k) = xmean + sigma * diagD .* arz(:,k);              % Eq. (1)
        else
          arx(:,k) = xmean + sigma * (BD * arz(:,k));                % Eq. (1)
        end
      else % re-evaluation solution with index > lambda
        if flgDiagonalOnly  
          arx(:,k) = arx(:,k-lambda) + (noiseEpsilon * sigma) * diagD .* randn(N,1);
        else
          arx(:,k) = arx(:,k-lambda) + (noiseEpsilon * sigma) * (BD * randn(N,1));
        end
      end
      
      % You may handle constraints here. You may either resample
      % arz(:,k) and/or multiply it with a factor between -1 and 1
      % (the latter will decrease the overall step size) and
      % recalculate arx accordingly. Do not change arx or arz in any
      % other way.
 
      if ~bnd.isactive
        arxvalid(:,k) = arx(:,k);
      else
        arxvalid(:,k) = xintobounds(arx(:,k), lbounds, ubounds);
      end
      % You may handle constraints here.  You may copy and alter
      % (columns of) arxvalid(:,k) only for the evaluation of the
      % fitness function. arx should not be changed.
      fitness.raw(k) = feval(fitfun, arxvalid(:,k), varargin{:}); 
      tries = tries + 1;
      if isnan(fitness.raw(k))
	countevalNaN = countevalNaN + 1;
      end
      if mod(tries, 100) == 0
	warning([num2str(tries) ...
                 ' NaN objective function values at evaluation ' ...
                 num2str(counteval)]);
      end
    end
    counteval = counteval + 1; % retries due to NaN are not counted
  end

  fitness.sel = fitness.raw; 

  % ----- handle boundaries -----
  if 1 < 3 && bnd.isactive
    % Get delta fitness values
    val = myprctile(fitness.raw, [25 75]);
    % more precise would be exp(mean(log(diagC)))
    val = (val(2) - val(1)) / N / mean(diagC) / sigma^2;
    %val = (myprctile(fitness.raw, 75) - myprctile(fitness.raw, 25)) ...
    %    / N / mean(diagC) / sigma^2;
    % Catch non-sensible values 
    if ~isfinite(val)
      warning('Non-finite fitness range');
      val = max(bnd.dfithist);  
    elseif val == 0 % happens if all points are out of bounds
      val = min(bnd.dfithist(bnd.dfithist>0));  % seems not to make sense, given all solutions are out of bounds
    elseif bnd.validfitval == 0 % flag that first sensible val was found
      bnd.dfithist = [];
      bnd.validfitval = 1;
    end

    % Store delta fitness values
    if length(bnd.dfithist) < 20+(3*N)/lambda
      bnd.dfithist = [bnd.dfithist val];
    else
      bnd.dfithist = [bnd.dfithist(2:end) val];
    end

    [tx ti]  = xintobounds(xmean, lbounds, ubounds);

    % Set initial weights
    if bnd.iniphase 
      if any(ti) 
        bnd.weights(find(bnd.isbounded)) = 2.0002 * median(bnd.dfithist);
	if bnd.flgscale == 0 % scale only initial weights then
	  dd = diagC; 
	  idx = find(bnd.isbounded); 
	  dd = dd(idx) / mean(dd); %  remove mean scaling
	  bnd.weights(idx) = bnd.weights(idx) ./ dd; 
	end
	if bnd.validfitval && countiter > 2
          bnd.iniphase = 0;
	end
      end
    end

    % Increase weights
    if  1 < 3 && any(ti) % any coordinate of xmean out of bounds
      % judge distance of xmean to boundary
      tx = xmean - tx;
      idx = (ti ~= 0 & abs(tx) > 3*max(1,sqrt(N)/mueff) ... 
	     * sigma*sqrt(diagC)) ;
      % only increase if xmean is moving away
      idx = idx & (sign(tx) == sign(xmean - xold));
      if ~isempty(idx) % increase
        % the factor became 1.2 instead of 1.1, because
        % changed from max to min in version 3.52
        bnd.weights(idx) = 1.2^(min(1, mueff/10/N)) * bnd.weights(idx); 
      end
    end

    % Calculate scaling biased to unity, product is one
    if bnd.flgscale ~= 0 
      bnd.scale = exp(0.9*(log(diagC)-mean(log(diagC)))); 
    end

    % Assigned penalized fitness
    bnd.arpenalty = (bnd.weights ./ bnd.scale)' * (arxvalid - arx).^2; 

    fitness.sel = fitness.raw + bnd.arpenalty;

  end % handle boundaries
  % ----- end handle boundaries -----
  
  % Sort by fitness 
  [fitness.raw, fitness.idx] = sort(fitness.raw); 
  [fitness.sel, fitness.idxsel] = sort(fitness.sel);  % minimization
  fitness.hist(2:end) = fitness.hist(1:end-1);    % record short history of
  fitness.hist(1) = fitness.raw(1);               % best fitness values
  if length(fitness.histbest) < 120+ceil(30*N/lambda) || ...
       (mod(countiter, 5) == 0  && length(fitness.histbest) < 2e4)  % 20 percent of 1e5 gen.
    fitness.histbest = [fitness.raw(1) fitness.histbest];          % best fitness values
    fitness.histmedian = [median(fitness.raw) fitness.histmedian]; % median fitness values
  else
    fitness.histbest(2:end) = fitness.histbest(1:end-1); 
    fitness.histmedian(2:end) = fitness.histmedian(1:end-1); 
    fitness.histbest(1) = fitness.raw(1);           % best fitness values
    fitness.histmedian(1) = median(fitness.raw);    % median fitness values
  end
  fitness.histsel(2:end) = fitness.histsel(1:end-1); % record short history of
  fitness.histsel(1) = fitness.sel(1);               % best sel fitness values

  % Calculate new xmean, this is selection and recombination 
  xold = xmean; % for speed up of Eq. (2) and (3)
  xmean = arx(:,fitness.idxsel(1:mu))*weights; 
  zmean = arz(:,fitness.idxsel(1:mu))*weights;%==D^-1*B'*(xmean-xold)/sigma
  if mu == 1
    fmean = fitness.sel(1);
  else
    fmean = NaN; % [] does not work in the latter assignment
    % fmean = feval(fitfun, xintobounds(xmean, lbounds, ubounds), varargin{:});
    % counteval = counteval + 1;
  end
  
  % Cumulation: update evolution paths
  ps = (1-cs)*ps + sqrt(cs*(2-cs)*mueff) * (B*zmean);          % Eq. (4)
  hsig = norm(ps)/sqrt(1-(1-cs)^(2*countiter))/chiN < 1.4 + 2/(N+1);
  if flg_future_setting
    hsig = sum(ps.^2) / (1-(1-cs)^(2*countiter)) / N < 2 + 4/(N+1); % just simplified
  end
%  hsig = norm(ps)/sqrt(1-(1-cs)^(2*countiter))/chiN < 1.4 + 2/(N+1);
%  hsig = norm(ps)/sqrt(1-(1-cs)^(2*countiter))/chiN < 1.5 + 1/(N-0.5);
%  hsig = norm(ps) < 1.5 * sqrt(N);
%  hsig = 1;

  pc = (1-cc)*pc ...
        + hsig*(sqrt(cc*(2-cc)*mueff)/sigma) * (xmean-xold);     % Eq. (2)
  if hsig == 0
    % disp([num2str(countiter) ' ' num2str(counteval) ' pc update stalled']);
  end

  % Adapt covariance matrix
  neg.ccov = 0;  % TODO: move parameter setting upwards at some point
  if ccov1 + ccovmu > 0                                                    % Eq. (3)
    if flgDiagonalOnly % internal linear(?) complexity
      diagC = (1-ccov1_sep-ccovmu_sep+(1-hsig)*ccov1_sep*cc*(2-cc)) * diagC ... % regard old matrix 
          + ccov1_sep * pc.^2 ...               % plus rank one update
          + ccovmu_sep ...                      % plus rank mu update
            * (diagC .* (arz(:,fitness.idxsel(1:mu)).^2 * weights));
%             * (repmat(diagC,1,mu) .* arz(:,fitness.idxsel(1:mu)).^2 * weights);
      diagD = sqrt(diagC); % replaces eig(C)
    else
      arpos = (arx(:,fitness.idxsel(1:mu))-repmat(xold,1,mu)) / sigma;
      % "active" CMA update: negative update, in case controlling pos. definiteness 
      if flgActiveCMA > 0
        % set parameters
        neg.mu = mu;  
        neg.mueff = mueff;
        if flgActiveCMA > 10  % flat weights with mu=lambda/2
          neg.mu = floor(lambda/2);  
          neg.mueff = neg.mu;
        end
        % neg.mu = ceil(min([N, lambda/4, mueff]));  neg.mueff = mu; % i.e. neg.mu <= N 
        % Parameter study: in 3-D lambda=50,100, 10-D lambda=200,400, 30-D lambda=1000,2000 a 
        % three times larger neg.ccov does not work. 
        %   increasing all ccov rates three times does work (probably because of the factor (1-ccovmu))
        %   in 30-D to looks fine

        neg.ccov = (1 - ccovmu) * 0.25 * neg.mueff / ((N+2)^1.5 + 2*neg.mueff);
        neg.minresidualvariance = 0.66;  % keep at least 0.66 in all directions, small popsize are most critical
        neg.alphaold = 0.5;  % where to make up for the variance loss, 0.5 means no idea what to do
                             % 1 is slightly more robust and gives a better "guaranty" for pos. def., 
                             % but does it make sense from the learning perspective for large ccovmu? 

        neg.ccovfinal = neg.ccov;

        % prepare vectors, compute negative updating matrix Cneg and checking matrix Ccheck
        arzneg = arz(:,fitness.idxsel(lambda:-1:lambda - neg.mu + 1));
        % i-th longest becomes i-th shortest
        % TODO: this is not in compliance with the paper Hansen&Ros2010, 
        %       where simply arnorms = arnorms(end:-1:1) ? 
        [arnorms idxnorms] = sort(sqrt(sum(arzneg.^2, 1))); 
        [ignore idxnorms] = sort(idxnorms);  % inverse index 
        arnormfacs = arnorms(end:-1:1) ./ arnorms; 
        % arnormfacs = arnorms(randperm(neg.mu)) ./ arnorms;
        arnorms = arnorms(end:-1:1); % for the record
        if flgActiveCMA < 20
          arzneg = arzneg .* repmat(arnormfacs(idxnorms), N, 1);  % E x*x' is N
          % arzneg = sqrt(N) * arzneg ./ repmat(sqrt(sum(arzneg.^2, 1)), N, 1);  % E x*x' is N
        end
        if flgActiveCMA < 10 && neg.mu == mu  % weighted sum
          if mod(flgActiveCMA, 10) == 1 % TODO: prevent this with a less tight but more efficient check (see below) 
            Ccheck = arzneg * diag(weights) * arzneg';  % in order to check the largest EV
          end
          artmp = BD * arzneg;
          Cneg = artmp * diag(weights) * artmp';
        else  % simple sum
          if mod(flgActiveCMA, 10) == 1
            Ccheck = (1/neg.mu) * arzneg*arzneg';  % in order to check largest EV
          end
          artmp = BD * arzneg;
          Cneg = (1/neg.mu) * artmp*artmp';

        end

        % check pos.def. and set learning rate neg.ccov accordingly, 
        % this check makes the original choice of neg.ccov extremly failsafe 
        % still assuming C == BD*BD', which is only approxim. correct 
        if mod(flgActiveCMA, 10) == 1 && 1 - neg.ccov * arnorms(idxnorms).^2 * weights < neg.minresidualvariance
          % TODO: the simple and cheap way would be to set
          %    fac = 1 - ccovmu - ccov1 OR 1 - mueff/lambda and
          %    neg.ccov = fac*(1 - neg.minresidualvariance) / (arnorms(idxnorms).^2 * weights)
          % this is the more sophisticated way: 
          % maxeigenval = eigs(arzneg * arzneg', 1, 'lm', eigsopts);  % not faster
          maxeigenval = max(eig(Ccheck));  % norm is much slower, because (norm()==max(svd())
          %disp([countiter log10([neg.ccov, maxeigenval, arnorms(idxnorms).^2 * weights, max(arnorms)^2]), ...
           %          neg.ccov * arnorms(idxnorms).^2 * weights])
          % pause
          % remove less than ??34*(1-cmu)%?? of variance in any direction
          %     1-ccovmu is the variance left from the old C
          neg.ccovfinal = min(neg.ccov, (1-ccovmu)*(1-neg.minresidualvariance)/maxeigenval); 
                                        % -ccov1 removed to avoid error message??
          if neg.ccovfinal < neg.ccov
            disp(['active CMA at iteration ' num2str(countiter) ...
                 ': max EV ==', num2str([maxeigenval, neg.ccov, neg.ccovfinal])]);
          end
        end
        % xmean = xold;  % the distribution does not degenerate!? 
        % update C
        C = (1-ccov1-ccovmu+neg.alphaold*neg.ccovfinal+(1-hsig)*ccov1*cc*(2-cc)) * C ... % regard old matrix 
            + ccov1 * pc*pc' ...     % plus rank one update
            + (ccovmu + (1-neg.alphaold)*neg.ccovfinal) ...  % plus rank mu update
              * arpos * (repmat(weights,1,N) .* arpos') ...
              - neg.ccovfinal * Cneg;                        % minus rank mu update
      else  % no active (negative) update
        C = (1-ccov1-ccovmu+(1-hsig)*ccov1*cc*(2-cc)) * C ... % regard old matrix 
            + ccov1 * pc*pc' ...     % plus rank one update
            + ccovmu ...             % plus rank mu update
              * arpos * (repmat(weights,1,N) .* arpos');
        % is now O(mu*N^2 + mu*N), was O(mu*N^2 + mu^2*N) when using diag(weights)
        %   for mu=30*N it is now 10 times faster, overall 3 times faster
      end
      diagC = diag(C);
    end
  end
  
  % the following is de-preciated and will be removed in future
  % better setting for cc makes this hack obsolete
  if 11 < 2 && ~flgscience  
    % remove momentum in ps, if ps is large and fitness is getting worse.
    % this should rarely happen. 
    % this might very well be counterproductive in dynamic environments
    if sum(ps.^2)/N > 1.5 + 10*(2/N)^.5 && ...
        fitness.histsel(1) > max(fitness.histsel(2:3))
      ps = ps * sqrt(N*(1+max(0,log(sum(ps.^2)/N))) / sum(ps.^2));
      if flgdisplay
        disp(['Momentum in ps removed at [niter neval]=' ...
              num2str([countiter counteval]) ']']);
      end
    end
  end

  % Adapt sigma
  if flg_future_setting  % according to a suggestion from Dirk Arnold (2000)
    % exp(1) is still not reasonably small enough
    sigma = sigma * exp(min(1, (sum(ps.^2)/N - 1)/2 * cs/damps));            % Eq. (5)
  else
    % exp(1) is still not reasonably small enough
    sigma = sigma * exp(min(1, (sqrt(sum(ps.^2))/chiN - 1) * cs/damps));             % Eq. (5)
  end
  % disp([countiter norm(ps)/chiN]);
  
  if 11 < 3   % testing with optimal step-size
      if countiter == 1
          disp('*********** sigma set to const * ||x|| ******************');
      end
      sigma = 0.04 * mueff * sqrt(sum(xmean.^2)) / N; % 20D,lam=1000:25e3
      sigma = 0.3 * mueff * sqrt(sum(xmean.^2)) / N; % 20D,lam=(40,1000):17e3
                                                     %      75e3 with def (1.5)
                                                     %      35e3 with damps=0.25
  end
  if 11 < 3 
          if countiter == 1
              disp('*********** xmean set to const ******************');
          end
      xmean = ones(N,1);
  end
  
  % Update B and D from C

  if ~flgDiagonalOnly && (ccov1+ccovmu+neg.ccov) > 0 && mod(countiter, 1/(ccov1+ccovmu+neg.ccov)/N/10) < 1
    C=triu(C)+triu(C,1)'; % enforce symmetry to prevent complex numbers
    [B,tmp] = eig(C);     % eigen decomposition, B==normalized eigenvectors
                          % effort: approx. 15*N matrix-vector multiplications
    diagD = diag(tmp); 

    if any(~isfinite(diagD))
      clear idx; % prevents error under octave 
      save(['tmp' opts.SaveFilename]);
      error(['function eig returned non-finited eigenvalues, cond(C)=' ...
	     num2str(cond(C)) ]);
    end
    if any(any(~isfinite(B)))
      clear idx; % prevents error under octave
      save(['tmp' opts.SaveFilename]);
      error(['function eig returned non-finited eigenvectors, cond(C)=' ...
	     num2str(cond(C)) ]);
    end

    % limit condition of C to 1e14 + 1
    if min(diagD) <= 0
	if stopOnWarnings
	  stopflag(end+1) = {'warnconditioncov'};
	else
	  warning(['Iteration ' num2str(countiter) ...
		   ': Eigenvalue (smaller) zero']);
	  diagD(diagD<0) = 0;
	  tmp = max(diagD)/1e14;
	  C = C + tmp*eye(N,N); diagD = diagD + tmp*ones(N,1); 
	end
    end
    if max(diagD) > 1e14*min(diagD) 
	if stopOnWarnings
	  stopflag(end+1) = {'warnconditioncov'};
	else
	  warning(['Iteration ' num2str(countiter) ': condition of C ' ...
		   'at upper limit' ]);
	  tmp = max(diagD)/1e14 - min(diagD);
	  C = C + tmp*eye(N,N); diagD = diagD + tmp*ones(N,1); 
	end
    end

    diagC = diag(C); 
    diagD = sqrt(diagD); % D contains standard deviations now
    % diagD = diagD / prod(diagD)^(1/N);  C = C / prod(diagD)^(2/N);
    BD = B.*repmat(diagD',N,1); % O(n^2)
  end % if mod

  % Align/rescale order of magnitude of scales of sigma and C for nicer output
  % not a very usual case
  if 1 < 2 && sigma > 1e10*max(diagD)
    fac = sigma / max(diagD);
    sigma = sigma/fac;
    pc = fac * pc;
    diagD = fac * diagD; 
    if ~flgDiagonalOnly
      C = fac^2 * C; % disp(fac);
      BD = B.*repmat(diagD',N,1); % O(n^2), but repmat might be inefficient todo?
    end
    diagC = fac^2 * diagC; 
  end

  if flgDiagonalOnly > 1 && countiter > flgDiagonalOnly 
    % full covariance matrix from now on 
    flgDiagonalOnly = 0; 
    B = eye(N,N);
    BD = diag(diagD);
    C = diag(diagC); % is better, because correlations are spurious anyway
  end

  % ----- numerical error management -----
  % Adjust maximal coordinate axis deviations
  if any(sigma*sqrt(diagC) > maxdx)
    sigma = min(maxdx ./ sqrt(diagC));
    %warning(['Iteration ' num2str(countiter) ': coordinate axis std ' ...
    %         'deviation at upper limit of ' num2str(maxdx)]);
    % stopflag(end+1) = {'maxcoorddev'};
  end
  % Adjust minimal coordinate axis deviations
  if any(sigma*sqrt(diagC) < mindx)
    sigma = max(mindx ./ sqrt(diagC)) * exp(0.05+cs/damps); 
    %warning(['Iteration ' num2str(countiter) ': coordinate axis std ' ...
    %         'deviation at lower limit of ' num2str(mindx)]);
    % stopflag(end+1) = {'mincoorddev'};;
  end
  % Adjust too low coordinate axis deviations
  if any(xmean == xmean + 0.2*sigma*sqrt(diagC)) 
    if stopOnWarnings
      stopflag(end+1) = {'warnnoeffectcoord'};
    else
      warning(['Iteration ' num2str(countiter) ': coordinate axis std ' ...
                'deviation too low' ]);
      if flgDiagonalOnly
        diagC = diagC + (ccov1_sep+ccovmu_sep) * (diagC .* ...
                                                  (xmean == xmean + 0.2*sigma*sqrt(diagC)));
      else
        C = C + (ccov1+ccovmu) * diag(diagC .* ...
                                      (xmean == xmean + 0.2*sigma*sqrt(diagC)));
      end
      sigma = sigma * exp(0.05+cs/damps); 
    end
  end
  % Adjust step size in case of (numerical) precision problem 
  if flgDiagonalOnly
    tmp = 0.1*sigma*diagD; 
  else
    tmp = 0.1*sigma*BD(:,1+floor(mod(countiter,N)));
  end
  if all(xmean == xmean + tmp)
    i = 1+floor(mod(countiter,N));
    if stopOnWarnings
	stopflag(end+1) = {'warnnoeffectaxis'};
    else
      warning(['Iteration ' num2str(countiter) ...
	       ': main axis standard deviation ' ...
	       num2str(sigma*diagD(i)) ' has no effect' ]);
	sigma = sigma * exp(0.2+cs/damps); 
    end
  end
  % Adjust step size in case of equal function values (flat fitness)
  % isequalfuncvalues = 0; 
  if fitness.sel(1) == fitness.sel(1+ceil(0.1+lambda/4))
    % isequalfuncvalues = 1; 
    if stopOnEqualFunctionValues
      arrEqualFunvals = [countiter arrEqualFunvals(1:end-1)];
      % stop if this happens in more than 33%
      if arrEqualFunvals(end) > countiter - 3 * length(arrEqualFunvals)
        stopflag(end+1) = {'equalfunvals'}; 
      end
    else
      if flgWarnOnEqualFunctionValues
        warning(['Iteration ' num2str(countiter) ...
		 ': equal function values f=' num2str(fitness.sel(1)) ...
		 ' at maximal main axis sigma ' ...
        num2str(sigma*max(diagD))]);
      end
      sigma = sigma * exp(0.2+cs/damps); 
    end
  end
  % Adjust step size in case of equal function values
  if countiter > 2 && myrange([fitness.hist fitness.sel(1)]) == 0  
    if stopOnWarnings
	stopflag(end+1) = {'warnequalfunvalhist'};
    else
      warning(['Iteration ' num2str(countiter) ...
	       ': equal function values in history at maximal main ' ...
	       'axis sigma ' num2str(sigma*max(diagD))]);
	sigma = sigma * exp(0.2+cs/damps); 
    end
  end
    
  % ----- end numerical error management -----
  
  % Keep overall best solution
  out.evals = counteval;
  out.solutions.evals = counteval;
  out.solutions.mean.x = xmean;
  out.solutions.mean.f = fmean;
  out.solutions.mean.evals = counteval;
  out.solutions.recentbest.x = arxvalid(:, fitness.idx(1));
  out.solutions.recentbest.f = fitness.raw(1);
  out.solutions.recentbest.evals = counteval + fitness.idx(1) - lambda;
  out.solutions.recentworst.x = arxvalid(:, fitness.idx(end));
  out.solutions.recentworst.f = fitness.raw(end);
  out.solutions.recentworst.evals = counteval + fitness.idx(end) - lambda;
  if fitness.hist(1) < out.solutions.bestever.f
    out.solutions.bestever.x = arxvalid(:, fitness.idx(1));
    out.solutions.bestever.f = fitness.hist(1);
    out.solutions.bestever.evals = counteval + fitness.idx(1) - lambda;
    bestever = out.solutions.bestever;
  end

  % Set stop flag
  if fitness.raw(1) <= stopFitness, stopflag(end+1) = {'fitness'}; end
  if counteval >= stopMaxFunEvals, stopflag(end+1) = {'maxfunevals'}; end
  if countiter >= stopMaxIter, stopflag(end+1) = {'maxiter'}; end
  if all(sigma*(max(abs(pc), sqrt(diagC))) < stopTolX) 
    stopflag(end+1) = {'tolx'};
  end
  if any(sigma*sqrt(diagC) > stopTolUpX) 
    stopflag(end+1) = {'tolupx'};
  end
  if sigma*max(diagD) == 0  % should never happen
    stopflag(end+1) = {'bug'};
  end
  if countiter > 2 && myrange([fitness.sel fitness.hist]) <= stopTolFun 
    stopflag(end+1) = {'tolfun'};
  end
  if countiter >= length(fitness.hist) && myrange(fitness.hist) <= stopTolHistFun 
    stopflag(end+1) = {'tolhistfun'};
  end
  l = floor(length(fitness.histbest)/3);
  if 1 < 2 && stopOnStagnation && ...  % leads sometimes early stop on ftablet, fcigtab
      countiter > N * (5+100/lambda) && ...  
      length(fitness.histbest) > 100 && ... 
      median(fitness.histmedian(1:l)) >= median(fitness.histmedian(end-l:end)) && ...
      median(fitness.histbest(1:l)) >= median(fitness.histbest(end-l:end))
    stopflag(end+1) = {'stagnation'};
  end

  if counteval >= stopFunEvals || countiter >= stopIter
    stopflag(end+1) = {'stoptoresume'};
    if length(stopflag) == 1 && flgsaving == 0
      error('To resume later the saving option needs to be set');
    end
  end
  % read stopping message from file signals.par 
  if flgreadsignals
    fid = fopen('./signals.par', 'rt');  % can be performance critical 
    while fid > 0
      strline = fgetl(fid); %fgets(fid, 300);
      if strline < 0 % fgets and fgetl returns -1 at end of file
        break;
      end
      % 'stop filename' sets stopflag to manual
      str = sscanf(strline, ' %s %s', 2);
      if strcmp(str, ['stop' opts.LogFilenamePrefix]) 
        stopflag(end+1) = {'manual'};
        break;
      end
      % 'skip filename run 3' skips a run, but not the last
      str = sscanf(strline, ' %s %s %s', 3);
      if strcmp(str, ['skip' opts.LogFilenamePrefix 'run'])
        i = strfind(strline, 'run');
        if irun == sscanf(strline(i+3:end), ' %d ', 1) && irun <= myeval(opts.Restarts)
          stopflag(end+1) = {'skipped'};
        end      
      end
    end % while, break 
    if fid > 0
      fclose(fid);
      clear fid; % prevents strange error under octave
    end
  end
  
  out.stopflag = stopflag;

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAVE STATE
stop = 0;
%exclude_stops = {'fitness' 'tolfun' 'tolhistfun' 'stagnation'};
if ~isempty(stopflag)
    stop = 1;
end;

arfitness = fitness.raw;
Xnew_sorted = arxvalid(:,fitness.idx);
invsqrtC = B * diag( diagD.^-1) * B';    % C^-1/2
newstate = cmaes_saveState(state, fitfun,xstart,insigma,inopts,varargin,cmaVersion,definput,defopts,flg_future_setting,nargin,input,opts,counteval, ...
 countevalNaN,irun, flgresume, xmean,N,numberofvariables,lambda0,popsize,lambda,lambda_last,stopFitness,stopMaxFunEvals, ... 
 stopMaxIter,stopFunEvals,stopIter, stopTolFun,stopTolHistFun,stopOnStagnation,stopOnWarnings,flgreadsignals,... 
 flgWarnOnEqualFunctionValues,flgEvalParallel,stopOnEqualFunctionValues,arrEqualFunvals, flgDiagonalOnly, flgActiveCMA, ... 
 noiseHandling,noiseMinMaxEvals,noiseAlphaEvals,noiseCallback,flgdisplay,flgplotting,verbosemodulo,flgscience,flgsaving, ... 
 strsaving,flgsavingfinal,savemodulo,savetime,time,maxdx,mindx,lbounds,ubounds,stopTolX,stopTolUpX,sigma,pc,diagD,diagC,B, ... 
 BD,C,fitness,bnd,out,startseed,chiN,countiter,outiter,filenameprefix,filenames, lambda_hist,mu,weights,mueff,cc,cs,ccov1, ... 
 ccovmu,ccov1_sep, ccovmu_sep,damps,noiseReevals,noiseAlpha,noiseEpsilon,noiseTheta,noisecum,noiseCutOff,arx, arxvalid,tries, ...
 noiseS,noiseSS,noiseN,xold,zmean,fmean,ps,neg, stopflag,noiseX,iterplotted,arfitness,Xnew_sorted,invsqrtC,stop);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAVE STATE


% ---------------------------------------------------------------  
% ---------------------------------------------------------------  
function flush
  if isoctave
    feval('fflush', stdout);
  end
  
% ---------------------------------------------------------------  
% ---------------------------------------------------------------  
% ----- replacements for statistic toolbox functions ------------
% ---------------------------------------------------------------  
% ---------------------------------------------------------------  
function res=myrange(x)
  res = max(x) - min(x);
  
  
% ---------------------------------------------------------------  
% ---------------------------------------------------------------  
function res = myprctile(inar, perc, idx)
%
% Computes the percentiles in vector perc from vector inar
% returns vector with length(res)==length(perc)
% idx: optional index-array indicating sorted order
%

N = length(inar);
flgtranspose = 0;

% sizes 
if size(perc,1) > 1
  perc = perc';
  flgtranspose = 1;
  if size(perc,1) > 1
    error('perc must not be a matrix');
  end
end
if size(inar, 1) > 1 && size(inar,2) > 1
  error('data inar must not be a matrix');
end
 
% sort inar
if nargin < 3 || isempty(idx)
  [sar idx] = sort(inar);
else
  sar = inar(idx);
end

res = [];
for p = perc
  if p <= 100*(0.5/N)
    res(end+1) = sar(1);
  elseif p >= 100*((N-0.5)/N)
    res(end+1) = sar(N);
  else
    % find largest index smaller than required percentile
    availablepercentiles = 100*((1:N)-0.5)/N;
    i = max(find(p > availablepercentiles));
    % interpolate linearly
    res(end+1) = sar(i) ...
	+ (sar(i+1)-sar(i))*(p - availablepercentiles(i)) ...
	/ (availablepercentiles(i+1) - availablepercentiles(i));

  end
end

if flgtranspose
  res = res';
end



% ---------------------------------------------------------------  
% --------------- TEST OBJECTIVE FUNCTIONS ----------------------  
% ---------------------------------------------------------------  

%%% Unimodal functions

function f=fjens1(x)
%
% use population size about 2*N
%
  f = sum((x>0) .* x.^1, 1);
  if any(any(x<0))
    idx = sum(x < 0, 1) > 0;
    f(idx) = 1e3;
%    f = f + 1e3 * sum(x<0, 1);
%    f = f + 10 * sum((x<0) .* x.^2, 1);
    f(idx) = f(idx) + 1e-3*abs(randn(1,sum(idx)));
%    f(idx) = NaN;
  end

function f=fsphere(x)
  f = sum(x.^2,1);

function f=fmax(x)
  f = max(abs(x), [], 1);

function f=fssphere(x)
  f=sqrt(sum(x.^2, 1));

%  lb = -0.512; ub = 512; 
%  xfeas = x; 
%  xfeas(x<lb) = lb;
%  xfeas(x>ub) = ub; 
%  f=sum(xfeas.^2, 1);
%  f = f + 1e-9 * sum((xfeas-x).^2); 
  
function f=fspherenoise(x, Nevals)
  if nargin < 2 || isempty(Nevals)
    Nevals = 1;
  end
  [N,popsi] = size(x);
%  x = x .* (1 +  0.3e-0 * randn(N, popsi)/(2*N)); % actuator noise
  fsum = 10.^(0*(0:N-1)/(N-1)) * x.^2; 
%  f = 0*rand(1,1) ...
%      + fsum ...
%      + fsum .* (2*randn(1,popsi) ./ randn(1,popsi).^0 / (2*N)) ...
%      + 1*fsum.^0.9 .* 2*randn(1,popsi) / (2*N); % 

%  f = fsum .* exp(0.1*randn(1,popsi));
  f = fsum .* (1 + (10/(N+10)/sqrt(Nevals))*randn(1,popsi));
%  f = fsum .* (1 + (0.1/N)*randn(1,popsi)./randn(1,popsi).^1);

  idx = rand(1,popsi) < 0.0;
  if sum(idx) > 0
    f(idx) = f(idx) + 1e3*exp(randn(1,sum(idx)));
  end
  
function f=fmixranks(x)
  N = size(x,1);
  f=(10.^(0*(0:(N-1))/(N-1))*x.^2).^0.5;
  if size(x, 2) > 1 % compute ranks, if it is a population 
    [ignore, idx] = sort(f);
    [ignore, ranks] = sort(idx);
    k = 9; % number of solutions randomly permuted, lambda/2-1
           % works still quite well (two time slower)
    for i = k+1:k-0:size(x,2)
      idx(i-k+(1:k)) = idx(i-k+randperm(k)); 
    end
    %disp([ranks' f'])
    [ignore, ranks] = sort(idx);
    %disp([ranks' f'])
    %pause
    f = ranks+1e-9*randn(1,1);
  end
  
function f = fsphereoneax(x)
  f = x(1)^2;
  f = mean(x)^2;
  
function f=frandsphere(x)
  N = size(x,1);
  idx = ceil(N*rand(7,1));
  f=sum(x(idx).^2);

function f=fspherelb0(x, M) % lbound at zero for 1:M needed
  if nargin < 2 M = 0; end
  N = size(x,1);
  % M active bounds, f_i = 1 for x = 0
  f = -M + sum((x(1:M) + 1).^2);
  f = f + sum(x(M+1:N).^2);
  
function f=fspherehull(x)
  % Patton, Dexter, Goodman, Punch
  % in -500..500
  % spherical ridge through zeros(N,1)
  % worst case start point seems x = 2*100*sqrt(N)
  % and small step size
  N = size(x,1);
  f = norm(x) + (norm(x-100*sqrt(N)) - 100*N)^2;
  
function f=fellilb0(x, idxM, scal) % lbound at zero for 1:M needed
  N = size(x,1);
  if nargin < 3 || isempty(scal)
    scal = 100;
  end
  scale=scal.^((0:N-1)/(N-1));
  if nargin < 2 || isempty(idxM)
    idxM = 1:N;
  end
  %scale(N) = 1e0;
  % M active bounds
  xopt = 0.1;
  x(idxM) = x(idxM) + xopt;
  f = scale.^2*x.^2;
  f = f - sum((xopt*scale(idxM)).^2); 
%  f = exp(f) - 1;
%  f = log10(f+1e-19) + 19;

  f = f + 1e-19;
  
function f=fcornersphere(x)
  w = ones(size(x,1));
  w(1) = 2.5; w(2)=2.5;
  idx = x < 0;
  f = sum(x(idx).^2);
  idx = x > 0;
  f = f + 2^2*sum(w(idx).*x(idx).^2);
  
function f=fsectorsphere(x, scal)
%
% This is deceptive for cumulative sigma control CSA in large dimension:
% The strategy (initially) diverges for N=50 and popsize = 150.  (Even
% for cs==1 this can be observed for larger settings of N and
% popsize.) The reason is obvious from the function topology. 
% Divergence can be avoided by setting boundaries or adding a
% penalty for large ||x||. Then, convergence can be observed again. 
% Conclusion: for popsize>N cumulative sigma control is not completely
% reasonable, but I do not know better alternatives. In particular:
% TPA takes longer to converge than CSA when the latter still works. 
%
  if nargin < 2 || isempty (scal)
    scal = 1e3;
  end
  f=sum(x.^2,1);
  idx = x<0;
  f = f + (scal^2 - 1) * sum((idx.*x).^2,1);
  if 11 < 3
    idxpen = find(f>1e9);
    if ~isempty(idxpen)
      f(idxpen) = f(idxpen) + 1e8*sum(x(:,idxpen).^2,1);
    end
  end
  
function f=fstepsphere(x, scal)
  if nargin < 2 || isempty (scal)
    scal = 1e0;
  end
  N = size(x,1);
  f=1e-11+sum(scal.^((0:N-1)/(N-1))*floor(x+0.5).^2);
  f=1e-11+sum(floor(scal.^((0:N-1)/(N-1))'.*x+0.5).^2);
%  f=1e-11+sum(floor(x+0.5).^2);

function f=fstep(x)
  % in -5.12..5.12 (bounded)
  N = size(x,1);
  f=1e-11+6*N+sum(floor(x));

function f=flnorm(x, scal, e)
if nargin < 2 || isempty(scal)
  scal = 1;
end
if nargin < 3 || isempty(e)
  e = 1;
end
if e==inf
  f = max(abs(x));
else
  N = size(x,1);
  scale = scal.^((0:N-1)/(N-1))';
  f=sum(abs(scale.*x).^e);
end

function f=fneumaier3(x) 
  % in -n^2..n^2
  % x^*-i = i(n+1-i)
  N = size(x,1);
%  f = N*(N+4)*(N-1)/6 + sum((x-1).^2) - sum(x(1:N-1).*x(2:N));
  f = sum((x-1).^2) - sum(x(1:N-1).*x(2:N));

function f = fmaxmindist(y)
  % y in [-1,1], y(1:2) is first point on a plane, y(3:4) second etc
  % points best
  %   5    1.4142
  %   8    1.03527618 
  %  10    0.842535997 
  %  20    0.5997   
  pop = size(y,2);
  N = size(y,1)/2;
  f = []; 
  for ipop = 1:pop
    if any(abs(y(:,ipop)) > 1)
      f(ipop) = NaN; 
    else
      x = reshape(y(:,ipop), [2, N]);
      f(ipop) = inf;
      for i = 1:N
        f(ipop) = min(f(ipop), min(sqrt(sum((x(:,[1:i-1 i+1:N]) - repmat(x(:,i), 1, N-1)).^2, 1))));
      end
    end
  end
  f = -f;

function f=fchangingsphere(x)
  N = size(x,1);
  global scale_G; global count_G; if isempty(count_G) count_G=-1; end
  count_G = count_G+1;
  if mod(count_G,10) == 0
    scale_G = 10.^(2*rand(1,N));
  end
  %disp(scale(1));
  f = scale_G*x.^2;
  
function f= flogsphere(x)
 f = 1-exp(-sum(x.^2));
  
function f= fexpsphere(x)
 f = exp(sum(x.^2)) - 1;
  
function f=fbaluja(x)
  % in [-0.16 0.16]
  y = x(1);
  for i = 2:length(x)
    y(i) = x(i) + y(i-1);
  end
  f = 1e5 - 1/(1e-5 + sum(abs(y)));

function f=fschwefel(x)
  f = 0;
  for i = 1:size(x,1),
    f = f+sum(x(1:i))^2;
  end

function f=fcigar(x, ar)
  if nargin < 2 || isempty(ar)
    ar = 1e3;
  end
  f = x(1,:).^2 + ar^2*sum(x(2:end,:).^2,1);
  
function f=fcigtab(x)
  f = x(1,:).^2 + 1e8*x(end,:).^2 + 1e4*sum(x(2:(end-1),:).^2, 1);
  
function f=ftablet(x)
  f = 1e6*x(1,:).^2 + sum(x(2:end,:).^2, 1);

function f=felli(x, lgscal, expon, expon2)
  % lgscal: log10(axis ratio)
  % expon: x_i^expon, sphere==2
  N = size(x,1); if N < 2 error('dimension must be greater one'); end

%  x = x - repmat(-0.5+(1:N)',1,size(x,2)); % optimum in 1:N
  if nargin < 2 || isempty(lgscal), lgscal = 3; end
  if nargin < 3 || isempty(expon), expon = 2; end
  if nargin < 4 || isempty(expon2), expon2 = 1; end

  f=((10^(lgscal*expon)).^((0:N-1)/(N-1)) * abs(x).^expon).^(1/expon2);
%  if rand(1,1) > 0.015
%    f = NaN;
%  end
%  f = f + randn(size(f));

function f=fellitest(x)
  beta = 0.9;
  N = size(x,1);
  f = (1e6.^((0:(N-1))/(N-1))).^beta * (x.^2).^beta; 

  
function f=fellii(x, scal)
  N = size(x,1); if N < 2 error('dimension must be greater one'); end
  if nargin < 2
    scal = 1;
  end
  f= (scal*(1:N)).^2 * (x).^2;

function f=fellirot(x)
  N = size(x,1);
  global ORTHOGONALCOORSYSTEM_G
  if isempty(ORTHOGONALCOORSYSTEM_G) ...
	|| length(ORTHOGONALCOORSYSTEM_G) < N ...
	|| isempty(ORTHOGONALCOORSYSTEM_G{N})
    coordinatesystem(N);
  end
  f = felli(ORTHOGONALCOORSYSTEM_G{N}*x);
  
function f=frot(x, fun, varargin)
  N = size(x,1);
  global ORTHOGONALCOORSYSTEM_G
  if isempty(ORTHOGONALCOORSYSTEM_G) ...
	|| length(ORTHOGONALCOORSYSTEM_G) < N ...
	|| isempty(ORTHOGONALCOORSYSTEM_G{N})
    coordinatesystem(N);
  end
  f = feval(fun, ORTHOGONALCOORSYSTEM_G{N}*x, varargin{:});
  
function coordinatesystem(N)
  if nargin < 1 || isempty(N)
    arN = 2:30;
  else
    arN = N;
  end
  global ORTHOGONALCOORSYSTEM_G
  ORTHOGONALCOORSYSTEM_G{1} = 1; 
  for N = arN
    ar = randn(N,N);
    for i = 1:N 
      for j = 1:i-1
	ar(:,i) = ar(:,i) - ar(:,i)'*ar(:,j) * ar(:,j);
      end
      ar(:,i) = ar(:,i) / norm(ar(:,i));
    end
    ORTHOGONALCOORSYSTEM_G{N} = ar; 
  end

function f=fplane(x)
  f=x(1);

function f=ftwoaxes(x)
  f = sum(x(1:floor(end/2),:).^2, 1) + 1e6*sum(x(floor(1+end/2):end,:).^2, 1);

function f=fparabR(x)
  f = -x(1,:) + 100*sum(x(2:end,:).^2,1);

function f=fsharpR(x)
  f = abs(-x(1, :)).^2 + 100 * sqrt(sum(x(2:end,:).^2, 1));
  
function f=frosen(x)
  if size(x,1) < 2 error('dimension must be greater one'); end
  N = size(x,1); 
  popsi = size(x,2); 
  f = 1e2*sum((x(1:end-1,:).^2 - x(2:end,:)).^2,1) + sum((x(1:end-1,:)-1).^2,1);
  % f = f + f^0.9 .* (2*randn(1,popsi) ./ randn(1,popsi).^0 / (2*N)); 

function f=frosenlin(x)
  if size(x,1) < 2 error('dimension must be greater one'); end

  x_org = x;
  x(x>30) = 30;
  x(x<-30) = -30;

  f = 1e2*sum(-(x(1:end-1,:).^2 - x(2:end,:)),1) + ...
      sum((x(1:end-1,:)-1).^2,1);

  f = f + sum((x-x_org).^2,1);
%  f(any(abs(x)>30,1)) = NaN; 

function f=frosenrot(x)
  N = size(x,1);
  global ORTHOGONALCOORSYSTEM_G
  if isempty(ORTHOGONALCOORSYSTEM_G) ...
	|| length(ORTHOGONALCOORSYSTEM_G) < N ...
	|| isempty(ORTHOGONALCOORSYSTEM_G{N})
    coordinatesystem(N);
  end
  f = frosen(ORTHOGONALCOORSYSTEM_G{N}*x);

function f=frosenmodif(x)
  f = 74 + 100*(x(2)-x(1)^2)^2 + (1-x(1))^2 ...
      - 400*exp(-sum((x+1).^2)/2/0.05);
  
function f=fschwefelrosen1(x)
  % in [-10 10] 
  f=sum((x.^2-x(1)).^2 + (x-1).^2);
  
function f=fschwefelrosen2(x)
  % in [-10 10] 
  f=sum((x(2:end).^2-x(1)).^2 + (x(2:end)-1).^2);

function f=fdiffpow(x)
  [N popsi] = size(x); if N < 2 error('dimension must be greater one'); end

  f = sum(abs(x).^repmat(2+10*(0:N-1)'/(N-1), 1, popsi), 1);
  f = sqrt(f); 

function f=fabsprod(x)
  f = sum(abs(x),1) + prod(abs(x),1);

function f=ffloor(x)
  f = sum(floor(x+0.5).^2,1); 

function f=fmaxx(x)
  f = max(abs(x), [], 1);

%%% Multimodal functions 

function f=fbirastrigin(x)
% todo: the volume needs to be a constant 
  N = size(x,1); 
  idx = (sum(x, 1) < 0.5*N); % global optimum
  f = zeros(1,size(x,2));
  f(idx) = 10*(N-sum(cos(2*pi*x(:,idx)),1)) + sum(x(:,idx).^2,1); 
  idx = ~idx;
  f(idx) = 0.1 + 10*(N-sum(cos(2*pi*(x(:,idx)-2)),1)) + sum((x(:,idx)-2).^2,1); 

function f=fackley(x)
  % -32.768..32.768
  % Adding a penalty outside the interval is recommended,  
  % because for large step sizes, fackley imposes like frand
  % 
  N = size(x,1); 
  f = 20-20*exp(-0.2*sqrt(sum(x.^2)/N)); 
  f = f + (exp(1) - exp(sum(cos(2*pi*x))/N));
  % add penalty outside the search interval
  f = f + sum((x(x>32.768)-32.768).^2) + sum((x(x<-32.768)+32.768).^2);
  
function f = fbohachevsky(x)
 % -15..15
  f = sum(x(1:end-1).^2 + 2 * x(2:end).^2 - 0.3 * cos(3*pi*x(1:end-1)) ...
	  - 0.4 * cos(4*pi*x(2:end)) + 0.7);
  
function f=fconcentric(x)
  % in  +-600
  s = sum(x.^2);
  f = s^0.25 * (sin(50*s^0.1)^2 + 1);

function f=fgriewank(x)
  % in [-600 600]
  [N, P] = size(x);
  f = 1 - prod(cos(x'./sqrt(1:N))) + sum(x.^2)/4e3;
  scale = repmat(sqrt(1:N)', 1, P);
  f = 1 - prod(cos(x./scale), 1) + sum(x.^2, 1)/4e3;
  % f = f + 1e4*sum(x(abs(x)>5).^2);
  % if sum(x(abs(x)>5).^2) > 0
  %   f = 1e4 * sum(x(abs(x)>5).^2) + 1e8 * sum(x(x>5)).^2;
  % end

function f=fgriewrosen(x)
% F13 or F8F2
  [N, P] = size(x);
  scale = repmat(sqrt(1:N)', 1, P);
  y = [x(2:end,:); x(1,:)];
  x = 100 * (x.^2 - y) + (x - 1).^2;  % Rosenbrock part
  f = 1 - prod(cos(x./scale), 1) + sum(x.^2, 1)/4e3;
  f = sum(1 - cos(x) + x.^2/4e3, 1);

function f=fspallpseudorastrigin(x, scal, skewfac, skewstart, amplitude)
% by default multi-modal about between -30 and 30
  if nargin < 5 || isempty(amplitude)
    amplitude = 40;
  end
  if nargin < 4 || isempty(skewstart)
    skewstart = 0;
  end
  if nargin < 3 || isempty(skewfac)
    skewfac = 1;
  end
  if nargin < 2 || isempty(scal)
    scal = 1;
  end
  N = size(x,1); 
  scale = 1;
  if N > 1
    scale=scal.^((0:N-1)'/(N-1)); 
  end
  % simple version: 
  % f = amplitude*(N - sum(cos(2*pi*(scale.*x)))) + sum((scale.*x).^2);

  % skew version: 
  y = repmat(scale, 1, size(x,2)) .* x;
  idx = find(x > skewstart);
  if ~isempty(idx)
    y(idx) =  skewfac*y(idx);
  end
  f = amplitude * (0*N-prod(cos((2*pi)^0*y),1)) + 0.05 * sum(y.^2,1) ...
      + randn(1,1);

function f=frastrigin(x, scal, skewfac, skewstart, amplitude)
% by default multi-modal about between -30 and 30
  if nargin < 5 || isempty(amplitude)
    amplitude = 10;
  end
  if nargin < 4 || isempty(skewstart)
    skewstart = 0;
  end
  if nargin < 3 || isempty(skewfac)
    skewfac = 1;
  end
  if nargin < 2 || isempty(scal)
    scal = 1;
  end
  N = size(x,1); 
  scale = 1;
  if N > 1
    scale=scal.^((0:N-1)'/(N-1)); 
  end
  % simple version: 
  % f = amplitude*(N - sum(cos(2*pi*(scale.*x)))) + sum((scale.*x).^2);

  % skew version: 
  y = repmat(scale, 1, size(x,2)) .* x;
  idx = find(x > skewstart);
  % idx = intersect(idx, 2:2:10); 
  if ~isempty(idx)
    y(idx) =  skewfac*y(idx);
  end
  f = amplitude * (N-sum(cos(2*pi*y),1)) + sum(y.^2,1);
  
function f=frastriginmax(x)
  N = size(x,1);
  f = (N/20)*807.06580387678 - (10 * (N-sum(cos(2*pi*x),1)) + sum(x.^2,1));
  f(any(abs(x) > 5.12)) = 1e2*N;

function f = fschaffer(x)
 % -100..100
  N = size(x,1);
  s = x(1:N-1,:).^2 + x(2:N,:).^2;
  f = sum(s.^0.25 .* (sin(50*s.^0.1).^2+1), 1);

function f=fschwefelmult(x)
  % -500..500
  % 
  N = size(x,1); 
  f = - sum(x.*sin(sqrt(abs(x))), 1);
  f = 418.9829*N - 1.27275661e-5*N - sum(x.*sin(sqrt(abs(x))), 1);
  % penalty term 
  f = f + 1e4*sum((abs(x)>500) .* (abs(x)-500).^2, 1);
  
function f=ftwomax(x)
  % Boundaries at +/-5
  N = size(x,1); 
  f = -abs(sum(x)) + 5*N;

function f=ftwomaxtwo(x)
  % Boundaries at +/-10
  N = size(x,1); 
  f = abs(sum(x));
  if f > 30
    f = f - 30;
  end
  f = -f;
  
function f=frand(x)
  f=1./(1-rand(1, size(x,2))) - 1;
