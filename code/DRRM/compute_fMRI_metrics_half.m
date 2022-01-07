%%		
% Description  -- function [NRMSE_raw, NRMSE_new, ...
%           R2_raw_mean,R2_new_mean, ...
%           AIC_raw, AIC_new, ...
%           pc_raw, pc_new] = compute_fMRI_metrics_half(fMRI, pRF, pRF0)
%		Compute metrics by fMRI.
%
% Parameter(s): 
%		fMRI[double array]  --  fMRI data.
%		pRF[double array]   --  prf model to calculate new R2 value.
%		pRF0[double array]  --  prf model to calculate R2 value.
% return:
%		NRMSE_raw[double]   --  Normalized RMSE of the pRF solution in fMRI	fitting.
%		NRMSE_new[double]   --  Normalized RMSE of the pRF0 solution in fMRI fitting
%		R2_raw_mean[double] --  Average Valiance explained for pRF solution
%		R2_new_mean[double] --  Average Valiance explained for pRF0 solution
%		AIC_raw[double]     --  AIC metric for pRF solution
%		AIC_new[double]     --  AIC metric for pRF0 solution
%		pc_raw[double]      --  Pearson correlation between fMRI signal and predicting for pRF solution
%		pc_new[double]      --  Pearson correlation between fMRI signal and predicting for pRF0 solution
%%
function [NRMSE_raw, NRMSE_new, ...
          R2_raw_mean,R2_new_mean, ...
          AIC_raw, AIC_new, ...
          pc_raw, pc_new] = compute_fMRI_metrics_half(fMRI, pRF, pRF0)
    % pRF = x y sigma gain n
    %% convert polar angel to pixel unit to pixel

    %% refit fMRI  
    stimulus =  load_hcp_stimuli();
    TR =1; hrf = getcanonicalhrf(TR,TR)';          % HRF that was used in the model
    degs =[3 3 3 3 3 3];  % vector of maximum polynomial degrees used in the model
    res = size(stimulus{1}(:,:,1));
    resmx = max(res);
    % Pre-compute cache for faster execution
    [~,xx,yy] = makegaussian2d(resmx,2,2,2,2);

    % Prepare the stimuli for use in the model
    stimulusPP = {};
    for p=1:length(stimulus)
      stimulusPP{p} = squish(stimulus{p},2)';  % this flattens the image so that the dimensionality is now frames x pixels
      stimulusPP{p} = stimulusPP{p}(151:300,:); % get later half
      stimulusPP{p} = [stimulusPP{p} p*ones(size(stimulusPP{p},1),1)];  % this adds a dummy column to indicate run breaks
    end

    % Construct projection matrices that fit and remove the polynomials.
    % Note that a separate projection matrix is constructed for each run.
    polymatrix = {};
    for p=1:length(degs)
      polymatrix{p} = projectionmatrix(constructpolynomialmatrix(150,0:degs(p)));
    end

    modelfun = @(pp,dd) conv2run(posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2))); 0]) .^ posrect(pp(5)),hrf,dd(:,prod(res)+1));


    R2array = [];
    rmsearray = [];
    errarray = []; 
    pcarray = [];
    for vx = 1 : size(pRF0,1)

        % we need vary power n
        nopt = 0.05;  

        parameters = [pRF0(vx,1:2)  pRF0(vx,6)*sqrt(nopt) pRF0(vx,3) sqrt(nopt)];
        [p,r,e, pc]=R2bypara(parameters,fMRI,polymatrix, degs, vx,modelfun, stimulusPP);
        R2array(vx) =p;
        rmsearray(vx)=r;
        errarray(vx,:)=e;
        pcarray(vx)=pc;
        vx;
    end


    %% Check new R2
    R2arraynew = []; 
    rmsearraynew = [];
    errarraynew = []; 
    pearsonnew = [];
    for vx = 1 : size(pRF,1)

        % we need vary power n
        nopt = 0.05;     
        parameters = [pRF(vx,1:2)  pRF(vx,6)*sqrt(nopt) pRF(vx,3) sqrt(nopt)];
        [p, r, e, pc]=R2bypara(parameters, fMRI,polymatrix, degs, vx,modelfun, stimulusPP);
        R2arraynew(vx) =p;
        rmsearraynew(vx)=r;
        errarraynew(vx,:)=e;
        pearsonnew(vx) = pc;
        vx;
    end

    AIC_raw = calAIC(errarray);
    AIC_new = calAIC(errarraynew);


    R2_new_mean= nanmean(R2arraynew);
    R2_raw_mean= nanmean(R2array);


    NRMSE_raw = nanmean(rmsearray);
    NRMSE_new =  nanmean(rmsearraynew);

    pc_raw =  nanmean(pcarray);
    pc_new =  nanmean(pearsonnew);

end


 

 
function [R2, rmse_res, err, pearcorr]=R2bypara(parameters,fMRI,polymatrix, degs, vx,modelfun, stimulusPP)    
    % For each run, collect the data and the model fit.  We project out polynomials
    % from both the data and the model fit.  This deals with the problem of
    % slow trends in the data.
    datats = {};
    modelts = {};
    N = size(fMRI,2);
    for p=1:length(degs)
        datats{p} =  polymatrix{p}*fMRI(vx,(p-1)*150+1:(p)*150)';
        modelts{p} = polymatrix{p}*modelfun(parameters,stimulusPP{p});
    end
    modelfit = cat(1,modelts{:});
     modelfit = (modelfit-mean(modelfit)) /  std(modelfit);
     if any(isnan(modelfit))
         modelfit = zeros(N,1);
     end
    fmridata = cat(1,datats{:});
     fmridata = (fmridata - mean(fmridata) )/ std(fmridata);
%     R2 = calccod(modelfit, fmridata);
    R2  = calccorrelation(modelfit, fmridata);
    pearcorr = corr(modelfit, fmridata);
    
    if (isnan(R2))
        R2 = nan;
    end
    
     if (isnan(pearcorr))
        pearcorr = nan;
    end
    if nargout > 1
        rmse_res= rmse( fmridata, modelfit);
        err = modelfit - fmridata;      
    end
end
 