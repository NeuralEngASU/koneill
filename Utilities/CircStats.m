% Extension of the CircStat package
function [circMean, circMed, circVar, vMParams, vMScale, RMSE ,vMCorr, vMR2, circStd, circSkew, circKurt] = CircStats(tmpDeltaPhi)

binEdge = [-pi:pi/100:pi];

[counts,centers] = hist(tmpDeltaPhi,binEdge);

% Mean
% [circMean(:,:,1), circMean(:,:,2), circMean(:,:,3)] = circ_mean(repmat(centers, 1, size(tmpDeltaPhi,3)), counts, 1);

% tic
[circMean(:,:,1), circMean(:,:,2), circMean(:,:,3)] = circ_mean(tmpDeltaPhi);
circMean = permute(circMean, [2,1,3]);
% deltaT = toc;
% fprintf('Mean time: %f\n', deltaT);

% Median
% tic
% circMed(:,:,1) = circ_median(tmpDeltaPhi);
% circMed = circMed(:);
circMed = -1*ones(size(circMean,1),1,1);
% deltaT = toc;
% fprintf('Median time: %f\n', deltaT);

% Varience
% tic
circVar(:,:,1) = circ_var(tmpDeltaPhi);
circVar = circVar(:);
% deltaT = toc;
% fprintf('Varience time: %f\n', deltaT);
% std
% [circStd(:,:,1), circStd(:,:,2)] = circ_std(tmpDeltaPhi);
% circStd = permute(circStd, [2,1,3]);
circStd = -1*ones(size(circVar,1),1,2);


% Skewness
% [circSkew(:,:,1), circSkew(:,:,2)] = circ_skewness(tmpDeltaPhi);
% circSkew = permute(circSkew, [2,1,3]);
circSkew = -1*ones(size(circVar,1),1,2);

% Kurtosis
% [circKurt(:,:,1), circKurt(:,:,2)] = circ_kurtosis(tmpDeltaPhi);
% circKurt = permute(circKurt, [2,1,3]);
circKurt = -1*ones(size(circVar,1),1,2);
% vMParams
sizeTime = size(tmpDeltaPhi,2);

vMParams = zeros(size(circVar,1),1,2);

tic
for ii = 1:sizeTime
    [vMParams(ii,:,1), vMParams(ii,:,2)] = circ_vmpar(centers, counts(:,ii));
    [vm, ~] = circ_vmpdf(centers, vMParams(ii,:,1), vMParams(ii,:,2));
    
    if ~sum(isnan(vm))
        options = optimoptions(@fminunc,'Algorithm','quasi-newton', 'Display', 'off');
        coeffs = fminunc(@(c) PLISqrErr(c,centers, counts(:,ii), centers, vm),[1;1], options);
        
        vMScale(ii,:,1) = coeffs(2);
        
        RMSE(ii,:,1) = sqrt(mean((counts(:,ii) - coeffs(2)*vm).^2));
        RMSE(ii,:,2) = sqrt(circ_mean((counts(:,ii) - coeffs(2)*vm).^2));
        
        %     [vMCorr(:,ii,1), vMCorr(:,ii,2)] = circ_corrcc(counts(:,ii),  coeffs(2)*vm);
        
        coreCoeffs = corrcoef(counts(:,ii),  coeffs(2)*vm);
        vMCorr(ii,:,1) = coreCoeffs(1,2);
    else
        vMScale(ii,:,1) = 0;
        RMSE(ii,:,1) = -1;
        RMSE(ii,:,2) = -1;
        vMCorr(ii,:,1) = 0;
    end %END IF 
        
%     if vMCorr(ii,:,1) < 0.5
%         disp('Stop here')
%         disp(ii)
%     end
    

end % END FOR
deltaT = toc;
fprintf('VM time: %f\n', deltaT);

vMR2(:,:,1) = vMCorr(:,:,1).^2;

end % END FUNCTION

% EOF