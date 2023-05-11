function [f, vals] = LL_ExpLearn(params, P)
%model implementing experiential learning only (learning token value only
%from experiencing outcomes associated with token)
%no learning of token reward magnitude (model treats reward as 1 and no
%reward as 0)
%additional parameter that boosts the probability to stay after a large (vs
%a small) reward
%edit 11/18/21: update magnitude value separately for each token and make
%boosting function of the magnitude difference; only update when reward is
%obtained; if not decay by half. Also decay by half the value of the
%unchosen token.
%edit 8/26/22: magnitude boosting mechanism outside the softmax

%description of P
%column 1: block number
%column 2: block ID
%column 3: trial number
%column 4: trial number per block
%column 5: common transition probability (associated with Box 1; 1:common, 0:rare)
%column 6: whether common (1) or rare (0) transition happens
%column 7: reward probability of orange token (1-blue token)
%column 8: goal token (1: orange, 2: blue)
%column 9: box chosen by partner (1: box 1 majority orange, 2: box 2 majority blue)
%column 10: whether partner's choice is correct (1) or not (0)
%column 11: obs token shown (1: orange, 2: blue)
%column 12: whether token is goal (1) or not (0)
%column 13: participant choice (1: left, 0: right)
%column 14: participant choice (1: majority orange box, 0: majority blue box)
%column 15: rt
%column 16: whether participant choice is correct (1) or not (0)
%column 17: token shown (1: orange, 2: blue)
%column 18: outcome
%column 19: missed trial (1)
%column 20: action chosen by partner (1: left, 0: right)
%column 21: reward if orange
%column 22: reward magnitude

%description of params
%params(1): softmax decision beta [0 10]
%params(2): learning rate experienced token values [0 1]
%params(3): boosting after large rewards

%transform parameters
params(1) = exp(params(1)); 
params(2) = 1/(1+exp(-params(2)));

nt = length(P(:,3));
choice   = P(:,14);
iscorr   = P(:,16);
token_sh = P(:,17);
outc     = P(:,18);
outc(outc>0) = 1;
outc_mag = P(:,18);

TokVals      = nan(nt,2); %value of tokens (orange vs blue) before update
TokVals_post = nan(nt,2); %value of tokens (orange vs blue) after update
eRPE         = nan(nt,1);
magnitude    = nan(nt,2); %magnitude of reward separately for each token

LH      = nan(nt,1); %likelihood values
Por     = nan(nt,1); %probability of choosing orange
ll_corr = nan(nt,1); %likelihood of being correct
for t=1:nt
      
    %initialize values
    if t==1
        TokVals(t,:) = [0.5 0.5];
        magnitude(t,:) = [0 0];
    else
        TokVals(t,:) = TokVals_post(t-1,:);
        magnitude(t,:) = magnitude(t-1,:);
    end
    
    val_diff = TokVals(t,1)-TokVals(t,2);
    mag_diff = magnitude(t,1)-magnitude(t,2);
    
    %calculate choice loglikelihood depending on participant's choice
    %if choice value is 1, use one part of likelihood contribution.
    Por(t) = 1./(1 + exp(-params(1)*val_diff - params(3)*mag_diff));
    if choice(t) == 1
        LH(t) = Por(t);
        if iscorr(t) == 1
            ll_corr(t) = Por(t);
        else
            ll_corr(t) = 1 - Por(t);
        end
    %if choice value is 0, use other part of likelihood contribution    
    elseif choice(t) == 0
        LH(t) = 1 - Por(t);    
        if iscorr(t) == 1
            ll_corr(t) = 1 - Por(t);
        else
            ll_corr(t) = Por(t);
        end
    end   
    
    %update token values accordingly
    if token_sh(t) == 1 %orange
        eRPE(t) = outc(t) - TokVals(t,1);
        TokVals_post(t,1) = TokVals(t,1) + params(2)*eRPE(t);
        TokVals_post(t,2) = 1 - TokVals_post(t,1);
    elseif token_sh(t) == 2 %blue
        eRPE(t) = outc(t) - TokVals(t,2);
        TokVals_post(t,2) = TokVals(t,2) + params(2)*eRPE(t);
        TokVals_post(t,1) = 1 - TokVals_post(t,2);
    else %no update
        TokVals_post(t,:) = TokVals(t,:);
    end
    
    if token_sh(t)>0
        if outc_mag(t)>0
            magnitude(t,token_sh(t)) = outc_mag(t)/100;
            magnitude(t,3-token_sh(t)) = 0.5*magnitude(t,3-token_sh(t));
        else
            magnitude(t,token_sh(t)) = 0.5*magnitude(t,token_sh(t));
        end
    end
    
end
f = nansum(log(LH + eps));
vals = [LH Por ll_corr TokVals magnitude eRPE];
end