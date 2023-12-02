function [f, vals] = LL_ExpLearn_mag(params, P)
%model implementing experiential learning only (learning token value only
%from experiencing outcomes associated with token)
%this version of the model assumes that subjects know the token values are
%not independent, i.e. P(reward blue) = 1 - P(reward orange)
%learning of reward magnitude

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
%params(1): softmax decision beta [0 +Inf]
%params(2): learning rate experienced token values [0 1]

%transform parameters
params(1) = exp(params(1)); 
params(2) = 1/(1+exp(-params(2)));

nt = length(P(:,3));
choice   = P(:,14);
iscorr   = P(:,16);
token_sh = P(:,17);
outc     = P(:,18);

TokVals      = nan(nt,2); %value of tokens (orange vs blue) before update
TokVals_post = nan(nt,2); %value of tokens (orange vs blue) after update
eRPE = nan(nt,1);

LH      = nan(nt,1); %likelihood values
Por     = nan(nt,1); %probability of choosing orange
ll_corr = nan(nt,1); %likelihood of being correct
for t=1:nt
      
    %initialize values
    if t==1
        TokVals(t,:) = [25 25];
    else
        TokVals(t,:) = TokVals_post(t-1,:);
    end
    
    %calculate choice loglikelihood depending on participant's choice
    %if choice value is 1, use one part of likelihood contribution.
    Por(t) = 1./(1 + exp(-params(1)*(TokVals(t,1)-TokVals(t,2))));
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
    TokVals_post(t,:) = TokVals(t,:);
    if token_sh(t) ~= 0
        eRPE(t) = outc(t) - TokVals(t,token_sh(t));
        TokVals_post(t,token_sh(t)) = TokVals(t,token_sh(t)) + params(2)*eRPE(t);
    end

end
f = nansum(log(LH + eps));
vals = [LH Por ll_corr TokVals eRPE];
end