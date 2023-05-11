function [f, vals] = LL_ObsLearn(params, P)
%model implementing observational learning only (learning token value from 
%observing partner's choices and inferring goal token)
%edit 10/2/2021 - made choice on first trial of each block random, as OL
%should not predict a clear preference on the first trial
%edit 11/11/2021 - added counterfactual update, assuming that when partner
%(and subject) updates the token probability associated with the chosen
%Box, they also update the unchosen Box accordingly, assuming that the more
%likely one Box is to give orange tokens, the more likely the other Box
%will give blue tokens.

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
%params(2): learning rate observed token values [0 1]

%transform parameters
params(1) = exp(params(1)); 
params(2) = 1/(1+exp(-params(2)));

nt = length(P(:,3));
token_sh_obs = P(:,11);
obs_ch = P(:,9);  
choice = P(:,14);
iscorr = P(:,16);

ActVals    = nan(nt,2); %value of partner's actions (probability that box 1 or box 2 lead to orange token)
TokVals    = nan(nt,2); %value of tokens (orange vs blue)
oSPE        = nan(nt,1); %obs state prediction error (=uncertainty of observational learning predictions)

LH = nan(nt,1); %likelihood values
Por     = nan(nt,1); %probability of choosing orange
ll_corr = nan(nt,1); %likelihood of being correct
for t=1:nt
      
    %initialize values
    if t==1 %action values initialized on trial 1
        ActVals_init = [0.5 0.5];
    else
        if P(t,1) ~= P(t-1,1) %action values also initialized at the beginning of each block (new boxes)
            ActVals_init = [0.5 0.5];
        else
            ActVals_init = ActVals(t-1,:);
        end
    end
        
    %learn association between partner's action and token
    if token_sh_obs(t)==1 %partner obtained orange token
        R_tok = 1;
    else %partner obtained blue token
        R_tok = 0;
    end
    if obs_ch(t)==1 %partner chose box 1 so update value of box 1 >ly and box 2 <ly
        oSPE(t) = R_tok - ActVals_init(1);
        ActVals(t,1) = ActVals_init(1) + params(2)*oSPE(t);
        ActVals(t,2) = ActVals_init(2) - params(2)*oSPE(t);
    else %partner chose box 2 so update value of box 2 >ly and box 2 <ly
        oSPE(t) = R_tok - ActVals_init(2);
        ActVals(t,2) = ActVals_init(2) + params(2)*oSPE(t);
        ActVals(t,1) = ActVals_init(1) - params(2)*oSPE(t);
    end
    
    %now infer the likelihood that partner's goal is orange vs blue given
    %their action - simplest inference is just to take the estimated transition
    %probability associated with the chosen action, and translate that into
    %the probability that orange is the goal
    if P(t,4) == 1 %random choice on the first trial of each block
        TokVals(t,:) = [0.5 0.5];
    else
        TokVals(t,:) = [ActVals(t,obs_ch(t)) 1-ActVals(t,obs_ch(t))];
    end
    
    val_diff = TokVals(t,1)-TokVals(t,2);
    
    %calculate choice loglikelihood depending on participant's choice
    %if choice value is 1, use one part of likelihood contribution.
    Por(t) = 1./(1 + exp(-params(1)*(val_diff)));
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

end
f = nansum(log(LH + eps));
vals = [LH Por ll_corr ActVals TokVals oSPE];
end