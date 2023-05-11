function [f, vals] = LL_FixArb(params, P)
%fixed weight arbitration model, learning of reward magnitude
%edit 10/2/2021 - made choice on first trial of each block random, as OL
%should not predict a clear preference on the first trial
%edit 10/25/2021 - added reward magnitude boost during EL
%edit 11/11/2021 - added counterfactual update for OL
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
%params(1): softmax OL decision beta [0 +Inf]
%params(2): softmax EL decision beta [0 +Inf]
%params(3): learning rate observed token values [0 1]
%params(4): learning rate experienced token values [0 1]
%params(5): weight of OL over EL [0 1]
%params(6): boosting after large rewards

%transform parameters
params(1) = exp(params(1)); 
params(2) = exp(params(2)); 
params(3) = 1/(1+exp(-params(3)));
params(4) = 1/(1+exp(-params(4)));
params(5) = 1/(1+exp(-params(5)));
mag_boost = params(6);

obs_ch       = P(:,9);
token_sh_obs = P(:,11);
choice       = P(:,14);
iscorr       = P(:,16);
token_sh     = P(:,17);
outc         = P(:,18);
outc(outc>0) = 1;
outc_mag     = P(:,18);

nt = length(P(:,3));

magnitude   = nan(nt,2); %magnitude of reward separately for each token
TokVals_obs = nan(nt,2); %value of tokens inferred from partner
TokVals_exp = nan(nt,2); %value of tokens inferred from outcomes
ActVals_obs = nan(nt,2); %value of partner's actions (probability that box 1 or box 2 lead to orange token)
oSPE        = nan(nt,1); %obs state prediction error (=uncertainty of observational learning predictions)
eRPE        = nan(nt,1); %reward prediction error (=uncertainty of experiential learning predictions)
Por_obs     = nan(nt,1); %probability of choosing orange according to OL
Por_exp     = nan(nt,1); %probability of choosing orange according to EL

LH      = nan(nt,1); %likelihood values
Por     = nan(nt,1); %probability of choosing orange
ll_corr = nan(nt,1); %likelihood of being correct
for t=1:nt
      
    %initialize values
    if t==1
        TokVals_exp_init = [0.5 0.5];
        magnitude(t,:) = [0 0];
    else
        TokVals_exp_init = TokVals_exp(t-1,:);
        magnitude(t,:) = magnitude(t-1,:);
    end
    if P(t,4)==1 %action values initialized at the beginning of each block (new boxes)
        ActVals_init = [0.5 0.5];
    else
        ActVals_init = ActVals_obs(t-1,:);
    end
        
    %learn association between partner's action and token
    if token_sh_obs(t)==1 %partner obtained orange token
        R_tok = 1;
    else %partner obtained blue token
        R_tok = 0;
    end
    if obs_ch(t)==1 %partner chose box 1 so update value of box 1 >ly and box 2 <ly
        oSPE(t) = R_tok - ActVals_init(1);
        ActVals_obs(t,1) = ActVals_init(1) + params(3)*oSPE(t);
        ActVals_obs(t,2) = ActVals_init(2) - params(3)*oSPE(t);
    else %partner chose box 2 so update value of box 2 >ly and box 2 <ly
        oSPE(t) = R_tok - ActVals_init(2);
        ActVals_obs(t,2) = ActVals_init(2) + params(3)*oSPE(t);
        ActVals_obs(t,1) = ActVals_init(1) - params(3)*oSPE(t);
    end
        
    %now infer the likelihood that partner's goal is orange vs blue given
    %their action - simplest inference is just to take the estimated transition
    %probability associated with the chosen action, and translate that into
    %the probability that orange is the goal
    if P(t,4) == 1 %random choice on the first trial of each block
        TokVals_obs(t,:) = [0.5 0.5];
    else
        TokVals_obs(t,:) = [ActVals_obs(t,obs_ch(t)) 1-ActVals_obs(t,obs_ch(t))];
    end
        
    % calculate choice probability according to OL
    val_diff_obs = TokVals_obs(t,1)-TokVals_obs(t,2);
    Por_obs(t) = 1./(1 + exp(-params(1)*(val_diff_obs)));
    
    % calculate choice probability according to EL (previous trial0
    val_diff_exp = TokVals_exp_init(1)-TokVals_exp_init(2);
    mag_diff = magnitude(t,1)-magnitude(t,2);
    Por_exp(t) = 1./(1 + exp(-params(2)*val_diff_exp - mag_boost*mag_diff));
    
    %calculate integrated choice probability
    Por(t) = params(5)*Por_obs(t) + (1-params(5))*Por_exp(t);

    %calculate choice loglikelihood depending on participant's choice
    %if choice value is 1, use one part of likelihood contribution.
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
        eRPE(t) = outc(t) - TokVals_exp_init(1);
        TokVals_exp(t,1) = TokVals_exp_init(1) + params(4)*eRPE(t);
        TokVals_exp(t,2) = 1 - TokVals_exp(t,1);
    elseif token_sh(t) == 2 %blue
        eRPE(t) = outc(t) - TokVals_exp_init(2);
        TokVals_exp(t,2) = TokVals_exp_init(2) + params(4)*eRPE(t);
        TokVals_exp(t,1) = 1 - TokVals_exp(t,2);
    else %no update
        TokVals_exp(t,:) = TokVals_exp_init;
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
vals = [LH Por ll_corr ActVals_obs TokVals_obs Por_obs TokVals_exp Por_exp magnitude oSPE eRPE];
end