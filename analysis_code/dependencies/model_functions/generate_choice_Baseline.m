function P_pred = generate_choice_Baseline(params, P)
%degenerate model implementing the 7 following biases:
%- hand bias (preference for left over right)
%- color bias (preference for orange over blue)
%- sticky color (bias towards choosing same token as previous)
%- action imitation (repeat left/right action that was just performed by partner)
%edit 8/26/2022 - remove softmax beta because of non-identifiability

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
%column 14: participant choice (1: orange token, 0: blue token)
%column 15: rt
%column 16: whether participant choice is correct (1) or not (0)
%column 17: token shown (1: orange, 2: blue)
%column 18: outcome
%column 19: missed trial (1)
%column 20: action chosen by partner (1: left, 0: right)
%column 21: reward if orange
%column 22: reward magnitude

%description of params
%params(1): color bias
%params(2): hand bias 
%params(3): sticky action
%params(4): action imitation

%transform parameters if needed
color_bias = params(1);
hand_bias = params(2);
sticky_action = params(3);
action_imitation = params(4);

%extract relevant variables from P matrix
nt = length(P(:,3));

action   = P(:,13); %participant's action (1: left, 0: right)
choice   = P(:,14); %participant's token choice (1: orange, 0: blue)
part_act = P(:,20); part_act(P(:,20)==0)=-1; %action performed by partner (1: left, -1: right)
goal_tok = P(:,8);

pred_act   = nan(nt,1); %predicted action (1: left, 0: right)
pred_ch    = nan(nt,1); %predicted choice
pred_token = nan(nt,1); %predicted token
iscorr     = nan(nt,1); %is generated choice correct
pred_outc  = nan(nt,1); %predicted outcome
pred_outc_mag = nan(nt,1); %predicted outcome magnitude

prev_act = nan(nt,1); %action performed on trial t-1
posOrTok = nan(nt,1); %position of orange token (1: left, -1 right)
LH      = nan(nt,1); %likelihood values
Por     = nan(nt,1); %probability of choosing orange
ll_corr = nan(nt,1); %likelihood of being correct

for t=1:nt
      
    %initialize values
    if t==1
        prev_act(t) = 0; 
    else
        if isnan(pred_act(t-1))
            prev_act(t) = 0;
        else
            if pred_act(t-1) == 0
                prev_act(t) = -1;
            else 
                prev_act(t) = 1;
            end
        end
    end
    
    %for the left-right parameters we need to compute whether the current
    %orange token is on the left (1) or on the right (-1)
    if ~isnan(choice(t))
        if action(t) == choice(t)
            posOrTok(t) = 1;
        else
            posOrTok(t) = -1;
        end
    else
        posOrTok(t) = 0;
    end  
    
    val_diff = color_bias + hand_bias*posOrTok(t) + ...
        sticky_action*prev_act(t)*posOrTok(t) + action_imitation*part_act(t)*posOrTok(t);
     
    %calculate choice loglikelihood depending on participant's choice
    %if choice value is 1, use one part of likelihood contribution.
    Por(t) = 1./(1 + exp(-val_diff));
    
    %generate choice, likelihood
    n=rand();
    if n<=Por(t)
        pred_ch(t) = 1;
        pred_token(t) = 1;
        LH(t) = Por(t);
        if goal_tok(t) == 1
            iscorr(t) = 1;
            ll_corr(t) = Por(t);
        else
            iscorr(t) = 0;
            ll_corr(t) = 1-Por(t);
        end
    else
        pred_ch(t) = 0;
        pred_token(t) = 2;
        LH(t) = 1 - Por(t);
        if goal_tok(t) == 2
            iscorr(t) = 1;
            ll_corr(t) = 1-Por(t);
        else
            iscorr(t) = 0;
            ll_corr(t) = Por(t);
        end
    end
    
    %determine left/right action
    if posOrTok(t) == 1 %orange is on the left, blue is on the right
        pred_act(t) = pred_ch(t);
    elseif posOrTok(t) == -1 %orange is on the right, blue is on the left
        pred_act(t) = 1-pred_ch(t);
    else
        pred_act(t) = nan;
    end
    
    %determine outcome
    if (pred_ch(t)==1 && P(t,21)==1) || (pred_ch(t)==0 && P(t,21)==0)
        pred_outc(t) = 1;
        pred_outc_mag(t) = P(t,22)/100; %scale magnitude
    else
        pred_outc(t) = 0;
        pred_outc_mag(t) = 0; 
    end
        
end
P_pred = [LH Por ll_corr pred_ch pred_token iscorr pred_outc pred_outc_mag pred_act];
end