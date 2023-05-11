// Actual trials
// Create timeline variable
var prac2Timeline_var = [];
for (var trial = 0; trial < numPract2Trials; trial++) {
    prac2Timeline_var.push({
        tv_trialID: prac2_trialID[trial],
        tv_leftObsStim: prac2_leftObsStim[trial],
        tv_rightObsStim: prac2_rightObsStim[trial],
        tv_leftObsChStim: prac2_leftObsChStim[trial],
        tv_rightObsChStim: prac2_rightObsChStim[trial],
        tv_obsTokStim: prac2_obsTokStim[trial],
        tv_leftPlayStim: prac2_leftPlayStim[trial],
        tv_rightPlayStim: prac2_rightPlayStim[trial],
        tv_play_LResp_OutStim: prac2_play_LResp_OutStim[trial],
        tv_play_RResp_OutStim: prac2_play_RResp_OutStim[trial],
        tv_goalToken: prac2_goalToken[trial],
        tv_tokenLeft: prac2_tokenLeft[trial],
        tv_tokenRight: prac2_tokenRight[trial],
        tv_leftOutVal: prac2_leftOutVal[trial],
        tv_rightOutVal: prac2_rightOutVal[trial],
        tv_rewMag: prac2_rewMag[trial],
        tv_corrResp: prac2_corrResp[trial]
    });
}

// Dummy ITI
// Create a 'dummy' ITI object of duration = 0. This will store task attributes
// for conditionals later on
var dummy_prac2 = {
    type: "html-keyboard-response",
    stimulus: "<head><style> body {background-color: black;} </style></head> <div></div>",
    choices: jsPsych.NO_KEYS,
    trial_duration: 0,
    response_ends_trial: false,
    data: {
        label: "dummy",
        trialID: jsPsych.timelineVariable("tv_trialID"),
        leftObsStim: jsPsych.timelineVariable("tv_leftObsStim"),
        rightObsStim: jsPsych.timelineVariable("tv_rightObsStim"),
        leftObsChStim: jsPsych.timelineVariable("tv_leftObsChStim"),
        rightObsChStim: jsPsych.timelineVariable("tv_rightObsChStim"),
        obsTokStim: jsPsych.timelineVariable("tv_obsTokStim"),
        leftPlayStim: jsPsych.timelineVariable("tv_leftPlayStim"),
        rightPlayStim: jsPsych.timelineVariable("tv_rightPlayStim"),
        play_LResp_OutStim: jsPsych.timelineVariable("tv_play_LResp_OutStim"),
        play_RResp_OutStim: jsPsych.timelineVariable("tv_play_RResp_OutStim"),
        goalToken: jsPsych.timelineVariable("tv_goalToken"),
        tokenLeft: jsPsych.timelineVariable("tv_tokenLeft"),
        tokenRight: jsPsych.timelineVariable("tv_tokenRight"),
        leftOutVal: jsPsych.timelineVariable("tv_leftOutVal"),
        rightOutVal: jsPsych.timelineVariable("tv_rightOutVal"),
        rewMag: jsPsych.timelineVariable("tv_rewMag"),
        corrResp: jsPsych.timelineVariable("tv_corrResp")
    }
};

var prac2_done = {
    type: 'html-keyboard-response',
    choices: [choice_keys.space],
    response_ends_trial: true,
    stimulus: '<body><div class="centered"><img src=' + instructDir + "instruction_final.png" + ' width="100%"></div></body>',
    data: {label: "prac2_done"},
    on_finish: function (data) {
      data.outcTotPracObs = jsPsych.data.get().filter({label: "playOut"}).select('currOutcome').sum();
    }
};

// ifnode - conditional to determine whether mini-timeline is run
// ifnode must come after the events called in the mini-timeline
var task_resp_ifNode = {
    timeline: [playOut],
    conditional_function: function () {
        let prev_data = jsPsych.data.get().last(1).values()[0];
        return prev_data.key_press == choice_keys.left || prev_data.key_press == choice_keys.right ;
    }
};

var task_noResp_ifNode = {
    timeline: [no_response],
    conditional_function: function (data) {
        let prev_data = jsPsych.data.get().last(1).values()[0];
        return prev_data.label == 'playOn';
    }
};

var prac2_done_ifNode = {
    timeline: [prac2_done],
    conditional_function: function (data) {
        return (jsPsych.data.get().filter({label: "observeOn"}).count() == prac2Timeline_var.length)
    }
};

// timeline for one trial, combine ifnodes and fixed events
var prac2_trialProcedure = {
    timeline: [dummy_prac2, iti, observeOn, observeCh, observeTok, playOn, task_resp_ifNode, task_noResp_ifNode, prac2_done_ifNode],
    timeline_variables: prac2Timeline_var,
};
