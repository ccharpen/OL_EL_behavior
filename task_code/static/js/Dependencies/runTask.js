// Actual trials
// Create timeline variable
var taskTimeline_var = [];
for (var trial = 0; trial < numTotalTrials; trial++) {
    taskTimeline_var.push({
        tv_blockID: blockID[trial],
        tv_trialID: trialID[trial],
        tv_leftObsStim: leftObsStim[trial],
        tv_rightObsStim: rightObsStim[trial],
        tv_leftObsChStim: leftObsChStim[trial],
        tv_rightObsChStim: rightObsChStim[trial],
        tv_obsTokStim: obsTokStim[trial],
        tv_leftPlayStim: leftPlayStim[trial],
        tv_rightPlayStim: rightPlayStim[trial],
        tv_play_LResp_OutStim: play_LResp_OutStim[trial],
        tv_play_RResp_OutStim: play_RResp_OutStim[trial],
        tv_goalToken: goalToken[trial],
        tv_tokenLeft: tokenLeft[trial],
        tv_tokenRight: tokenRight[trial],
        tv_leftOutVal: leftOutVal[trial],
        tv_rightOutVal: rightOutVal[trial],
        tv_rewMag: rewMag[trial],
        tv_corrResp: corrResp[trial],
        tv_isBreak: isBreak[trial]
    });
}

// Fullscreen
var go_full = {
    type: "fullscreen",
    fullscreen_mode: true,
    data: {
        label: "go_full",
    },
};

// Dummy ITI
// Create a 'dummy' ITI object of duration = 0. This will store task attributes
// for conditionals later on
var dummy = {
    type: "html-keyboard-response",
    stimulus: "<head><style> body {background-color: black;} </style></head> <div></div>",
    choices: jsPsych.NO_KEYS,
    trial_duration: 0,
    response_ends_trial: false,
    data: {
        label: "dummy",
        blockID: jsPsych.timelineVariable("tv_blockID"),
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
        corrResp: jsPsych.timelineVariable("tv_corrResp"),
        isBreak: jsPsych.timelineVariable("tv_isBreak")
    },
};

var observeOn = {
    type: "html-keyboard-response",
    choices: jsPsych.NO_KEYS,
    trial_duration: obsOnDur,
    stimulus:
        function () {
            let prev_data = jsPsych.data.get().filter({label: "dummy"}).last(1).values()[0];
            let leftHTML = "<img src='" + prev_data.leftObsStim + "' height='100%' />";
            let rightHTML = "<img src='" + prev_data.rightObsStim + "' height='100%' />";
            return get2boxesHTMLstr(leftHTML, rightHTML, "OBSERVE")
        },
    response_ends_trial: false,
    data: {label: "observeOn"},
    on_start: function () {
        console.log("Observe On")
    }
};

var observeCh = {
    type: "html-keyboard-response",
    choices: jsPsych.NO_KEYS,
    trial_duration: obsChDur,
    stimulus:
        function () {
            let prev_data = jsPsych.data.get().filter({label: "dummy"}).last(1).values()[0];
            let leftHTML = "<img src='" + prev_data.leftObsChStim + "' height='100%' />";
            let rightHTML = "<img src='" + prev_data.rightObsChStim + "' height='100%' />";
            return get2boxesHTMLstr(leftHTML, rightHTML, "OBSERVE")
        },
    response_ends_trial: false,
    data: {label: "observeCh"},
    on_start: function () {
        console.log("Observe Choice")
    }
};

var observeTok = {
    type: "html-keyboard-response",
    choices: jsPsych.NO_KEYS,
    trial_duration: obsTokDur,
    stimulus:
        function () {
            let prev_data = jsPsych.data.get().filter({label: "dummy"}).last(1).values()[0];
            let htmlstr = '<head><style> body {background-color: black;} </style></head>' +
                "<img src='" + prev_data.obsTokStim + "' class='boxopen'/>";
            return htmlstr
        },
    response_ends_trial: false,
    data: {label: "observeTok"},
    on_start: function () {
        console.log("Observe Token")
    }
};

var playOn = {
    type: "html-keyboard-response",
    choices: allowable_keys,
    trial_duration: cueDur,
    stimulus:
        function () {
            let prev_data = jsPsych.data.get().filter({label: "dummy"}).last(1).values()[0];
            let leftHTML ="<img src='" + prev_data.leftPlayStim + "' height='100%' />";
            let rightHTML = "<img src='" + prev_data.rightPlayStim + "' height='100%' />";
            return get2tokensHTMLstr(leftHTML, rightHTML, "CHOOSE");
        },
    response_ends_trial: true,
    data: {label: "playOn"},
    on_start: function () {
        console.log("Play On")
    },
    on_finish: function(data) {
        data.respKey = jsPsych.pluginAPI.convertKeyCodeToKeyCharacter(data.key_press);
        let cr = jsPsych.data.get().filter({label: "dummy"}).last(1).values()[0].corrResp;
        if ((data.respKey == "leftarrow" && cr == "left") || (data.respKey == "rightarrow" && cr == "right")) {
            data.isCorr = 1;
        } else {
            data.isCorr = 0;
        };
    }
};

var playOut = {
    type: "html-keyboard-response",
    choices: jsPsych.NO_KEYS,
    trial_duration: outDur,
    stimulus:
        function () {
            let prev_data = jsPsych.data.get().filter({label: "dummy"}).last(1).values()[0];
            let key_press = jsPsych.data.get().last(1).values()[0].key_press;
            if (key_press == choice_keys.left) {
                outTok = prev_data.play_LResp_OutStim;
                outVal = prev_data.leftOutVal;
            } else if (key_press == choice_keys.right) {
                outTok = prev_data.play_RResp_OutStim;
                outVal = prev_data.rightOutVal;
            } else console.log("Error.. key press: " + key_press);
            if (outVal == 0) {
                outText = "0 point";
            } else {
                outText = "+" + outVal.toString() + " points";
            }
            let out_image = "<img src='" + outTok + "' height='100%' />";
            return getOutcomeHTMLstr(out_image, outText);
        },
    response_ends_trial: false,
    data: {label: "playOut"},
    on_start: function () {
        console.log("Play Outcome")
    },
    on_finish: function(data){
      let prev_data = jsPsych.data.get().filter({label: "dummy"}).last(1).values()[0];
      let key_press = jsPsych.data.get().last(2).values()[0].key_press;
      if (key_press == choice_keys.left) {
          data.choice = 1;
          data.currOutcome = prev_data.leftOutVal;
          data.tokenShown = prev_data.tokenLeft;
      } else if (key_press == choice_keys.right) {
          data.choice = 0;
          data.currOutcome = prev_data.rightOutVal;
          data.tokenShown = prev_data.tokenRight;
      } else {
          data.currOutcome = 0;
          data.tokenShown = "none";
      };
      data.isGoal = (prev_data.goalToken == data.tokenShown) ? 1 : 0;
      if (prev_data.blockID) { //main task trials
        fieldname = 'maintask_bl' + prev_data.blockID + '.trial';
        block = prev_data.blockID;
        isBreak = prev_data.isBreak;
      } else { //trials from practice2
        fieldname = 'prac2.trial';
        block = 1;
        isBreak = 0;
      };
    }
};

var no_response = {
    type: 'html-keyboard-response',
    choices: jsPsych.NO_KEYS,
    stimulus: function () {
        var html = '<head><style> body {background-color: black;} </style></head>' +
            "<div class=\"centered whiteText\"><p>" +
            "Missed response! <br /><br /> Make your choice within 3 seconds." +
            "</p></div>"
        return html;
    },
    trial_duration: outDur + outDur,
    response_ends_trial: false,
    data: {label: "no_response"},
    on_finish: function(data){
      let prev_data = jsPsych.data.get().filter({label: "dummy"}).last(1).values()[0];
      if (prev_data.blockID) { //main task trials
        fieldname = 'maintask_bl' + prev_data.blockID + '.trial';
        block = prev_data.blockID;
        isBreak = prev_data.isBreak;
      } else { //trials from practice2
        fieldname = 'prac2.trial';
        block = 1;
        isBreak = 0;
      };
    }
};

var task_break = {
    type: "html-keyboard-response",
    choices: [choice_keys.space],
    stimulus: function () {
        let outcomes = jsPsych.data.get().filter({label: "playOut"}).select('currOutcome');
        let out_prac_exp = jsPsych.data.get().filter({label: "prac_done"}).select('pracTotOutc');
        let out_prac_obs = jsPsych.data.get().filter({label: "prac2_done"}).select('outcTotPracObs');
        let curr_out = outcomes.sum() + out_prac_exp.sum() - out_prac_obs.sum();
        htmlstr = "<head><style> body {background-color: black;} </style></head>" + //'<head><link rel="alternate stylesheet" href="/static/css/custom.css" type="text/css" /></head>' +
            '<body><div class="centered whiteText"><p>' +
            '  Time for a short break.<br /><br />' +
            '  You\'ve earned a total of ' + curr_out + ' points so far.<br /><br />' +
            '  Press space when you\'re ready to continue' +
            '</p></div></body>';
        return htmlstr;
    },
    response_ends_trial: true,
    trial_duration: null,
    data: {label: "task_break"},
    on_finish: function(data) {
        let prev_data = jsPsych.data.get().filter({label: "dummy"}).last(1).values()[0];
        let outcomes = jsPsych.data.get().filter({label: "playOut"}).select('currOutcome');
        let out_prac_exp = jsPsych.data.get().filter({label: "prac_done"}).select('pracTotOutc');
        let out_prac_obs = jsPsych.data.get().filter({label: "prac2_done"}).select('outcTotPracObs');
        let curr_out = outcomes.sum() + out_prac_exp.sum() - out_prac_obs.sum();
    }
};

var task_done = {
    type: 'html-keyboard-response',
    choices: [choice_keys.space],
    response_ends_trial: true,
    stimulus: function () {
        let outcomes = jsPsych.data.get().filter({label: "playOut"}).select('currOutcome');
        let out_prac_exp = jsPsych.data.get().filter({label: "prac_done"}).select('pracTotOutc');
        let out_prac_obs = jsPsych.data.get().filter({label: "prac2_done"}).select('outcTotPracObs');
        let score = outcomes.sum() + out_prac_exp.sum() - out_prac_obs.sum();
        let htmlstr = "<head><style> body {background-color: black;} </style></head>" + //'<head><link rel="alternate stylesheet" href="/static/css/custom.css" type="text/css" /></head>' +
            '<body><div class="centered whiteText"><p>' +
            '  You\'re done.<br /><br />' +
            '  You earned a total of ' + score + ' points.<br /><br />' +
            '  Press space to continue.\n' +
            '</p></div></body>';
        return htmlstr
    },
    data: {
        label: "task_done",
        taskname: "Box Opening",
    },
    on_finish: function(data) {
        let outcomes = jsPsych.data.get().filter({label: "playOut"}).select('currOutcome');
        let out_prac_exp = jsPsych.data.get().filter({label: "prac_done"}).select('pracTotOutc');
        let out_prac_obs = jsPsych.data.get().filter({label: "prac2_done"}).select('outcTotPracObs');
        let score = outcomes.sum() + out_prac_exp.sum()  - out_prac_obs.sum();
        data.finalPoints = score;
    }
};

var iti = {
    type: "html-keyboard-response",
    trial_duration: itiDur,
    stimulus: function () {
        let htmlstr = "<head><style> body {background-color: black;} </style></head>" +
            '<body><div class="centered whiteText largeFont"><p>' +
            '  +  ' +
            '</p></div></body>';
        return htmlstr
    },
    choices: jsPsych.NO_KEYS,
    response_ends_trial: false,
    data: {label: "iti"},
    on_finish: function (data) {
        data.trial_timeStamp = t - t_start - d_start;
        t = new Date().getTime();
        data.time = t;
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

var task_break_ifNode = {
    timeline: [task_break],
    conditional_function: function () {
        var prev_data = jsPsych.data.get().filter({label: "dummy"}).last(1).values()[0];
        return (prev_data.isBreak == 1 && prev_data.blockID != 8);
    }
};

var task_done_ifNode = {
    timeline: [task_done],
    conditional_function: function (data) {
        return (jsPsych.data.get().filter({label: "observeOn"}).count() == numTotalTrials + numPract2Trials)
    }
};

// timeline for one trial, combine ifnodes and fixed events
var task_trialProcedure = {
    timeline: [dummy, iti, observeOn, observeCh, observeTok, playOn, task_resp_ifNode, task_noResp_ifNode, task_break_ifNode, task_done_ifNode],
    timeline_variables: taskTimeline_var,
};
