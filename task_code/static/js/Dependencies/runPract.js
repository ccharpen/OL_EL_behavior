// Instrumental trials
// Create timeline variable
var pracTimeline_var = [];
for (var trial = 0; trial < numPractTrials; trial++) {
    pracTimeline_var.push({
        tv_trialID: prac_trialID[trial],
        tv_leftStim: prac_leftStim[trial],
        tv_rightStim: prac_rightStim[trial],
        tv_leftStim_LResp: prac_leftStim_LResp[trial],
        tv_rightStim_LResp: prac_rightStim_LResp[trial],
        tv_leftStim_RResp: prac_leftStim_RResp[trial],
        tv_rightStim_RResp: prac_rightStim_RResp[trial],
        tv_TokStim_LResp: prac_TokStim_LResp[trial],
        tv_TokStim_RResp: prac_TokStim_RResp[trial],
        tv_OutStim_LResp: prac_OutStim_LResp[trial],
        tv_OutStim_RResp: prac_OutStim_RResp[trial],
        tv_TokValStim: prac_TokValStim[trial],
        tv_goalToken: prac_goalToken[trial],
        tv_tokenIfLeft: prac_tokenIfLeft[trial],
        tv_tokenIfRight: prac_tokenIfRight[trial],
        tv_leftOutVal: prac_leftOutVal[trial],
        tv_rightOutVal: prac_rightOutVal[trial],
        tv_rewMag: prac_rewMag[trial],
        tv_corrResp: prac_corrResp[trial],
    });
}

// Dummy ITI
// Create a 'dummy' ITI object of duration = 0. This will store task attributes
// for conditionals later on
var prac_dummy = {
    type: "html-keyboard-response",
    stimulus: "<head><style> body {background-color: black;} </style></head> <div></div>",
    choices: jsPsych.NO_KEYS,
    trial_duration: 0,
    response_ends_trial: false,
    data: {
        label: "prac_dummy",
        prac_trialID: jsPsych.timelineVariable("tv_trialID"),
        prac_leftStim: jsPsych.timelineVariable("tv_leftStim"),
        prac_rightStim: jsPsych.timelineVariable("tv_rightStim"),
        prac_leftStim_LResp: jsPsych.timelineVariable("tv_leftStim_LResp"),
        prac_rightStim_LResp: jsPsych.timelineVariable("tv_rightStim_LResp"),
        prac_leftStim_RResp: jsPsych.timelineVariable("tv_leftStim_RResp"),
        prac_rightStim_RResp: jsPsych.timelineVariable("tv_rightStim_RResp"),
        prac_TokStim_LResp: jsPsych.timelineVariable("tv_TokStim_LResp"),
        prac_TokStim_RResp: jsPsych.timelineVariable("tv_TokStim_RResp"),
        prac_OutStim_LResp: jsPsych.timelineVariable("tv_OutStim_LResp"),
        prac_OutStim_RResp: jsPsych.timelineVariable("tv_OutStim_RResp"),
        prac_TokValStim: jsPsych.timelineVariable("tv_TokValStim"),
        prac_goalToken: jsPsych.timelineVariable("tv_goalToken"),
        prac_tokenIfLeft: jsPsych.timelineVariable("tv_tokenIfLeft"),
        prac_tokenIfRight: jsPsych.timelineVariable("tv_tokenIfRight"),
        prac_leftOutVal: jsPsych.timelineVariable("tv_leftOutVal"),
        prac_rightOutVal: jsPsych.timelineVariable("tv_rightOutVal"),
        prac_rewMag: jsPsych.timelineVariable("tv_rewMag"),
        prac_corrResp: jsPsych.timelineVariable("tv_corrResp")
    },
};

var prac_tokVal = {
    type: "html-keyboard-response",
    choices: jsPsych.NO_KEYS,
    trial_duration: 2000,
    stimulus:
        function () {
            let prev_data = jsPsych.data.get().filter({label: "prac_dummy"}).last(1).values()[0];
            let htmlstr = '<head><style> body {background-color: black;} </style></head>' +
                "<img src='" + prev_data.prac_TokValStim + "' height='100%' />";
            return htmlstr
        },
    response_ends_trial: false,
    data: {label: "prac_tokVal"},
    on_start: function () {
        console.log("Practice Token Values")
    }
};

var prac_playOn = {
    type: "html-keyboard-response",
    choices: allowable_keys,
    trial_duration: cueDur,
    stimulus:
        function () {
            let prev_data = jsPsych.data.get().filter({label: "prac_dummy"}).last(1).values()[0];
            let leftHTML ="<img src='" + prev_data.prac_leftStim + "' height='100%' />";
            let rightHTML = "<img src='" + prev_data.prac_rightStim + "' height='100%' />";
            return get2boxesHTMLstr(leftHTML, rightHTML, "CHOOSE");
        },
    response_ends_trial: true,
    data: {label: "prac_playOn"},
    on_start: function () {
        console.log("Practice Play On")
    },
    on_finish: function(data) {
        data.respKey = jsPsych.pluginAPI.convertKeyCodeToKeyCharacter(data.key_press);
        let cr = jsPsych.data.get().filter({label: "prac_dummy"}).last(1).values()[0].prac_corrResp;
        if ((data.respKey == "leftarrow" && cr == "left") || (data.respKey == "rightarrow" && cr == "right")) {
            data.isCorr = 1;
        } else {
            data.isCorr = 0;
        };
    }
};

var prac_playCh = {
    type: "html-keyboard-response",
    choices: jsPsych.NO_KEYS,
    trial_duration: respDur,
    stimulus:
        function () {
            let prev_data = jsPsych.data.get().filter({label: "prac_dummy"}).last(1).values()[0];
            let key_press = jsPsych.data.get().last(1).values()[0].key_press;
            if (key_press == choice_keys.left) {
                leftBox = prev_data.prac_leftStim_LResp;
                rightBox = prev_data.prac_rightStim_LResp;
            } else if (key_press == choice_keys.right) {
                leftBox = prev_data.prac_leftStim_RResp;
                rightBox = prev_data.prac_rightStim_RResp;
            } else console.log("Error.. key press: " + key_press);
            let leftHTML ="<img src='" + leftBox + "' height='100%' />";
            let rightHTML = "<img src='" + rightBox + "' height='100%' />";
            return get2boxesHTMLstr(leftHTML, rightHTML, "");
        },
    response_ends_trial: false,
    data: {label: "prac_playCh"},
    on_start: function () {
        console.log("Practice Play Choice")
    }
};

var prac_playTok = {
    type: "html-keyboard-response",
    choices: jsPsych.NO_KEYS,
    trial_duration: tokDur,
    stimulus:
        function () {
            let prev_data = jsPsych.data.get().filter({label: "prac_dummy"}).last(1).values()[0];
            let key_press = jsPsych.data.get().last(2).values()[0].key_press;
            if (key_press == choice_keys.left) {
                tokBox = prev_data.prac_TokStim_LResp;
            } else if (key_press == choice_keys.right) {
                tokBox = prev_data.prac_TokStim_RResp;
            } else console.log("Error.. key press: " + key_press);
            let htmlstr = '<head><style> body {background-color: black;} </style></head>' +
                "<img src='" + tokBox + "' class='boxopen' />";
            return htmlstr
        },
    response_ends_trial: false,
    data: {label: "prac_playTok"},
    on_start: function () {
        console.log("Practice Play Token")
    }
};

var prac_playOut = {
    type: "html-keyboard-response",
    choices: jsPsych.NO_KEYS,
    trial_duration: outDur,
    stimulus:
        function () {
            let prev_data = jsPsych.data.get().filter({label: "prac_dummy"}).last(1).values()[0];
            let key_press = jsPsych.data.get().last(3).values()[0].key_press;
            if (key_press == choice_keys.left) {
                outTok = prev_data.prac_OutStim_LResp;
                outVal = prev_data.prac_leftOutVal;
            } else if (key_press == choice_keys.right) {
                outTok = prev_data.prac_OutStim_RResp;
                outVal = prev_data.prac_rightOutVal;
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
    data: {label: "prac_playOut"},
    on_start: function () {
        console.log("Practice Play Outcome")
    },
    on_finish: function(data){
      let prev_data = jsPsych.data.get().filter({label: "prac_dummy"}).last(1).values()[0];
      let key_press = jsPsych.data.get().last(4).values()[0].key_press;
      if (key_press == choice_keys.left) {
          data.choice = 1;
          data.currOutcome = prev_data.prac_leftOutVal;
          data.tokenShown = prev_data.prac_tokenIfLeft;
      } else if (key_press == choice_keys.right) {
          data.choice = 0;
          data.currOutcome = prev_data.prac_rightOutVal;
          data.tokenShown = prev_data.prac_tokenIfRight;
      } else {
          data.currOutcome = 0;
          data.tokenShown = "none";
      };
      data.isGoal = (prev_data.prac_goalToken == data.tokenShown) ? 1 : 0;
    }
};

var prac_no_response = {
    type: 'html-keyboard-response',
    choices: jsPsych.NO_KEYS,
    stimulus: function () {
        var html = '<head><style> body {background-color: black;} </style></head>' +
            "<div class=\"centered whiteText\"><p>" +
            "Missed response! <br /><br /> Make your choice within 3 seconds." +
            "</p></div>"
        return html;
    },
    trial_duration: respDur + tokDur + outDur,
    response_ends_trial: false,
    data: {label: "prac_no_response"},
};

var prac_done = {
    type: 'html-keyboard-response',
    choices: [choice_keys.space],
    response_ends_trial: true,
    stimulus: function () {
        let outcome_tot = jsPsych.data.get().filter({label: "prac_playOut"}).select('currOutcome').sum();
        if (jsPsych.data.get().filter({label: "prac_redo"}).count() == 1) { // subtract outcomes from the first practice
           outcomes = outcome_tot - jsPsych.data.get().filter({label: "prac_redo"}).select('outcTotPrac1').sum();
        } else {
           outcomes = outcome_tot;
        };
        htmlstr = "<head><style> body {background-color: black;} </style></head>" + //'<head><link rel="alternate stylesheet" href="/static/css/custom.css" type="text/css" /></head>' +
            '<body><div class="centered whiteText"><p>' +
            '  You\'re done with this practice.<br /><br />'+
            '  You\'ve earned ' + outcomes + ' points so far.<br /><br />' +
            '  Press space to continue.<br /><br />'+
            '</p></div></body>';
        return htmlstr;
    },
    data: {label: "prac_done"},
    on_finish: function (data) {
      let outcome_tot = jsPsych.data.get().filter({label: "prac_playOut"}).select('currOutcome').sum();
      if (jsPsych.data.get().filter({label: "prac_redo"}).count() == 1) { // subtract outcomes from the first practice
         outcomes = outcome_tot - jsPsych.data.get().filter({label: "prac_redo"}).select('outcTotPrac1').sum();
      } else {
         outcomes = outcome_tot;
      };
      data.pracTotOutc = outcomes;
      data.redo = 0;
    }
};

var prac_redo = {
    type: 'html-keyboard-response',
    choices: [choice_keys.space],
    response_ends_trial: true,
    stimulus: function () {
        htmlstr = "<head><style> body {background-color: black;} </style></head>" + //'<head><link rel="alternate stylesheet" href="/static/css/custom.css" type="text/css" /></head>' +
            '<body><div class="centered whiteText"><p>' +
            '  Your accuracy on these practice trials was too low.<br /><br />'+
            '  Please press space to try again.<br /><br />'+
            '</p></div></body>';
        return htmlstr;
    },
    data: {label: "prac_redo"},
    on_finish: function (data) {
      data.redo = 1;
      data.outcTotPrac1 = jsPsych.data.get().filter({label: "prac_playOut"}).select('currOutcome').sum();
    }
};

var last_instruct_page = {
    type: "html-keyboard-response",
    // dummy stimulus and key choices are reset in on_start()
    stimulus: '<body><div class="centered"><img src=' + instructDir + 'instr_10.png' + ' alt="instruction" width="100%"></div></body>',
    choices: [choice_keys.space],
    response_ends_trial: true,
    data: {label: 'last_instruct_page'},
};

// ifnode - conditional to determine whether mini-timeline is run
// ifnode must come after the events called in the mini-timeline
var prac_resp_ifNode = {
    timeline: [prac_playCh, prac_playTok, prac_playOut],
    conditional_function: function () {
        let prev_data = jsPsych.data.get().last(1).values()[0];
        return prev_data.key_press == choice_keys.left || prev_data.key_press == choice_keys.right ;
    }
};

var prac_noResp_ifNode = {
    timeline: [prac_no_response],
    conditional_function: function (data) {
        let prev_data = jsPsych.data.get().last(1).values()[0];
        return prev_data.label == 'prac_playOn';
    }
};

var prac_done_ifNode = {
    timeline: [prac_done],
    conditional_function: function (data) {
        let accuracy = jsPsych.data.get().filter({label: "prac_playOn"}).select('isCorr');
        let totalNT = jsPsych.data.get().filter({label: "prac_dummy"}).count();
        return totalNT == pracTimeline_var.length && accuracy.mean()>0.65;
    }
};

var prac_redo_ifNode = {
    timeline: [prac_redo, last_instruct_page],
    conditional_function: function (data) {
        let accuracy = jsPsych.data.get().filter({label: "prac_playOn"}).select('isCorr');
        let totalNT = jsPsych.data.get().filter({label: "prac_dummy"}).count();
        return totalNT == pracTimeline_var.length && accuracy.mean()<0.65;
    }
};

var prac_done2_ifNode = {
    timeline: [prac_done],
    conditional_function: function (data) {
        let totalNT = jsPsych.data.get().filter({label: "prac_dummy"}).count();
        return totalNT == 2*pracTimeline_var.length;
    }
};

// timeline for one trial, combine ifnodes and fixed events
var prac_trialProcedure = {
    timeline: [prac_dummy, iti, prac_tokVal, prac_playOn, prac_resp_ifNode, prac_noResp_ifNode, prac_done_ifNode, prac_redo_ifNode],
    timeline_variables: pracTimeline_var,
};

var prac_trialProcedure_redo = {
    timeline: [prac_dummy, iti, prac_tokVal, prac_playOn, prac_resp_ifNode, prac_noResp_ifNode, prac_done2_ifNode],
    timeline_variables: pracTimeline_var,
    conditional_function: function (data) {
        return jsPsych.data.get().filter({label: "prac_redo"}).count() == 1; // && accuracy.mean()<0.65;
    }
};
