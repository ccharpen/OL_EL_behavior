var pagenum = 0;    // holds current page number
var pagenum_2 = 0;    // holds current page number
var instr_index = 0;
var instr2_index = 0;
var choice = [];    // holds keys available from current page

var initInstruct = {
    type: "html-keyboard-response",
    stimulus: "<body><div></div></body>",
    choices: jsPsych.NO_KEYS,
    trial_duration: 0,
    response_ends_trial: false,
    data: {
        label: "initInstruct",
        page_index: pagenum,
        pageStep: 1,    // holds 1 if stepping to next page, -1 if stepping to previous
        prac_designNumber: prac_designNo,
        designNumber: designNo
    }
};

// event to display instruction images
var instruct_page = {
    type: "html-keyboard-response",
    // dummy stimulus and key choices are reset in on_start()
    stimulus: "<body><div></div></body>",
    choices: jsPsych.NO_KEYS,
    trial_duration: null,
    response_ends_trial: true,
    data: {label: 'instruction'},
    on_start: function (trial) {
        var prev_data = jsPsych.data.get().last(1).values()[0];
        pagenum = prev_data.page_index + prev_data.pageStep;
        trial.stimulus = '<body><div class="centered"><img src=' + instruct_pages[pagenum-1] + ' alt="instruction" width="100%"></div></body>';
        choice = choiceKey_array[pagenum-1];
        trial.choices = choice;     // set keys available from current page
    },
    on_finish: function (data) {
        data.page_index = pagenum;
        // set pageStep, which will be added to page_index to compute the next pagenum
        if (data.key_press == choice[0]) {
            data.pageStep = 1
        } else if  (data.key_press == choice[1]) {
            data.pageStep = -1
        } else if (data.key_press == choice[2] || data.key_press == choice[3]) {
            data.pageStep = 0
        } else throw new Error("invalid key press: " + data.key_press);
        instr_index = instr_index + 1;
    }
};

//second set of instructions, after practice
var initInstruct_2 = {
    type: "html-keyboard-response",
    stimulus: "<body><div></div></body>",
    choices: jsPsych.NO_KEYS,
    trial_duration: 0,
    response_ends_trial: false,
    data: {
        label: "initInstruct_2",
        page_index: pagenum_2,
        pageStep: 1,    // holds 1 if stepping to next page, -1 if stepping to previous
    }
};

var instruct_page_2 = {
    type: "html-keyboard-response",
    // dummy stimulus and key choices are reset in on_start()
    stimulus: "<body><div></div></body>",
    choices: jsPsych.NO_KEYS,
    trial_duration: null,
    response_ends_trial: true,
    data: {label: 'instruction_2'},
    on_start: function (trial) {
        var prev_data = jsPsych.data.get().last(1).values()[0];
        pagenum_2 = prev_data.page_index + prev_data.pageStep;
        trial.stimulus = '<body><div class="centered"><img src=' + instruct_pages_2[pagenum_2-1] + ' alt="instruction" width="100%"></div></body>';
        choice = choiceKey_array_2[pagenum_2-1];
        trial.choices = choice;     // set keys available from current page
    },
    on_finish: function (data) {
        data.page_index = pagenum_2;
        // set pageStep, which will be added to page_index to compute the next pagenum
        if (data.key_press == choice[0]) {
            data.pageStep = 1
        } else if  (data.key_press == choice[1]) {
            data.pageStep = -1
        } else if (data.key_press == choice[2] || data.key_press == choice[3]) {
            data.pageStep = 0
        } else throw new Error("invalid key press: " + data.key_press)
        instr2_index = instr2_index + 1;
    }
};

// timeline for one trial, combine ifnodes and fixed events
var instruct_trialLoop = {
    timeline: [instruct_page],
    // use a loop function to play instruct_page event with different instruction pages
    // until all instruction pages have been shown
    loop_function: function (data) {
        var prev_data = data.values()[0];
        if (prev_data.pageStep < 0 || prev_data.page_index < Ninstruct) {
            return true;
        }
        else {
            console.log("out of pages.");
            return false;
        }
    }
};

var instruct_trialLoop_2 = {
    timeline: [instruct_page_2],
    // use a loop function to play instruct_page event with different instruction pages
    // until all instruction pages have been shown
    loop_function: function (data) {
        var prev_data = data.values()[0];
        if (prev_data.pageStep < 0 || prev_data.page_index < Ninstruct_2) {
            return true;
        }
        else {
            console.log("out of pages.");
            return false;
        }
    }
};

var instruct_trialProcedure = {
    timeline: [initInstruct, instruct_trialLoop],
};

var instruct_trialProcedure_2 = {
    timeline: [initInstruct_2, instruct_trialLoop_2],
};
