// End page slide
var endTaskSlide = {
    type: 'html-keyboard-response',
    stimulus: '<div>End of the demo! You can now close your browser.<br><br>' +
    '</div>',
    choices: jsPsych.NO_KEYS,
    data: {label: "all_done"}
};

// timeline for one trial, combine ifnodes and fixed events
var endTask_procedure = {
    timeline: [endTaskSlide],
};
