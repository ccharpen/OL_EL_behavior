// Initialize the timeline
var timeline = [];

// Set up trial event durations
var obsOnDur = 2000;
var obsChDur = 500;
var obsTokDur = 1000;
var cueDur = 3000;
var respDur = 500;
var tokDur = 1000;
var outDur = 1000;
var itiDur = 1500;

// Set up clock
var d_start = new Date();
var t_start = 0;

// Specify image directories
var imageDir = "static/images/";
var instructDir = imageDir + "instructions/";
var stimDir = imageDir + "stimuli/";
// Define paths to general images
var iti = imageDir + "fix.png";

var tokenDict = {
    "1": "orange",
    "2": "blue",
};

// look up table for available keys
var choice_keys = {
    "left": jsPsych.pluginAPI.convertKeyCharacterToKeyCode('leftarrow'),
    "right": jsPsych.pluginAPI.convertKeyCharacterToKeyCode('rightarrow'),
    "f": jsPsych.pluginAPI.convertKeyCharacterToKeyCode('f'),
    "h": jsPsych.pluginAPI.convertKeyCharacterToKeyCode('h'),
    "k": jsPsych.pluginAPI.convertKeyCharacterToKeyCode('k'),
    "z": jsPsych.pluginAPI.convertKeyCharacterToKeyCode('z'),
    "b": jsPsych.pluginAPI.convertKeyCharacterToKeyCode('b'),
    "space": jsPsych.pluginAPI.convertKeyCharacterToKeyCode('space'),
};
var allowable_keys = [choice_keys.left, choice_keys.right];
var left_allowableKey = choice_keys.left;
var right_allowableKey = choice_keys.right;

// STORE DESIGN ID FOR EACH PARTICIPANT
// Load in design from csv file
var numDesigns = designArray.main.length;
var designNo = getRandomIntInclusive(0, numDesigns - 1)
var designMat = designArray.main[designNo];

// These variables will change depending on the task
var blockID = Array(designMat.length).fill(0);
var trialID = Array(designMat.length).fill(0);
var leftObsStim = [];
var rightObsStim = [];
var leftObsChStim = [];
var rightObsChStim = [];
var obsTokStim = [];
var leftPlayStim = [];
var rightPlayStim = [];
var play_LResp_OutStim = [];
var play_RResp_OutStim = [];
var goalToken = Array(designMat.length).fill(0);
var tokenLeft = Array(designMat.length).fill(0);
var tokenRight = Array(designMat.length).fill(0);
var leftOutVal = Array(designMat.length).fill(0);
var rightOutVal = Array(designMat.length).fill(0);
var rewMag = Array(designMat.length).fill(0);
var corrResp = Array(designMat.length).fill(0);
var isBreak = Array(designMat.length).fill(0);

// Convert from JSON object to set of arrays
for (t = 0; t < designMat.length; t++) {
    blockID[t] = Number(designMat[t].blNb);
    trialID[t] = Number(designMat[t].trNbBl);
    goalToken[t] = designMat[t].goalToken;
    tokenLeft[t] = designMat[t].tokenLeft;
    tokenRight[t] = designMat[t].tokenRight;
    leftOutVal[t] = Number(designMat[t].outcomeIfLeft);
    rightOutVal[t] = Number(designMat[t].outcomeIfRight);
    rewMag[t] = Number(designMat[t].rewMag);
    corrResp[t] =  designMat[t].corrResp;
    isBreak[t] = (Number(designMat[t].trNbBl) == 20) ? 1 : 0; //break every 30 trials

    leftObsStim.push(stimDir + designMat[t].leftBoxObs + ".png");
    rightObsStim.push(stimDir + designMat[t].rightBoxObs + ".png");
    if ((Number(designMat[t].chBoxObs) == 1) && (designMat[t].posBoxObs1 == "left") || (Number(designMat[t].chBoxObs) == 2) && (designMat[t].posBoxObs1 == "right")) {
        //left box chosen by partner
        leftObsChStim.push(stimDir + designMat[t].leftBoxObs + "_ch.png");
        rightObsChStim.push(stimDir + designMat[t].rightBoxObs + ".png");
        obsTokStim.push(stimDir + designMat[t].leftBoxObs + "_" + designMat[t].tokObs + ".png");
    } else {
        //right box chosen by partner
        leftObsChStim.push(stimDir + designMat[t].leftBoxObs + ".png");
        rightObsChStim.push(stimDir + designMat[t].rightBoxObs + "_ch.png");
        obsTokStim.push(stimDir + designMat[t].rightBoxObs + "_" + designMat[t].tokObs + ".png");
    };

    leftPlayStim.push(stimDir + "token_" + designMat[t].tokenLeft + ".png");
    rightPlayStim.push(stimDir + "token_" + designMat[t].tokenRight + ".png");
    play_LResp_OutStim.push(stimDir + "token_" + designMat[t].tokenLeft + ".png");
    play_RResp_OutStim.push(stimDir + "token_" + designMat[t].tokenRight + ".png");
}

var preloadSet = new Set();
leftObsStim.forEach(img => preloadSet.add(img));
rightObsStim.forEach(img => preloadSet.add(img));
leftObsChStim.forEach(img => preloadSet.add(img));
rightObsChStim.forEach(img => preloadSet.add(img));
obsTokStim.forEach(img => preloadSet.add(img));
leftPlayStim.forEach(img => preloadSet.add(img));
rightPlayStim.forEach(img => preloadSet.add(img));
play_LResp_OutStim.forEach(img => preloadSet.add(img));
play_RResp_OutStim.forEach(img => preloadSet.add(img));


// Initialize list of preload images
console.log("Preloading task stim...");
preloadImages = [...preloadSet];
console.log("done");

// Initialize task parameters
var numBlocks = blockID.filter(onlyUnique).length;
var numBlockTrials = Math.max.apply(Math, trialID);
var numTotalTrials = designMat.length;

// Functions to display multiple stimuli on screen
function get2boxesHTMLstr(leftHTML, rightHTML, string) {
    let html = '<head><style> body {background-color: black;} </style></head>' +
        "<div id='mainContainer'>" +
        "    <div id='boxsContainer'>" +
        "        <div class='box'>" + leftHTML + "</div>" +
        "        <div class='box'>" + rightHTML + "</div>" +
        "    </div>" +
        "    <div id='instrContainer' class='centered whiteText largeFont'>" + string + "</div>" +
        "</div>";
    return html;
}

function get2tokensHTMLstr(leftHTML, rightHTML, string) {
    let html = '<head><style> body {background-color: black;} </style></head>' +
        "<div id='mainContainer'>" +
        "    <div id='boxsContainer'>" +
        "        <div class='token'>" + leftHTML + "</div>" +
        "        <div class='token'>" + rightHTML + "</div>" +
        "    </div>" +
        "    <div id='instrContainer' class='centered whiteText largeFont'>" + string + "</div>" +
        "</div>";
    return html;
}

function getOutcomeHTMLstr(middle, string) {
    let html = '<head><style> body {background-color: black;} </style></head>' +
        "<div id='mainContainer'>" +
        "    <div id='tokenContainer'>" + middle + "</div>" +
        "    <div id='outcContainer' class='centered whiteText largeFont'>" + string + "</div>" +
        "</div>";
    return html;
}
