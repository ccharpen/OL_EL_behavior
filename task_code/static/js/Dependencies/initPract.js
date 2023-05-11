// load design
var prac_numDesigns = designArray.practice.length;
var prac_designNo = getRandomIntInclusive(0, prac_numDesigns - 1);
var prac_designMat = designArray.practice[prac_designNo];

// These variables will change depending on the task
var prac_trialID = Array(prac_designMat.length).fill(0);
var prac_leftStim = [];
var prac_rightStim = [];
var prac_leftStim_LResp = [];
var prac_rightStim_LResp = [];
var prac_leftStim_RResp = [];
var prac_rightStim_RResp = [];
var prac_TokStim_LResp = [];
var prac_TokStim_RResp = [];
var prac_OutStim_LResp = [];
var prac_OutStim_RResp = [];
var prac_TokValStim = [];
var prac_goalToken = Array(prac_designMat.length).fill(0);
var prac_tokenIfLeft = Array(prac_designMat.length).fill(0);
var prac_tokenIfRight = Array(prac_designMat.length).fill(0);
var prac_leftOutVal = Array(prac_designMat.length).fill(0);
var prac_rightOutVal = Array(prac_designMat.length).fill(0);
var prac_rewMag = Array(prac_designMat.length).fill(0);
var prac_corrResp = Array(prac_designMat.length).fill(0);

// Convert from JSON object to set of arrays
for (t = 0; t < prac_designMat.length; t++) {
    prac_trialID[t] = Number(prac_designMat[t].trNb);
    prac_goalToken[t] = prac_designMat[t].goalToken;
    prac_tokenIfLeft[t] = prac_designMat[t].tokenIfLeft;
    prac_tokenIfRight[t] = prac_designMat[t].tokenIfRight;
    prac_leftOutVal[t] = Number(prac_designMat[t].outcomeIfLeft);
    prac_rightOutVal[t] = Number(prac_designMat[t].outcomeIfRight);
    prac_rewMag[t] = Number(prac_designMat[t].rewMag);
    prac_corrResp[t] =  prac_designMat[t].corrResp;
    prac_leftStim.push(stimDir + prac_designMat[t].leftBox + ".png");
    prac_rightStim.push(stimDir + prac_designMat[t].rightBox + ".png");
    prac_leftStim_LResp.push(stimDir + prac_designMat[t].leftBox + "_ch.png");
    prac_rightStim_LResp.push(stimDir + prac_designMat[t].rightBox + ".png");
    prac_leftStim_RResp.push(stimDir + prac_designMat[t].leftBox + ".png");
    prac_rightStim_RResp.push(stimDir + prac_designMat[t].rightBox + "_ch.png");
    prac_TokStim_LResp.push(stimDir + prac_designMat[t].leftBox + "_" + prac_designMat[t].tokenIfLeft + ".png");
    prac_TokStim_RResp.push(stimDir + prac_designMat[t].rightBox + "_" + prac_designMat[t].tokenIfRight + ".png");
    prac_OutStim_LResp.push(stimDir + "token_" + prac_designMat[t].tokenIfLeft + ".png");
    prac_OutStim_RResp.push(stimDir + "token_" + prac_designMat[t].tokenIfRight + ".png");
    switch (prac_designMat[t].rewPorange) {
        case "0.8":
            prac_TokValStim.push(stimDir + "tokval_o80_b20.png");
            break;
        case "0.6":
            prac_TokValStim.push(stimDir + "tokval_o60_b40.png");
            break;
        case "0.4":
            prac_TokValStim.push(stimDir + "tokval_o40_b60.png");
            break;
        case "0.2":
            prac_TokValStim.push(stimDir + "tokval_o20_b80.png");
            break;
    };
}

prac_leftStim.forEach(img => preloadSet.add(img));
prac_rightStim.forEach(img => preloadSet.add(img));
prac_leftStim_LResp.forEach(img => preloadSet.add(img));
prac_rightStim_LResp.forEach(img => preloadSet.add(img));
prac_leftStim_RResp.forEach(img => preloadSet.add(img));
prac_rightStim_RResp.forEach(img => preloadSet.add(img));
prac_TokStim_LResp.forEach(img => preloadSet.add(img));
prac_TokStim_RResp.forEach(img => preloadSet.add(img));
prac_OutStim_LResp.forEach(img => preloadSet.add(img));
prac_OutStim_RResp.forEach(img => preloadSet.add(img));
prac_TokValStim.forEach(img => preloadSet.add(img));

// Initialize list of preload images
console.log("Preloading practice stim...");
preloadImages = [...preloadSet];
console.log("done");

var numPractTrials = prac_designMat.length;
