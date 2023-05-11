// load design
var prac2_numDesigns = designArray.practice_obs.length;
var prac2_designNo = getRandomIntInclusive(0, prac2_numDesigns - 1);
var prac2_designMat = designArray.practice_obs[prac2_designNo];

// These variables will change depending on the task
var prac2_trialID = Array(prac2_designMat.length).fill(0);
var prac2_leftObsStim = [];
var prac2_rightObsStim = [];
var prac2_leftObsChStim = [];
var prac2_rightObsChStim = [];
var prac2_obsTokStim = [];
var prac2_leftPlayStim = [];
var prac2_rightPlayStim = [];
var prac2_play_LResp_OutStim = [];
var prac2_play_RResp_OutStim = [];
var prac2_goalToken = Array(prac2_designMat.length).fill(0);
var prac2_tokenLeft = Array(prac2_designMat.length).fill(0);
var prac2_tokenRight = Array(prac2_designMat.length).fill(0);
var prac2_leftOutVal = Array(prac2_designMat.length).fill(0);
var prac2_rightOutVal = Array(prac2_designMat.length).fill(0);
var prac2_rewMag = Array(prac2_designMat.length).fill(0);
var prac2_corrResp = Array(prac2_designMat.length).fill(0);

// Convert from JSON object to set of arrays
for (t = 0; t < prac2_designMat.length; t++) {
    prac2_trialID[t] = Number(prac2_designMat[t].trNb);
    prac2_goalToken[t] = prac2_designMat[t].goalToken;
    prac2_tokenLeft[t] = prac2_designMat[t].tokenLeft;
    prac2_tokenRight[t] = prac2_designMat[t].tokenRight;
    prac2_leftOutVal[t] = Number(prac2_designMat[t].outcomeIfLeft);
    prac2_rightOutVal[t] = Number(prac2_designMat[t].outcomeIfRight);
    prac2_rewMag[t] = Number(prac2_designMat[t].rewMag);
    prac2_corrResp[t] =  prac2_designMat[t].corrResp;

    prac2_leftObsStim.push(stimDir + prac2_designMat[t].leftBoxObs + ".png");
    prac2_rightObsStim.push(stimDir + prac2_designMat[t].rightBoxObs + ".png");
    if ((Number(prac2_designMat[t].chBoxObs) == 1) && (prac2_designMat[t].posBoxObs1 == "left") || (Number(prac2_designMat[t].chBoxObs) == 2) && (prac2_designMat[t].posBoxObs1 == "right")) {
        //left box chosen by partner
        prac2_leftObsChStim.push(stimDir + prac2_designMat[t].leftBoxObs + "_ch.png");
        prac2_rightObsChStim.push(stimDir + prac2_designMat[t].rightBoxObs + ".png");
        prac2_obsTokStim.push(stimDir + prac2_designMat[t].leftBoxObs + "_" + prac2_designMat[t].tokObs + ".png");
    } else {
        //right box chosen by partner
        prac2_leftObsChStim.push(stimDir + prac2_designMat[t].leftBoxObs + ".png");
        prac2_rightObsChStim.push(stimDir + prac2_designMat[t].rightBoxObs + "_ch.png");
        prac2_obsTokStim.push(stimDir + prac2_designMat[t].rightBoxObs + "_" + prac2_designMat[t].tokObs + ".png");
    };

    prac2_leftPlayStim.push(stimDir + "token_" + prac2_designMat[t].tokenLeft + ".png");
    prac2_rightPlayStim.push(stimDir + "token_" + prac2_designMat[t].tokenRight + ".png");
    prac2_play_LResp_OutStim.push(stimDir + "token_" + prac2_designMat[t].tokenLeft + ".png");
    prac2_play_RResp_OutStim.push(stimDir + "token_" + prac2_designMat[t].tokenRight + ".png");

};

prac2_leftObsStim.forEach(img => preloadSet.add(img));
prac2_rightObsStim.forEach(img => preloadSet.add(img));
prac2_leftObsChStim.forEach(img => preloadSet.add(img));
prac2_rightObsChStim.forEach(img => preloadSet.add(img));
prac2_obsTokStim.forEach(img => preloadSet.add(img));
prac2_leftPlayStim.forEach(img => preloadSet.add(img));
prac2_rightPlayStim.forEach(img => preloadSet.add(img));
prac2_play_LResp_OutStim.forEach(img => preloadSet.add(img));
prac2_play_RResp_OutStim.forEach(img => preloadSet.add(img));

// List of preload images
console.log("Preloading obs practice stim...");
preloadImages = [...preloadSet];
console.log("done");

var numPract2Trials = prac2_designMat.length;
