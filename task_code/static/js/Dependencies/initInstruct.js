// Preload instruction slides
var instruct_pages = []; //for first set of instructions
var Ninstruct = 10;
var instruct_pages_2 = []; //for second set of instructions
var Ninstruct_2 = 17;

let prev_next = [choice_keys.space, choice_keys.b, null, null]; //default keys to move through instructions

// pages from which key(s) other than space/"b" are used for flipping slides
//for first set of instructions
var exception_pages = {
    1: [choice_keys.space, null, null, null],
    4: [choice_keys.left, choice_keys.b, choice_keys.space, null],
};
var choiceKey_array = [];   // array holding valid keys for each instruction page
for (s = 0; s < Ninstruct; s++) {
    let slideNo = s + 1 + "";
    // add options from exception_pages for instructions with questions
    if (! Object.keys(exception_pages).includes(slideNo)) choiceKey_array.push(prev_next);
    else choiceKey_array.push(exception_pages[slideNo]);
    if (slideNo < 10) {slideNo = "0"+slideNo}   // padding zero
    instruct_pages.push(instructDir + 'instr_' + slideNo + '.png');
}

//for second set of instructions
var exception_pages_2 = {
    9: [choice_keys.h, choice_keys.b, choice_keys.f, choice_keys.k],
    11: [choice_keys.f, choice_keys.b, choice_keys.h, choice_keys.k],
    13: [choice_keys.k, choice_keys.b, choice_keys.f, choice_keys.h],
    15: [choice_keys.f, choice_keys.b, choice_keys.h, choice_keys.k],
};
var choiceKey_array_2 = [];   // array holding valid keys for each instruction page
for (s = 0; s < Ninstruct_2; s++) {
    let slideNo = s + 1 + "";
    // add options from exception_pages for instructions with questions
    if (! Object.keys(exception_pages_2).includes(slideNo)) choiceKey_array_2.push(prev_next);
    else choiceKey_array_2.push(exception_pages_2[slideNo]);
    if (slideNo < 10) {slideNo = "0"+slideNo}   // padding zero
    instruct_pages_2.push(instructDir + 'instr2_' + slideNo + '.png');
}

console.log("Preloading instruction pages...");
preloadImages.push(...instruct_pages);
preloadImages.push(...instruct_pages_2);
console.log("done");
