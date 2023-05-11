/*
 * Main task script
 * Gets consent
 * Runs the task
 */

var check_consent = function (elem) {
  if (document.getElementById('consent_yes').checked == true && document.getElementById('consent_no').checked == false) {

      //run all the steps for the task
      timeline.push(go_full);
      timeline.push(instruct_trialProcedure);
      timeline.push(prac_trialProcedure); //experiential practice
      timeline.push(prac_trialProcedure_redo); //repeat instructions and practice if accuracy is too low
      timeline.push(instruct_trialProcedure_2); //second set of instructions for main task
      timeline.push(prac2_trialProcedure); //practice of main task (4 trials only, mostly for timing)
      timeline.push(task_trialProcedure);
      timeline.push(endTask_procedure);

      jsPsych.init({
          timeline: timeline,
          preload_images: unique(preloadImages),
          on_finish: function () {
              //jsPsych.data.displayData(); //can be commented out later on
              console.log('The experiment is over')
          },
      });
  }
  else {
    alert("Please tick Yes and click Continue to participate in the demo.");
    return false;
  }
};

function getQueryVariable(variable) {
    var query = window.location.search.substring(1);
    var vars = query.split("&");
    for (var i=0;i<vars.length;i++) {
        var pair = vars[i].split("=");
        if (pair[0] == variable) {
          return pair[1];
        }
    }
    return(false);
}

document.getElementById('header_title').innerHTML = "Welcome! Please click continue to proceed.";
document.getElementById('consent').innerHTML = "<p>\n" +
    "This demo will take you through the instructions, practice and 10 trials of the task.\n" +
    "</p>\n" +
    "<hr/>\n" +
    "<h4>Do you wish to continue?</h4>\n" +
    "\n" +
    "<label class=\"container\">Yes\n" +
        "<input type=\"checkbox\" id=\"consent_yes\">\n" +
    "</label>\n" +
    "\n" +
    "<label class=\"container\">No\n" +
        "<input type=\"checkbox\" id=\"consent_no\">\n" +
    "</label>\n" +
    "<p>\n" +
    "<h3> It will take a moment to load the next page. Please be patient.</h3>\n" +
    "</p>\n" +
    "<br><br>\n" +
    "<button type=\"button\" id=\"start\" class=\"submit_button\">continue</button>\n" +
    "<br><br>";

/*******************
 * Run Task
 ******************/
document.getElementById("start").onclick = check_consent;

if (/Android|webOS|iPhone|iPad|iPod|BlackBerry|IEMobile|Opera Mini/i.test(navigator.userAgent)) {
    alert("Sorry, this experiment does not work on mobile devices");
    document.getElementById('consent').innerHTML = "";
}
