// Function to round to a specified number of decimal places
function round(value, decimals) {
  return Number(Math.round(value+'e'+decimals)+'e-'+decimals);
}

// Funciton to round to the nearest multiple of 5
function round5(x) {
    return Math.round(x/5)*5;
}

// Function go draw from a uniform distribution of specified bounds, inclusive
function getRandomIntInclusive(min, max) {
  min = Math.ceil(min);
  max = Math.floor(max);
  return Math.floor(Math.random() * (max - min + 1)) + min; //The maximum is inclusive and the minimum is inclusive
}

// Function to create a random array of length 'length', drawn from a
// uniform distribution with min-max specified by 'params',
// filtered by a boolean vector. Depends on the 'round5' function above
function randomArray(min, max, length, outBool) {
  if (outBool.length != length) {
    throw "outBool length does match length parameter!";
  }
  var randArray = Array(length).fill(0);
  for (i = 0; i < length; i++) {
    if (outBool[i] == 0) {
        randArray[i] = "NA";
    } else {
        randArray[i] = round5(getRandomIntInclusive(min, max));
    }
  }
  return randArray
}

// Linspace function
function linspace(a,b,n) {
    if(typeof n === "undefined") n = Math.max(Math.round(b-a)+1,1);
    if(n<2) { return n===1?[a]:[]; }
    var i,ret = Array(n);
    n--;
    for(i=n;i>=0;i--) { ret[i] = (i*b+(n-i)*a)/n; }
    return ret;
}

// Range function, depends on the 'round' function above
// Inclusive of the stop parameter (unlike python range)
function arange(start, stop, step, decimals){
      step = step || 1;
      var arr = [];
      for (var i=start;i<stop+step;i+=step){
         arr.push(round(i,decimals));
      }
      return arr;
};

// Boolean coin toss functions
function coinFlip() {
    return (Math.floor(Math.random() * 2) == 0);
}

// Find index of matching string in arrays
function indexOfAll(array, searchItem) {
  var i = array.indexOf(searchItem),
      indexes = [];
  while (i !== -1) {
    indexes.push(i);
    i = array.indexOf(searchItem, ++i);
  }
  return indexes;
}

// Shuffle multiple arrays (of same length) using the same randomization
var isArray = Array.isArray || function(value) {
  return {}.toString.call(value) !== "[object Array]"
};
function jointShuffle() {
  var arrLength = 0;
  var argsLength = arguments.length;
  var rnd, tmp;
  for (var index = 0; index < argsLength; index += 1) {
    if (!isArray(arguments[index])) {
      throw new TypeError("Argument is not an array.");
    }
    if (index === 0) {
      arrLength = arguments[0].length;
    }
    if (arrLength !== arguments[index].length) {
      throw new RangeError("Array lengths do not match.");
    }
  }
  while (arrLength) {
    rnd = Math.floor(Math.random() * arrLength);
    arrLength -= 1;
    for (argsIndex = 0; argsIndex < argsLength; argsIndex += 1) {
      tmp = arguments[argsIndex][arrLength];
      arguments[argsIndex][arrLength] = arguments[argsIndex][rnd];
      arguments[argsIndex][rnd] = tmp;
    }
  }
}

// Compute factorial
function factorial(num) {
  if (num < 0)
        return -1;
  else if (num == 0)
      return 1;
  else {
      return (num * factorial(num - 1));
  }
}

// Get only unique values
function onlyUnique(value, index, self) {
    return self.indexOf(value) === index;
}


// Get sum of array
const arrSum = arr => arr.reduce((a,b) => a + b, 0);

// Tile an array with # repetitions
const makeRepeated = (arr, repeats) =>
  [].concat(...Array.from({ length: repeats }, () => arr));


// Get only unique values of array
function unique(array) {
    var filteredArray = array.filter(function(item, pos){
        return array.indexOf(item) == pos;
    });
    return filteredArray
};
