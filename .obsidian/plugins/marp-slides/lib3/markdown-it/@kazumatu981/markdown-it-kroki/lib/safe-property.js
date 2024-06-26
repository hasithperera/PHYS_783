'use strict';

function safeProperty(test, name, type, defaultValue) {
    if (test == null || test == undefined) return defaultValue;
    if (typeof test[name] !== type) return defaultValue;
    if (typeof test[name] === "string" && test[name] === '') return defaultValue;
    return test[name];
}

function safeChoice(test, candidates, defaultValue) {
    return candidates.includes(test) ?
        test : defaultValue;
}

module.exports = { safeProperty, safeChoice };
