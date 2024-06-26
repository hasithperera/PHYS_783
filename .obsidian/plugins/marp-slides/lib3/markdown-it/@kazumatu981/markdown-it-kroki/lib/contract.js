'use strict';

module.exports = {
    /**
     * contract `test` to be non-empty string.
     * @param {string} test test string
     * @param {string} msg message on exception
     */
    toNonEmptyString: function (test, msg) {
        if (typeof test !== 'string') throw new Error(msg);
        if (test === ''
            || test === null
            || test === undefined) throw new Error(msg);
    },
    /**
     * contract `test` to be true.
     * @param {boolean} test test boolean. 
     * @param {sting} msg massage on excetion.
     */
    toTrue: function (test, msg) {
        if (typeof test !== 'boolean') throw new Error(msg);
        if (!test) throw new Error(msg);
    },
    toBeUrlString: function (test, msg) {
        this.toNonEmptyString(test, msg);
        try {
            new URL(test);
        } catch {
            throw new Error(msg);
        }
    },
    toBeClassName: function (test, msg) {
        if (!/^[A-Za-z0-9]+(-?[A-Za-z0-9]+)*$/.exec(test)) {
            throw new Error(msg);
        }
    }
};
