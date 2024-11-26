import test from './tap-esm.js'

test('Console log', async function (t) {
    console.log("Log message");
    console.debug("Debug message");
    console.warn("Warning message");
    console.error("Error message");
	t.end();
});

test('Syntax error', async function (t) {
    let empty = null;
	t.arrayEqual(empty[0], null, "Array index 0 of null, should crash");
	t.end();
});