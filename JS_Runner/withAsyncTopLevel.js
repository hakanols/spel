function setPRNG(fn) {
  let x = new Uint8Array(8);
  fn(x, 8);
  console.log(x)
};

function cleanup(arr) {
  for (var i = 0; i < arr.length; i++) arr[i] = 0;
}

(async function() {
  // Initialize PRNG if environment provides CSPRNG.
  // If not, methods calling randombytes will throw.
  if (typeof crypto !== 'undefined' && crypto.getRandomValues) {
    // Browsers.
    const QUOTA = 65536;
    setPRNG(function(x, n) {
      let i, v = new Uint8Array(n);
      for (i = 0; i < n; i += QUOTA) {
        crypto.getRandomValues(v.subarray(i, i + Math.min(n - i, QUOTA)));
      }
      for (i = 0; i < n; i++) x[i] = v[i];
      cleanup(v);
    });
  } else {
    // Node.js.
    const { randomBytes } = await import('crypto');
    if (randomBytes) {
      setPRNG(function(x, n) {
        let i, v = randomBytes(n);
        for (i = 0; i < n; i++) x[i] = v[i];
        cleanup(v);
      });
    }
  }
})();