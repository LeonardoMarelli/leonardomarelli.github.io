// HMWK 3 - Caesar cipher, toy elliptic-curve encryption, and key-space attacks
// Plain JavaScript ES module. The modern-scheme demo uses WebCrypto (P-256, HKDF, AES-256-GCM),
// which needs a secure context: an https:// page or localhost.

/* =========================================================================
   1. Reference data
   ========================================================================= */

export const SAMPLE_EN =
  'Cryptography has a long history. For centuries, military commanders and diplomats relied on secret ' +
  'writing to protect their plans, and most early methods were simple substitutions. The Caesar cipher is ' +
  'one of the oldest and best known of them. Each letter is moved a fixed number of places along the alphabet, ' +
  'so a message that reads attack at dawn becomes something else entirely. Anyone who knows the shift can read ' +
  'it at once, yet anyone who only sees the output can try the few possible shifts and check which one gives ' +
  'text that looks like the language. Letter frequencies make this check easy, because in ordinary English the ' +
  'letter E appears far more often than Z, and that pattern survives the shift unchanged.';

export const SAMPLE_IT =
  'La crittografia ha una storia lunga. Per secoli i comandanti militari e i diplomatici si sono affidati ' +
  'alla scrittura segreta per proteggere i propri piani, e i metodi più antichi erano semplici sostituzioni. ' +
  'Il cifrario di Cesare è uno dei più antichi e noti. Ogni lettera viene spostata di un numero fisso di ' +
  'posizioni nell\'alfabeto, quindi un messaggio qualsiasi diventa qualcosa di completamente diverso. Chi conosce ' +
  'lo spostamento può leggerlo subito, ma chi vede solo il risultato può provare i pochi spostamenti possibili ' +
  'e verificare quale produce un testo simile alla lingua. Le frequenze delle lettere rendono questo controllo ' +
  'facile, perché nell\'italiano comune la E compare molto più spesso della Z, e questa regolarità resta ' +
  'invariata dopo lo spostamento.';

// Reference letter frequencies in percent, A to Z.
// English: Oxford University Mathematical Institute table.
export const ENGLISH_FREQ = [
  8.167, 1.492, 2.782, 4.253, 12.702, 2.228, 2.015, 6.094, 6.966, 0.153, 0.772, 4.025, 2.406,
  6.749, 7.507, 1.929, 0.095, 5.987, 6.327, 9.056, 2.758, 0.978, 2.360, 0.150, 1.974, 0.074,
];
// Italian: sttmedia.com table. J, K, W, X, Y are not Italian letters; they get 0.05 instead of 0
// so that the chi-square statistic stays defined.
export const ITALIAN_FREQ = [
  10.85, 1.05, 4.30, 3.39, 11.49, 1.01, 1.65, 1.43, 10.18, 0.05, 0.05, 5.70, 2.87,
  7.02, 9.97, 2.96, 0.45, 6.19, 5.48, 6.97, 3.16, 1.75, 0.05, 0.05, 0.05, 0.85,
];

export const CRIT_BINS_DF15 = 25.0;      // chi-square, 16 equal bins, 15 df, 5% level
export const CRIT_LETTERS_DF25 = 37.65;  // chi-square, 26 letters, 25 df, 5% level
export const CRIT_BYTES_DF255 = 293.3;   // chi-square, 256 byte values, 255 df, 5% level

/* =========================================================================
   2. Caesar cipher: encryption, decryption, statistical attack
   ========================================================================= */

// Shift every letter by k places (k may be negative). Non-letters are left unchanged.
export function shiftText(text, k) {
  const s = ((k % 26) + 26) % 26;
  return text.replace(/[A-Za-z]/g, ch => {
    const base = ch <= 'Z' ? 65 : 97;
    return String.fromCharCode(((ch.charCodeAt(0) - base + s) % 26) + base);
  });
}

// C = (P + k) mod 26 and P = (C - k) mod 26, applied to letters only.
export const caesarEncrypt = (text, k) => shiftText(text, k);
export const caesarDecrypt = (text, k) => shiftText(text, -k);

// Keep only A-Z, upper case. Used by the statistical tests.
export function lettersOnly(text) {
  return text.toUpperCase().replace(/[^A-Z]/g, '');
}

export function letterCounts(text) {
  const counts = Array(26).fill(0);
  for (const ch of text.toUpperCase()) {
    const i = ch.charCodeAt(0) - 65;
    if (i >= 0 && i < 26) counts[i]++;
  }
  return counts;
}

// Chi-square: sum over letters of (observed - expected)^2 / expected. Lower = closer to the language.
export function chiSquareLetters(counts, freq) {
  const total = counts.reduce((a, b) => a + b, 0);
  if (total === 0) return Infinity;
  const sumFreq = freq.reduce((a, b) => a + b, 0);
  let chi = 0;
  for (let i = 0; i < 26; i++) {
    const expected = (total * freq[i]) / sumFreq;
    chi += (counts[i] - expected) ** 2 / expected;
  }
  return chi;
}

// Unigram log-likelihood of the counts under the reference distribution. Higher = more plausible.
export function logLikelihoodLetters(counts, freq) {
  const sumFreq = freq.reduce((a, b) => a + b, 0);
  let ll = 0;
  for (let i = 0; i < 26; i++) ll += counts[i] * Math.log(freq[i] / sumFreq);
  return ll;
}

// Exhaustive search over the 26 keys. Candidate k decrypts by shifting back k places.
// Sorted by chi-square, lowest (most plausible) first.
export function rankShifts(cipher, freq) {
  const rows = [];
  for (let k = 0; k < 26; k++) {
    const plain = caesarDecrypt(cipher, k);
    const counts = letterCounts(plain);
    rows.push({
      key: k,
      chi: chiSquareLetters(counts, freq),
      ll: logLikelihoodLetters(counts, freq),
      plain,
    });
  }
  return rows.sort((a, b) => a.chi - b.chi);
}

// Best candidate according to the log-likelihood score (instead of chi-square).
export const bestByLL = ranked => ranked.reduce((a, b) => (b.ll > a.ll ? b : a));

// Small deterministic pseudo-random generator, so that experiments can be repeated.
export function prng(seed) {
  let s = seed >>> 0;
  return () => {
    s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
    return s / 4294967296;
  };
}

// Attack experiment on texts of different lengths.
// For each length L: take `trials` windows of L letters from the corpus (random offsets),
// encrypt each with a random key 1..25, and check whether the attack recovers the key.
export function lengthExperiment(corpus, freq, lengths, trials, seed = 2024) {
  const text = lettersOnly(corpus);
  const rnd = prng(seed);
  return lengths.map(L => {
    if (L > text.length) throw new Error(`corpus has only ${text.length} letters`);
    let okChi = 0, okLL = 0, rankSum = 0;
    for (let t = 0; t < trials; t++) {
      const off = Math.floor(rnd() * (text.length - L + 1));
      const key = 1 + Math.floor(rnd() * 25);
      const cipher = caesarEncrypt(text.slice(off, off + L), key);
      const ranked = rankShifts(cipher, freq);
      if (ranked[0].key === key) okChi++;
      if (bestByLL(ranked).key === key) okLL++;
      rankSum += ranked.findIndex(r => r.key === key) + 1;
    }
    return {
      L,
      successChi: okChi / trials,
      successLL: okLL / trials,
      meanRank: rankSum / trials,
    };
  });
}

/* =========================================================================
   3. Toy elliptic-curve scheme (same equation shape as secp256k1: a = 0, b = 7)
   ElGamal-style encryption over a small prime field.
   ========================================================================= */

export const TOY_CURVE = { p: 1019, a: 0, b: 7 };  // y^2 = x^3 + 7 (mod 1019)
export const KEY_MAX = 1000;                       // restricted private-key range: d = 1 .. KEY_MAX
export const SYMBOLS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ '; // 27 symbols, space included

const INF = null;                                  // point at infinity O
const keyOf = P => (P === INF ? 'O' : `${P.x},${P.y}`);
const mod = (n, p) => ((n % p) + p) % p;

export function modInv(a, p) {
  let r0 = mod(a, p), r1 = p;
  let s0 = 1, s1 = 0;
  while (r1 !== 0) {
    const q = Math.floor(r0 / r1);
    [r0, r1] = [r1, r0 - q * r1];
    [s0, s1] = [s1, s0 - q * s1];
  }
  if (r0 !== 1) throw new Error(`${a} has no inverse mod ${p}`);
  return mod(s0, p);
}

export function pointAdd(P, Q, c) {
  if (P === INF) return Q;
  if (Q === INF) return P;
  if (P.x === Q.x && mod(P.y + Q.y, c.p) === 0) return INF;
  let lambda;
  if (P.x === Q.x) {
    lambda = mod((3 * P.x * P.x + c.a) * modInv(2 * P.y, c.p), c.p);
  } else {
    lambda = mod((Q.y - P.y) * modInv(Q.x - P.x, c.p), c.p);
  }
  const x = mod(lambda * lambda - P.x - Q.x, c.p);
  const y = mod(lambda * (P.x - x) - P.y, c.p);
  return { x, y };
}

export const pointNeg = (P, c) => (P === INF ? INF : { x: P.x, y: mod(-P.y, c.p) });

// Double-and-add: kP in O(log k) group operations.
export function scalarMult(k, P, c) {
  let result = INF;
  let addend = P;
  let n = k;
  while (n > 0) {
    if (n % 2 === 1) result = pointAdd(result, addend, c);
    addend = pointAdd(addend, addend, c);
    n = Math.floor(n / 2);
  }
  return result;
}

export function orderOf(P, c) {
  let Q = P;
  let n = 1;
  while (Q !== INF) {
    Q = pointAdd(Q, P, c);
    n++;
  }
  return n;
}

// Base point: among the first affine points, take the one of largest order (as in hw1.js).
export function findBasePoint(c, candidates = 40) {
  const pts = [];
  for (let x = 0; x < c.p && pts.length < candidates; x++) {
    const rhs = mod(x ** 3 + c.a * x + c.b, c.p);
    for (let y = 0; y < c.p; y++) {
      if (mod(y * y, c.p) === rhs) pts.push({ x, y });
    }
  }
  return pts.reduce((best, P) => (orderOf(P, c) > orderOf(best, c) ? P : best), pts[0]);
}

// Public parameters of the toy scheme.
export function setupToy() {
  const c = TOY_CURVE;
  const G = findBasePoint(c);
  const n = orderOf(G, c);
  if (n < SYMBOLS.length + 1) throw new Error('base point order too small for the message alphabet');
  return { c, G, n, keyMax: Math.min(KEY_MAX, n - 1) };
}

// Message encoding: symbol i is mapped to the point M_i = (i+1)G.
// The table allows the decoder to map a recovered point back to its symbol.
export function symbolTable(G, c) {
  const points = [];
  const dec = new Map();
  let P = G;
  for (let i = 0; i < SYMBOLS.length; i++) {
    points.push(P);
    dec.set(keyOf(P), i);
    P = pointAdd(P, G, c);
  }
  return { points, dec };
}

// Key generation: private d in 1..keyMax, public Q = dG.
export function toyKeyGen(toy, rnd = Math.random) {
  const d = 1 + Math.floor(rnd() * toy.keyMax);
  const Q = scalarMult(d, toy.G, toy.c);
  return { d, Q };
}

// Encryption: for each symbol, pick a fresh r, send C1 = rG and C2 = M + rQ.
// The returned `r` is for the demo table only; it is not part of the ciphertext.
export function toyEncrypt(text, Q, toy, table, rnd = Math.random) {
  const blocks = [];
  for (const ch of text.toUpperCase()) {
    const i = SYMBOLS.indexOf(ch);
    if (i < 0) continue;                       // characters outside the alphabet are dropped
    const r = 1 + Math.floor(rnd() * toy.keyMax);
    const C1 = scalarMult(r, toy.G, toy.c);
    const S = scalarMult(r, Q, toy.c);         // shared point rQ
    const C2 = pointAdd(table.points[i], S, toy.c);
    blocks.push({ ch, r, C1, C2 });
  }
  return blocks;
}

// Decryption: S = d*C1 = rdG = r*Q, then M = C2 - S, then look M up in the symbol table.
export function toyDecrypt(blocks, d, toy, table) {
  return blocks.map(b => {
    const S = scalarMult(d, b.C1, toy.c);
    const M = pointAdd(b.C2, pointNeg(S, toy.c), toy.c);
    const i = table.dec.get(keyOf(M));
    return i === undefined ? '?' : SYMBOLS[i];
  }).join('');
}

/* ---------- Attacks on the restricted key space ---------- */

// Exhaustive search: compute 1G, 2G, ... until kG = Q. Uses only the public key Q, G and the curve.
export function attackExhaustive(Q, toy) {
  const t0 = performance.now();
  let P = INF;
  let steps = 0;
  for (let k = 1; k <= toy.keyMax; k++) {
    P = pointAdd(P, toy.G, toy.c);
    steps++;
    if (keyOf(P) === keyOf(Q)) {
      return { d: k, steps, ms: performance.now() - t0 };
    }
  }
  return { d: null, steps, ms: performance.now() - t0 };
}

// Baby-step giant-step: with m = ceil(sqrt(keyMax+1)), write d = i*m + j.
// Store jG for j = 0..m-1 in a table, then walk Q - i*mG and look for a match. About 2*sqrt(range) steps.
export function attackBSGS(Q, toy) {
  const t0 = performance.now();
  const c = toy.c, G = toy.G;
  const m = Math.ceil(Math.sqrt(toy.keyMax + 1));
  const baby = new Map();
  let P = INF;
  for (let j = 0; j < m; j++) {
    baby.set(keyOf(P), j);
    P = pointAdd(P, G, c);
  }
  const negMG = pointNeg(scalarMult(m, G, c), c);
  let T = Q;
  for (let i = 0; i <= m; i++) {
    if (baby.has(keyOf(T))) {
      const d = i * m + baby.get(keyOf(T));
      if (d >= 1 && d <= toy.keyMax && keyOf(scalarMult(d, G, c)) === keyOf(Q)) {
        return { d, steps: m + i + 1, ms: performance.now() - t0 };
      }
    }
    T = pointAdd(T, negMG, c);
  }
  return { d: null, steps: m + m + 1, ms: performance.now() - t0 };
}

// Feasibility of the same attacks for a key space of `bits` bits (ops/s is an assumption).
export function feasibilityRows(opsPerSecond = 1e12) {
  const secondsPerYear = 3.156e7;
  return [10, 20, 32, 64, 128, 256].map(bits => {
    const exhaustive = 2 ** bits;
    const bsgs = 2 ** (bits / 2);
    return {
      bits,
      exhaustive,
      bsgs,
      yearsBSGS: bsgs / opsPerSecond / secondsPerYear,
    };
  });
}

/* =========================================================================
   4. Modern scheme for comparison: ECIES-style with WebCrypto
   ECDH on P-256, HKDF-SHA-256, AES-256-GCM.
   ========================================================================= */

const enc = new TextEncoder();
const dec = new TextDecoder();
const ECDH_PARAMS = { name: 'ECDH', namedCurve: 'P-256' };
const KDF_INFO = enc.encode('hw3-ecies-aes-256-gcm');

async function deriveAesKey(sharedBits) {
  const base = await crypto.subtle.importKey('raw', sharedBits, 'HKDF', false, ['deriveKey']);
  return crypto.subtle.deriveKey(
    { name: 'HKDF', hash: 'SHA-256', salt: new Uint8Array(0), info: KDF_INFO },
    base,
    { name: 'AES-GCM', length: 256 },
    false,
    ['encrypt', 'decrypt'],
  );
}

export function generateRecipientKeys() {
  return crypto.subtle.generateKey(ECDH_PARAMS, true, ['deriveBits']);
}

export async function eccEncrypt(message, recipientPublic) {
  const ephemeral = await crypto.subtle.generateKey(ECDH_PARAMS, false, ['deriveBits']);
  const shared = await crypto.subtle.deriveBits(
    { name: 'ECDH', public: recipientPublic }, ephemeral.privateKey, 256);
  const aesKey = await deriveAesKey(shared);
  const iv = crypto.getRandomValues(new Uint8Array(12));
  const ct = new Uint8Array(await crypto.subtle.encrypt({ name: 'AES-GCM', iv }, aesKey, enc.encode(message)));
  const ephPub = new Uint8Array(await crypto.subtle.exportKey('raw', ephemeral.publicKey));
  return { ephPub, iv, ct };
}

export async function eccDecrypt(packet, recipientPrivate) {
  const ephPublic = await crypto.subtle.importKey('raw', packet.ephPub, ECDH_PARAMS, false, []);
  const shared = await crypto.subtle.deriveBits({ name: 'ECDH', public: ephPublic }, recipientPrivate, 256);
  const aesKey = await deriveAesKey(shared);
  const pt = await crypto.subtle.decrypt({ name: 'AES-GCM', iv: packet.iv }, aesKey, packet.ct);
  return dec.decode(pt);
}

/* ---------- byte statistics ---------- */

export function byteCounts(bytes) {
  const counts = Array(256).fill(0);
  for (const b of bytes) counts[b]++;
  return counts;
}

export function chiSquareUniform(counts) {
  const total = counts.reduce((a, b) => a + b, 0);
  const expected = total / counts.length;
  return counts.reduce((s, c) => s + (c - expected) ** 2 / expected, 0);
}

/* =========================================================================
   5. Page (runs only in the browser)
   ========================================================================= */

/* ---------- direct comparison: Caesar vs toy scheme on the same letters ---------- */

// Attacker's view of the toy ciphertext: the x-coordinate of each C2, reduced mod 26 and read as a letter.
export function toyCiphertextLetters(blocks) {
  return blocks.map(b => String.fromCharCode(65 + (b.C2.x % 26))).join('');
}

// Chi-square of values split into equal-width bins over [0, range), against uniform.
export function chiSquareBins(values, nbins, range) {
  const counts = Array(nbins).fill(0);
  for (const v of values) counts[Math.floor((v * nbins) / range)]++;
  const expected = values.length / nbins;
  return counts.reduce((s, c) => s + (c - expected) ** 2 / expected, 0);
}

// Statistics of the Caesar ciphertext of the letters of `text` (key k).
export function caesarStatistics(text, k) {
  const letters = lettersOnly(text);
  const cipher = caesarEncrypt(letters, k);
  const ranked = rankShifts(cipher, ENGLISH_FREQ);
  return {
    n: letters.length,
    unshiftedChi: chiSquareLetters(letterCounts(cipher), ENGLISH_FREQ),
    attackChi: ranked[0].chi,
    attackKey: ranked[0].key,
  };
}

// Statistics of the toy ciphertext of the same letters: a fresh r for every letter.
export function toyStatistics(text, Q, toy, table, seed = 99) {
  const letters = lettersOnly(text);
  const blocks = toyEncrypt(letters, Q, toy, table, prng(seed));
  return {
    n: blocks.length,
    letterChi: rankShifts(toyCiphertextLetters(blocks), ENGLISH_FREQ)[0].chi,
    xChi: chiSquareBins(blocks.map(b => b.C2.x), 16, toy.c.p),
    yChi: chiSquareBins(blocks.map(b => b.C2.y), 16, toy.c.p),
  };
}

/* ---------- page helpers for the new sections ---------- */

function renderToyStatistics() {
  const cs = caesarStatistics(SAMPLE_EN, 7);
  const ts = toyStatistics(SAMPLE_EN, toyKeys.Q, toy, toyTable);
  const wrap = $('toy-stats');
  if (!wrap) return;
  wrap.textContent = '';
  wrap.appendChild(table(['Ciphertext of the same letters', 'Test', 'Value', 'Critical (5%)', 'Verdict'], [
    ['Caesar (key 7)', 'chi-square of the ciphertext, not decrypted', cs.unshiftedChi.toFixed(1),
      CRIT_LETTERS_DF25, verdict(cs.unshiftedChi < CRIT_LETTERS_DF25)],
    ['Caesar (key 7)', `chi-square after best key (recovered key ${cs.attackKey})`, cs.attackChi.toFixed(1),
      CRIT_LETTERS_DF25, verdict(cs.attackChi < CRIT_LETTERS_DF25)],
    ['Toy scheme', 'x-coordinate of C2 mod 26, best shift', ts.letterChi.toFixed(1),
      CRIT_LETTERS_DF25, verdict(ts.letterChi < CRIT_LETTERS_DF25)],
    ['Toy scheme', 'x-coordinates of C2 in 16 bins vs uniform', ts.xChi.toFixed(1),
      CRIT_BINS_DF15, ts.xChi < CRIT_BINS_DF15 ? 'consistent with uniform' : 'not uniform'],
    ['Toy scheme', 'y-coordinates of C2 in 16 bins vs uniform', ts.yChi.toFixed(1),
      CRIT_BINS_DF15, ts.yChi < CRIT_BINS_DF15 ? 'consistent with uniform' : 'not uniform'],
  ]));
  wrap.appendChild(el('p',
    `Both ciphertexts encrypt the same ${cs.n} letters. The toy scheme produced ${ts.n} blocks, one per letter. ` +
    'The coordinate values are computed with a fixed seed, so the numbers are the same on each page load.'));
}

function renderItalianLengths() {
  const lengths = [10, 20, 50, 100, 200, 300, 400];
  const results = lengthExperiment(SAMPLE_IT, ITALIAN_FREQ, lengths, 300);
  const rows = results.map(r => [
    r.L,
    300,
    pctOf(r.successChi * 300, 300) + '%',
    pctOf(r.successLL * 300, 300) + '%',
    r.meanRank.toFixed(2),
  ]);
  const out = $('caesar-length-it');
  if (!out) return;
  out.textContent = '';
  out.appendChild(table(
    ['Ciphertext letters', 'Trials', 'Key found (chi-square)', 'Key found (log-likelihood)', 'Mean rank of true key'],
    rows));
}

const $ = id => document.getElementById(id);
const LETTERS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

function el(tag, text, cls) {
  const node = document.createElement(tag);
  if (text !== undefined) node.textContent = text;
  if (cls) node.className = cls;
  return node;
}

function table(headers, rows, highlight = -1) {
  const t = el('table');
  const head = el('tr');
  headers.forEach(h => head.appendChild(el('th', h)));
  t.appendChild(head);
  rows.forEach((r, i) => {
    const tr = el('tr');
    if (i === highlight) tr.style.fontWeight = 'bold';
    r.forEach(cell => tr.appendChild(el('td', String(cell))));
    t.appendChild(tr);
  });
  return t;
}

const pctOf = (count, total) => (total ? ((100 * count) / total).toFixed(1) : '0.0');
const toBase64 = bytes => btoa(Array.from(bytes, b => String.fromCharCode(b)).join(''));
const toHex = (bytes, n) => Array.from(bytes.slice(0, n), b => b.toString(16).padStart(2, '0')).join('');
const fmtPoint = P => (P === INF ? 'O' : `(${P.x}, ${P.y})`);
const verdict = ok => (ok ? 'looks like the language' : 'does not look like the language');

function guard(fn) {
  try {
    fn();
  } catch (err) {
    console.error(err);
    const box = $('toy-setup');
    if (box) box.textContent = `Error: ${err.message}`;
  }
}

// ---- Caesar: reference table, encryption and decryption ----

function renderReference() {
  const wrap = $('reference');
  if (!wrap) return;
  wrap.textContent = '';
  const rows = [...LETTERS].map((L, i) => [L, ENGLISH_FREQ[i].toFixed(2), ITALIAN_FREQ[i].toFixed(2)]);
  wrap.appendChild(table(['Letter', 'English %', 'Italian %'], rows));
}

function runCaesarEncrypt() {
  const k = Number($('caesar-key').value);
  if (!Number.isInteger(k) || k < 1 || k > 25) {
    $('caesar-cipher-out').textContent = 'Please enter an integer key from 1 to 25.';
    return;
  }
  const cipher = caesarEncrypt($('caesar-plain').value, k);
  $('caesar-cipher-out').textContent = cipher;
  $('caesar-decrypt-in').value = cipher;
  $('caesar-cipher-in').value = cipher;
}

function runCaesarDecrypt() {
  const k = Number($('caesar-key').value);
  const plain = caesarDecrypt($('caesar-decrypt-in').value, k);
  $('caesar-decrypt-out').textContent = plain;
  $('caesar-roundtrip').textContent =
    `Same key, decrypted text equals the original plaintext: ${plain === $('caesar-plain').value}`;
}

// ---- Caesar: statistical attack ----

function renderCaesarAttack() {
  const cipher = $('caesar-cipher-in').value;
  const freq = $('caesar-lang').value === 'italian' ? ITALIAN_FREQ : ENGLISH_FREQ;
  const ranked = rankShifts(cipher, freq);
  const bestChi = ranked[0];
  const bestLL = bestByLL(ranked);

  const out = $('caesar-result');
  out.textContent = '';
  out.appendChild(el('p',
    `Best key by chi-square: ${bestChi.key} (chi-square ${bestChi.chi.toFixed(1)}; 5% critical value ${CRIT_LETTERS_DF25}). ` +
    `Best key by log-likelihood: ${bestLL.key}.`));
  out.appendChild(el('pre', bestChi.plain, 'output'));

  const cc = letterCounts(cipher), dc = letterCounts(bestChi.plain);
  const cTotal = cc.reduce((a, b) => a + b, 0), dTotal = dc.reduce((a, b) => a + b, 0);
  const sumRef = freq.reduce((a, b) => a + b, 0);
  const rows = [...LETTERS].map((L, i) => [
    L,
    ((100 * freq[i]) / sumRef).toFixed(2),
    pctOf(cc[i], cTotal),
    pctOf(dc[i], dTotal),
  ]);
  out.appendChild(el('h3', 'Letter frequencies: reference, ciphertext and decoded text (%)'));
  out.appendChild(table(['Letter', 'Reference', 'Ciphertext', 'Decoded with best key'], rows));

  out.appendChild(el('h3', 'All 26 keys, ranked by chi-square'));
  const rankRows = ranked.map(r => [r.key, r.chi.toFixed(1), r.ll.toFixed(1), r.plain.slice(0, 50) + '...']);
  out.appendChild(table(['Key', 'Chi-square', 'Log-likelihood', 'Decoded start'], rankRows, 0));
  return bestChi;
}

// ---- Caesar: effect of text length ----

function renderLengthExperiment() {
  const lengths = [10, 20, 50, 100, 200, 300, 400];
  const results = lengthExperiment(SAMPLE_EN, ENGLISH_FREQ, lengths, 300);
  const rows = results.map(r => [
    r.L,
    300,
    pctOf(r.successChi * 300, 300) + '%',
    pctOf(r.successLL * 300, 300) + '%',
    r.meanRank.toFixed(2),
  ]);
  const out = $('caesar-length-out');
  out.textContent = '';
  out.appendChild(table(
    ['Ciphertext letters', 'Trials', 'Key found (chi-square)', 'Key found (log-likelihood)', 'Mean rank of true key'],
    rows));
  out.appendChild(el('p',
    'Each trial uses a random window of the English sample and a random key 1..25. ' +
    'Windows for long lengths overlap heavily, because the sample has only about 600 letters.'));
}

// ---- Comparison of Caesar ciphertext with the language ----

function renderCompare(best) {
  const cipher = $('caesar-cipher-in').value;
  const freq = $('caesar-lang').value === 'italian' ? ITALIAN_FREQ : ENGLISH_FREQ;
  const raw = chiSquareLetters(letterCounts(cipher), freq);
  const wrap = $('compare');
  if (!wrap) return;
  wrap.textContent = '';
  wrap.appendChild(table(['Text', 'Chi-square vs language', 'Verdict (5%)'], [
    ['Caesar ciphertext, not decrypted', raw.toFixed(1), verdict(raw < CRIT_LETTERS_DF25)],
    ['Best candidate after decryption', best.chi.toFixed(1), verdict(best.chi < CRIT_LETTERS_DF25)],
  ]));
  wrap.appendChild(el('p', 'The ECIES ciphertext statistics are in the modern-scheme section above.'));
}

// ---- Toy elliptic-curve scheme ----

let toy = null, toyTable = null, toyKeys = null, toyBlocks = null;

function renderToySetup() {
  toy = setupToy();
  toyTable = symbolTable(toy.G, toy.c);
  toyKeys = toyKeyGen(toy);
  const wrap = $('toy-setup');
  wrap.textContent = '';
  wrap.appendChild(table(['Parameter', 'Value'], [
    ['Curve', `y^2 = x^3 + ${toy.c.a}x + ${toy.c.b}  (mod ${toy.c.p})`],
    ['Base point G', fmtPoint(toy.G)],
    ['Order of G', toy.n],
    ['Private key range', `d = 1 .. ${toy.keyMax}`],
    ['Private key d (owner only, shown for the demo)', toyKeys.d],
    ['Public key Q = dG', fmtPoint(toyKeys.Q)],
  ]));
}

function runToyEncrypt() {
  toyBlocks = toyEncrypt($('toy-plain').value, toyKeys.Q, toy, toyTable);
  const wrap = $('toy-cipher');
  wrap.textContent = '';
  wrap.appendChild(table(
    ['Block', 'C1 = rG', 'C2 = M + rQ'],
    toyBlocks.map((b, i) => [i + 1, fmtPoint(b.C1), fmtPoint(b.C2)])));
  $('toy-decrypt-out').textContent =
    `Decrypted with the private key d: ${toyDecrypt(toyBlocks, toyKeys.d, toy, toyTable)}`;
}

// The attacker sees only Q, G, the curve and the blocks (C1, C2). Not d, not r.
function runToyAttack() {
  if (!toyBlocks) return;
  const out = $('toy-attack-out');
  out.textContent = '';
  const ex = attackExhaustive(toyKeys.Q, toy);
  const bs = attackBSGS(toyKeys.Q, toy);
  out.appendChild(table(['Attack', 'Recovered d', 'Group operations', 'Time'], [
    ['Exhaustive search', ex.d === null ? 'not found' : ex.d, ex.steps, ex.ms.toFixed(2) + ' ms'],
    ['Baby-step giant-step', bs.d === null ? 'not found' : bs.d, bs.steps, bs.ms.toFixed(2) + ' ms'],
  ]));
  const d = ex.d !== null ? ex.d : bs.d;
  if (d !== null) {
    out.appendChild(el('p',
      `Recovered d = ${d}. Decrypting the intercepted blocks with it gives: ${toyDecrypt(toyBlocks, d, toy, toyTable)}`));
  }
}

function renderFeasibility() {
  const rows = feasibilityRows().map(r => [
    r.bits,
    r.exhaustive.toExponential(2),
    r.bsgs.toExponential(2),
    r.yearsBSGS.toExponential(2),
  ]);
  const wrap = $('feasibility');
  wrap.textContent = '';
  wrap.appendChild(table(
    ['Key bits', 'Exhaustive operations (2^b)', 'BSGS operations (2^(b/2))', 'Years for BSGS at 10^12 operations/s'],
    rows));
}

// ---- Modern scheme (WebCrypto) ----

let recipient = null;
let lastPacket = null;

async function eccSetup() {
  recipient = await generateRecipientKeys();
  $('ecc-plain').value = SAMPLE_EN;
  await eccEncryptShow(SAMPLE_EN);
  await eccStats();
}

async function eccEncryptShow(text) {
  lastPacket = await eccEncrypt(text, recipient.publicKey);
  const out = $('ecc-out');
  out.textContent = '';
  out.appendChild(table(['Field', 'Value (truncated)', 'Length (bytes)'], [
    ['Ephemeral public key', toHex(lastPacket.ephPub, 16) + '...', lastPacket.ephPub.length],
    ['IV (nonce)', toHex(lastPacket.iv, 12), lastPacket.iv.length],
    ['Ciphertext + GCM tag', toBase64(lastPacket.ct).slice(0, 48) + '...', lastPacket.ct.length],
  ]));
  $('ecc-decrypt-out').textContent = '';
}

async function eccDecryptShow() {
  if (!lastPacket) return;
  const text = await eccDecrypt(lastPacket, recipient.privateKey);
  const original = $('ecc-plain').value;
  $('ecc-decrypt-out').textContent =
    `Decrypted text: ${text}\nMatches the original message: ${text === original}`;
}

// Statistics on a longer ciphertext (the sample repeated 4 times), so that each byte value
// has enough expected counts for the chi-square test to be meaningful.
async function eccStats() {
  const long = SAMPLE_EN.repeat(4);
  const p = await eccEncrypt(long, recipient.publicKey);
  const byteChi = chiSquareUniform(byteCounts(p.ct));
  const b64 = toBase64(p.ct).replace(/[^A-Za-z]/g, '');
  const ranked = rankShifts(b64, ENGLISH_FREQ);
  const wrap = $('ecc-stats');
  wrap.textContent = '';
  wrap.appendChild(table(['Test', 'Chi-square', 'Critical (5%)', 'Verdict'], [
    ['ECIES ciphertext bytes vs uniform random (256 values)', byteChi.toFixed(1), CRIT_BYTES_DF255,
      byteChi < CRIT_BYTES_DF255 ? 'consistent with random bytes' : 'not random'],
    ['ECIES ciphertext letters (base64) vs English, best shift', ranked[0].chi.toFixed(1), CRIT_LETTERS_DF25,
      ranked[0].chi < CRIT_LETTERS_DF25 ? 'looks like English' : 'not English'],
  ]));
  wrap.appendChild(el('p',
    `Ciphertext length: ${p.ct.length} bytes, so each byte value is expected about ` +
    `${(p.ct.length / 256).toFixed(1)} times.`));
}

function showEccError(err) {
  const out = $('ecc-out');
  if (out) out.textContent = `WebCrypto error: ${err.message}. The page must run over https:// or localhost.`;
}

function init() {
  if (!$('caesar-plain')) return;

  // Caesar
  $('caesar-plain').value = SAMPLE_EN;
  $('caesar-key').value = '7';
  $('caesar-encrypt').addEventListener('click', runCaesarEncrypt);
  $('caesar-decrypt').addEventListener('click', runCaesarDecrypt);
  $('caesar-analyse').addEventListener('click', () => renderCompare(renderCaesarAttack()));
  $('caesar-lang').addEventListener('change', () => renderCompare(renderCaesarAttack()));
  $('caesar-length').addEventListener('click', () => guard(renderLengthExperiment));
  runCaesarEncrypt();
  renderReference();
  renderCompare(renderCaesarAttack());
  renderLengthExperiment();
  renderItalianLengths();

  // Toy elliptic-curve scheme
  $('toy-plain').value = 'ATTACK AT DAWN';
  $('toy-encrypt').addEventListener('click', () => guard(runToyEncrypt));
  $('toy-attack').addEventListener('click', () => guard(runToyAttack));
  guard(() => {
    renderToySetup();
    renderFeasibility();
    runToyEncrypt();
    renderToyStatistics();
  });

  // Modern scheme
  $('ecc-encrypt').addEventListener('click', () => eccEncryptShow($('ecc-plain').value).catch(showEccError));
  $('ecc-decrypt').addEventListener('click', () => eccDecryptShow().catch(showEccError));
  eccSetup().catch(showEccError);
}

if (typeof document !== 'undefined') {
  document.addEventListener('DOMContentLoaded', init);
}