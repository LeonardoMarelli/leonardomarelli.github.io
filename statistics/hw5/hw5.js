// HMWK 5 - Welford's recurrence: sampling mean and variance of penetration scores
// Scores X in {0,...,10} ~ Binomial(10, 0.3). Pure JavaScript, seeded PRNG.

/* ---------- random numbers ---------- */

export function mulberry32(seed) {
  let a = seed >>> 0;
  return function () {
    a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/* ---------- score distribution: X ~ Binomial(10, 0.3) on {0,...,10} ---------- */

const N_SCORE = 10;
const P_SCORE = 0.3;

function binom(n, k) {
  let c = 1;
  for (let i = 1; i <= k; i++) c = (c * (n - k + i)) / i;
  return c;
}

const pmf = Array.from({ length: N_SCORE + 1 }, (_, k) =>
  binom(N_SCORE, k) * P_SCORE ** k * (1 - P_SCORE) ** (N_SCORE - k));

// Cumulative sums of the pmf, for inverse-transform sampling.
const cdf = [];
pmf.reduce((s, q) => { cdf.push(s + q); return s + q; }, 0);

// Maps u in [0,1) to a score in {0,...,10}.
export function sampleScore(u) {
  let k = 0;
  while (k < N_SCORE && u > cdf[k]) k++;
  return k;
}

// Theoretical mean, variance and fourth central moment of the score distribution.
export function theory() {
  const mu = pmf.reduce((s, q, k) => s + k * q, 0);
  const sigma2 = pmf.reduce((s, q, k) => s + (k - mu) ** 2 * q, 0);
  const mu4 = pmf.reduce((s, q, k) => s + (k - mu) ** 4 * q, 0);
  return { mu, sigma2, mu4 };
}

/* ---------- Welford's online update ---------- */

export const newState = () => ({ n: 0, mean: 0, S: 0 });

// Adds one score x to the running state. Core logic: the three lines marked below.
export function welford(st, x) {
  st.n += 1;
  const d = x - st.mean;       // δ = x_n − M_{n−1}
  st.mean += d / st.n;         // M_n = M_{n−1} + δ / n
  st.S += d * (x - st.mean);   // SS_n = SS_{n−1} + δ(x_n − M_n)
  return st;
}

// Sampling variance s² = SS_n / (n − 1).
export const sampleVariance = st => st.S / (st.n - 1);

/* ---------- simulation ---------- */

// m hackers, each attempting n systems: returns each hacker's sampling mean and sampling variance.
export function simulate({ m, n, seed }) {
  const rng = mulberry32(seed);
  const means = new Float64Array(m);
  const vars = new Float64Array(m);
  for (let j = 0; j < m; j++) {
    const st = newState();
    for (let i = 0; i < n; i++) welford(st, sampleScore(rng()));
    means[j] = st.mean;
    vars[j] = sampleVariance(st);
  }
  return { means, vars };
}

/* ---------- statistics ---------- */

export function mean(a) {
  let s = 0;
  for (const x of a) s += x;
  return s / a.length;
}

export function variance(a) {
  const mu = mean(a);
  let s = 0;
  for (const x of a) s += (x - mu) ** 2;
  return s / (a.length - 1);
}

// Welford's SS against the two-pass formula on one trajectory: should be ~1e-12.
export function check(n, seed) {
  const rng = mulberry32(seed);
  const xs = Array.from({ length: n }, () => sampleScore(rng()));
  const st = newState();
  xs.forEach(x => welford(st, x));
  const mu = mean(xs);
  const ss = xs.reduce((s, x) => s + (x - mu) ** 2, 0);
  return Math.abs(st.S - ss);
}

/* ---------- the experiment ---------- */

// For each n: averages across hackers of the sampling mean and sampling variance,
// and the variability of each across hackers, against theory.
export function experiment({ m, ns, seed }) {
  const th = theory();
  return ns.map(n => {
    const { means, vars } = simulate({ m, n, seed: seed + n });
    return {
      n,
      avgMean: mean(means), varMean: variance(means), varMeanTh: th.sigma2 / n,
      avgVar: mean(vars), varVar: variance(vars),
      varVarTh: (th.mu4 - th.sigma2 ** 2 * (n - 3) / (n - 1)) / n,
    };
  });
}

/* =====================================================================
   Page (browser only)
   ===================================================================== */

const $ = id => document.getElementById(id);

export function el(tag, text, cls) {
  const node = document.createElement(tag);
  if (text !== undefined) node.textContent = text;
  if (cls) node.className = cls;
  return node;
}

// Wide tables scroll horizontally inside their box; cells never wrap.
export function table(headers, rows) {
  const wrap = el('div');
  wrap.style.overflowX = 'auto';
  const t = el('table');
  t.style.width = 'auto';
  t.style.minWidth = '100%';
  t.style.fontSize = '0.8rem';
  const head = el('tr');
  // Headers are static strings written in this file, so they may contain <sub> and <sup>.
  headers.forEach(h => {
    const th = el('th');
    th.innerHTML = h;
    head.appendChild(th);
  });
  t.appendChild(head);
  rows.forEach(r => {
    const tr = el('tr');
    r.forEach(cell => tr.appendChild(el('td', String(cell))));
    t.appendChild(tr);
  });
  t.querySelectorAll('th, td').forEach(c => {
    c.style.whiteSpace = 'nowrap';
    c.style.padding = '6px 8px';
  });
  wrap.appendChild(t);
  return wrap;
}

const f3 = x => (Math.abs(x) < 1e-3 && x !== 0 ? x.toExponential(2) : x.toFixed(3));
const sci = x => (x === 0 ? '0' : x.toExponential(2));

export function run() {
  const m = Number($('m').value);
  const seed = Number($('seed').value) || 1;
  const err = $('error');
  err.textContent = '';
  if (!Number.isInteger(m) || m < 100 || m > 20000) {
    err.textContent = 'Use an integer m between 100 and 20000.';
    return;
  }
  const th = theory();
  const rows = experiment({ m, ns: [5, 10, 50, 200, 1000], seed });

  const out = $('results');
  out.textContent = '';
  const div = el('div', undefined, 'output');
  div.appendChild(el('h3', `μ = ${th.mu.toFixed(2)}, σ² = ${th.sigma2.toFixed(2)}, m = ${m} hackers`));
  div.appendChild(table(
    ['n', 'avg of x̄<sub>j</sub>', 'μ', 'Var(x̄<sub>j</sub>)', 'σ²/n',
     'avg of s²<sub>j</sub>', 'σ²', 'Var(s²<sub>j</sub>)', 'Var(s²) theory'],
    rows.map(r => [r.n, f3(r.avgMean), f3(th.mu), sci(r.varMean), sci(r.varMeanTh),
      f3(r.avgVar), f3(th.sigma2), sci(r.varVar), sci(r.varVarTh)]),
  ));
  out.appendChild(div);

  $('check').textContent =
    `Welford check on 1000 scores: |SS_Welford − SS_two-pass| = ${check(1000, seed).toExponential(2)}`;
}

export function init() {
  if (!$('run')) return;
  $('run').addEventListener('click', run);
  run();
}

if (typeof document !== 'undefined') {
  document.addEventListener('DOMContentLoaded', init);
}