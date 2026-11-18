// HMWK 7 - Statistical independence: contingency table, expected frequencies, chi-square test
// Pure JavaScript, no libraries. Seeded PRNG, so results are reproducible.

/* ---------- datasets: lists of observed pairs (x_k, f_k) ---------- */

export function mulberry32(seed) {
  let a = seed >>> 0;
  return function () {
    a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

// Real data: Titanic (2201 people, R's built-in Titanic table). X = travel class, F = outcome.
const TITANIC_COUNTS = {
  'Crew': { Survived: 212, Died: 673 },
  '1st': { Survived: 202, Died: 123 },
  '2nd': { Survived: 118, Died: 167 },
  '3rd': { Survived: 178, Died: 528 },
};

export function titanicPairs() {
  const pairs = [];
  for (const [x, row] of Object.entries(TITANIC_COUNTS))
    for (const [f, c] of Object.entries(row))
      for (let i = 0; i < c; i++) pairs.push([x, f]);
  return pairs;
}

// Simulated data where X and F are independent by construction (n = 600).
export function independentPairs(seed = 7) {
  const rng = mulberry32(seed);
  const pick = (labels, probs) => {
    const u = rng();
    let s = 0;
    for (let i = 0; i < labels.length; i++) { s += probs[i]; if (u < s) return labels[i]; }
    return labels[labels.length - 1];
  };
  return Array.from({ length: 600 }, () => [
    pick(['Linux', 'Windows', 'macOS'], [0.5, 0.3, 0.2]),
    pick(['Breached', 'Not breached'], [0.35, 0.65]),
  ]);
}

const DATASETS = {
  titanic: { name: 'Titanic: travel class (X) vs outcome (F)', xName: 'Class', fName: 'Outcome', pairs: titanicPairs },
  indep: { name: 'Simulated: operating system (X) vs breach (F), independent by construction',
    xName: 'OS', fName: 'Breach', pairs: () => independentPairs(7) },
};

/* ---------- contingency table ---------- */

// Builds the table n_ij from raw pairs. Rows = distinct x (first-seen order), columns = distinct f.
export function contingency(pairs) {
  const xs = [...new Set(pairs.map(p => p[0]))];
  const fs = [...new Set(pairs.map(p => p[1]))];
  const n = xs.map(() => fs.map(() => 0));
  for (const [x, f] of pairs) n[xs.indexOf(x)][fs.indexOf(f)]++;
  return { xs, fs, n };
}

// Marginal counts, marginal / conditional / joint relative frequencies.
export function frequencies({ n }) {
  const R = n.length, C = n[0].length;
  const rowSum = n.map(r => r.reduce((a, b) => a + b, 0));
  const colSum = Array.from({ length: C }, (_, j) => n.reduce((s, r) => s + r[j], 0));
  const N = rowSum.reduce((a, b) => a + b, 0);
  return {
    N, rowSum, colSum,
    joint: n.map(r => r.map(v => v / N)),                               // f_ij = n_ij / N
    fx: rowSum.map(v => v / N),                                         // f_i.
    ff: colSum.map(v => v / N),                                         // f_.j
    fGivenX: n.map((r, i) => r.map(v => v / rowSum[i])),                // f(F=j | X=i)
    xGivenF: n.map((r, i) => r.map((v, j) => v / colSum[j])),           // f(X=i | F=j)
  };
}

/* ---------- chi-square test ---------- */

// Expected counts under independence: E_ij = n_i. * n_.j / N.
export function expected({ rowSum, colSum, N }) {
  return rowSum.map(r => colSum.map(c => (r * c) / N));
}

export function chiSquare(n, E) {
  let chi2 = 0;
  const contrib = n.map((r, i) => r.map((o, j) => {
    const c = (o - E[i][j]) ** 2 / E[i][j];
    chi2 += c;
    return c;
  }));
  const df = (n.length - 1) * (n[0].length - 1);
  return { chi2, df, contrib };
}

// log Gamma (Lanczos approximation, g = 7).
function lgamma(z) {
  const g = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
    -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
  if (z < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * z)) - lgamma(1 - z);
  z -= 1;
  let x = g[0];
  for (let i = 1; i < 9; i++) x += g[i] / (z + i);
  const t = z + 7.5;
  return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x);
}

// Regularised upper incomplete gamma Q(a, x): series for x < a+1, continued fraction otherwise.
export function gammaQ(a, x) {
  if (x <= 0) return 1;
  const lg = a * Math.log(x) - x - lgamma(a);
  if (x < a + 1) {
    let term = 1 / a, sum = term;
    for (let n = 1; n < 1000; n++) {
      term *= x / (a + n);
      sum += term;
      if (Math.abs(term) < Math.abs(sum) * 1e-16) break;
    }
    return 1 - sum * Math.exp(lg);
  }
  const tiny = 1e-300;
  let b = x + 1 - a, c = 1 / tiny, d = 1 / b, h = d;
  for (let i = 1; i < 1000; i++) {
    const an = -i * (i - a);
    b += 2;
    d = an * d + b; if (Math.abs(d) < tiny) d = tiny;
    c = b + an / c; if (Math.abs(c) < tiny) c = tiny;
    d = 1 / d;
    const del = d * c;
    h *= del;
    if (Math.abs(del - 1) < 1e-16) break;
  }
  return Math.exp(lg) * h;
}

// p-value = P(chi2_df >= observed) = Q(df/2, chi2/2).
export const pValue = (chi2, df) => gammaQ(df / 2, chi2 / 2);

// Cramér's V: effect size in [0,1], independent of the sample size.
export const cramersV = (chi2, N, R, C) => Math.sqrt(chi2 / (N * (Math.min(R, C) - 1)));

export function analyse(pairs) {
  const T = contingency(pairs);
  const F = frequencies(T);
  const E = expected(F);
  const X = chiSquare(T.n, E);
  return { ...T, ...F, E, ...X,
    p: pValue(X.chi2, X.df),
    V: cramersV(X.chi2, F.N, T.n.length, T.n[0].length),
    minE: Math.min(...E.flat()) };
}

/* =====================================================================
   Page (browser only)
   ===================================================================== */

const $ = id => document.getElementById(id);

function el(tag, text, cls) {
  const node = document.createElement(tag);
  if (text !== undefined) node.textContent = text;
  if (cls) node.className = cls;
  return node;
}

// Wide tables scroll horizontally inside their box; cells never wrap.
function table(headers, rows) {
  const wrap = el('div');
  wrap.style.overflowX = 'auto';
  const t = el('table');
  t.style.fontSize = '0.8rem';
  const head = el('tr');
  headers.forEach(h => { const th = el('th'); th.innerHTML = h; head.appendChild(th); });
  t.appendChild(head);
  rows.forEach(r => {
    const tr = el('tr');
    r.forEach(c => tr.appendChild(el('td', String(c))));
    t.appendChild(tr);
  });
  t.querySelectorAll('th, td').forEach(c => { c.style.whiteSpace = 'nowrap'; c.style.padding = '6px 8px'; });
  wrap.appendChild(t);
  return wrap;
}

const f3 = x => x.toFixed(3);
const f2 = x => x.toFixed(2);
const pFmt = p => (p < 1e-4 ? p.toExponential(2) : p.toFixed(4));

export function render(key, out) {
  const D = DATASETS[key];
  const A = analyse(D.pairs());
  const { xs, fs, n } = A;
  const sec = el('div', undefined, 'output');
  sec.appendChild(el('h3', `${D.name} (N = ${A.N})`));

  sec.appendChild(el('h3', 'A. Contingency table and marginal counts'));
  sec.appendChild(table([`${D.xName} \\ ${D.fName}`, ...fs, 'Total n<sub>i·</sub>'],
    [...xs.map((x, i) => [x, ...n[i], A.rowSum[i]]), ['Total n<sub>·j</sub>', ...A.colSum, A.N]]));

  sec.appendChild(el('h3', 'B. Joint and marginal relative frequencies'));
  sec.appendChild(table([`${D.xName} \\ ${D.fName}`, ...fs.map(f => `f(${f})`), 'f<sub>i·</sub>'],
    [...xs.map((x, i) => [x, ...A.joint[i].map(f3), f3(A.fx[i])]), ['f<sub>·j</sub>', ...A.ff.map(f3), '1.000']]));

  sec.appendChild(el('h3', `C. Conditional frequencies of ${D.fName} given ${D.xName}`));
  sec.appendChild(table([D.xName, ...fs.map(f => `f(${f} | x)`), 'Sum'],
    [...xs.map((x, i) => [x, ...A.fGivenX[i].map(f3), f3(A.fGivenX[i].reduce((a, b) => a + b, 0))]),
      ['Marginal f<sub>·j</sub>', ...A.ff.map(f3), '1.000']]));
  sec.appendChild(el('p', `Under independence every row of this table equals the marginal row: f(${D.fName} | ${D.xName}) = f(${D.fName}).`));

  sec.appendChild(el('h3', `D. Conditional frequencies of ${D.xName} given ${D.fName}`));
  sec.appendChild(table([D.xName, ...fs.map(f => `f(x | ${f})`)],
    xs.map((x, i) => [x, ...A.xGivenF[i].map(f3)]).concat([['Sum', ...fs.map((_, j) => f3(A.xGivenF.reduce((s, r) => s + r[j], 0)))]])));
  sec.appendChild(el('p', `Under independence each column equals the marginal f(x): ${xs.map((x, i) => `${x} ${f3(A.fx[i])}`).join(', ')}.`));

  sec.appendChild(el('h3', 'E. Expected counts E<sub>ij</sub> = n<sub>i·</sub> n<sub>·j</sub> / N'));
  sec.appendChild(table([`${D.xName} \\ ${D.fName}`, ...fs], xs.map((x, i) => [x, ...A.E[i].map(f2)])));

  sec.appendChild(el('h3', 'F. Contributions (O − E)² / E'));
  sec.appendChild(table([`${D.xName} \\ ${D.fName}`, ...fs], xs.map((x, i) => [x, ...A.contrib[i].map(f2)])));

  sec.appendChild(el('h3', 'G. Test result'));
  sec.appendChild(table(['χ²', 'df = (r−1)(c−1)', 'p-value', "Cramér's V", 'min E<sub>ij</sub>'],
    [[f2(A.chi2), A.df, pFmt(A.p), f3(A.V), f2(A.minE)]]));

  const reject = A.p < 0.05;
  const strength = A.V < 0.1 ? 'negligible' : A.V < 0.3 ? 'weak to moderate' : 'strong';
  sec.appendChild(el('p',
    `Decision at α = 0.05: ${reject ? 'reject' : 'do not reject'} H₀ (independence). ` +
    (reject
      ? `The observed counts are far from those expected under independence, so ${D.xName} and ${D.fName} are statistically dependent. Association strength (Cramér's V = ${f3(A.V)}): ${strength}.`
      : `The data are compatible with independence. This is not proof of independence, only lack of evidence against it. Cramér's V = ${f3(A.V)}.`) +
    (A.minE >= 5 ? ' All expected counts are at least 5, so the χ² approximation is acceptable.'
      : ' Some expected counts are below 5: the χ² approximation is questionable.')));

  out.appendChild(sec);
}

export function run() {
  const out = $('results');
  out.textContent = '';
  const key = $('dataset').value;
  (key === 'both' ? Object.keys(DATASETS) : [key]).forEach(k => render(k, out));
}

export function init() {
  if (!$('run')) return;
  $('run').addEventListener('click', run);   // tables are generated on click
}

if (typeof document !== 'undefined') {
  document.addEventListener('DOMContentLoaded', init);
}