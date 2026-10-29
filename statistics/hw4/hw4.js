// HMWK 4 - Penetration score as a sum of i.i.d. increments (Bernoulli and Rademacher)
// Model: n systems per attacker, m attackers (independent trajectories).
// Pure JavaScript, no libraries. Seeded PRNG, so results are reproducible.

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

// One increment X_i. Bernoulli: {0,1} with P(1)=p. Rademacher: {-1,+1} with P(+1)=P(-1)=1/2.
export function increment(model, p, u) {
  if (model === 'bernoulli') return u < p ? 1 : 0;
  return u < 0.5 ? -1 : 1;
}

// Theoretical mean and variance of one increment.
export function moments(model, p) {
  return model === 'bernoulli' ? { mu: p, sigma2: p * (1 - p) } : { mu: 0, sigma2: 1 };
}

/* ---------- simulation ---------- */

// m independent trajectories of length max(ns). Returns, for each checkpoint n in ns,
// an array with S_n of every trajectory.
export function simulate({ model, p, m, ns, seed }) {
  const rng = mulberry32(seed);
  const nMax = Math.max(...ns);
  const out = ns.map(() => new Float64Array(m));
  for (let j = 0; j < m; j++) {
    let S = 0, c = 0;
    for (let i = 1; i <= nMax; i++) {
      S += increment(model, p, rng());
      if (i === ns[c]) { out[c][j] = S; c++; }
    }
  }
  return out;
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

// Standard normal CDF via the Abramowitz-Stegun erf approximation (error < 1.5e-7).
export function normCdf(x) {
  const z = x / Math.SQRT2;
  const sign = z < 0 ? -1 : 1;
  const t = 1 / (1 + 0.3275911 * Math.abs(z));
  const poly = ((((1.061405429 * t - 1.453152027) * t + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t;
  const erf = sign * (1 - poly * Math.exp(-z * z));
  return 0.5 * (1 + erf);
}

export function normPdf(x) {
  return Math.exp(-x * x / 2) / Math.sqrt(2 * Math.PI);
}

// Kolmogorov distance between the empirical CDF and the standard normal CDF.
export function ksDistance(samples) {
  const xs = Float64Array.from(samples).sort();
  const n = xs.length;
  let d = 0;
  for (let i = 0; i < n; i++) {
    const F = normCdf(xs[i]);
    d = Math.max(d, (i + 1) / n - F, F - i / n);
  }
  return d;
}

/* ---------- the experiment ---------- */

// For each n: the three quantities, their mean and variance against theory,
// the LLN deviation probability (eps = 0.05) against Chebyshev's bound, and the normal fit.
export function experiment({ model, p, m, ns, seed }) {
  const { mu, sigma2 } = moments(model, p);
  const S = simulate({ model, p, m, ns, seed });
  const eps = 0.05;

  return ns.map((n, c) => {
    const total = S[c];
    const avg = Array.from(total, v => v / n);
    const norm = Array.from(total, v => v / Math.sqrt(n));
    // Standardised for the CLT check: (S_n - n*mu) / sqrt(n*sigma2) -> N(0,1)
    const std = Array.from(total, v => (v - n * mu) / Math.sqrt(n * sigma2));

    const devCount = avg.filter(v => Math.abs(v - mu) > eps).length;
    return {
      n,
      meanTotal: mean(total), meanTotalTh: n * mu,
      varTotal: variance(total), varTotalTh: n * sigma2,
      meanAvg: mean(avg), meanAvgTh: mu,
      varAvg: variance(avg), varAvgTh: sigma2 / n,
      meanNorm: mean(norm), meanNormTh: Math.sqrt(n) * mu,
      varNorm: variance(norm), varNormTh: sigma2,
      pDev: devCount / m,
      chebyshev: Math.min(1, sigma2 / (n * eps * eps)),
      ks: ksDistance(std),
    };
  });
}

// Aggregate score of m attackers, each running n steps: A = sum_j S_n^(j).
// Samples `samples` aggregates. Checks Var(A) = m n sigma2 and Var(A/(mn)) = sigma2/(mn).
export function aggregate({ model, p, n, mList, samples, seed }) {
  const { mu, sigma2 } = moments(model, p);
  const rng = mulberry32(seed);
  return mList.map(m => {
    const A = new Float64Array(samples);
    for (let s = 0; s < samples; s++) {
      let tot = 0;
      for (let j = 0; j < m; j++) {
        for (let i = 0; i < n; i++) tot += increment(model, p, rng());
      }
      A[s] = tot;
    }
    const perStep = Array.from(A, v => v / (m * n));
    return {
      m,
      meanA: mean(A), meanATh: m * n * mu,
      varA: variance(A), varATh: m * n * sigma2,
      varPerStep: variance(perStep), varPerStepTh: sigma2 / (m * n),
    };
  });
}

// Histogram counts over [lo, hi] with `bins` equal bins.
export function histogram(samples, bins, lo, hi) {
  const counts = Array(bins).fill(0);
  const w = (hi - lo) / bins;
  for (const x of samples) {
    const i = Math.floor((x - lo) / w);
    if (i >= 0 && i < bins) counts[i]++;
  }
  return counts;
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

export function table(headers, rows) {
    const t = el('table');

    const thead = el('thead');
    const head = el('tr');

    headers.forEach(h => {
        head.appendChild(el('th', h));
    });

    thead.appendChild(head);
    t.appendChild(thead);

    const tbody = el('tbody');

    rows.forEach(row => {
        const tr = el('tr');

        // Gestisce anche un'eventuale riga racchiusa
        // accidentalmente in un ulteriore array.
        const cells =
            Array.isArray(row) &&
            row.length === 1 &&
            Array.isArray(row[0])
                ? row[0]
                : row;

        cells.forEach(cell => {
            tr.appendChild(
                el('td', cell == null ? '' : String(cell))
            );
        });

        tbody.appendChild(tr);
    });

    t.appendChild(tbody);

    return t;
}

const f3 = x => (Math.abs(x) < 1e-3 && x !== 0 ? x.toExponential(2) : x.toFixed(3));
const sci = x => (x === 0 ? '0' : x.toExponential(2));

export function renderModel(title, rows) {
  const wrap = el('div', undefined, 'output');
  wrap.appendChild(el('h3', title));
  wrap.appendChild(table(
    [
    'n',
    'E[Sn]', 'E[Sn] th.',
    'Var[Sn]', 'Var[Sn] th.',
    'E[Sn/n]', 'E[Sn/n] th.',
    'Var[Sn/n]', 'Var[Sn/n] th.',
    'P(deviation > 0.05)', 'Chebyshev bound',
    'E[Sn/√n]', 'E[Sn/√n] th.',
    'Var[Sn/√n]', 'Var[Sn/√n] th.',
    'KS distance'
],
    rows.map(r => [
      [
    r.n,

    f3(r.meanTotal),
    f3(r.meanTotalTh),

    f3(r.varTotal),
    f3(r.varTotalTh),

    f3(r.meanAvg),
    f3(r.meanAvgTh),

    sci(r.varAvg),
    sci(r.varAvgTh),

    f3(r.pDev),
    f3(r.chebyshev),

    f3(r.meanNorm),
    f3(r.meanNormTh),

    f3(r.varNorm),
    f3(r.varNormTh),

    f3(r.ks)
]
    ]),
  ));
  return wrap;
}

// SVG histogram of the standardised total score at the largest n, with the N(0,1) density.
export function histogramSvg(samples, title) {
  const NS = 'http://www.w3.org/2000/svg';
  const W = 640, H = 260, pad = 36, bins = 40, lo = -4, hi = 4;
  const counts = histogram(samples, bins, lo, hi);
  const w = (hi - lo) / bins;
  const maxDensity = Math.max(...counts.map(c => c / (samples.length * w)), normPdf(0));
  const sx = x => pad + ((x - lo) / (hi - lo)) * (W - 2 * pad);
  const sy = d => H - pad - (d / maxDensity) * (H - 2 * pad);

  const svg = document.createElementNS(NS, 'svg');
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  svg.setAttribute('width', '100%');
  svg.setAttribute('role', 'img');
  svg.setAttribute('aria-label', title);
  const add = (tag, attrs) => {
    const node = document.createElementNS(NS, tag);
    Object.entries(attrs).forEach(([k, v]) => node.setAttribute(k, String(v)));
    svg.appendChild(node);
  };
  counts.forEach((c, i) => {
    const d = c / (samples.length * w);
    add('rect', { x: sx(lo + i * w), y: sy(d), width: sx(lo + (i + 1) * w) - sx(lo + i * w) - 1,
      height: H - pad - sy(d), fill: '#547f98', opacity: 0.6 });
  });
  const pts = [];
  for (let x = lo; x <= hi + 1e-9; x += 0.02) pts.push(`${sx(x).toFixed(1)},${sy(normPdf(x)).toFixed(1)}`);
  add('polyline', { points: pts.join(' '), fill: 'none', stroke: '#202b33', 'stroke-width': 2 });
  add('line', { x1: pad, y1: H - pad, x2: W - pad, y2: H - pad, stroke: 'currentColor' });
  return svg;
}

export function run() {
  const model = $('model').value;
  const p = Number($('p').value);
  const m = Number($('m').value);
  const seed = Number($('seed').value) || 1;
  const err = $('error');
  err.textContent = '';
  if (!(p > 0 && p < 1) || !Number.isInteger(m) || m < 100 || m > 20000) {
    err.textContent = 'Use 0 < p < 1 and an integer m between 100 and 20000.';
    return;
  }
  const ns = [10, 100, 1000, 10000];
  const models = model === 'both' ? ['bernoulli', 'rademacher'] : [model];

  const out = $('results');
  out.textContent = '';
  models.forEach(md => {
    const rows = experiment({ model: md, p, m, ns, seed });
    out.appendChild(renderModel(
      md === 'bernoulli' ? `Bernoulli, p = ${p}` : 'Rademacher, P(±1) = 1/2', rows));

    // Histogram of the standardised total at the largest n (CLT check)
    const sim = simulate({ model: md, p, m, ns: [ns[ns.length - 1]], seed: seed + 1 })[0];
    const { mu, sigma2 } = moments(md, p);
    const n = ns[ns.length - 1];
    const std = Array.from(sim, v => (v - n * mu) / Math.sqrt(n * sigma2));
    const hw = el('div', undefined, 'output');
    hw.appendChild(el('h3', `Standardised total at n = ${n} vs N(0,1)`));
    hw.appendChild(histogramSvg(std, `Histogram of standardised total, ${md}`));
    out.appendChild(hw);
  });

  // Aggregate over attackers, at a fixed system size
  const aggN = 100;
  const aggModels = models;
  const aggWrap = $('aggregate');
  aggWrap.textContent = '';
  aggModels.forEach(md => {
    const rows = aggregate({ model: md, p, n: aggN, mList: [1, 10, 100], samples: 2000, seed: seed + 2 });
    const div = el('div', undefined, 'output');
    div.appendChild(el('h3', `${md === 'bernoulli' ? `Bernoulli, p = ${p}` : 'Rademacher'}: n = ${aggN} systems per attacker`));
    div.appendChild(table(
      ['m attackers', 'E[A]', 'E[A] th.', 'Var[A]', 'Var[A] th.', 'Var[A/(mn)]', 'Var[A/(mn)] th.'],
      rows.map(r => [r.m, f3(r.meanA), f3(r.meanATh), f3(r.varA), f3(r.varATh),
        sci(r.varPerStep), sci(r.varPerStepTh)]),
    ));
    aggWrap.appendChild(div);
  });
}

export function init() {
  if (!$('run')) return;
  $('run').addEventListener('click', run);
  run();
}

if (typeof document !== 'undefined') {
  document.addEventListener('DOMContentLoaded', init);
}