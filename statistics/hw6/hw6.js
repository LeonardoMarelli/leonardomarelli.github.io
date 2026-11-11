// HMWK 6 - Ex.1 Riemann vs Lebesgue integration; Ex.2 numerical standard normal CDF
// Pure JavaScript, no libraries.

const PI = Math.PI;
const SQRT2 = Math.SQRT2;
const $ = id => document.getElementById(id);

/* =====================================================================
   Ex.1: integrals
   ===================================================================== */

// Regular functions: f evaluated in floating point; exact integral known.
const FUNCS = [
  { name: 'sin x on [0, π]', a: 0, b: PI, f: Math.sin, exact: 2 },
  { name: 'x² on [0, 1]', a: 0, b: 1, f: x => x * x, exact: 1 / 3 },
  { name: '1 if x<½, 3 if x≥½ on [0, 1]', a: 0, b: 1, f: x => (x < 0.5 ? 1 : 3), exact: 2 },
  { name: '|x − ½| on [0, 1]', a: 0, b: 1, f: x => Math.abs(x - 0.5), exact: 0.25 },
  { name: '1/√x on (0, 1]', a: 0, b: 1, f: x => 1 / Math.sqrt(x), exact: 2, cap: 100 },
  { name: '1 if x<½, 0 if x≥½ on [0, 1]', a: 0, b: 1, f: x => (x < 0.5 ? 1 : 0), exact: 0.5 },
  { name: 'sin x / x on [0, π]', a: 0, b: PI, f: x => (x === 0 ? 1 : Math.sin(x) / x), exact: 1.8519370520 },
];

// Pathological functions: no float evaluator. Exact oracles instead, because every float
// is a dyadic rational and sampling cannot see irrationals.
// oracle(i, N): [sup, inf] of f on I_i = [i/N, (i+1)/N].   levelMeasure(t): μ{f > t}, t > 0.

// Thomae: sup on I_i is 1/q, q = smallest denominator of a fraction p/q in I_i (q ≤ N always works).
function thomaeSup(i, N) {
  for (let q = 1; q <= N; q++) {
    const p = Math.ceil((i * q) / N);
    if (p * N <= (i + 1) * q) return 1 / q;
  }
}

const PATHO = [
  { name: 'Dirichlet', a: 0, b: 1, exact: 0, top: 1,
    oracle: () => [1, 0],               // sup = 1, inf = 0 on every interval
    levelMeasure: () => 0 },            // {f > t} = rationals, measure 0 for t > 0
  { name: 'Thomae', a: 0, b: 1, exact: 0, top: 1,
    oracle: (i, N) => [thomaeSup(i, N), 0],
    levelMeasure: () => 0 },            // {f > t} is finite for t > 0
];

// Riemann: N equal subintervals. Darboux sums U, L and the midpoint sum.
// Regular functions: sup/inf estimated from M+1 samples per subinterval.
function riemann(F, N, M = 20) {
  const h = (F.b - F.a) / N;
  let U = 0, L = 0, mid = 0;
  for (let i = 0; i < N; i++) {
    const x0 = F.a + i * h;
    let hi, lo;
    if (F.oracle) {
      [hi, lo] = F.oracle(i, N);
    } else {
      hi = -Infinity; lo = Infinity;
      for (let j = 0; j <= M; j++) {
        const v = F.f(x0 + (j * h) / M);
        if (Number.isFinite(v)) { if (v > hi) hi = v; if (v < lo) lo = v; }
      }
      mid += F.f(x0 + h / 2) * h;
    }
    U += hi * h;
    L += lo * h;
  }
  return { U, L, mid: F.oracle ? NaN : mid };
}

// Lebesgue (layer cake): ∫f = ∫₀^T μ{f > t} dt, midpoint rule in t with step delta.
// μ is estimated on NL midpoints of the domain: sort the values once, binary-search each level.
function lebesgue(F, NL, delta = 1e-3) {
  let measure, T;
  if (F.levelMeasure) {
    measure = F.levelMeasure;
    T = F.top;
  } else {
    const h = (F.b - F.a) / NL;
    const vals = Float64Array.from({ length: NL }, (_, i) => F.f(F.a + (i + 0.5) * h)).sort();
    measure = t => {
      let lo = 0, hi = NL;
      while (lo < hi) {
        const m = (lo + hi) >> 1;
        if (vals[m] > t) hi = m; else lo = m + 1;
      }
      return (NL - lo) / NL;
    };
    T = F.cap ?? vals[NL - 1];
  }
  const K = Math.ceil(T / delta);
  let s = 0;
  for (let k = 0; k < K; k++) s += measure((k + 0.5) * delta);
  return s * delta;
}

/* ---------- Ex.2 helpers ---------- */

const phi = x => Math.exp(-x * x / 2) / Math.sqrt(2 * PI);

// Reference values of Φ(x): Abramowitz & Stegun Table 26.1; DLMF chapter 7.
const REF = [[-3, 0.0013498980], [-2, 0.0227501319], [-1, 0.1586552539], [-0.5, 0.3085375387],
  [-1.96, 0.0249978951], [0, 0.5], [0.5, 0.6914624613], [1, 0.8413447461], [1.5, 0.9331927987],
  [1.96, 0.9750021049], [2, 0.9772498681], [2.5, 0.9937903347], [3, 0.9986501020]];
const REF_2 = 0.9772498681;

// 1. Trapezoid rule on [0,|x|]:  Φ(x) = ½ + sign(x) ∫₀^|x| φ
function cdfTrap(x, N = 1000) {
  const z = Math.abs(x), h = z / N;
  let s = (phi(0) + phi(z)) / 2;
  for (let i = 1; i < N; i++) s += phi(i * h);
  return 0.5 + Math.sign(x) * s * h;
}

// 2. Simpson's rule on [0,|x|] (N even)
function cdfSimpson(x, N = 200) {
  const z = Math.abs(x), h = z / N;
  let s = phi(0) + phi(z);
  for (let i = 1; i < N; i++) s += phi(i * h) * (i % 2 ? 4 : 2);
  return 0.5 + Math.sign(x) * s * h / 3;
}

// 3. Power series: Φ(x) = ½ + φ(x) Σ x^(2n+1)/(1·3·…·(2n+1)). Terms share the sign of x: no cancellation.
function cdfSeries(x) {
  let t = x, s = x;
  for (let n = 1; n < 1000; n++) {
    t *= x * x / (2 * n + 1);
    s += t;
    if (Math.abs(t) <= 1e-17 * Math.abs(s)) break;
  }
  return 0.5 + phi(x) * s;
}

// erf via Kummer's series: erf z = (2/√π) e^(−z²) Σ 2^n z^(2n+1)/(1·3·…·(2n+1)). All terms positive.
function erfSeries(z) {
  let t = z, s = z;
  for (let n = 1; n < 1000; n++) {
    t *= 2 * z * z / (2 * n + 1);
    s += t;
    if (Math.abs(t) <= 1e-17 * Math.abs(s)) break;
  }
  return (2 / Math.sqrt(PI)) * Math.exp(-z * z) * s;
}

// erfc via Laplace's continued fraction (z > 0), fixed depth D. Fast and accurate for z ≥ 2.
function erfcCF(z, D = 100) {
  let t = z;
  for (let k = D; k >= 1; k--) t = z + (k / 2) / t;
  return Math.exp(-z * z) / (Math.sqrt(PI) * t);
}

// 4. erf hybrid: Φ(x) = ½(1 + erf(x/√2)); in the tails, Φ(x) = ½ erfc(|x|/√2) for x < 0.
function cdfHybrid(x) {
  const z = Math.abs(x) / SQRT2;
  if (z < 2) return 0.5 * (1 + Math.sign(x) * erfSeries(z));
  const tail = 0.5 * erfcCF(z);     // P(X ≤ −|x|)
  return x > 0 ? 1 - tail : tail;
}

// 5. Abramowitz–Stegun 26.2.17: polynomial approximation, |error| < 7.5e-8.
function cdfAS(x) {
  const z = Math.abs(x);
  const t = 1 / (1 + 0.2316419 * z);
  const poly = t * (0.319381530 + t * (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429))));
  const Q = phi(z) * poly;          // P(X > |x|)
  return x >= 0 ? 1 - Q : Q;
}

const METHODS = [
  { name: 'Trapezoid (N=1000)', fn: x => cdfTrap(x, 1000) },
  { name: 'Simpson (N=200)', fn: x => cdfSimpson(x, 200) },
  { name: 'Power series', fn: cdfSeries },
  { name: 'erf hybrid', fn: cdfHybrid },
  { name: 'A&S 26.2.17', fn: cdfAS },
];

/* ---------- output helpers ---------- */

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
  // Headers are static strings written in this file, so they may contain <sub> and <sup>.
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

const f4 = x => (Number.isFinite(x) ? x.toFixed(4) : '–');
const e2 = x => (x === 0 ? '0' : x.toExponential(2));
const timed = fn => { const t0 = performance.now(); const v = fn(); return { v, ms: performance.now() - t0 }; };

/* ---------- runners ---------- */

function runEx1(out) {
  // A. Summary at fixed resolution: Riemann N = 1000, Lebesgue NL = 10^5
  const rowsA = FUNCS.map(F => {
    const r = timed(() => riemann(F, 1000));
    const l = timed(() => lebesgue(F, 100000));
    return [F.name, F.exact.toFixed(6), f4(r.v.mid), f4(r.v.U - r.v.L), f4(l.v),
      e2(Math.abs(r.v.mid - F.exact)), e2(Math.abs(l.v - F.exact)),
      r.ms.toFixed(2), l.ms.toFixed(2)];
  });
  out.appendChild(el('h3', 'A. Summary (Riemann: N = 1000 midpoint; Lebesgue: NL = 10⁵)'));
  out.appendChild(table(
    ['Function', 'Exact', 'Riemann mid', 'U<sub>N</sub> − L<sub>N</sub>', 'Lebesgue',
     '|err| Riemann', '|err| Lebesgue', 'ms Riemann', 'ms Lebesgue'],
    rowsA));

  // B. Convergence: N = NL, |error| Riemann / |error| Lebesgue
  const NS = [10, 100, 1000, 10000];
  const rowsB = FUNCS.map(F => [F.name, ...NS.map(N => {
    const r = riemann(F, N).mid;
    const l = lebesgue(F, N);
    return `${e2(Math.abs(r - F.exact))} / ${e2(Math.abs(l - F.exact))}`;
  })]);
  out.appendChild(el('h3', 'B. Convergence: |error| Riemann / |error| Lebesgue, with N = NL'));
  out.appendChild(table(['Function', ...NS.map(N => `N = ${N}`)], rowsB));

  // C. Pathological functions: exact oracles
  const rowsC = [];
  PATHO.forEach(P => [10, 100, 1000].forEach(N => {
    const r = riemann(P, N);
    rowsC.push([P.name, N, f4(r.U), f4(r.L), f4(r.U - r.L), f4(lebesgue(P, 0))]);
  }));
  out.appendChild(el('h3', 'C. Pathological functions (exact oracles)'));
  out.appendChild(table(
    ['Function', 'N', 'U<sub>N</sub>', 'L<sub>N</sub>', 'U<sub>N</sub> − L<sub>N</sub>', 'Lebesgue'],
    rowsC));
}

function runEx2(out) {
  // A. Errors at the reference points
  const rowsA = REF.map(([x, ref]) => [x, ref.toFixed(10),
    ...METHODS.map(M => e2(Math.abs(M.fn(x) - ref)))]);
  out.appendChild(el('h3', 'A. |error| at the reference points'));
  out.appendChild(table(['x', 'Φ(x) reference', ...METHODS.map(M => M.name)], rowsA));

  // B. Max error and cost (ms for 10^4 evaluations on [−8, 8])
  const GRID = Array.from({ length: 10000 }, (_, i) => -8 + (16 * i) / 9999);
  const rowsB = METHODS.map(M => {
    const maxErr = Math.max(...REF.map(([x, ref]) => Math.abs(M.fn(x) - ref)));
    const t = timed(() => { let s = 0; for (const x of GRID) s += M.fn(x); return s; });
    return [M.name, e2(maxErr), t.ms.toFixed(2)];
  });
  out.appendChild(el('h3', 'B. Accuracy and cost (cost = ms for 10⁴ evaluations on [−8, 8])'));
  out.appendChild(table(['Method', 'max |error|', 'ms / 10⁴ evals'], rowsB));

  // C. Convergence of the quadrature rules at x = 2
  const rowsC = [10, 100, 1000].map(N => [N,
    e2(Math.abs(cdfTrap(2, N) - REF_2)), e2(Math.abs(cdfSimpson(2, N) - REF_2))]);
  out.appendChild(el('h3', 'C. Convergence at x = 2: |error| vs N'));
  out.appendChild(table(['N', 'Trapezoid', 'Simpson'], rowsC));
}

export function runAll() {
  const ex1 = $('ex1'), ex2 = $('ex2');
  ex1.textContent = '';
  ex2.textContent = '';
  runEx1(ex1);
  runEx2(ex2);
}

export function init() {
  if (!$('run')) return;
  $('run').addEventListener('click', runAll);
  runAll();
}

if (typeof document !== 'undefined') {
  document.addEventListener('DOMContentLoaded', init);
}