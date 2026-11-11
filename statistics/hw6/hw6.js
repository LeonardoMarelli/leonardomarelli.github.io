// HMWK 6 - Ex.1 Riemann vs Lebesgue integration; Ex.2 numerical standard normal CDF
// Pure JavaScript, no libraries.

const PI = Math.PI;
const SQRT2 = Math.SQRT2;
const $ = id => document.getElementById(id);

/* =====================================================================
   Ex.1: integrals
   ===================================================================== */

// Table A budget: Riemann with N_A intervals and M_A+1 samples each, Lebesgue with the same number of samples.
const N_A = 1000, M_A = 20;
const NL_A = N_A * (M_A + 1);

// Layer-cake step δ = 1/K with K ≈ 10√N. With δ = 1/K, the jump levels t = 1 and t = 3 of the
// step functions lie on grid boundaries.
const deltaFor = N => 1 / Math.round(10 * Math.sqrt(N));

// Each function: f (float evaluator), exact integral, top = sup of f (upper limit of the t-integral).
// Optional: oracle(i, N) = [sup, inf] on I_i (exact Darboux values), levelMeasure(t) = exact μ{f > t},
// tail(T) = ∫_T^∞ μ{f > t} dt (analytic correction for truncation).
const FUNCS = [
  { name: 'sin x on [0, π]', a: 0, b: PI, f: Math.sin, exact: 2, top: 1 },
  { name: 'x² on [0, 1]', a: 0, b: 1, f: x => x * x, exact: 1 / 3, top: 1 },
  { name: '1 if x<½, 3 if x≥½ on [0, 1]', a: 0, b: 1, f: x => (x < 0.5 ? 1 : 3), exact: 2, top: 3 },
  { name: '|x − ½| on [0, 1]', a: 0, b: 1, f: x => Math.abs(x - 0.5), exact: 0.25, top: 0.5 },
  {
    name: '1/√x on (0, 1]', a: 0, b: 1, f: x => 1 / Math.sqrt(x), exact: 2, top: 100,
    // Sup on the first interval is +∞; on [i/N, (i+1)/N] the inf is f((i+1)/N) because f is decreasing.
    oracle: (i, N) => [i === 0 ? Infinity : Math.sqrt(N / i), Math.sqrt(N / (i + 1))],
    // μ{f > t} = min(1, 1/t²) exactly, and ∫_T^∞ 1/t² dt = 1/T.
    levelMeasure: t => (t <= 1 ? 1 : 1 / (t * t)),
    tail: T => 1 / T,
  },
  { name: '1 if x<½, 0 if x≥½ on [0, 1]', a: 0, b: 1, f: x => (x < 0.5 ? 1 : 0), exact: 0.5, top: 1 },
  { name: 'sin x / x on [0, π]', a: 0, b: PI, f: x => (x === 0 ? 1 : Math.sin(x) / x), exact: 1.851937052, top: 1 },
];

// Pathological functions: no float evaluator. Exact oracles, because every float is a dyadic rational.
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
    levelMeasure: () => 0 },            // {f > t} is finite for t > 0, measure 0
];

// Riemann: N equal subintervals. Darboux sums U, L and the midpoint sum.
// Sampled functions: sup/inf from M+1 equally spaced samples; the centre sample (M even) is reused for the midpoint.
// Oracle functions: sup/inf from F.oracle; the midpoint is evaluated separately if F.f exists.
function riemann(F, N, M = 20) {
  const h = (F.b - F.a) / N;
  let U = 0, L = 0, mid = 0, evals = 0;
  for (let i = 0; i < N; i++) {
    const x0 = F.a + i * h;
    let hi, lo;
    if (F.oracle) {
      [hi, lo] = F.oracle(i, N);
      if (F.f) { mid += F.f(x0 + h / 2) * h; evals++; }
    } else {
      hi = -Infinity; lo = Infinity;
      for (let j = 0; j <= M; j++) {
        const v = F.f(x0 + (j * h) / M);
        evals++;
        if (v > hi) hi = v;
        if (v < lo) lo = v;
        if (j === M / 2) mid += v * h;
      }
    }
    U += hi * h;
    L += lo * h;
  }
  return { U, L, mid: (F.f || !F.oracle) ? mid : NaN, evals };
}

// Lebesgue (layer cake): ∫f = ∫₀^T μ{f > t} dt, midpoint rule in t with step delta.
// Exact level measure when available; otherwise μ ≈ (b − a) · #{f > t} / NL with sorted samples.
function lebesgue(F, NL, delta = 1e-3) {
  let measure, tail = 0;
  if (F.levelMeasure) {
    measure = F.levelMeasure;               // exact μ{f > t}
    if (F.tail) tail = F.tail(F.top);       // analytic ∫_T^∞ μ{f > t} dt
  } else {
    const h = (F.b - F.a) / NL;
    const vals = Float64Array.from({ length: NL }, (_, i) => F.f(F.a + (i + 0.5) * h)).sort();
    const len = F.b - F.a;
    measure = t => {
      let lo = 0, hi = NL;
      while (lo < hi) {
        const m = (lo + hi) >> 1;
        if (vals[m] > t) hi = m; else lo = m + 1;
      }
      return (len * (NL - lo)) / NL;        // μ ≈ (b − a) · fraction of samples with f > t
    };
  }
  const K = Math.round(F.top / delta);
  let s = 0;
  for (let k = 0; k < K; k++) s += measure((k + 0.5) * delta);
  return s * delta + tail;
}

/* =====================================================================
   Ex.2: normal CDF
   ===================================================================== */

const phi = x => Math.exp(-x * x / 2) / Math.sqrt(2 * PI);

// Reference values of Φ(x). Standard tabulated values; see Reference values in the page for sources.
const REF = [[-3, 0.0013498980], [-2, 0.0227501319], [-1, 0.1586552539], [-0.5, 0.3085375387],
  [-1.96, 0.0249978951], [0, 0.5], [0.5, 0.6914624613], [1, 0.8413447461], [1.5, 0.9331927987],
  [1.96, 0.9750021049], [2, 0.9772498681], [2.5, 0.9937903347], [3, 0.9986501020]];
const REF_2 = 0.9772498681;

// Upper tail Q(x) = P(Z > x), reference values for Table D.
const REF_Q = [[3, 1.349898032e-3], [5, 2.866515719e-7], [7, 1.279812544e-12]];

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

// 3. Power series: Φ(x) = ½ + φ(x) Σ x^(2n+1)/(1·3·…·(2n+1)).
// Accurate in absolute terms for moderate |x|. For large negative x, ½ + φS subtracts nearly equal
// quantities, so the relative accuracy of the left tail is lost. Not used for large |x| (overflow).
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

// erfc for z ≥ 0: series for z < 2 (1 − erf), continued fraction otherwise.
function erfcAny(z) {
  return z < 2 ? 1 - erfSeries(z) : erfcCF(z);
}

// Survival function Q(x) = P(Z > x) = ½ erfc(x/√2). Keeps relative accuracy in the upper tail.
function survival(x) {
  const z = x / SQRT2;
  return z >= 0 ? 0.5 * erfcAny(z) : 1 - 0.5 * erfcAny(-z);
}

// 4. erf hybrid: Φ(x) = 1 − Q(x) for x ≥ 0 (cancellation near 1 is unavoidable here, see Table D),
// and Φ(x) = Q(|x|) for x < 0 (no cancellation).
function cdfHybrid(x) {
  return x >= 0 ? 1 - survival(x) : survival(-x);
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

/* =====================================================================
   Output and benchmark helpers
   ===================================================================== */

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

// Number formatting: finite values with 4 decimals, ±∞ shown as symbols, NaN as a dash.
const f4 = x => (x === Infinity ? '+∞' : x === -Infinity ? '−∞' : Number.isFinite(x) ? x.toFixed(4) : '–');
const e2 = x => (x === 0 ? '0' : x.toExponential(2));

let SINK = 0;                       // keeps the benchmarked results alive
const median = a => { const s = [...a].sort((x, y) => x - y); return s[s.length >> 1]; };

// One warm-up run, then the median of `reps` timed runs. Times depend on the browser and the machine.
function bench(fn, reps = 7) {
  const v0 = fn();
  if (typeof v0 === 'number') SINK += v0;
  const ts = [];
  for (let r = 0; r < reps; r++) {
    const t0 = performance.now();
    const v = fn();
    ts.push(performance.now() - t0);
    if (typeof v === 'number') SINK += v;
  }
  return median(ts);
}

/* =====================================================================
   Runners
   ===================================================================== */

function runEx1(out) {
  // A. Summary at equal cost: Riemann N_A with M_A+1 samples each, Lebesgue NL_A samples.
  const rowsA = FUNCS.map(F => {
    const r = riemann(F, N_A, M_A);
    const l = lebesgue(F, NL_A);
    const msR = bench(() => riemann(F, N_A, M_A).mid);
    const msL = bench(() => lebesgue(F, NL_A));
    const evalsL = F.levelMeasure ? 0 : NL_A;
    return [F.name, F.exact.toFixed(6), f4(r.mid), f4(r.U - r.L), f4(l),
      e2(Math.abs(r.mid - F.exact)), e2(Math.abs(l - F.exact)),
      `${r.evals} / ${evalsL}`, msR.toFixed(2), msL.toFixed(2)];
  });
  out.appendChild(el('h3', `A. Summary at equal cost (Riemann: N = ${N_A}, M = ${M_A}; Lebesgue: NL = ${NL_A})`));
  out.appendChild(table(
    ['Function', 'Exact', 'Riemann mid', 'U<sub>N</sub> − L<sub>N</sub>', 'Lebesgue',
     '|err| Riemann', '|err| Lebesgue', 'f-evals R / L', 'ms R (median)', 'ms L (median)'],
    rowsA));

  // B. Convergence: N = NL, Lebesgue step δ shown in the header.
  const NS = [10, 100, 1000, 10000];
  const rowsB = FUNCS.map(F => [F.name, ...NS.map(N => {
    const r = riemann(F, N).mid;
    const l = lebesgue(F, N, deltaFor(N));
    return `${e2(Math.abs(r - F.exact))} / ${e2(Math.abs(l - F.exact))}`;
  })]);
  out.appendChild(el('h3', 'B. Convergence: |error| Riemann / |error| Lebesgue, with N = NL and δ shown'));
  out.appendChild(table(
    ['Function', ...NS.map(N => `N = NL = ${N}, δ = ${deltaFor(N).toPrecision(3)}`)],
    rowsB));
  out.appendChild(el('p',
    'For 1/√x the Lebesgue measure is exact, so NL does not enter; its error comes only from δ and the midpoint rule in t.'));

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
  const GRID = Array.from({ length: 10000 }, (_, i) => -8 + (16 * i) / 9999);

  // A. Errors at the reference points
  const rowsA = REF.map(([x, ref]) => [x, ref.toFixed(10),
    ...METHODS.map(M => e2(Math.abs(M.fn(x) - ref)))]);
  out.appendChild(el('h3', 'A. |error| at the reference points'));
  out.appendChild(table(['x', 'Φ(x) reference', ...METHODS.map(M => M.name)], rowsA));

  // B. Maximum error on the reference points, and median cost for the same 10⁴ points on [−8, 8].
  const rowsB = METHODS.map(M => {
    const maxErr = Math.max(...REF.map(([x, ref]) => Math.abs(M.fn(x) - ref)));
    const ms = bench(() => { let s = 0; for (const x of GRID) s += M.fn(x); return s; }, 5);
    return [M.name, e2(maxErr), ms.toFixed(2)];
  });
  out.appendChild(el('h3', 'B. Accuracy and cost'));
  out.appendChild(table(['Method', 'Maximum error on reference points', 'ms / 10⁴ evals (median)'], rowsB));
  out.appendChild(el('p',
    'Timings depend on the browser, the machine and the current load. Every method is evaluated at the same 10⁴ points.'));

  // C. Convergence of the quadrature rules at x = 2
  const rowsC = [10, 100, 1000].map(N => [N,
    e2(Math.abs(cdfTrap(2, N) - REF_2)), e2(Math.abs(cdfSimpson(2, N) - REF_2))]);
  out.appendChild(el('h3', 'C. Convergence at x = 2: |error| vs N'));
  out.appendChild(table(['N', 'Trapezoid', 'Simpson'], rowsC));

  // D. Upper tail: 1 − Φ(x) against the survival function Q(x)
  const rowsD = REF_Q.map(([x, q]) => {
    const naive = 1 - cdfHybrid(x);
    const direct = survival(x);
    return [x, q.toExponential(8),
      e2(Math.abs(naive - q) / q), e2(Math.abs(direct - q) / q)];
  });
  out.appendChild(el('h3', 'D. Upper tail Q(x) = P(Z > x): relative error of 1 − Φ(x) vs direct survival'));
  out.appendChild(table(['x', 'Q(x) reference', 'rel. err 1 − Φ(x)', 'rel. err survival(x)'], rowsD));

  // E. Symmetry check Φ(x) + Φ(−x) = 1 on the grid, for the erf hybrid.
  let worst = 0;
  for (const x of GRID) worst = Math.max(worst, Math.abs(cdfHybrid(x) + cdfHybrid(-x) - 1));
  out.appendChild(el('p', `Symmetry check for the erf hybrid on [−8, 8]: max |Φ(x) + Φ(−x) − 1| = ${e2(worst)}`));
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
}

if (typeof document !== 'undefined') {
  document.addEventListener('DOMContentLoaded', init);
}