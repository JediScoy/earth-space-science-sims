/* Milankovitch Cycles Explorer — vanilla JS, no dependencies. */
(() => {
'use strict';

/* ------------------------------------------------------------------ *
 * Constants & helpers
 * ------------------------------------------------------------------ */
const D = Math.PI / 180, TWO_PI = Math.PI * 2;
const T0 = -400000, T1 = 100000, STEP = 1000;           // years from today
const N = (T1 - T0) / STEP + 1;
const S0 = 1361;                                          // solar constant, W/m²
const COL = { e: '#b28cff', o: '#3fe0a5', p: '#ffbf47', ins: '#ff7a85', ice: '#6cc4ff' };

const $ = (s, r = document) => r.querySelector(s);
const $$ = (s, r = document) => [...r.querySelectorAll(s)];
const clamp = (v, a, b) => Math.min(b, Math.max(a, v));
const mod = (a, n) => ((a % n) + n) % n;

/* Reference ice-age pattern (SCHEMATIC): approximate dates of glacial terminations,
 * years from today. Replace with real data (e.g. EPICA / LR04) if desired. */
const REFERENCE = { terminations: [-533000, -424000, -337000, -248000, -130000, -14000] };

/* ------------------------------------------------------------------ *
 * Orbital elements — simplified sums of the main periodic terms
 * ------------------------------------------------------------------ */
const orbit = {
  ecc(t) {
    const x = TWO_PI * t;
    const e = 0.0275
      + 0.0135 * Math.cos(x / 413000 + 2.520)
      + 0.0090 * Math.cos(x / 95000 + 0.456)
      + 0.0065 * Math.cos(x / 124000 + 3.973)
      + 0.0045 * Math.cos(x / 99000 + 2.414);
    return clamp(e, 0.004, 0.058);
  },
  obl(t) {                                                // degrees
    const x = TWO_PI * t;
    return 23.15 + 1.0 * Math.cos(x / 41000 + 1.292) + 0.216 * Math.cos(x / 53000 + 1.493) + 0.139 * Math.cos(x / 29000 + 4.660);
  },
  varpi(t) {                                              // longitude of perihelion, degrees from moving vernal equinox
    return mod(283 + 360 * t / 20900, 360);
  }
};
const NOW = { e: orbit.ecc(0), o: orbit.obl(0), w: orbit.varpi(0) };

/** Daily-mean insolation (W/m²) at latitude lat for solar longitude lam (90° = June solstice). */
function insolation(e, eps, varpiDeg, latDeg, lamDeg = 90) {
  const lam = lamDeg * D, phi = clamp(latDeg, -89.5, 89.5) * D;
  const delta = Math.asin(Math.sin(eps * D) * Math.sin(lam));
  const rho = (1 - e * e) / (1 + e * Math.cos(lam - varpiDeg * D));
  const c = -Math.tan(phi) * Math.tan(delta);
  const h0 = c >= 1 ? 0 : c <= -1 ? Math.PI : Math.acos(c);
  return S0 / Math.PI / (rho * rho) * (h0 * Math.sin(phi) * Math.sin(delta) + Math.cos(phi) * Math.cos(delta) * Math.sin(h0));
}

/* ------------------------------------------------------------------ *
 * State
 * ------------------------------------------------------------------ */
const state = {
  t: 0, lam: 90, exag: 8, lat: 65,
  sc: { e: true, o: true, p: true },
  playing: false, spin: false, speed: 30, iceSrc: 'model'
};

const series = { E: null, O: null, W: null, P: null, Qall: null, Q: null, Vmodel: null, Vref: null };
let norm = { mean: 0, std: 1 };

function eff(t) {                                           // effective elements given the toggles
  return {
    e: state.sc.e ? orbit.ecc(t) : NOW.e,
    o: state.sc.o ? orbit.obl(t) : NOW.o,
    w: state.sc.p ? orbit.varpi(t) : NOW.w
  };
}
function qAt(t, sc, lat) {
  const e = sc.e ? orbit.ecc(t) : NOW.e, o = sc.o ? orbit.obl(t) : NOW.o, w = sc.p ? orbit.varpi(t) : NOW.w;
  return insolation(e, o, w, lat);
}

/* ------------------------------------------------------------------ *
 * Ice models
 * ------------------------------------------------------------------ */
function refIce(t) {                                        // schematic saw-tooth, 0..1 (null in the future)
  if (t > 0) return null;
  const B = REFERENCE.terminations;
  if (t >= B[B.length - 1]) return 0.04;
  let i = 0;
  while (i < B.length - 2 && t >= B[i + 1]) i++;
  const a = B[i], b = B[i + 1], u = clamp((t - a) / (b - a), 0, 1);
  if (u < 0.08) return 0.04 + 0.1 * u / 0.08 * 0;           // interglacial
  const grow = 0.04 + 0.96 * Math.pow((Math.min(u, 0.92) - 0.08) / 0.84, 0.85);
  return u < 0.92 ? grow : grow * (1 - (u - 0.92) / 0.08);  // fast deglaciation
}

/* Toy relaxation-oscillator ice model (in the spirit of Paillard 1998): ice grows slowly while northern
 * summers are not too warm; once the ice sheet is large, a warm-summer peak can trigger fast melting. */
const ICE = { TAU_G: 47000, TAU_M: 6400, K: 2.1, F_GROW: 0.98, V_HI: 0.73, F_TERM: 0.14, V_LO: 0.045 };
function iceModel(sc, lat) {
  const out = new Float32Array(N);
  let V = 0.3, melting = false;
  for (let t = -1000000; t <= T1; t += STEP) {
    const F = (qAt(t, sc, lat) - norm.mean) / norm.std;           // + = warm summer
    if (!melting && V > ICE.V_HI && F > ICE.F_TERM) melting = true;
    else if (melting && V < ICE.V_LO) melting = false;
    if (melting) V -= V * (1 - Math.exp(-STEP / ICE.TAU_M));
    else {
      const eff = 1 / (1 + Math.exp(ICE.K * (F - ICE.F_GROW)));   // growth throttled by warm summers
      V += (1 - V) * eff * (1 - Math.exp(-STEP / ICE.TAU_G));
    }
    if (t >= T0) out[(t - T0) / STEP] = V;
  }
  return out;
}

function computeStatic() {
  series.E = new Float32Array(N); series.O = new Float32Array(N); series.W = new Float32Array(N); series.P = new Float32Array(N);
  for (let i = 0; i < N; i++) {
    const t = T0 + i * STEP;
    series.E[i] = orbit.ecc(t); series.O[i] = orbit.obl(t); series.W[i] = orbit.varpi(t);
    series.P[i] = series.E[i] * Math.sin(series.W[i] * D);
  }
  series.Vref = new Float32Array(N);
  for (let i = 0; i < N; i++) { const v = refIce(T0 + i * STEP); series.Vref[i] = v == null ? NaN : v; }
}

function computeScenario() {
  const all = { e: true, o: true, p: true };
  series.Qall = new Float32Array(N);
  let s = 0, s2 = 0;
  for (let i = 0; i < N; i++) { const q = qAt(T0 + i * STEP, all, state.lat); series.Qall[i] = q; s += q; s2 += q * q; }
  norm.mean = s / N; norm.std = Math.sqrt(s2 / N - norm.mean ** 2) || 1;
  series.Q = new Float32Array(N);
  for (let i = 0; i < N; i++) series.Q[i] = qAt(T0 + i * STEP, state.sc, state.lat);
  series.Vmodel = iceModel(state.sc, state.lat);
}

const idxOf = t => clamp((t - T0) / STEP, 0, N - 1);
function sample(arr, t) {
  const f = idxOf(t), i = Math.floor(f), j = Math.min(N - 1, i + 1), k = f - i;
  return arr[i] * (1 - k) + arr[j] * k;
}
function iceAt(t) {
  if (state.iceSrc === 'ref') { const v = refIce(t); return v == null ? sample(series.Vmodel, t) : v; }
  return sample(series.Vmodel, t);
}

/* ------------------------------------------------------------------ *
 * Canvas utilities
 * ------------------------------------------------------------------ */
function fit(canvas) {
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const w = canvas.clientWidth, h = canvas.clientHeight;
  if (canvas.width !== Math.round(w * dpr) || canvas.height !== Math.round(h * dpr)) {
    canvas.width = Math.round(w * dpr); canvas.height = Math.round(h * dpr);
  }
  const ctx = canvas.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return { ctx, w, h };
}

const fmtYears = t => {
  if (t === 0) return 'Today';
  const k = Math.abs(t) / 1000;
  return t < 0 ? `${k.toLocaleString()},000 years ago`.replace(/^0,000/, '0') : `In ${k.toLocaleString()},000 years`;
};
const SEASONS = ['Spring equinox', 'Summer solstice', 'Autumn equinox', 'Winter solstice'];
function seasonOf(lam) {                                    // northern-hemisphere season for solar longitude
  const l = mod(lam, 360);
  return ['Spring', 'Summer', 'Autumn', 'Winter'][Math.floor(l / 90)];
}
const perihelionSeason = w => seasonOf(w);

/* ------------------------------------------------------------------ *
 * Orbit view
 * ------------------------------------------------------------------ */
const orbitCanvas = $('#orbit');
function drawOrbit() {
  const { ctx, w, h } = fit(orbitCanvas);
  ctx.clearRect(0, 0, w, h);
  const { e, o, w: vp } = eff(state.t);
  const ev = Math.min(0.85, e * state.exag);
  const rOf = lam => (1 - ev * ev) / (1 + ev * Math.cos((lam - vp) * D));
  const posOf = lam => { const L = (lam + 180) * D, r = rOf(lam); return [r * Math.cos(L), -r * Math.sin(L)]; };

  // fit the orbit in the canvas
  let minx = 1e9, maxx = -1e9, miny = 1e9, maxy = -1e9;
  for (let l = 0; l < 360; l += 3) { const [x, y] = posOf(l); minx = Math.min(minx, x); maxx = Math.max(maxx, x); miny = Math.min(miny, y); maxy = Math.max(maxy, y); }
  const pad = 46, span = Math.max(maxx - minx, maxy - miny, 1.2);
  const sc = (Math.min(w, h) - pad * 2) / span;
  const cx = w / 2 - ((minx + maxx) / 2) * sc, cy = h / 2 - ((miny + maxy) / 2) * sc;
  const P = lam => { const [x, y] = posOf(lam); return [cx + x * sc, cy + y * sc]; };

  // reference circle (e = 0) for contrast
  ctx.strokeStyle = 'rgba(255,255,255,.12)'; ctx.setLineDash([4, 5]); ctx.lineWidth = 1;
  ctx.beginPath(); ctx.arc(cx, cy, sc, 0, TWO_PI); ctx.stroke(); ctx.setLineDash([]);

  // orbit ellipse
  ctx.beginPath();
  for (let l = 0; l <= 360; l += 2) { const [x, y] = P(l); l ? ctx.lineTo(x, y) : ctx.moveTo(x, y); }
  ctx.closePath(); ctx.strokeStyle = 'rgba(160,185,255,.75)'; ctx.lineWidth = 2; ctx.stroke();

  // sun glow
  const g = ctx.createRadialGradient(cx, cy, 2, cx, cy, 34);
  g.addColorStop(0, 'rgba(255,235,150,1)'); g.addColorStop(.35, 'rgba(255,190,70,.55)'); g.addColorStop(1, 'rgba(255,170,40,0)');
  ctx.fillStyle = g; ctx.beginPath(); ctx.arc(cx, cy, 34, 0, TWO_PI); ctx.fill();
  ctx.fillStyle = '#ffd966'; ctx.beginPath(); ctx.arc(cx, cy, 9, 0, TWO_PI); ctx.fill();

  // perihelion / aphelion markers
  const [px, py] = P(vp), [ax, ay] = P(vp + 180);
  ctx.strokeStyle = 'rgba(255,191,71,.55)'; ctx.lineWidth = 1; ctx.setLineDash([3, 4]);
  ctx.beginPath(); ctx.moveTo(ax, ay); ctx.lineTo(px, py); ctx.stroke(); ctx.setLineDash([]);
  ctx.font = '600 11px system-ui,sans-serif'; ctx.textAlign = 'center';
  const lab = (txt, x, y, col) => {
    const dx = x - cx, dy = y - cy, d = Math.hypot(dx, dy) || 1;
    ctx.fillStyle = col; ctx.beginPath(); ctx.arc(x, y, 4, 0, TWO_PI); ctx.fill();
    ctx.fillText(txt, x + dx / d * 24, y + dy / d * 24 + 4);
  };
  lab('Perihelion', px, py, COL.p);
  lab('Aphelion', ax, ay, '#8e99b8');

  // season ticks
  [[0, 'Mar'], [90, 'Jun'], [180, 'Sep'], [270, 'Dec']].forEach(([l, t]) => {
    const [x, y] = P(l);
    ctx.fillStyle = 'rgba(255,255,255,.55)'; ctx.beginPath(); ctx.arc(x, y, 2.5, 0, TWO_PI); ctx.fill();
    const dx = x - cx, dy = y - cy, d = Math.hypot(dx, dy) || 1;
    ctx.font = '11px system-ui,sans-serif'; ctx.fillStyle = 'rgba(200,210,240,.7)';
    ctx.fillText(t, x - dx / d * 16, y - dy / d * 16 + 4);
  });

  // Earth
  const [ex, ey] = P(state.lam);
  const sdx = cx - ex, sdy = cy - ey, sd = Math.hypot(sdx, sdy) || 1;
  const R = 13;
  const eg = ctx.createRadialGradient(ex + sdx / sd * 5, ey + sdy / sd * 5, 1, ex, ey, R);
  eg.addColorStop(0, '#9fd0ff'); eg.addColorStop(.6, '#2c6bd1'); eg.addColorStop(1, '#0c1d45');
  ctx.fillStyle = eg; ctx.beginPath(); ctx.arc(ex, ey, R, 0, TWO_PI); ctx.fill();
  ctx.strokeStyle = 'rgba(255,255,255,.35)'; ctx.lineWidth = 1; ctx.stroke();
  // pole marker leans toward where the Sun is at the June solstice (screen up in this view)
  const lean = Math.sin(o * D) * R * 1.0;
  ctx.fillStyle = '#fff'; ctx.beginPath(); ctx.arc(ex, ey - lean, 2.6, 0, TWO_PI); ctx.fill();
  ctx.strokeStyle = 'rgba(255,255,255,.7)'; ctx.beginPath(); ctx.moveTo(ex, ey); ctx.lineTo(ex, ey - lean - 9); ctx.stroke();

  // eccentricity caption
  ctx.textAlign = 'left'; ctx.fillStyle = 'rgba(200,210,240,.8)'; ctx.font = '12px system-ui,sans-serif';
  ctx.fillText(`e = ${e.toFixed(3)}`, 12, 20);
  ctx.textAlign = 'right'; ctx.fillText('dashed circle = zero eccentricity', w - 12, 20);
}

/* ------------------------------------------------------------------ *
 * Tilt view
 * ------------------------------------------------------------------ */
const tiltCanvas = $('#tilt');
function drawTilt() {
  const { ctx, w, h } = fit(tiltCanvas);
  ctx.clearRect(0, 0, w, h);
  const { o } = eff(state.t);
  const cx = w * 0.56, cy = h * 0.5, R = Math.min(h * 0.36, w * 0.28);
  const eps = o * D;

  // sunlight arrows from the left
  ctx.strokeStyle = 'rgba(255,205,90,.55)'; ctx.fillStyle = 'rgba(255,205,90,.8)'; ctx.lineWidth = 2;
  for (let i = -3; i <= 3; i++) {
    const y = cy + i * R * 0.33, x0 = 14, x1 = cx - R * 1.35 + Math.abs(i) * 4;
    ctx.beginPath(); ctx.moveTo(x0, y); ctx.lineTo(x1, y); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(x1 + 7, y); ctx.lineTo(x1 - 2, y - 4); ctx.lineTo(x1 - 2, y + 4); ctx.fill();
  }
  ctx.fillStyle = 'rgba(255,205,90,.9)'; ctx.font = '600 11px system-ui,sans-serif'; ctx.textAlign = 'left';
  ctx.fillText('Sunlight (June solstice)', 12, 18);

  // Earth disk with day-side shading
  const g = ctx.createLinearGradient(cx - R, 0, cx + R, 0);
  g.addColorStop(0, '#3f86e6'); g.addColorStop(.55, '#1d4a99'); g.addColorStop(1, '#0a1a3d');
  ctx.fillStyle = g; ctx.beginPath(); ctx.arc(cx, cy, R, 0, TWO_PI); ctx.fill();

  ctx.save(); ctx.translate(cx, cy);
  // axis tilts toward the Sun (left) at the June solstice
  const ax = -Math.sin(eps), ay = -Math.cos(eps);
  const latRing = (lat, col, wd, dash) => {
    const p = lat * D, rx = Math.cos(p) * R, off = Math.sin(p) * R;
    ctx.strokeStyle = col; ctx.lineWidth = wd; ctx.setLineDash(dash || []);
    ctx.beginPath(); ctx.ellipse(ax * off, ay * off, rx, rx * 0.28, -eps, 0, TWO_PI); ctx.stroke(); ctx.setLineDash([]);
  };
  latRing(0, 'rgba(255,255,255,.35)', 1);
  latRing(state.lat, COL.ins, 2.5);
  latRing(90 - o, 'rgba(160,230,255,.9)', 1.5, [4, 3]);
  // axis
  ctx.strokeStyle = '#fff'; ctx.lineWidth = 2;
  ctx.beginPath(); ctx.moveTo(-ax * R * 1.35, -ay * R * 1.35); ctx.lineTo(ax * R * 1.35, ay * R * 1.35); ctx.stroke();
  ctx.fillStyle = '#fff'; ctx.font = '700 11px system-ui,sans-serif'; ctx.textAlign = 'center';
  ctx.fillText('N', ax * R * 1.5, ay * R * 1.5 + 4);
  // vertical reference + angle arc
  ctx.strokeStyle = 'rgba(255,255,255,.35)'; ctx.setLineDash([3, 4]); ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(0, 0); ctx.lineTo(0, -R * 1.35); ctx.stroke(); ctx.setLineDash([]);
  ctx.strokeStyle = COL.o; ctx.lineWidth = 2;
  ctx.beginPath(); ctx.arc(0, 0, R * 0.55, -Math.PI / 2 - eps, -Math.PI / 2); ctx.stroke();
  ctx.fillStyle = COL.o; ctx.textAlign = 'left';
  ctx.fillText(`${o.toFixed(1)}°`, 6, -R * 0.62);
  ctx.restore();

  ctx.font = '11px system-ui,sans-serif'; ctx.textAlign = 'left';
  ctx.fillStyle = COL.ins; ctx.fillText(`● ${state.lat}°N`, w - 78, h - 28);
  ctx.fillStyle = 'rgba(160,230,255,.95)'; ctx.fillText(`┄ Arctic Circle ${(90 - o).toFixed(1)}°N`, w - 150, h - 10);
  $('#tiltNote').textContent = `Axis tilt ${o.toFixed(2)}° (range ≈ 22.1°–24.5°). The Arctic Circle sits at ${(90 - o).toFixed(1)}°N; more tilt means more summer sunlight at high latitudes.`;
}

/* ------------------------------------------------------------------ *
 * Timeline lanes
 * ------------------------------------------------------------------ */
const lanesEl = $('#lanes');
const lanes = [
  { id: 'ecc', key: 'e', title: 'Eccentricity', sub: 'orbit shape · ~100 & 413 kyr', color: COL.e, data: () => series.E, range: [0, 0.062], fmt: v => v.toFixed(3) },
  { id: 'obl', key: 'o', title: 'Tilt', sub: 'axis angle · ~41 kyr', color: COL.o, data: () => series.O, range: [21.8, 24.8], fmt: v => v.toFixed(2) + '°' },
  { id: 'pre', key: 'p', title: 'Precession index', sub: 'e · sin ϖ · ~19–23 kyr', color: COL.p, data: () => series.P, range: [-0.062, 0.062], fmt: v => v.toFixed(3) },
  { id: 'ins', title: 'Summer sunlight', sub: 'June-solstice, daily mean', color: COL.ins, data: () => series.Q, ghost: () => series.Qall, range: [400, 570], fmt: v => v.toFixed(0) + ' W/m²', lat: true },
  { id: 'ice', title: 'Ice volume', sub: 'model vs. ice-age pattern', color: COL.ice, data: () => series.Vmodel, ghost: () => series.Vref, ghostDash: true, range: [0, 1.05], fmt: v => Math.round(v * 100) + '%', big: true }
];

function buildLanes() {
  lanes.forEach(L => {
    const el = document.createElement('div');
    el.className = 'lane' + (L.big ? ' big' : '');
    el.innerHTML = `
      <div class="lane-head">
        <span class="dot" style="background:${L.color}"></span>
        <h4>${L.id === 'ins' ? 'Summer sunlight at <span class="latv">' + state.lat + '</span>°N' : L.title}</h4>
        <span class="sub">${L.sub}</span>
        ${L.lat ? '<label class="lat-ctl">Latitude <input type="range" id="lat" min="30" max="85" step="1" value="' + state.lat + '" aria-label="Latitude"></label>' : ''}
        <span class="val" data-v style="color:${L.color}"></span>
      </div>
      <div class="chartbox"><canvas class="base"></canvas><canvas class="over"></canvas></div>
      ${L.id === 'ins' ? '<div class="legend"><span style="--c:' + COL.ins + '">Selected cycles</span><span class="dash" style="--c:#8e99b8">All three cycles</span></div>' : ''}
      ${L.id === 'ice' ? '<div class="legend"><span style="--c:' + COL.ice + '">Model (from sunlight above)</span><span class="dash" style="--c:#e8ecf8">Ice-age pattern (schematic)</span></div>' : ''}`;
    lanesEl.appendChild(el);
    L.el = el; L.base = $('.base', el); L.over = $('.over', el); L.box = $('.chartbox', el); L.val = $('[data-v]', el);
    const pick = ev => {
      const r = L.box.getBoundingClientRect();
      const x = clamp((ev.clientX - r.left - PADX) / (r.width - PADX * 2), 0, 1);
      setTime(Math.round((T0 + x * (T1 - T0)) / STEP) * STEP, true);
    };
    let down = false;
    L.box.addEventListener('pointerdown', ev => { down = true; L.box.setPointerCapture(ev.pointerId); pick(ev); });
    L.box.addEventListener('pointermove', ev => { if (down) pick(ev); });
    L.box.addEventListener('pointerup', () => { down = false; });
    L.box.addEventListener('pointercancel', () => { down = false; });
  });
  $('#lat').addEventListener('input', ev => {
    state.lat = +ev.target.value;
    $$('.latv').forEach(s => s.textContent = state.lat);
    scenarioChanged();
  });
}

const PADX = 10, PADY = 8;
const xOf = (t, w) => PADX + (t - T0) / (T1 - T0) * (w - PADX * 2);

function drawLaneBase(L) {
  const { ctx, w, h } = fit(L.base);
  ctx.clearRect(0, 0, w, h);
  const [lo, hi] = L.range;
  const yOf = v => h - PADY - 12 - (clamp(v, lo, hi) - lo) / (hi - lo) * (h - PADY * 2 - 12);
  // grid
  ctx.font = '10px system-ui,sans-serif'; ctx.textAlign = 'center';
  for (let t = T0; t <= T1; t += 100000) {
    const x = xOf(t, w);
    ctx.strokeStyle = t === 0 ? 'rgba(255,255,255,.28)' : 'rgba(255,255,255,.07)'; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(x, PADY); ctx.lineTo(x, h - 12); ctx.stroke();
    ctx.fillStyle = 'rgba(160,172,205,.8)';
    ctx.fillText(t === 0 ? 'now' : `${t / 1000}k`, clamp(x, 14, w - 14), h - 1);
  }
  const path = (arr, close) => {
    ctx.beginPath(); let started = false, lastX = 0;
    for (let i = 0; i < N; i++) {
      const v = arr[i]; if (Number.isNaN(v)) { started = false; continue; }
      const x = xOf(T0 + i * STEP, w), y = yOf(v);
      if (!started) { ctx.moveTo(x, y); started = true; } else ctx.lineTo(x, y);
      lastX = x;
    }
    return lastX;
  };
  if (L.ghost) {
    const gh = L.ghost();
    if (gh) {
      ctx.strokeStyle = L.ghostDash ? 'rgba(232,236,248,.85)' : 'rgba(142,153,184,.7)'; ctx.lineWidth = 1.4;
      ctx.setLineDash([5, 4]); path(gh); ctx.stroke(); ctx.setLineDash([]);
    }
  }
  const arr = L.data();
  ctx.lineJoin = 'round';
  const lastX = path(arr);
  ctx.strokeStyle = L.color; ctx.lineWidth = 2; ctx.stroke();
  ctx.lineTo(lastX, h - 12); ctx.lineTo(xOf(T0, w), h - 12); ctx.closePath();
  const g = ctx.createLinearGradient(0, 0, 0, h);
  g.addColorStop(0, L.color + '44'); g.addColorStop(1, L.color + '00');
  ctx.fillStyle = g; ctx.fill();
  L.yOf = yOf;
}

function drawLaneOver(L) {
  const { ctx, w, h } = fit(L.over);
  ctx.clearRect(0, 0, w, h);
  const x = xOf(state.t, w);
  ctx.strokeStyle = '#fff'; ctx.lineWidth = 1.5;
  ctx.beginPath(); ctx.moveTo(x, PADY - 2); ctx.lineTo(x, h - 12); ctx.stroke();
  const v = sample(L.data(), state.t);
  if (L.yOf) {
    ctx.fillStyle = L.color; ctx.strokeStyle = '#fff'; ctx.lineWidth = 2;
    ctx.beginPath(); ctx.arc(x, L.yOf(v), 4.5, 0, TWO_PI); ctx.fill(); ctx.stroke();
  }
  L.val.textContent = L.fmt(v);
}

function drawLanesBase() {
  lanes.forEach(L => { L.el.classList.toggle('off', L.key ? !state.sc[L.key] : false); drawLaneBase(L); });
  lanes.forEach(drawLaneOver);
}
const drawLanesOver = () => lanes.forEach(drawLaneOver);

/* ------------------------------------------------------------------ *
 * Glacier map
 * ------------------------------------------------------------------ */
const mapCanvas = $('#map');
// coarse outlines [lon, lat]
const LAND = {
  na: [[-165,64],[-168,66],[-162,70],[-156,71],[-141,70],[-128,70],[-115,68],[-100,68],[-90,69],[-82,68],[-86,64],[-94,60],[-92,57],[-85,55],[-82,52],[-79,52],[-78,58],[-77,62],[-70,61],[-65,60],[-61,56],[-56,52],[-60,47],[-66,45],[-70,42],[-74,40],[-76,35],[-81,31],[-80,26],[-82,25],[-83,29],[-90,29],[-94,29],[-97,26],[-97,22],[-105,21],[-110,23],[-112,30],[-117,32],[-120,34],[-124,40],[-124,47],[-130,54],[-135,58],[-145,60],[-152,59],[-158,57],[-165,60]],
  gl: [[-73,78],[-60,82],[-30,83],[-20,80],[-20,72],[-30,68],[-42,60],[-50,63],[-55,68],[-58,75]],
  eu: [[-9,37],[-9,43],[-1,46],[-4,48],[2,51],[8,54],[10,57],[11,59],[5,59],[5,62],[14,67],[20,70],[28,71],[40,67],[44,68],[60,69],[70,73],[80,73],[100,77],[113,74],[130,71],[150,71],[170,70],[180,68],[178,65],[170,60],[163,58],[160,52],[156,57],[142,59],[137,54],[141,48],[135,43],[129,41],[127,38],[126,35],[122,38],[121,40],[118,38],[122,31],[120,25],[110,21],[100,13],[97,17],[92,22],[88,22],[80,15],[77,10],[73,17],[70,22],[67,25],[57,25],[56,27],[50,30],[48,29],[52,25],[57,22],[44,13],[40,20],[35,28],[32,31],[10,34],[-5,35]],
  uk: [[-5,50],[1,51],[2,53],[-2,56],[-3,58],[-6,58],[-5,54],[-3,53]],
  is: [[-24,65],[-14,66],[-14,64],[-22,63.5]]
};
const bump = (lon, c, wd) => { let d = Math.abs(mod(lon - c + 180, 360) - 180); return Math.exp(-((d / wd) ** 2)); };
const reach = lon => 12 + 40 * bump(lon, -95, 32) + 14 * bump(lon, -128, 14) + 32 * bump(lon, 15, 24) + 11 * bump(lon, 90, 50);

function drawMap() {
  const { ctx, w, h } = fit(mapCanvas);
  ctx.clearRect(0, 0, w, h);
  const R = Math.min(w, h) / 2 - 6, cx = w / 2, cy = h / 2;
  const proj = (lon, lat) => { const r = (90 - lat) / 60 * R, a = lon * D; return [cx + r * Math.sin(a), cy - r * Math.cos(a)]; };
  const V = Math.min(1.1, iceAt(state.t) / 0.8);   // 1 = full glacial maximum

  ctx.save(); ctx.beginPath(); ctx.arc(cx, cy, R, 0, TWO_PI); ctx.clip();
  const og = ctx.createRadialGradient(cx, cy, 0, cx, cy, R);
  og.addColorStop(0, '#133a63'); og.addColorStop(1, '#0b2a4d');
  ctx.fillStyle = og; ctx.fillRect(0, 0, w, h);

  const poly = (pts, fill, stroke) => {
    ctx.beginPath(); pts.forEach(([lo, la], i) => { const [x, y] = proj(lo, la); i ? ctx.lineTo(x, y) : ctx.moveTo(x, y); });
    ctx.closePath(); if (fill) { ctx.fillStyle = fill; ctx.fill(); } if (stroke) { ctx.strokeStyle = stroke; ctx.lineWidth = 1; ctx.stroke(); }
  };
  // Lower sea level exposes shelves: thicken land slightly (cosmetic) — skip, keep clean.
  ['eu', 'na', 'uk', 'is'].forEach(k => poly(LAND[k], '#3b5a3d', 'rgba(255,255,255,.18)'));
  poly(LAND.gl, '#eaf4ff', 'rgba(255,255,255,.5)');

  // ice sheets + sea ice
  const iceRing = (fn, fill, stroke) => {
    ctx.beginPath();
    for (let lon = -180; lon <= 180; lon += 2) { const [x, y] = proj(lon, 90 - fn(lon)); lon === -180 ? ctx.moveTo(x, y) : ctx.lineTo(x, y); }
    ctx.closePath(); ctx.fillStyle = fill; ctx.fill(); if (stroke) { ctx.strokeStyle = stroke; ctx.lineWidth = 1.5; ctx.stroke(); }
  };
  const seaC = 13 + 9 * V;
  iceRing(() => seaC, 'rgba(210,235,255,.28)');
  if (V > 0.02) iceRing(lon => Math.max(seaC * 0.8, V * reach(lon)), 'rgba(235,246,255,.82)', 'rgba(255,255,255,.95)');

  // graticule
  ctx.strokeStyle = 'rgba(255,255,255,.14)'; ctx.lineWidth = 1; ctx.fillStyle = 'rgba(255,255,255,.55)'; ctx.font = '10px system-ui,sans-serif'; ctx.textAlign = 'left';
  [75, 60, 45].forEach(la => { ctx.beginPath(); ctx.arc(cx, cy, (90 - la) / 60 * R, 0, TWO_PI); ctx.stroke(); const [x, y] = proj(180, la); ctx.fillText(la + '°N', x + 4, y - 3); });
  ctx.restore();
  ctx.strokeStyle = 'rgba(255,255,255,.35)'; ctx.beginPath(); ctx.arc(cx, cy, R, 0, TWO_PI); ctx.stroke();

  // Great Lakes marker
  const gl = { lon: -85, lat: 44 };
  const [mx, my] = proj(gl.lon, gl.lat);
  const covered = 90 - V * reach(gl.lon) <= gl.lat;
  ctx.fillStyle = covered ? '#ff7a85' : '#ffd966'; ctx.strokeStyle = '#000'; ctx.lineWidth = 1;
  ctx.beginPath(); ctx.arc(mx, my, 5, 0, TWO_PI); ctx.fill(); ctx.stroke();
  ctx.font = '600 11px system-ui,sans-serif'; ctx.fillStyle = '#fff'; ctx.textAlign = 'left';
  ctx.fillText('Great Lakes', mx + 9, my + 4);
  ctx.textAlign = 'center'; ctx.fillStyle = 'rgba(255,255,255,.7)'; ctx.font = '10px system-ui,sans-serif';
  ctx.fillText('N. America', proj(-100, 58)[0], proj(-100, 58)[1]);
  ctx.fillText('Eurasia', proj(80, 58)[0], proj(80, 58)[1]);
  ctx.fillText('Greenland', proj(-40, 74)[0], proj(-40, 74)[1] + 4);

  // stats
  const sea = -Math.round(120 * Math.min(V, 1));
  const marginLaur = Math.round(90 - V * reach(-95));
  const marginFen = Math.round(90 - V * reach(15));
  $('#glStats').innerHTML = `
    <div class="stat"><span class="k">Ice volume</span><span class="v">${Math.round(Math.min(V, 1) * 100)}%</span><span class="s">of glacial maximum</span></div>
    <div class="stat"><span class="k">Sea level</span><span class="v">${sea === 0 ? '≈ today' : sea + ' m'}</span><span class="s">water locked in ice</span></div>
    <div class="stat"><span class="k">N. American ice edge</span><span class="v">${V < 0.05 ? '—' : '~' + marginLaur + '°N'}</span><span class="s">Laurentide sheet</span></div>
    <div class="stat"><span class="k">European ice edge</span><span class="v">${V < 0.05 ? '—' : '~' + marginFen + '°N'}</span><span class="s">Scandinavian sheet</span></div>`;
  const co = $('#glCallout');
  co.className = 'callout' + (covered ? ' ice' : '');
  co.textContent = covered
    ? `❄ The Great Lakes region (~44°N) would be buried under ice sheet here. The Great Lakes themselves only formed as the ice melted back.`
    : V > 0.5 ? 'The ice edge is advancing south, but the Great Lakes region is still just ice-free in this scenario.'
    : 'The Great Lakes region (~44°N) is ice-free. Only Greenland and the Arctic hold ice at this moment.';
}

/* ------------------------------------------------------------------ *
 * Readouts
 * ------------------------------------------------------------------ */
function drawReadouts() {
  const { e, o, w } = eff(state.t);
  const q = qAt(state.t, state.sc, state.lat);
  const q0 = qAt(0, { e: true, o: true, p: true }, state.lat);
  const dq = q - q0;
  $('#readouts').innerHTML = `
    <div class="stat"><span class="k">Perihelion</span><span class="v">${perihelionSeason(w)}</span><span class="s">ϖ = ${w.toFixed(0)}°</span></div>
    <div class="stat"><span class="k">Eccentricity</span><span class="v" style="color:${COL.e}">${e.toFixed(3)}</span><span class="s">today 0.017</span></div>
    <div class="stat"><span class="k">Tilt</span><span class="v" style="color:${COL.o}">${o.toFixed(2)}°</span><span class="s">today 23.44°</span></div>
    <div class="stat"><span class="k">${state.lat}°N June sun</span><span class="v" style="color:${COL.ins}">${q.toFixed(0)}</span><span class="s">${dq >= 0 ? '+' : '−'}${Math.abs(dq).toFixed(0)} W/m² vs today</span></div>`;
}

/* ------------------------------------------------------------------ *
 * Update orchestration
 * ------------------------------------------------------------------ */
const timeInput = $('#time');
function setTime(t, fromUser) {
  state.t = clamp(t, T0, T1);
  if (fromUser) { state.playing = false; $('#playBtn').textContent = '▶'; }
  refresh();
}
function refresh() {
  timeInput.value = state.t / 1000;
  $('#timeLabel').textContent = fmtYears(state.t);
  drawLanesOver(); drawOrbit(); drawTilt(); drawMap(); drawReadouts();
  $$('.chip[data-jump]').forEach(b => b.classList.toggle('active', +b.dataset.jump === state.t));
}
function scenarioChanged() {
  computeScenario();
  $$('.tog').forEach(b => { const on = state.sc[b.dataset.k]; b.classList.toggle('on', on); b.setAttribute('aria-pressed', on); });
  drawLanesBase(); refresh();
}

/* ------------------------------------------------------------------ *
 * Wiring
 * ------------------------------------------------------------------ */
function init() {
  computeStatic();
  buildLanes();
  computeScenario();

  timeInput.min = T0 / 1000; timeInput.max = T1 / 1000;
  timeInput.addEventListener('input', () => setTime(+timeInput.value * 1000, true));

  $('#playBtn').addEventListener('click', () => {
    state.playing = !state.playing;
    if (state.playing && state.t >= T1) state.t = T0;
    $('#playBtn').textContent = state.playing ? '❚❚' : '▶';
  });
  $('#speed').addEventListener('change', ev => state.speed = +ev.target.value);
  $$('.chip[data-jump]').forEach(b => b.addEventListener('click', () => setTime(+b.dataset.jump, true)));

  $$('.tog').forEach(b => b.addEventListener('click', () => { state.sc[b.dataset.k] = !state.sc[b.dataset.k]; scenarioChanged(); }));
  $$('[data-only]').forEach(b => b.addEventListener('click', () => {
    const k = b.dataset.only;
    state.sc = { e: k === 'all' || k === 'e', o: k === 'all' || k === 'o', p: k === 'all' || k === 'p' };
    scenarioChanged();
  }));

  $('#exag').addEventListener('input', ev => { state.exag = +ev.target.value; $('#exagVal').textContent = state.exag; drawOrbit(); });
  const seasonEl = $('#season');
  const seasonLabel = () => {
    const l = state.lam, names = ['Spring equinox', 'Summer solstice', 'Autumn equinox', 'Winter solstice'];
    const near = names.findIndex((_, i) => Math.abs(mod(l - i * 90 + 180, 360) - 180) < 6);
    $('#seasonName').textContent = near >= 0 ? names[near] : seasonOf(l);
  };
  seasonEl.addEventListener('input', () => { state.lam = +seasonEl.value; seasonLabel(); drawOrbit(); });
  $('#spinBtn').addEventListener('click', () => {
    state.spin = !state.spin;
    $('#spinBtn').textContent = state.spin ? '❚❚ Pause' : '▶ Spin through a year';
  });

  $$('.seg-btn').forEach(b => b.addEventListener('click', () => {
    state.iceSrc = b.dataset.src;
    $$('.seg-btn').forEach(x => x.classList.toggle('on', x === b)); drawMap();
  }));

  // tooltips
  const tips = {
    orbit: 'Looking down on the orbit from above the North Pole. Perihelion is where Earth is closest to the Sun. The orbit\'s direction rotates over ~21,000 years (precession), and its shape stretches and relaxes (eccentricity).',
    tilt: 'Earth\'s axis tilts about 23°. Over ~41,000 years that angle slowly swings between ~22° and ~24.5°. The red ring is the latitude whose summer sunlight we track.'
  };
  const tip = $('#tip');
  $$('.info').forEach(b => {
    const show = () => { const r = b.getBoundingClientRect(); tip.textContent = tips[b.dataset.tip]; tip.hidden = false; tip.style.left = clamp(r.left, 8, innerWidth - 300) + 'px'; tip.style.top = (r.bottom + 8) + 'px'; };
    b.addEventListener('mouseenter', show); b.addEventListener('focus', show); b.addEventListener('click', show);
    b.addEventListener('mouseleave', () => tip.hidden = true); b.addEventListener('blur', () => tip.hidden = true);
  });

  let rz; window.addEventListener('resize', () => { clearTimeout(rz); rz = setTimeout(() => { drawLanesBase(); refresh(); }, 80); });

  seasonLabel();
  drawLanesBase(); refresh();

  let last = performance.now();
  (function frame(ts) {
    const dt = Math.min(0.1, (ts - last) / 1000); last = ts;
    if (state.spin) { state.lam = mod(state.lam + dt * 70, 360); seasonEl.value = Math.round(state.lam) % 360; seasonLabel(); }
    if (state.playing) {
      state.t += dt * state.speed * 1000;
      if (state.t >= T1) { state.t = T1; state.playing = false; $('#playBtn').textContent = '▶'; }
      refresh();
    } else if (state.spin) { drawOrbit(); }
    requestAnimationFrame(frame);
  })(last);
}

window.MILANKOVITCH = { state, orbit, series, insolation, refIce, ICE, computeScenario, drawLanesBase, refresh };
init();
})();
