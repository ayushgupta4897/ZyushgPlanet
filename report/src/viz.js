/* ZyushgPlanet report graphics. Everything is drawn from report/data.json (embedded as window.ZP_DATA).
   Planet and star renderings are illustrations; every chart is to scale from pipeline outputs. */
(function () {
  "use strict";
  const C = { void: "#050506", ink: "#0b0c0e", panel: "#101114", line: "#25282d", steel: "#737a83", silver: "#b8bdc4", chrome: "#e3e5e8", white: "#f7f8fa" };
  const reduce = typeof matchMedia === "function" && matchMedia("(prefers-reduced-motion: reduce)").matches;
  const STILL = /[?&]still\b/.test(location.search) || window.ZP_STILL === true;
  const NS = "http://www.w3.org/2000/svg";

  /* ---------- noise ---------- */
  function hash(x, y) { let h = Math.sin(x * 127.1 + y * 311.7) * 43758.5453; return h - Math.floor(h); }
  function vnoise(x, y) {
    const xi = Math.floor(x), yi = Math.floor(y), xf = x - xi, yf = y - yi;
    const u = xf * xf * (3 - 2 * xf), v = yf * yf * (3 - 2 * yf);
    const a = hash(xi, yi), b = hash(xi + 1, yi), c = hash(xi, yi + 1), d = hash(xi + 1, yi + 1);
    return a + (b - a) * u + (c - a) * v + (a - b - c + d) * u * v;
  }
  function fbm(x, y) { let s = 0, a = 0.5, f = 1; for (let i = 0; i < 4; i++) { s += a * vnoise(x * f, y * f); f *= 2; a *= 0.5; } return s; }

  /* ---------- canvas sizing ---------- */
  function fit(canvas) {
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const r = canvas.getBoundingClientRect();
    const w = Math.max(1, Math.round(r.width)), h = Math.max(1, Math.round(r.height));
    if (canvas.width !== Math.round(w * dpr) || canvas.height !== Math.round(h * dpr)) { canvas.width = Math.round(w * dpr); canvas.height = Math.round(h * dpr); }
    const ctx = canvas.getContext("2d"); ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    return { ctx, w, h };
  }

  /* ---------- halftone dotted planet (illustration) ---------- */
  function drawPlanet(ctx, cx, cy, R, t, o) {
    o = Object.assign({ spacing: Math.max(2.6, R / 46), tilt: -0.38, light: [-0.62, -0.42, 0.66], color: "232,235,240", bands: 9, rim: true, night: 0.07 }, o || {});
    const s = o.spacing, L = o.light, Ln = Math.hypot(L[0], L[1], L[2]);
    const lx = L[0] / Ln, ly = L[1] / Ln, lz = L[2] / Ln;
    const ct = Math.cos(o.tilt), st = Math.sin(o.tilt);
    const LEV = 9, paths = Array.from({ length: LEV }, () => new Path2D());
    const rot = t * 0.12;
    for (let y = -R; y <= R; y += s) {
      const row = Math.round(y / s);
      const off = (row & 1) ? s / 2 : 0;
      for (let x = -R + off; x <= R; x += s) {
        const nx0 = x / R, ny0 = y / R, rr = nx0 * nx0 + ny0 * ny0;
        if (rr > 1) continue;
        const nz = Math.sqrt(1 - rr);
        const nx = nx0 * ct - ny0 * st, ny = nx0 * st + ny0 * ct;
        const lat = Math.asin(Math.max(-1, Math.min(1, ny)));
        const lon = Math.atan2(nx, nz) + rot;
        const band = 0.5 + 0.28 * Math.sin(lat * o.bands + 1.6 * fbm(lon * 1.3 + 3, lat * 3)) + 0.22 * (fbm(lon * 3.2, lat * 9) - 0.5);
        const lam = Math.max(0, nx0 * lx + ny0 * ly + nz * lz);
        let v = Math.pow(lam, 0.85) * (0.45 + 0.55 * band);
        if (o.rim) v += 0.22 * Math.pow(1 - nz, 3) * Math.pow(lam + 0.05, 0.3);
        v = Math.max(o.night, Math.min(1, v));
        const lev = Math.min(LEV - 1, Math.floor(v * LEV));
        const r = s * 0.5 * Math.sqrt((lev + 0.6) / LEV);
        paths[lev].moveTo(cx + x + r, cy + y); paths[lev].arc(cx + x, cy + y, r, 0, 6.2832);
      }
    }
    for (let i = 0; i < LEV; i++) { ctx.fillStyle = `rgba(${o.color},${(0.18 + 0.82 * i / (LEV - 1)).toFixed(3)})`; ctx.fill(paths[i]); }
    // thin atmosphere halo of dots
    const halo = new Path2D(), n = Math.round(R * 1.2);
    for (let i = 0; i < n; i++) {
      const a = (i / n) * 6.2832, d = R * (1.035 + 0.03 * hash(i, 7));
      const facing = Math.max(0, -(Math.cos(a) * lx + Math.sin(a) * ly) * -1);
      if (facing < 0.05) continue;
      const r = 0.6 + 1.0 * facing;
      halo.moveTo(cx + Math.cos(a) * d + r, cy + Math.sin(a) * d); halo.arc(cx + Math.cos(a) * d, cy + Math.sin(a) * d, r, 0, 6.2832);
    }
    ctx.fillStyle = `rgba(${o.color},0.35)`; ctx.fill(halo);
  }

  /* ---------- limb-darkened dotted star (u from the fitted model) ---------- */
  function drawStar(ctx, cx, cy, R, u, o) {
    o = Object.assign({ spacing: 7, occ: null, t: 0, clip: null, color: "247,248,250" }, o || {});
    const s = o.spacing, LEV = 8, paths = Array.from({ length: LEV }, () => new Path2D());
    const x0 = o.clip ? o.clip[0] : cx - R, x1 = o.clip ? o.clip[1] : cx + R, y0 = o.clip ? o.clip[2] : cy - R, y1 = o.clip ? o.clip[3] : cy + R;
    for (let y = Math.max(cy - R, y0); y <= Math.min(cy + R, y1); y += s) {
      const row = Math.round(y / s), off = (row & 1) ? s / 2 : 0;
      for (let x = Math.max(cx - R, x0) + off; x <= Math.min(cx + R, x1); x += s) {
        const dx = (x - cx) / R, dy = (y - cy) / R, rr = dx * dx + dy * dy;
        if (rr > 1) continue;
        if (o.occ) { const ox = x - o.occ[0], oy = y - o.occ[1]; if (ox * ox + oy * oy < o.occ[2] * o.occ[2]) continue; }
        const mu = Math.sqrt(1 - rr), m1 = 1 - mu;
        let I = 1 - u[0] * m1 - u[1] * m1 * m1;
        I *= 0.95 + 0.05 * vnoise(x * 0.11 + o.t * 0.15, y * 0.11);
        const lev = Math.max(0, Math.min(LEV - 1, Math.floor(I * LEV)));
        const r = s * 0.48 * Math.sqrt((lev + 1) / LEV);
        paths[lev].moveTo(x + r, y); paths[lev].arc(x, y, r, 0, 6.2832);
      }
    }
    for (let i = 0; i < LEV; i++) { ctx.fillStyle = `rgba(${o.color},${(0.15 + 0.85 * i / (LEV - 1)).toFixed(3)})`; ctx.fill(paths[i]); }
  }

  /* ---------- drifting dot field ---------- */
  function dotField(canvas, o) {
    o = Object.assign({ density: 0.00022, speed: 0.012 }, o || {});
    let pts = [], w = 0, h = 0, ctx;
    function init() { const f = fit(canvas); ctx = f.ctx; w = f.w; h = f.h; const n = Math.round(w * h * o.density); pts = Array.from({ length: n }, (_, i) => ({ x: hash(i, 1) * w, y: hash(i, 2) * h, z: 0.25 + 0.75 * hash(i, 3), p: hash(i, 4) * 6.28 })); }
    function frame(t) {
      ctx.clearRect(0, 0, w, h);
      for (const p of pts) {
        const x = (p.x + t * o.speed * 60 * p.z) % w, tw = 0.45 + 0.55 * Math.sin(p.p + t * 0.9 * p.z);
        ctx.fillStyle = `rgba(214,218,224,${(0.12 + 0.5 * p.z * tw).toFixed(3)})`;
        ctx.fillRect(x, p.y, p.z * 1.4, p.z * 1.4);
      }
    }
    init(); window.addEventListener("resize", init);
    return frame;
  }

  /* ---------- animation loop with visibility ---------- */
  const loops = [];
  function loop(el, fn) { loops.push({ el, fn, vis: true }); fn(0); }
  function start() {
    if (STILL) return;
    if ("IntersectionObserver" in window) {
      const io = new IntersectionObserver(es => es.forEach(e => { const L = loops.find(l => l.el === e.target); if (L) L.vis = e.isIntersecting; }), { rootMargin: "100px" });
      loops.forEach(l => io.observe(l.el));
    }
    if (reduce) return;
    const t0 = performance.now();
    function tick(now) { const t = (now - t0) / 1000; for (const l of loops) if (l.vis) l.fn(t); requestAnimationFrame(tick); }
    requestAnimationFrame(tick);
  }

  /* ---------- SVG helpers ---------- */
  function el(tag, attrs, parent, text) {
    const e = document.createElementNS(NS, tag);
    for (const k in attrs || {}) e.setAttribute(k, attrs[k]);
    if (text !== undefined) e.textContent = text;
    if (parent) parent.appendChild(e);
    return e;
  }
  function svgRoot(host, W, H, label) {
    host.innerHTML = "";
    const s = el("svg", { viewBox: `0 0 ${W} ${H}`, width: "100%", role: "img", "aria-label": label || "", preserveAspectRatio: "xMidYMid meet" }, host);
    s.style.display = "block"; s.style.maxWidth = "100%"; s.style.height = "auto";
    return s;
  }
  const lin = (d0, d1, r0, r1) => v => r0 + (v - d0) / (d1 - d0) * (r1 - r0);
  const logs = (d0, d1, r0, r1) => v => r0 + (Math.log10(v) - Math.log10(d0)) / (Math.log10(d1) - Math.log10(d0)) * (r1 - r0);
  function ticks(a, b, n) {
    const span = b - a, step0 = span / n, mag = Math.pow(10, Math.floor(Math.log10(step0))), r = step0 / mag;
    const step = (r < 1.5 ? 1 : r < 3 ? 2 : r < 7 ? 5 : 10) * mag;
    const out = []; for (let v = Math.ceil(a / step) * step; v <= b + 1e-9; v += step) out.push(+v.toFixed(10)); return out;
  }
  const fmt = (v, d) => v.toLocaleString("en-US", { minimumFractionDigits: d || 0, maximumFractionDigits: d || 0 });
  function rearth(t, prefix) { el("tspan", {}, t, prefix + " R"); el("tspan", { "baseline-shift": "sub", "font-size": "70%" }, t, "E"); return t; }
  function axes(s, box, xs, ys, xt, yt, o) {
    o = o || {};
    const g = el("g", { class: "axes" }, s);
    for (const v of yt) { const y = ys(v); el("line", { x1: box.x0, x2: box.x1, y1: y, y2: y, class: "grid" }, g); el("text", { x: box.x0 - 10, y: y + 4, "text-anchor": "end", class: "tick" }, g, o.yfmt ? o.yfmt(v) : fmt(v)); }
    for (const v of xt) { const x = xs(v); el("line", { x1: x, x2: x, y1: box.y1, y2: box.y1 + 6, class: "axis" }, g); el("text", { x, y: box.y1 + 22, "text-anchor": "middle", class: "tick" }, g, o.xfmt ? o.xfmt(v) : fmt(v, o.xd || 0)); }
    el("line", { x1: box.x0, x2: box.x1, y1: box.y1, y2: box.y1, class: "axis" }, g);
    if (o.xlabel) el("text", { x: (box.x0 + box.x1) / 2, y: box.y1 + 46, "text-anchor": "middle", class: "label" }, g, o.xlabel);
    if (o.ylabel) el("text", { x: 0, y: 0, transform: `translate(${box.x0 - 58},${(box.y0 + box.y1) / 2}) rotate(-90)`, "text-anchor": "middle", class: "label" }, g, o.ylabel);
    return g;
  }

  /* ---------- charts ---------- */
  function lightCurve(host, c, o) {
    o = Object.assign({ W: 900, H: 460, pts: true, title: "" }, o || {});
    const W = o.W, H = o.H, box = { x0: 96, x1: W - 24, y0: 24, y1: H - 70 };
    const s = svgRoot(host, W, H, `Phase-folded light curve of ${c.fname}`);
    const xr = Math.max(...c.lc_binned.map(p => Math.abs(p[0])));
    const d = c.depth_ppm, ylo = -d * 2.1, yhi = d * 1.0;
    const xs = lin(-xr, xr, box.x0, box.x1), ys = lin(yhi, ylo, box.y0, box.y1);
    axes(s, box, xs, ys, ticks(-xr, xr, 6), ticks(ylo, yhi, 5), { xlabel: "hours from mid-transit", ylabel: "flux change (ppm)", xd: 0 });
    if (o.pts) { const g = el("g", { class: "raw" }, s); for (const p of c.lc_points) { const y = ys(p[1]); if (y < box.y0 || y > box.y1) continue; el("circle", { cx: xs(p[0]).toFixed(1), cy: y.toFixed(1), r: 1.15 }, g); } }
    const m = c.lc_model.map(p => `${xs(p[0]).toFixed(1)},${ys(p[1]).toFixed(1)}`).join(" ");
    el("polyline", { points: m, class: "model" }, s);
    const gb = el("g", { class: "bin" }, s);
    for (const p of c.lc_binned) { el("line", { x1: xs(p[0]), x2: xs(p[0]), y1: ys(p[1] - p[2]), y2: ys(p[1] + p[2]), class: "err" }, gb); el("circle", { cx: xs(p[0]), cy: ys(p[1]), r: 3.4 }, gb); }
    el("line", { x1: box.x0, x2: box.x1, y1: ys(-d), y2: ys(-d), class: "ref" }, s);
    el("text", { x: box.x1 - 4, y: ys(-d) - 8, "text-anchor": "end", class: "ann" }, s, `depth ${fmt(d)} ppm`);
    return s;
  }

  function periodogram(host, c, o) {
    o = Object.assign({ W: 620, H: 330 }, o || {});
    const W = o.W, H = o.H, box = { x0: 70, x1: W - 20, y0: 34, y1: H - 62 };
    const s = svgRoot(host, W, H, `Box-least-squares periodogram, real and inverted flux, ${c.fname}`);
    const ymax = Math.ceil(Math.max(...c.periodogram_real.map(p => p[1] || 0)) / 5) * 5 + 2;
    const xs = logs(1, 50, box.x0, box.x1), ys = lin(-2, ymax, box.y1, box.y0);
    axes(s, box, xs, ys, [1, 2, 5, 10, 20, 50], ticks(0, ymax, 4), { xlabel: "trial period (days, log scale)", ylabel: "SDE", xfmt: v => v });
    const path = arr => arr.filter(p => p[0] && p[1] !== null).map((p, i) => `${i ? "L" : "M"}${xs(p[0]).toFixed(1)},${ys(Math.max(-2, p[1])).toFixed(1)}`).join("");
    el("path", { d: path(c.periodogram_inverted), class: "pg-inv" }, s);
    el("path", { d: path(c.periodogram_real), class: "pg-real" }, s);
    const pk = c.periodogram_real.reduce((a, b) => (b[1] > a[1] ? b : a));
    el("circle", { cx: xs(pk[0]), cy: ys(pk[1]), r: 5, class: "pk" }, s);
    el("text", { x: xs(pk[0]) + 10, y: ys(pk[1]) + 4, class: "ann" }, s, `${pk[0].toFixed(3)} d · SDE ${pk[1].toFixed(1)}`);
    const ik = c.periodogram_inverted.reduce((a, b) => (b[1] > a[1] ? b : a));
    el("line", { x1: box.x0 + 4, x2: box.x0 + 26, y1: box.y0 - 18, y2: box.y0 - 18, class: "pg-inv" }, s);
    el("text", { x: box.x0 + 32, y: box.y0 - 14, class: "ann dim" }, s, `inverted flux (2 − f): highest peak SDE ${ik[1].toFixed(1)}`);
    return s;
  }

  function eventStrip(host, c, o) {
    o = Object.assign({ W: 620, H: 330 }, o || {});
    const W = o.W, H = o.H, box = { x0: 70, x1: W - 20, y0: 20, y1: H - 62 };
    const s = svgRoot(host, W, H, `Depth of every individual transit, ${c.fname}`);
    const ev = c.events, t0 = 120, t1 = 1600;
    const vals = ev.map(e => e[2]), hi = Math.max(c.depth_ppm * 3.2, ...vals.map(v => Math.min(v, c.depth_ppm * 3.2))), lo = Math.min(-c.depth_ppm * 1.2, 0);
    const xs = lin(t0, t1, box.x0, box.x1), ys = lin(lo, hi, box.y1, box.y0);
    axes(s, box, xs, ys, [200, 500, 800, 1100, 1400], ticks(lo, hi, 4), { xlabel: "time (BJD − 2,454,833, days)", ylabel: "depth (ppm)" });
    el("text", { x: box.x1, y: ys(c.depth_ppm) - 8, "text-anchor": "end", class: "ann" }, s, `mean ${fmt(c.depth_ppm)} ppm · ${c.events.length} transits`);
    el("line", { x1: box.x0, x2: box.x1, y1: ys(c.depth_ppm), y2: ys(c.depth_ppm), class: "ref" }, s);
    el("line", { x1: box.x0, x2: box.x1, y1: ys(0), y2: ys(0), class: "zero" }, s);
    const g = el("g", { class: "ev" }, s);
    for (const e of ev) { const y = ys(Math.max(lo, Math.min(hi, e[2]))); el("line", { x1: xs(e[1]), x2: xs(e[1]), y1: ys(Math.max(lo, e[2] - e[3])), y2: ys(Math.min(hi, e[2] + e[3])), class: "err" }, g); el("circle", { cx: xs(e[1]), cy: y, r: 2.6 }, g); }
    return s;
  }

  function fppChart(host, D, o) {
    o = Object.assign({ W: 900, H: 380, big: false }, o || {});
    const W = o.W, H = o.H, box = { x0: 210, x1: W - 40, y0: 26, y1: H - 64 };
    const s = svgRoot(host, W, H, "False-positive probability before and after Gaia companion constraints");
    const names = Object.keys(D.candidates), xmax = 4.5;
    const xs = lin(0, xmax, box.x0, box.x1), rowH = (box.y1 - box.y0) / names.length;
    axes(s, box, xs, () => 0, [0, 1, 2, 3, 4], [], { xlabel: "false-positive probability (%)", xfmt: v => v + " %" });
    el("line", { x1: xs(1.5), x2: xs(1.5), y1: box.y0 - 10, y2: box.y1, class: "thresh" }, s);
    el("text", { x: xs(1.5) + 8, y: box.y0 - 2, class: "ann" }, s, "1.5 % validation threshold");
    names.forEach((n, i) => {
      const c = D.candidates[n], y = box.y0 + rowH * (i + 0.5);
      el("text", { x: box.x0 - 18, y: y - 4, "text-anchor": "end", class: "rowlab" }, s, "KIC " + c.kic);
      rearth(el("text", { x: box.x0 - 18, y: y + 16, "text-anchor": "end", class: "tick" }, s), `${c.period_d.toFixed(2)} d · ${c.rp_earth.median.toFixed(2)}`);
      const b = c.fpp.baseline.mean * 100, m = c.fpp.molusc.mean * 100;
      el("line", { x1: xs(m), x2: xs(b), y1: y, y2: y, class: "arrow" }, s);
      for (const v of c.fpp.molusc.seeds) el("circle", { cx: xs(v * 100), cy: y + 14, r: 2.2, class: "seed" }, s);
      el("circle", { cx: xs(b), cy: y, r: 7, class: "before" }, s);
      el("circle", { cx: xs(m), cy: y, r: 9, class: "after" }, s);
      el("text", { x: xs(b) + 14, y: y + 5, class: "ann dim" }, s, `${b.toFixed(1)} %`);
      el("text", { x: xs(m) - 14, y: y - 12, "text-anchor": "end", class: "val" }, s, `${m.toFixed(2)} %`);
    });
    return s;
  }

  function densityChart(host, D, o) {
    o = Object.assign({ W: 900, H: 300 }, o || {});
    const W = o.W, H = o.H, box = { x0: 210, x1: W - 40, y0: 20, y1: H - 62 };
    const s = svgRoot(host, W, H, "Stellar density from the transit shape versus Gaia");
    const names = Object.keys(D.candidates), xs = logs(0.2, 10, box.x0, box.x1), rowH = (box.y1 - box.y0) / names.length;
    axes(s, box, xs, () => 0, [0.2, 0.5, 1, 2, 5, 10], [], { xlabel: "stellar density (g/cm³, log scale)", xfmt: v => v });
    names.forEach((n, i) => {
      const c = D.candidates[n], f = c.fit.rho_star_gcc, y = box.y0 + rowH * (i + 0.5);
      el("text", { x: box.x0 - 18, y: y + 5, "text-anchor": "end", class: "rowlab" }, s, "KIC " + c.kic);
      el("line", { x1: xs(Math.max(0.2, f.median - f.minus)), x2: xs(Math.min(10, f.median + f.plus)), y1: y - 8, y2: y - 8, class: "ci" }, s);
      el("circle", { cx: xs(f.median), cy: y - 8, r: 5, class: "after" }, s);
      el("rect", { x: xs(c.density.gaia) - 1.5, y: y + 2, width: 3, height: 14, class: "gaia" }, s);
      el("text", { x: xs(c.density.gaia) + 8, y: y + 14, class: "ann dim" }, s, `Gaia ${c.density.gaia.toFixed(2)}`);
      el("text", { x: xs(Math.min(10, f.median + f.plus)) + 8, y: y - 4, class: "ann" }, s, `${c.density.tension_sigma > 0 ? "+" : ""}${c.density.tension_sigma.toFixed(1)} σ`);
    });
    return s;
  }

  function contrastChart(host, D, o) {
    o = Object.assign({ W: 620, H: 360 }, o || {});
    const W = o.W, H = o.H, box = { x0: 70, x1: W - 30, y0: 56, y1: H - 62 };
    const s = svgRoot(host, W, H, "Gaia DR3 empirical contrast sensitivity");
    const xs = lin(0, 4, box.x0, box.x1), ys = lin(0, 8.5, box.y0, box.y1);
    axes(s, box, xs, ys, [0, 1, 2, 3, 4], [0, 2, 4, 6, 8], { xlabel: "separation from the star (arcsec)", ylabel: "Δ G (mag)", xfmt: v => v + "″" });
    const names = Object.keys(D.candidates), dash = ["", "6 5", "2 5"];
    names.forEach((n, i) => {
      const cc = D.candidates[n].gaia_contrast_curve; if (!cc) return;
      const d = cc.map((p, k) => `${k ? "L" : "M"}${xs(p[0]).toFixed(1)},${ys(p[1]).toFixed(1)}`).join("");
      el("path", { d, class: "cc", "stroke-dasharray": dash[i] }, s);
      const lx = box.x0 + 10 + i * 175;
      el("line", { x1: lx, x2: lx + 26, y1: 18, y2: 18, class: "cc", "stroke-dasharray": dash[i] }, s);
      el("text", { x: lx + 32, y: 22, class: "ann dim" }, s, "KIC " + D.candidates[n].kic);
    });
    el("text", { x: xs(2.1), y: ys(1.6), class: "ann" }, s, "brighter companions here");
    el("text", { x: xs(2.1), y: ys(1.6) + 16, class: "ann" }, s, "would be in Gaia DR3");
    return s;
  }

  function pixelMap(host, m, o) {
    o = Object.assign({ cell: 46, mode: "snr", label: "" }, o || {});
    const H = m.image.length, W = m.image[0].length, cs = o.cell, pad = 2;
    const s = svgRoot(host, W * cs + pad * 2, H * cs + pad * 2, o.label);
    const g = el("g", {}, s);
    let lo = Infinity, hi = -Infinity;
    for (const r of m.image) for (const v of r) if (v > 0) { lo = Math.min(lo, Math.log10(v)); hi = Math.max(hi, Math.log10(v)); }
    for (let y = 0; y < H; y++) for (let x = 0; x < W; x++) {
      const yy = H - 1 - y; let fill;
      if (o.mode === "image") { const v = m.image[y][x]; const t = v > 0 ? (Math.log10(v) - lo) / (hi - lo + 1e-9) : 0; const c = Math.round(18 + 220 * Math.pow(t, 0.9)); fill = `rgb(${c},${c + 2},${c + 5})`; }
      else { const v = Math.max(-5, Math.min(6, m.snr[y][x])); const t = (v + 5) / 11; const c = Math.round(10 + 240 * Math.pow(t, 1.6)); fill = `rgb(${c},${c + 1},${c + 4})`; }
      el("rect", { x: pad + x * cs + 1, y: pad + yy * cs + 1, width: cs - 2, height: cs - 2, fill, rx: 2 }, g);
      if (o.mode === "snr" && o.values) el("text", { x: pad + x * cs + cs / 2, y: pad + yy * cs + cs / 2 + 4, "text-anchor": "middle", class: "cellv" }, g, m.snr[y][x].toFixed(1));
    }
    const P = (p) => [pad + (p[0] + 0.5) * cs, pad + (H - 1 - p[1] + 0.5) * cs];
    for (const nb of m.neighbours) { if (Math.hypot(nb[0] - m.target[0], nb[1] - m.target[1]) < 0.3) continue; const [x, y] = P(nb); if (x < -cs || y < -cs || x > W * cs + cs || y > H * cs + cs) continue; el("circle", { cx: x, cy: y, r: cs * 0.22, class: "nbmark" }, s); }
    const [tx, ty] = P(m.target);
    el("circle", { cx: tx, cy: ty, r: cs * 0.42, class: "tmark" }, s);
    el("line", { x1: tx - cs * 0.7, x2: tx - cs * 0.5, y1: ty, y2: ty, class: "tmark" }, s); el("line", { x1: tx + cs * 0.5, x2: tx + cs * 0.7, y1: ty, y2: ty, class: "tmark" }, s);
    el("line", { x1: tx, x2: tx, y1: ty - cs * 0.7, y2: ty - cs * 0.5, class: "tmark" }, s); el("line", { x1: tx, x2: tx, y1: ty + cs * 0.5, y2: ty + cs * 0.7, class: "tmark" }, s);
    return s;
  }

  /* ---------- funnel dot field ---------- */
  function funnel(canvas, D, o) {
    o = Object.assign({ stage: null }, o || {});
    const F = D.funnel;
    const stages = [
      { n: F.targets, label: "stars searched" },
      { n: F.first_screen_real, label: "pass the first screen" },
      { n: F.path_stage2, label: "pass stage-2 vetting" },
      { n: F.path_pixel_on_target, label: "lose light on the target's own pixels" },
      { n: F.path_redteam_survives, label: "survive re-derivation and red team" },
      { n: F.candidates, label: "candidates with FPP below 1 %" }];
    const N = F.targets, names = Object.keys(D.candidates).map(k => "KIC " + D.candidates[k].kic);
    let ctx, w, h, cols, rows, sp, rank;
    function init() {
      const f = fit(canvas); ctx = f.ctx; w = f.w; h = f.h; sp = Math.sqrt((w * h) / (N * 1.05)); cols = Math.floor(w / sp); rows = Math.ceil(N / cols);
      const spots = [[0.27, 0.36], [0.56, 0.62], [0.78, 0.3]].map(([fx, fy]) => Math.floor(fy * rows) * cols + Math.floor(fx * cols));
      const order = Array.from({ length: N }, (_, i) => i).filter(i => !spots.includes(i)).sort((a, b) => hash(a, 11) - hash(b, 11));
      rank = new Int32Array(N); spots.forEach((v, i) => (rank[v] = i)); order.forEach((v, i) => (rank[v] = i + spots.length));
    }
    init(); window.addEventListener("resize", init);
    function draw(t) {
      const cyc = o.stage !== null ? o.stage + 0.999 : ((t / 2.2) % (stages.length + 1.5));
      const si = Math.min(stages.length - 1, Math.floor(cyc)), fr = Math.min(1, (cyc - si) * 1.6);
      ctx.clearRect(0, 0, w, h);
      const keep = stages[si].n, prev = si ? stages[si - 1].n : N;
      const p1 = new Path2D(), p2 = new Path2D(), p3 = new Path2D();
      for (let i = 0; i < N; i++) {
        const x = (i % cols) * sp + sp / 2, y = Math.floor(i / cols) * sp + sp / 2, r = rank[i];
        if (r < keep) { const rr = si === stages.length - 1 ? sp * 1.1 : sp * 0.36; p1.moveTo(x + rr, y); p1.arc(x, y, rr, 0, 6.2832); }
        else if (r < prev) { const rr = sp * 0.26 * (1 - fr) + sp * 0.14; p2.moveTo(x + rr, y); p2.arc(x, y, rr, 0, 6.2832); }
        else { const rr = sp * 0.13; p3.moveTo(x + rr, y); p3.arc(x, y, rr, 0, 6.2832); }
      }
      ctx.fillStyle = "rgba(140,147,156,0.6)"; ctx.fill(p3);
      ctx.fillStyle = `rgba(184,189,196,${(0.6 * (1 - fr) + 0.2).toFixed(3)})`; ctx.fill(p2);
      ctx.fillStyle = "rgba(247,248,250,0.95)"; ctx.fill(p1);
      if (si === stages.length - 1) {
        ctx.font = "500 12px 'IBM Plex Mono', ui-monospace, monospace"; ctx.fillStyle = "rgba(227,229,232,0.95)"; ctx.strokeStyle = "rgba(247,248,250,0.7)";
        for (let i = 0; i < N; i++) if (rank[i] < 3) { const x = (i % cols) * sp + sp / 2, y = Math.floor(i / cols) * sp + sp / 2;
          ctx.beginPath(); ctx.arc(x, y, sp * 3, 0, 6.2832); ctx.setLineDash([2, 4]); ctx.stroke(); ctx.setLineDash([]);
          ctx.fillText(names[rank[i]], x + sp * 3.6, y + 4); }
      }
      if (o.onStage) o.onStage(si);
    }
    return { draw, stages };
  }

  /* ---------- transit scene: star limb + planet + model light curve ---------- */
  function transitScene(canvas, c, o) {
    o = Object.assign({ phase: null }, o || {});
    const u = c.limb_darkening_u, model = c.lc_model, tmin = model[0][0], tmax = model[model.length - 1][0];
    let ctx, w, h;
    function init() { const f = fit(canvas); ctx = f.ctx; w = f.w; h = f.h; }
    init(); window.addEventListener("resize", init);
    return function (t) {
      ctx.clearRect(0, 0, w, h);
      const R = Math.max(w, h) * 0.9, cx = w * 0.5, cy = h * 0.58 + R * 0.62;
      const ph = o.phase !== null ? o.phase : ((t * 0.07) % 1.25) - 0.1;
      const tt = tmin + Math.max(0, Math.min(1, ph)) * (tmax - tmin);
      const chordY = cy - R * 0.83, px = w * 0.06 + Math.max(0, Math.min(1, ph)) * w * 0.88, pr = Math.max(14, w * 0.03);
      drawStar(ctx, cx, cy, R, u, { spacing: Math.max(6, w / 150), occ: [px, chordY, pr * 1.06], t, clip: [0, w, 0, h * 0.62] });
      drawPlanet(ctx, px, chordY, pr, t * 2, { spacing: Math.max(1.8, pr / 9), light: [0.0, 0.8, -0.6], night: 0.03, color: "40,42,46" });
      ctx.strokeStyle = "rgba(247,248,250,0.55)"; ctx.lineWidth = 1; ctx.beginPath(); ctx.arc(px, chordY, pr * 1.03, 0, 6.2832); ctx.stroke();
      const g0 = h * 0.7, g1 = h * 0.94, x0 = w * 0.06, x1 = w * 0.94, dmin = Math.min(...model.map(p => p[1]));
      const ys = v => g0 + (v / dmin) * (g1 - g0) * 0.9 + (g1 - g0) * 0.05;
      ctx.fillStyle = "rgba(115,122,131,0.6)";
      for (let i = 0; i < model.length; i++) { const x = x0 + (model[i][0] - tmin) / (tmax - tmin) * (x1 - x0); ctx.fillRect(x - 0.7, ys(model[i][1]) - 0.7, 1.4, 1.4); }
      ctx.fillStyle = "rgba(247,248,250,1)";
      for (let i = 0; i < model.length; i++) { if (model[i][0] > tt) break; const x = x0 + (model[i][0] - tmin) / (tmax - tmin) * (x1 - x0); ctx.beginPath(); ctx.arc(x, ys(model[i][1]), 2.1, 0, 6.2832); ctx.fill(); }
      ctx.fillStyle = "rgba(184,189,196,0.9)"; ctx.font = "500 11px 'IBM Plex Mono', ui-monospace, monospace";
      ctx.fillText(`t = ${tt >= 0 ? "+" : "−"}${Math.abs(tt).toFixed(2)} h`, x0, g1 + 16);
      const cur = model.find(p => p[0] >= tt) || model[model.length - 1];
      ctx.textAlign = "right"; ctx.fillText(`${cur[1] <= 0 ? "−" : "+"}${fmt(Math.abs(cur[1]))} ppm`, x1, g1 + 16); ctx.textAlign = "left";
    };
  }

  window.ZP = { C, fit, drawPlanet, drawStar, dotField, loop, start, lightCurve, periodogram, eventStrip, fppChart, densityChart, contrastChart, pixelMap, funnel, transitScene, fmt, STILL, reduce, hash };
})();
