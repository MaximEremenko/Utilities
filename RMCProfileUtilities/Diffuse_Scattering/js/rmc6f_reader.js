"use strict";

(function (global) {
  const PI2 = 2 * Math.PI;

  class RMC6fParser {
    constructor() {
      this.header = [];
      this.rows = [];
      this.super = [1, 1, 1];
      this.cellDeg = [NaN, NaN, NaN, NaN, NaN, NaN];
      this.cellRad = [NaN, NaN, NaN, NaN, NaN, NaN];
      this.cols = 0;
      this.names = [];
    }

    parse(txt) {
      const lines = String(txt || "")
        .split(/\r?\n/)
        .map((line) => line.replace(/\x00/g, ""));
      const start = this.findHeader(lines);
      this.header = lines.slice(0, start);
      this.rows = this.parseRows(lines.slice(start));
      this.getSuper();
      this.getCell();
    }

    findHeader(lines) {
      for (let i = 0; i < Math.min(lines.length, 151); i++) {
        if (lines[i].includes("Atoms:")) return i + 1;
      }
      return 0;
    }

    parseRows(lines) {
      const out = [];
      for (const raw of lines) {
        const line = raw.trim();
        if (!line) continue;
        const p = line.split(/\s+/);
        if (!this.cols) {
          this.cols = p.length;
          if (this.cols === 10) {
            this.names = [
              "atomNumber",
              "element",
              "id",
              "x",
              "y",
              "z",
              "refNumber",
              "cellRefNumX",
              "cellRefNumY",
              "cellRefNumZ",
            ];
          } else if (this.cols === 9) {
            this.names = [
              "atomNumber",
              "element",
              "x",
              "y",
              "z",
              "refNumber",
              "cellRefNumX",
              "cellRefNumY",
              "cellRefNumZ",
            ];
          } else {
            throw new Error(`Unsupported RMC6f: ${this.cols} columns`);
          }
        }
        const row = {};
        for (let i = 0; i < this.cols; i++) {
          const name = this.names[i] || `c${i}`;
          row[name] = [
            "atomNumber",
            "x",
            "y",
            "z",
            "refNumber",
            "cellRefNumX",
            "cellRefNumY",
            "cellRefNumZ",
          ].includes(name)
            ? Number(p[i])
            : p[i];
        }
        out.push(row);
      }
      return out;
    }

    getSuper() {
      const line = this.header.find((v) => v.includes("Supercell"));
      if (!line) {
        this.super = [1, 1, 1];
        return;
      }
      const nums = line
        .split(/\s+/)
        .map((v) => parseInt(v, 10))
        .filter(Number.isInteger);
      this.super = nums.length >= 3 ? nums.slice(0, 3) : [1, 1, 1];
    }

    getCell() {
      const line =
        this.header.find((v) => v.includes("Cell (Ang/deg)")) ||
        this.header.find((v) => /^\s*Cell\b/i.test(v));
      if (!line) return;
      const values = line
        .split(/\s+/)
        .map(Number)
        .filter(Number.isFinite);
      if (values.length >= 6) {
        this.cellDeg = values.slice(0, 6);
        this.cellRad = [
          values[0],
          values[1],
          values[2],
          (values[3] * Math.PI) / 180,
          (values[4] * Math.PI) / 180,
          (values[5] * Math.PI) / 180,
        ];
      }
    }
  }

  function cross(a, b) {
    return [
      a[1] * b[2] - a[2] * b[1],
      a[2] * b[0] - a[0] * b[2],
      a[0] * b[1] - a[1] * b[0],
    ];
  }

  function dot(a, b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  }

  function col3(m, j) {
    return [m[0][j], m[1][j], m[2][j]];
  }

  function inv3(m) {
    const a = m[0][0];
    const b = m[0][1];
    const c = m[0][2];
    const d = m[1][0];
    const e = m[1][1];
    const f = m[1][2];
    const g = m[2][0];
    const h = m[2][1];
    const i = m[2][2];
    const A = e * i - f * h;
    const B = -(d * i - f * g);
    const C = d * h - e * g;
    const D = -(b * i - c * h);
    const E = a * i - c * g;
    const F = -(a * h - b * g);
    const G = b * f - c * e;
    const H = -(a * f - c * d);
    const I = a * e - b * d;
    const det = a * A + b * B + c * C;
    if (!Number.isFinite(det) || Math.abs(det) < 1e-20) {
      throw new Error("Singular matrix");
    }
    const s = 1 / det;
    return [
      [A * s, D * s, G * s],
      [B * s, E * s, H * s],
      [C * s, F * s, I * s],
    ];
  }

  function mul(v, m) {
    return [
      v[0] * m[0][0] + v[1] * m[1][0] + v[2] * m[2][0],
      v[0] * m[0][1] + v[1] * m[1][1] + v[2] * m[2][1],
      v[0] * m[0][2] + v[1] * m[1][2] + v[2] * m[2][2],
    ];
  }

  function scaleCols(m, s) {
    return [
      [m[0][0] * s[0], m[0][1] * s[1], m[0][2] * s[2]],
      [m[1][0] * s[0], m[1][1] * s[1], m[1][2] * s[2]],
      [m[2][0] * s[0], m[2][1] * s[1], m[2][2] * s[2]],
    ];
  }

  function scaleMat(m, s) {
    return [
      [m[0][0] * s, m[0][1] * s, m[0][2] * s],
      [m[1][0] * s, m[1][1] * s, m[1][2] * s],
      [m[2][0] * s, m[2][1] * s, m[2][2] * s],
    ];
  }

  function cell2vec(a, b, c, al, be, ga) {
    const v = [
      [0, 0, 0],
      [0, 0, 0],
      [0, 0, 0],
    ];
    v[2][0] = a * Math.cos(be);
    v[1][1] = b * Math.sin(al);
    v[2][1] = b * Math.cos(al);
    v[2][2] = c;
    v[1][0] = (a * b * Math.cos(ga) - v[2][0] * v[2][1]) / v[1][1];
    v[0][0] = Math.sqrt(
      Math.max(0, a * a - v[1][0] * v[1][0] - v[2][0] * v[2][0])
    );
    return [
      [v[0][0], v[1][0], v[2][0]],
      [v[0][1], v[1][1], v[2][1]],
      [v[0][2], v[1][2], v[2][2]],
    ];
  }

  function vec2space(v) {
    const av = col3(v, 0);
    const bv = col3(v, 1);
    const cv = col3(v, 2);
    const den = dot(av, cross(bv, cv));
    if (!Number.isFinite(den) || Math.abs(den) < 1e-20) {
      throw new Error("Invalid reciprocal basis");
    }
    const ap = cross(bv, cv).map((x) => x / den);
    const bp = cross(cv, av).map((x) => x / den);
    const cp = cross(av, bv).map((x) => x / den);
    const Bp = [ap, bp, cp];
    return { B: inv3(Bp), Bp };
  }

  function wrapDelta(x) {
    if (x < -0.5) return x + 1;
    if (x > 0.5) return x - 1;
    return x;
  }

  function parseRmc6f(name, text) {
    const r = new RMC6fParser();
    r.parse(text);
    if (!r.rows.length) throw new Error("No atom rows found");
    if (!r.cellRad.every(Number.isFinite)) throw new Error("Cell parameters missing");

    const [sx, sy, sz] = r.super;
    const [a, b, c, al, be, ga] = r.cellRad;
    const vpc = cell2vec(a / sx, b / sy, c / sz, al, be, ga);
    const sp = vec2space(vpc);
    const Bscaled = scaleCols(sp.B, [sx, sy, sz]);
    const Bq = scaleMat(sp.Bp, PI2);
    const n = r.rows.length;
    const x = new Float64Array(n);
    const y = new Float64Array(n);
    const z = new Float64Array(n);
    const xa = new Float64Array(n);
    const ya = new Float64Array(n);
    const za = new Float64Array(n);
    const dx = new Float64Array(n);
    const dy = new Float64Array(n);
    const dz = new Float64Array(n);
    const elements = new Array(n);

    for (let i = 0; i < n; i++) {
      const row = r.rows[i];
      const fx = Number(row.cellRefNumX) / sx;
      const fy = Number(row.cellRefNumY) / sy;
      const fz = Number(row.cellRefNumZ) / sz;
      const d0 = wrapDelta(Number(row.x) - fx);
      const d1 = wrapDelta(Number(row.y) - fy);
      const d2 = wrapDelta(Number(row.z) - fz);
      const pA = mul([fx, fy, fz], Bscaled);
      const pD = mul([d0, d1, d2], Bscaled);

      xa[i] = pA[0];
      ya[i] = pA[1];
      za[i] = pA[2];
      dx[i] = pD[0];
      dy[i] = pD[1];
      dz[i] = pD[2];
      x[i] = pA[0] + pD[0];
      y[i] = pA[1] + pD[1];
      z[i] = pA[2] + pD[2];
      elements[i] = String(row.element || "").trim();
    }

    return {
      file: String(name || ""),
      atoms: n,
      super: [sx, sy, sz],
      cellDeg: r.cellDeg.slice(0, 6),
      Bp: sp.Bp,
      Bq,
      x,
      y,
      z,
      xa,
      ya,
      za,
      dx,
      dy,
      dz,
      elements,
    };
  }

  global.RMC6fReader = Object.freeze({
    parse: parseRmc6f,
    Parser: RMC6fParser,
  });
})(typeof window !== "undefined" ? window : globalThis);
