"use strict";

(function (global) {
  const PI2 = 2 * Math.PI;
  const ELEMENT_DATA = (
    "H:1.008 He:4.0026 Li:6.94 Be:9.0122 B:10.81 C:12.011 N:14.007 O:15.999 " +
    "F:18.998 Ne:20.180 Na:22.990 Mg:24.305 Al:26.982 Si:28.085 P:30.974 " +
    "S:32.06 Cl:35.45 Ar:39.948 K:39.098 Ca:40.078 Sc:44.956 Ti:47.867 " +
    "V:50.942 Cr:51.996 Mn:54.938 Fe:55.845 Co:58.933 Ni:58.693 Cu:63.546 " +
    "Zn:65.38 Ga:69.723 Ge:72.630 As:74.922 Se:78.971 Br:79.904 Kr:83.798 " +
    "Rb:85.468 Sr:87.62 Y:88.906 Zr:91.224 Nb:92.906 Mo:95.95 Tc:98 " +
    "Ru:101.07 Rh:102.906 Pd:106.42 Ag:107.868 Cd:112.414 In:114.818 " +
    "Sn:118.710 Sb:121.760 Te:127.60 I:126.904 Xe:131.293 Cs:132.905 " +
    "Ba:137.327 La:138.905 Ce:140.116 Pr:140.908 Nd:144.242 Pm:145 " +
    "Sm:150.36 Eu:151.964 Gd:157.25 Tb:158.925 Dy:162.500 Ho:164.930 " +
    "Er:167.259 Tm:168.934 Yb:173.045 Lu:174.967 Hf:178.49 Ta:180.948 " +
    "W:183.84 Re:186.207 Os:190.23 Ir:192.217 Pt:195.084 Au:196.967 " +
    "Hg:200.592 Tl:204.38 Pb:207.2 Bi:208.980 Po:209 At:210 Rn:222 " +
    "Fr:223 Ra:226 Ac:227 Th:232.038 Pa:231.036 U:238.029"
  ).split(/\s+/).map(function (item) {
    const p = item.split(":");
    return { symbol: p[0], mass: Number(p[1]) };
  });
  const ELEMENTS = new Set(ELEMENT_DATA.map(function (item) { return item.symbol; }));

  // [atom-type column, first coordinate column], zero based.
  const STYLE_FIELDS = Object.freeze({
    atomic: [1, 2], charge: [1, 3], angle: [2, 3], bond: [2, 3],
    molecular: [2, 3], full: [2, 4], body: [1, 4], "bpm/sphere": [2, 5],
    dielectric: [1, 3], dipole: [1, 3], dpd: [1, 3], edpd: [1, 4],
    electron: [1, 5], ellipsoid: [1, 4], line: [2, 5], mdpd: [1, 3],
    peri: [1, 4], rheo: [1, 4], "rheo/thermal": [1, 5], smd: [1, 10],
    sph: [1, 5], sphere: [1, 4], spin: [1, 2], tdpd: [1, 2],
    template: [1, 5], tri: [2, 5], hybrid: [1, 2]
  });

  function clean(line) { return String(line || "").replace(/\x00/g, ""); }
  function body(line) { return clean(line).split("#", 1)[0].trim(); }
  function parts(line) {
    const raw = clean(line), at = raw.indexOf("#");
    return { body: (at < 0 ? raw : raw.slice(0, at)).trim(), comment: at < 0 ? "" : raw.slice(at + 1).trim() };
  }
  function number(value, label) {
    const n = Number(value);
    if (!Number.isFinite(n)) throw new Error("LAMMPS data: invalid " + label);
    return n;
  }
  function element(value) {
    const raw = String(value || "").trim();
    if (!/^[A-Za-z]{1,2}$/.test(raw)) return "";
    const symbol = raw[0].toUpperCase() + raw.slice(1).toLowerCase();
    return ELEMENTS.has(symbol) ? symbol : "";
  }
  function inferElement(mass) {
    if (!(mass > 0)) return "";
    let best = null;
    for (const item of ELEMENT_DATA) {
      const error = Math.abs(item.mass - mass);
      if (!best || error < best.error) best = { symbol: item.symbol, error: error };
    }
    return best && best.error <= Math.max(0.25, mass * 0.015) ? best.symbol : "";
  }
  function findSection(lines, name) {
    const wanted = name.toLowerCase();
    for (let i = 0; i < lines.length; i++) {
      const p = parts(lines[i]);
      if (p.body.toLowerCase() === wanted) return { index: i, comment: p.comment };
    }
    return null;
  }
  function sectionRows(lines, section, count, name) {
    if (!section) return [];
    const out = [];
    for (let i = section.index + 1; i < lines.length && out.length < count; i++) {
      const p = parts(lines[i]);
      if (p.body) out.push({ body: p.body, comment: p.comment, line: i + 1 });
    }
    if (out.length !== count) throw new Error("LAMMPS data: " + name + " has " + out.length + " rows, expected " + count);
    return out;
  }
  function headerCount(lines, suffix, required) {
    const target = suffix.toLowerCase();
    for (const line of lines) {
      const p = body(line).split(/\s+/);
      if (p.length >= 2 && p.slice(1).join(" ").toLowerCase() === target && /^\d+$/.test(p[0])) return Number(p[0]);
    }
    if (required) throw new Error("LAMMPS data: missing " + suffix + " count");
    return 0;
  }
  function header3(lines, keyword) {
    const target = keyword.toLowerCase().split(/\s+/);
    for (const line of lines) {
      const p = body(line).split(/\s+/);
      if (p.length === 3 + target.length && p.slice(3).join(" ").toLowerCase() === target.join(" ")) {
        return p.slice(0, 3).map(Number);
      }
    }
    return null;
  }
  function bounds(lines, axis) {
    const target = axis + "lo " + axis + "hi";
    for (const line of lines) {
      const p = body(line).split(/\s+/);
      if (p.length === 4 && p.slice(2).join(" ").toLowerCase() === target) return p.slice(0, 2).map(Number);
    }
    return null;
  }
  function units(lines, warnings) {
    const m = /\bunits\s*=\s*([A-Za-z]+)/i.exec(clean(lines[0]));
    const style = m ? m[1].toLowerCase() : "";
    const factors = { metal: 1, real: 1, electron: 0.529177210903, nano: 10, micro: 1e4, si: 1e10, cgs: 1e8 };
    if (!style) {
      warnings.push("LAMMPS unit style is not stored; assuming coordinates are in Angstrom.");
      return { style: "unknown", factor: 1, atomicMasses: true };
    }
    if (style === "lj") throw new Error("LAMMPS data: units=lj has no absolute Angstrom scale");
    if (!Object.prototype.hasOwnProperty.call(factors, style)) throw new Error("LAMMPS data: unsupported unit style " + style);
    return { style: style, factor: factors[style], atomicMasses: style === "metal" || style === "real" };
  }
  function dot(a, b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
  function cross(a, b) { return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]; }
  function norm(a) { return Math.hypot(a[0], a[1], a[2]); }
  function angle(a, b) { return Math.acos(Math.max(-1, Math.min(1, dot(a, b) / (norm(a) * norm(b))))) * 180 / Math.PI; }
  function box(lines, factor) {
    const av = header3(lines, "avec"), bv = header3(lines, "bvec"), cv = header3(lines, "cvec"), ao = header3(lines, "abc origin");
    let A, B, C, origin, kind;
    if (av || bv || cv || ao) {
      if (!av || !bv || !cv || !ao) throw new Error("LAMMPS data: incomplete general triclinic box");
      A = av; B = bv; C = cv; origin = ao; kind = "general triclinic";
    } else {
      const xb = bounds(lines, "x"), yb = bounds(lines, "y"), zb = bounds(lines, "z");
      if (!xb || !yb || !zb) throw new Error("LAMMPS data: incomplete box bounds");
      const tilt = header3(lines, "xy xz yz") || [0, 0, 0];
      A = [xb[1] - xb[0], 0, 0]; B = [tilt[0], yb[1] - yb[0], 0]; C = [tilt[1], tilt[2], zb[1] - zb[0]];
      origin = [xb[0], yb[0], zb[0]]; kind = tilt.some(function (v) { return v !== 0; }) ? "restricted triclinic" : "orthogonal";
    }
    A = A.map(function (v) { return number(v, "box vector") * factor; });
    B = B.map(function (v) { return number(v, "box vector") * factor; });
    C = C.map(function (v) { return number(v, "box vector") * factor; });
    origin = origin.map(function (v) { return number(v, "box origin") * factor; });
    if (!(dot(A, cross(B, C)) > 1e-20)) throw new Error("LAMMPS data: box vectors must form a right-handed 3-D cell");
    return { A: A, B: B, C: C, origin: origin, kind: kind };
  }
  function fractional(point, cell) {
    const r = [point[0] - cell.origin[0], point[1] - cell.origin[1], point[2] - cell.origin[2]], v = dot(cell.A, cross(cell.B, cell.C));
    return [dot(r, cross(cell.B, cell.C)) / v, dot(r, cross(cell.C, cell.A)) / v, dot(r, cross(cell.A, cell.B)) / v];
  }
  function wrappedPoint(point, cell) {
    const f = fractional(point, cell).map(function (v) { const w = v - Math.floor(v); return Math.abs(w - 1) < 1e-12 || Math.abs(w) < 1e-12 ? 0 : w; });
    return [
      cell.A[0] * f[0] + cell.B[0] * f[1] + cell.C[0] * f[2],
      cell.A[1] * f[0] + cell.B[1] * f[1] + cell.C[1] * f[2],
      cell.A[2] * f[0] + cell.B[2] * f[1] + cell.C[2] * f[2]
    ];
  }
  function reciprocal(cell) {
    const v = dot(cell.A, cross(cell.B, cell.C));
    return [cross(cell.B, cell.C).map(function (x) { return x / v; }), cross(cell.C, cell.A).map(function (x) { return x / v; }), cross(cell.A, cell.B).map(function (x) { return x / v; })];
  }
  function atomStyle(section, rows, warnings) {
    const explicit = String(section.comment || "").split(/\s+/)[0].toLowerCase();
    if (explicit) {
      if (!STYLE_FIELDS[explicit]) throw new Error("LAMMPS data: unsupported atom style " + explicit);
      return explicit;
    }
    const count = rows[0].body.split(/\s+/).length;
    if (count === 5 || count === 8) return "atomic";
    if (count === 6 || count === 9) { warnings.push("Atoms section has no style tag; interpreting six-column rows as atom_style charge."); return "charge"; }
    if (count === 7 || count === 10) { warnings.push("Atoms section has no style tag; interpreting seven-column rows as atom_style full."); return "full"; }
    throw new Error("LAMMPS data: Atoms section needs a supported # atom-style tag");
  }
  function typeMaps(lines, count, unitInfo, warnings) {
    const labels = new Map(), masses = new Map(), resolved = new Map();
    const ls = findSection(lines, "Atom Type Labels"), ms = findSection(lines, "Masses");
    for (const row of sectionRows(lines, ls, ls ? count : 0, "Atom Type Labels")) { const p = row.body.split(/\s+/); labels.set(String(p[0]), String(p[1] || "")); }
    for (const row of sectionRows(lines, ms, ms ? count : 0, "Masses")) {
      const p = row.body.split(/\s+/), key = String(p[0]);
      masses.set(key, { mass: number(p[1], "mass for atom type " + key), comment: element(row.comment.split(/\s+/)[0]) });
    }
    for (let i = 1; i <= count; i++) {
      const key = String(i), m = masses.get(key), found = element(labels.get(key)) || (m && m.comment) || (m && unitInfo.atomicMasses ? inferElement(m.mass) : "");
      if (found) resolved.set(key, found);
    }
    if (!resolved.size) warnings.push("No chemical elements could be resolved from type labels or Masses.");
    return { labels: labels, resolved: resolved };
  }
  function resolveType(token, maps) {
    const key = String(token), label = maps.labels.get(key) || key;
    return element(label) || maps.resolved.get(key) || "Type" + key;
  }
  function parse(name, text) {
    const lines = String(text || "").split(/\r?\n/).map(clean), warnings = [];
    const atomCount = headerCount(lines, "atoms", true), typeCount = headerCount(lines, "atom types", true);
    if (!(atomCount > 0) || !(typeCount > 0)) throw new Error("LAMMPS data: invalid atom counts");
    const unitInfo = units(lines, warnings), cell = box(lines, unitInfo.factor), section = findSection(lines, "Atoms");
    if (!section) throw new Error("LAMMPS data: missing Atoms section");
    const rows = sectionRows(lines, section, atomCount, "Atoms"), style = atomStyle(section, rows, warnings), fields = STYLE_FIELDS[style], maps = typeMaps(lines, typeCount, unitInfo, warnings);
    const atoms = rows.map(function (row) {
      const p = row.body.split(/\s+/);
      if (p.length < fields[1] + 3) throw new Error("LAMMPS data: short atom row at line " + row.line);
      const id = Number(p[0]);
      if (!Number.isInteger(id) || id < 0) throw new Error("LAMMPS data: invalid atom ID at line " + row.line);
      const point = [p[fields[1]], p[fields[1] + 1], p[fields[1] + 2]].map(function (v) { return number(v, "coordinate at line " + row.line) * unitInfo.factor; });
      return { id: id, element: resolveType(p[fields[0]], maps), point: wrappedPoint(point, cell) };
    }).sort(function (a, b) { return a.id - b.id; });
    const n = atoms.length, x = new Float64Array(n), y = new Float64Array(n), z = new Float64Array(n), xa = new Float64Array(n), ya = new Float64Array(n), za = new Float64Array(n), dx = new Float64Array(n), dy = new Float64Array(n), dz = new Float64Array(n), elements = new Array(n);
    for (let i = 0; i < n; i++) { x[i] = xa[i] = atoms[i].point[0]; y[i] = ya[i] = atoms[i].point[1]; z[i] = za[i] = atoms[i].point[2]; elements[i] = atoms[i].element; }
    const Bp = reciprocal(cell), Bq = Bp.map(function (row) { return row.map(function (v) { return v * PI2; }); });
    return { file: String(name || ""), atoms: n, super: [1, 1, 1], hasSupercell: false, cellDeg: [norm(cell.A), norm(cell.B), norm(cell.C), angle(cell.B, cell.C), angle(cell.A, cell.C), angle(cell.A, cell.B)], Bp: Bp, Bq: Bq, x: x, y: y, z: z, xa: xa, ya: ya, za: za, dx: dx, dy: dy, dz: dz, elements: elements, sourceFormat: "LAMMPS data (" + style + ", " + cell.kind + ")", sourceWarnings: warnings, lammps: { atomStyle: style, boxKind: cell.kind, units: unitInfo.style } };
  }

  global.LammpsDataReader = Object.freeze({ parse: parse });
})(typeof window !== "undefined" ? window : globalThis);
