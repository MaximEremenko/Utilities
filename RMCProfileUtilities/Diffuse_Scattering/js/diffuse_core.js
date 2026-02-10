"use strict";

(function (global) {
  const NEUTRON = {
    va: 0,
    h: -0.3739,
    d: 0.6671,
    he: 0.326,
    li: -0.19,
    "7l": -0.222,
    be: 0.779,
    b: 0.53,
    c: 0.6646,
    n: 0.936,
    o: 0.5803,
    f: 0.5654,
    ne: 0.4566,
    na: 0.363,
    mg: 0.5375,
    al: 0.3449,
    si: 0.41491,
    p: 0.513,
    s: 0.2847,
    cl: 0.9577,
    ar: 0.1909,
    k: 0.367,
    ca: 0.47,
    sc: 1.229,
    ti: -0.3438,
    v: -0.03824,
    cr: 0.3635,
    mn: -0.373,
    fe: 0.945,
    co: 0.249,
    ni: 1.03,
    cu: 0.7718,
    zn: 0.568,
    ga: 0.7288,
    ge: 0.8185,
    as: 0.658,
    se: 0.797,
    br: 0.6795,
    kr: 0.781,
    rb: 0.709,
    sr: 0.702,
    y: 0.775,
    zr: 0.716,
    nb: 0.7054,
    mo: 0.6715,
    tc: 0.68,
    ru: 0.703,
    rh: 0.588,
    pd: 0.591,
    ag: 0.5922,
    cd: 0.487,
    in: 0.4065,
    sn: 0.6225,
    sb: 0.557,
    te: 0.58,
    i: 0.528,
    xe: 0.492,
    cs: 0.542,
    ba: 0.507,
    la: 0.824,
    ce: 0.484,
    pr: 0.458,
    nd: 0.769,
    pm: 1.26,
    sm: 0.08,
    eu: 0.722,
    gd: 0.65,
    tb: 0.738,
    dy: 1.69,
    ho: 0.801,
    er: 0.779,
    tm: 0.707,
    yb: 1.243,
    lu: 0.721,
    hf: 0.777,
    ta: 0.691,
    w: 0.486,
    re: 0.92,
    os: 1.07,
    ir: 1.06,
    pt: 0.96,
    au: 0.763,
    hg: 1.2692,
    tl: 0.8776,
    pb: 0.9405,
    bi: 0.8532,
    po: 0,
    at: 0,
    rn: 0,
    fr: 0,
    ra: 1,
    ac: 0,
    th: 1.031,
    pa: 0.91,
    u: 0.8417,
    np: 1.055,
    pu: 0,
    am: 0.83,
    cm: 0,
  };

  function mul(v, m) {
    return [
      v[0] * m[0][0] + v[1] * m[1][0] + v[2] * m[2][0],
      v[0] * m[0][1] + v[1] * m[1][1] + v[2] * m[2][1],
      v[0] * m[0][2] + v[1] * m[1][2] + v[2] * m[2][2],
    ];
  }

  function interleave3(x, y, z) {
    const n = x.length;
    const out = new Float64Array(n * 3);
    for (let i = 0; i < n; i++) {
      out[i * 3] = x[i];
      out[i * 3 + 1] = y[i];
      out[i * 3 + 2] = z[i];
    }
    return out;
  }

  function complexReal(a) {
    const out = new Float64Array(a.length * 2);
    for (let i = 0; i < a.length; i++) out[i * 2] = a[i];
    return out;
  }

  function complexOnes(atomCount) {
    const out = new Float64Array(Math.max(0, atomCount) * 2);
    for (let i = 0; i < atomCount; i++) out[i * 2] = 1;
    return out;
  }

  function targetsChunk(h, k, l, Bq, start, count) {
    const nh = h.length;
    const nk = k.length;
    const nl = l.length;
    const total = nh * nk * nl;
    if (start < 0 || count < 0 || start + count > total) {
      throw new Error(
        `Target chunk out of range: start=${start}, count=${count}, total=${total}`
      );
    }
    const t = new Float64Array(count * 3);
    const plane = nk * nl;
    for (let q = 0; q < count; q++) {
      const linear = start + q;
      const ih = Math.floor(linear / plane);
      const rem = linear - ih * plane;
      const ik = Math.floor(rem / nl);
      const il = rem - ik * nl;
      const v = mul([h[ih], k[ik], l[il]], Bq);
      t[q * 3] = v[0];
      t[q * 3 + 1] = v[1];
      t[q * 3 + 2] = v[2];
    }
    return t;
  }

  function attachNeutronCoefficients(parsed, fallback = 10) {
    if (!parsed || !Array.isArray(parsed.elements)) {
      throw new Error("Invalid parsed RMC structure");
    }
    const n = parsed.elements.length;
    const fca = new Float64Array(n);
    const present = new Set();
    const unknown = new Set();
    for (let i = 0; i < n; i++) {
      const elem = String(parsed.elements[i] || "").trim();
      const key = elem.toLowerCase();
      present.add(elem);
      if (Object.prototype.hasOwnProperty.call(NEUTRON, key)) {
        fca[i] = NEUTRON[key];
      } else {
        fca[i] = fallback;
        unknown.add(elem);
      }
    }
    return {
      ...parsed,
      fca,
      present: [...present].sort(),
      unknown: [...unknown].sort(),
    };
  }

  function accumulateIntensityChunk(
    out,
    start,
    count,
    q,
    qa,
    qd,
    mInv,
    subtractAverage
  ) {
    let min = Infinity;
    let max = -Infinity;
    for (let i = 0; i < count; i++) {
      const ar = q[i * 2];
      const ai = q[i * 2 + 1];
      let br = 0;
      let bi = 0;
      if (subtractAverage) {
        const avr = qa[i * 2];
        const avi = qa[i * 2 + 1];
        const dr = qd[i * 2];
        const di = qd[i * 2 + 1];
        br = (avr * dr - avi * di) * mInv;
        bi = (avr * di + avi * dr) * mInv;
      }
      const rr = ar - br;
      const ii = ai - bi;
      const v = rr * rr + ii * ii;
      out[start + i] = v;
      if (v < min) min = v;
      if (v > max) max = v;
    }
    return { min, max };
  }

  global.DiffuseCore = Object.freeze({
    NEUTRON,
    interleave3,
    complexReal,
    complexOnes,
    targetsChunk,
    attachNeutronCoefficients,
    accumulateIntensityChunk,
  });
})(typeof window !== "undefined" ? window : globalThis);
