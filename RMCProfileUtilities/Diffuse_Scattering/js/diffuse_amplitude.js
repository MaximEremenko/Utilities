"use strict";

(function (global) {
  function needCore() {
    if (!global.DiffuseCore) {
      throw new Error("DiffuseCore missing. Load ./js/diffuse_core.js first.");
    }
    return global.DiffuseCore;
  }

  function needCoeffSource() {
    if (!global.ScatteringCoeffSource) {
      throw new Error(
        "ScatteringCoeffSource missing. Load ./js/scattering_coeff_sources.js first."
      );
    }
    return global.ScatteringCoeffSource;
  }

  function nowMs() {
    if (typeof performance !== "undefined" && performance.now) {
      return performance.now();
    }
    return Date.now();
  }

  function normElem(sym) {
    return String(sym || "").trim().toLowerCase();
  }

  function normalizeScatteringConfig(cfg) {
    const src = cfg || {};
    const type = String(src.type || "neutron").trim().toLowerCase();
    const rawTableNum = Number.isFinite(Number(src.tableNum))
      ? Math.floor(Number(src.tableNum))
      : 1;
    const out = {
      type: type === "xray" || type === "electron" ? type : "neutron",
      model: "fast",
      tableNum: Math.max(1, Math.min(5, rawTableNum)),
      defaultValence: Number.isFinite(Number(src.defaultValence))
        ? Math.trunc(Number(src.defaultValence))
        : 0,
      valenceMapText: String(src.valenceMapText || ""),
    };
    if (out.type === "neutron") {
      const model = String(src.model || "fast").trim().toLowerCase();
      out.model = model === "grouped_exact" ? "grouped_exact" : "fast";
    } else if (out.type === "xray") {
      const model = String(src.model || "waasmaier").trim().toLowerCase();
      out.model = model === "table" ? "table" : "waasmaier";
      out.tableNum = Math.max(0, Math.min(5, rawTableNum));
    } else {
      const model = String(src.model || "neutral_table").trim().toLowerCase();
      out.model = model === "ion_peng" ? "ion_peng" : "neutral_table";
    }
    return out;
  }

  function parseValenceMap(text) {
    const map = Object.create(null);
    const raw = String(text || "").trim();
    if (!raw) return map;
    const parts = raw.split(/[;,]+/);
    for (const part of parts) {
      const s = part.trim();
      if (!s) continue;
      const m = /^([A-Za-z0-9]+)\s*:\s*([+-]?\d+)$/.exec(s);
      if (!m) continue;
      const key = normElem(m[1]);
      const val = Number.parseInt(m[2], 10);
      if (!Number.isInteger(val)) continue;
      map[key] = val;
    }
    return map;
  }

  function packTripletsByIndex(a, b, c, indices) {
    const n = indices.length;
    const out = new Float64Array(n * 3);
    for (let i = 0; i < n; i++) {
      const ix = indices[i];
      const j = i * 3;
      out[j] = a[ix];
      out[j + 1] = b[ix];
      out[j + 2] = c[ix];
    }
    return out;
  }

  function groupAtomsByElement(parsed) {
    const groups = new Map();
    for (let i = 0; i < parsed.atoms; i++) {
      const raw = String(parsed.elements[i] || "").trim();
      const key = normElem(raw);
      if (!groups.has(key)) {
        groups.set(key, {
          key,
          element: raw || key,
          indices: [],
        });
      }
      groups.get(key).indices.push(i);
    }
    return Array.from(groups.values());
  }

  function makeQMagnitudes(targetsPacked) {
    const n = Math.floor(targetsPacked.length / 3);
    const out = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      const j = i * 3;
      const qx = targetsPacked[j];
      const qy = targetsPacked[j + 1];
      const qz = targetsPacked[j + 2];
      out[i] = Math.hypot(qx, qy, qz);
    }
    return out;
  }

  function fillConstant(out, v) {
    for (let i = 0; i < out.length; i++) out[i] = v;
    return out;
  }

  function accumulateScaledComplex(dst, src, scale) {
    const n = scale.length;
    for (let i = 0; i < n; i++) {
      const s = scale[i];
      const j = i * 2;
      dst[j] += src[j] * s;
      dst[j + 1] += src[j + 1] * s;
    }
  }

  function describeConfig(cfg) {
    if (cfg.type === "neutron") {
      return cfg.model === "grouped_exact"
        ? "Neutron grouped exact"
        : "Neutron fast";
    }
    if (cfg.type === "xray") {
      return cfg.model === "table"
        ? `X-ray (table ${cfg.tableNum})`
        : "X-ray (Waasmaier)";
    }
    return cfg.model === "ion_peng"
      ? "Electron ion (Peng)"
      : `Electron neutral (table ${cfg.tableNum})`;
  }

  function buildGroupedEvaluator(cfg) {
    if (cfg.type === "neutron") {
      return function evalGroup(group, qMag) {
        return fillConstant(new Float64Array(qMag.length), group.neutronB);
      };
    }
    const coeff = needCoeffSource();
    if (cfg.type === "xray") {
      const xrayModel =
        cfg.model === "table"
          ? { model: "table", tableNum: cfg.tableNum }
          : { model: "waasmaier" };
      return function evalGroup(group, qMag) {
        const out = new Float64Array(qMag.length);
        for (let i = 0; i < qMag.length; i++) {
          out[i] = coeff.evalXray(group.key, qMag[i], xrayModel);
        }
        return out;
      };
    }
    if (cfg.model === "ion_peng") {
      return function evalGroup(group, qMag) {
        const out = new Float64Array(qMag.length);
        for (let i = 0; i < qMag.length; i++) {
          out[i] = coeff.evalElectronIonPeng(group.key, group.valence, qMag[i]);
        }
        return out;
      };
    }
    return function evalGroup(group, qMag) {
      const out = new Float64Array(qMag.length);
      for (let i = 0; i < qMag.length; i++) {
        out[i] = coeff.evalElectronNeutral(group.key, qMag[i], cfg.tableNum);
      }
      return out;
    };
  }

  function validateGroups(groups, evalGroup) {
    const errs = [];
    const qProbe = new Float64Array([0]);
    for (const g of groups) {
      try {
        evalGroup(g, qProbe);
      } catch (e) {
        errs.push(`${g.element}: ${e && e.message ? e.message : String(e)}`);
      }
    }
    if (errs.length) {
      throw new Error(`Unsupported scattering coefficients: ${errs.join(" | ")}`);
    }
  }

  function buildProfile(parsed, scatteringCfg) {
    const core = needCore();
    const cfg = normalizeScatteringConfig(scatteringCfg);
    const profile = {
      config: cfg,
      description: describeConfig(cfg),
      atomCount: parsed.atoms,
      warnings: [],
    };

    const tPrep = nowMs();
    profile.src = core.interleave3(parsed.x, parsed.y, parsed.z);
    profile.srcA = core.interleave3(parsed.xa, parsed.ya, parsed.za);
    profile.srcD = core.interleave3(parsed.dx, parsed.dy, parsed.dz);
    profile.onesAll = core.complexOnes(parsed.atoms);

    if (cfg.type === "neutron" && cfg.model === "fast") {
      profile.engine = "neutron_fast";
      profile.cAtoms = core.complexReal(parsed.fca);
      profile.prepareMs = nowMs() - tPrep;
      return profile;
    }

    const valenceMap = parseValenceMap(cfg.valenceMapText);
    const incSet = parsed.includedElements || null;
    const rawGroups = groupAtomsByElement(parsed)
      .filter((g) => !incSet || incSet.has(g.element) || incSet.has(g.key));
    const groups = rawGroups.map((g) => {
      const gi = g.indices[0];
      const valence = Number.isInteger(valenceMap[g.key])
        ? valenceMap[g.key]
        : cfg.defaultValence;
      return {
        key: g.key,
        element: g.element,
        count: g.indices.length,
        valence,
        neutronB: Number(parsed.fca[gi]),
        src: packTripletsByIndex(parsed.x, parsed.y, parsed.z, g.indices),
        srcD: packTripletsByIndex(parsed.dx, parsed.dy, parsed.dz, g.indices),
        ones: core.complexOnes(g.indices.length),
      };
    });

    const evalGroup = buildGroupedEvaluator(cfg);
    validateGroups(groups, evalGroup);

    profile.engine = "grouped_exact";
    profile.groups = groups;
    profile.evalGroupCoeffs = evalGroup;
    profile.prepareMs = nowMs() - tPrep;
    return profile;
  }

  async function computeIntensity(args) {
    const core = needCore();
    const parsed = args.parsed;
    const h = args.h;
    const k = args.k;
    const l = args.l;
    const Bq = args.Bq;
    const backend = args.backend;
    const opts = args.opts || {};
    const sub = !!args.sub;
    const runType3 = args.runType3;
    const onStatus = args.onStatus;
    const onStageTiming = args.onStageTiming;
    if (typeof runType3 !== "function") {
      throw new Error("computeIntensity requires runType3 callback");
    }

    const profile = args.profile || buildProfile(parsed, args.scattering);
    const grid = h.length * k.length * l.length;
    const chunkSize = Math.max(
      1,
      Math.floor(Number.isFinite(args.chunkSize) ? args.chunkSize : grid)
    );
    const totalChunks = Math.max(1, Math.ceil(grid / chunkSize));
    const I = new Float64Array(grid);
    const timings = {
      prepare: profile.prepareMs || 0,
      a: 0,
      aavg: 0,
      adelta: 0,
      finalize: 0,
    };
    const mInv = 1 / Math.max(1, parsed.atoms);
    let min = Infinity;
    let max = -Infinity;

    for (let chunk = 0, start = 0; start < grid; chunk++, start += chunkSize) {
      const count = Math.min(chunkSize, grid - start);
      const trg = core.targetsChunk(h, k, l, Bq, start, count);
      const chunkTag = totalChunks > 1 ? ` (${chunk + 1}/${totalChunks})` : "";
      let q = null;
      let qa = null;
      let qd = null;

      if (profile.engine === "neutron_fast") {
        if (typeof onStatus === "function") onStatus(`Computing A(hkl)${chunkTag} ...`);
        const tA = nowMs();
        q = (
          await runType3(
            {
              dim: 3,
              isign: 1,
              sourcesPacked: profile.src,
              targetsPacked: trg,
              strengths: profile.cAtoms,
            },
            opts,
            backend,
            onStageTiming
          )
        ).out;
        timings.a += nowMs() - tA;

        if (sub) {
          if (typeof onStatus === "function") {
            onStatus(`Computing Aavg(hkl)${chunkTag} ...`);
          }
          const tAa = nowMs();
          qa = (
            await runType3(
              {
                dim: 3,
                isign: 1,
                sourcesPacked: profile.srcA,
                targetsPacked: trg,
                strengths: profile.onesAll,
              },
              opts,
              backend,
              onStageTiming
            )
          ).out;
          timings.aavg += nowMs() - tAa;

          if (typeof onStatus === "function") {
            onStatus(`Computing Adelta(hkl)${chunkTag} ...`);
          }
          const tAd = nowMs();
          qd = (
            await runType3(
              {
                dim: 3,
                isign: 1,
                sourcesPacked: profile.srcD,
                targetsPacked: trg,
                strengths: profile.cAtoms,
              },
              opts,
              backend,
              onStageTiming
            )
          ).out;
          timings.adelta += nowMs() - tAd;
        }
      } else {
        const qMag = makeQMagnitudes(trg);
        q = new Float64Array(count * 2);
        qd = sub ? new Float64Array(count * 2) : null;

        if (sub) {
          if (typeof onStatus === "function") {
            onStatus(`Computing Aavg(hkl)${chunkTag} ...`);
          }
          const tAa = nowMs();
          qa = (
            await runType3(
              {
                dim: 3,
                isign: 1,
                sourcesPacked: profile.srcA,
                targetsPacked: trg,
                strengths: profile.onesAll,
              },
              opts,
              backend,
              onStageTiming
            )
          ).out;
          timings.aavg += nowMs() - tAa;
        }

        for (let gi = 0; gi < profile.groups.length; gi++) {
          const g = profile.groups[gi];
          const gTag =
            profile.groups.length > 1
              ? ` [${gi + 1}/${profile.groups.length} ${g.element}]`
              : "";
          const coeff = profile.evalGroupCoeffs(g, qMag);

          if (typeof onStatus === "function") {
            onStatus(`Computing A(hkl)${chunkTag}${gTag} ...`);
          }
          const tA = nowMs();
          const qG = (
            await runType3(
              {
                dim: 3,
                isign: 1,
                sourcesPacked: g.src,
                targetsPacked: trg,
                strengths: g.ones,
              },
              opts,
              backend,
              onStageTiming
            )
          ).out;
          timings.a += nowMs() - tA;
          accumulateScaledComplex(q, qG, coeff);

          if (sub) {
            if (typeof onStatus === "function") {
              onStatus(`Computing Adelta(hkl)${chunkTag}${gTag} ...`);
            }
            const tAd = nowMs();
            const qdG = (
              await runType3(
                {
                  dim: 3,
                  isign: 1,
                  sourcesPacked: g.srcD,
                  targetsPacked: trg,
                  strengths: g.ones,
                },
                opts,
                backend,
                onStageTiming
              )
            ).out;
            timings.adelta += nowMs() - tAd;
            accumulateScaledComplex(qd, qdG, coeff);
          }
        }
      }

      if (typeof onStatus === "function") {
        onStatus(`Finalizing intensity${chunkTag} ...`);
      }
      const tFin = nowMs();
      const mmChunk = core.accumulateIntensityChunk(
        I,
        start,
        count,
        q,
        qa,
        qd,
        mInv,
        sub
      );
      timings.finalize += nowMs() - tFin;
      if (mmChunk.min < min) min = mmChunk.min;
      if (mmChunk.max > max) max = mmChunk.max;
    }

    return {
      I,
      min,
      max,
      timings,
      profile,
      chunkSize,
      totalChunks,
    };
  }

  global.DiffuseAmplitude = Object.freeze({
    normalizeScatteringConfig,
    parseValenceMap,
    buildProfile,
    computeIntensity,
  });
})(typeof window !== "undefined" ? window : globalThis);
