/*
 * Diffuse-scattering data format converter core.
 *
 * Converts between:
 *   - Unified data format HDF5 (/scattering/data and /entry/data layouts,
 *     as written by RMCProfile's unified_hdf5_io.f90 and DiffuseCode)
 *   - Yell 1.0 HDF5 (flat /data, /lower_limits, /step_sizes, /unit_cell)
 *   - RMCProfile old 3d-diffuse text format ("npoints nsec" header followed
 *     by "i j k qx qy qz intensity" rows with Q in cartesian 1/Angstrom)
 *
 * The internal model mirrors rmc_unified_data_t:
 *   dims    [nh, nk, nl]
 *   corner  hkl of pixel (1,1,1)
 *   vectors vectors[axis][component] hkl increment per pixel step
 *   values  Float64Array, h fastest: idx = (il*nk + ik)*nh + ih
 *   cellLengths/cellAngles: cell stored with the data (may be unit metric)
 *
 * Axis-order conventions (verified against files produced by the Fortran
 * writers and by DISCUS):
 *   /scattering/data/data     C dims [nl,nk,nh]  (flat: h fastest)
 *   /entry/data/data_values   C dims [nh,nk,nl]  (flat: l fastest)
 *   Yell /data                C dims [nh,nk,nl]  (flat: l fastest)
 *   /scattering/data/step_vectors      flat[axis*3 + comp]
 *   /entry/data/data_increment_vector  flat[comp*3 + axis]
 *
 * Scatty VTK (STRUCTURED_POINTS, ASCII): ORIGIN/SPACING are cartesian Q in
 * 1/Angstrom with the 2*pi convention (verified against Scatty's paired
 * *_list.txt hkl output); point values are x fastest, z slowest, which
 * matches the internal h-fastest layout directly.
 */
(function (root, factory) {
    if (typeof module === 'object' && module.exports) {
        module.exports = factory();
    } else {
        root.Converter = factory();
    }
}(typeof self !== 'undefined' ? self : this, function () {
    'use strict';

    const DEG = Math.PI / 180.0;
    const UNIT_METRIC_TOL = 1e-6;

    // ----------------------------------------------------------------- cell math

    // Mirrors cell_to_lattice: columns of A are the lattice vectors.
    function cellToLattice(lengths, angles) {
        const [a, b, c] = lengths;
        const [al, be, ga] = angles.map(x => x * DEG);
        if (!(a > 0 && b > 0 && c > 0)) throw new Error('cell lengths must be positive');
        const ca = Math.cos(al), cb = Math.cos(be), cg = Math.cos(ga), sg = Math.sin(ga);
        if (Math.abs(sg) < 1e-12) throw new Error('gamma angle gives a singular lattice');
        const A = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
        A[0][0] = a;
        A[0][1] = b * cg; A[1][1] = b * sg;
        A[0][2] = c * cb;
        A[1][2] = c * (ca - cb * cg) / sg;
        let z2 = c * c - A[0][2] * A[0][2] - A[1][2] * A[1][2];
        if (z2 < 0 && Math.abs(z2) < 1e-10) z2 = 0;
        if (z2 < 0) throw new Error('cell angles do not produce a real lattice');
        A[2][2] = Math.sqrt(z2);
        return A;
    }

    function latticeToCell(A) {
        const col = j => [A[0][j], A[1][j], A[2][j]];
        const norm = v => Math.hypot(v[0], v[1], v[2]);
        const dot = (u, v) => u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
        const a = norm(col(0)), b = norm(col(1)), c = norm(col(2));
        const clamp = x => Math.max(-1, Math.min(1, x));
        return {
            lengths: [a, b, c],
            angles: [
                Math.acos(clamp(dot(col(1), col(2)) / (b * c))) / DEG,
                Math.acos(clamp(dot(col(0), col(2)) / (a * c))) / DEG,
                Math.acos(clamp(dot(col(0), col(1)) / (a * b))) / DEG,
            ],
        };
    }

    function cross(u, v) {
        return [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
    }

    // Columns of B are the reciprocal basis vectors (without 2*pi).
    function reciprocalBasis(A) {
        const col = j => [A[0][j], A[1][j], A[2][j]];
        const a1 = col(0), a2 = col(1), a3 = col(2);
        const vol = a1[0] * cross(a2, a3)[0] + a1[1] * cross(a2, a3)[1] + a1[2] * cross(a2, a3)[2];
        if (Math.abs(vol) < 1e-12) throw new Error('cell gives a singular reciprocal basis');
        const b1 = cross(a2, a3).map(x => x / vol);
        const b2 = cross(a3, a1).map(x => x / vol);
        const b3 = cross(a1, a2).map(x => x / vol);
        return [[b1[0], b2[0], b3[0]], [b1[1], b2[1], b3[1]], [b1[2], b2[2], b3[2]]];
    }

    function hklToQ(B, hkl) {
        const t = 2 * Math.PI;
        return [
            t * (B[0][0] * hkl[0] + B[0][1] * hkl[1] + B[0][2] * hkl[2]),
            t * (B[1][0] * hkl[0] + B[1][1] * hkl[1] + B[1][2] * hkl[2]),
            t * (B[2][0] * hkl[0] + B[2][1] * hkl[1] + B[2][2] * hkl[2]),
        ];
    }

    // hkl = A^T q / 2pi  (mirrors cart_q_to_hkl)
    function qToHkl(A, q) {
        const t = 2 * Math.PI;
        return [
            (A[0][0] * q[0] + A[1][0] * q[1] + A[2][0] * q[2]) / t,
            (A[0][1] * q[0] + A[1][1] * q[1] + A[2][1] * q[2]) / t,
            (A[0][2] * q[0] + A[1][2] * q[1] + A[2][2] * q[2]) / t,
        ];
    }

    function isUnitMetric(lengths, angles) {
        return lengths.every(x => Math.abs(x - 1) <= UNIT_METRIC_TOL) &&
               angles.every(x => Math.abs(x - 90) <= UNIT_METRIC_TOL);
    }

    // ----------------------------------------------------------------- structure

    // Parse an .rmc6f header; returns parent cell = supercell / dimensions.
    function parseRmc6f(text) {
        let unitCells = null, cell = null, latticeRows = null;
        const lines = text.split(/\r?\n/);
        for (let n = 0; n < lines.length; n++) {
            const line = lines[n];
            if (/^Supercell dimensions:/i.test(line)) {
                unitCells = line.split(':')[1].trim().split(/\s+/).map(Number);
            } else if (/^Cell \(Ang\/deg\):/i.test(line)) {
                cell = line.split(':')[1].trim().split(/\s+/).map(Number);
            } else if (/^Lattice vectors \(Ang\):/i.test(line)) {
                latticeRows = [lines[n + 1], lines[n + 2], lines[n + 3]]
                    .map(s => s.trim().split(/\s+/).map(Number));
            } else if (/^Atoms/i.test(line)) {
                break;
            }
        }
        if (!cell && latticeRows) {
            // Each rmc6f line holds one lattice vector: column i of A.
            const A = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
            for (let i = 0; i < 3; i++) for (let r = 0; r < 3; r++) A[r][i] = latticeRows[i][r];
            const derived = latticeToCell(A);
            cell = derived.lengths.concat(derived.angles);
        }
        if (!cell) throw new Error('rmc6f: no "Cell (Ang/deg)" or "Lattice vectors" found');
        if (!unitCells || unitCells.some(x => !(x > 0))) {
            throw new Error('rmc6f: no valid "Supercell dimensions" found');
        }
        return {
            lengths: [cell[0] / unitCells[0], cell[1] / unitCells[1], cell[2] / unitCells[2]],
            angles: [cell[3], cell[4], cell[5]],
            supercell: unitCells,
        };
    }

    // Parent cell from a unified structure HDF5 file (h5wasm File object).
    function readUnifiedStructure(f) {
        const need = p => {
            const d = f.get(p);
            if (!d) throw new Error('unified structure: missing ' + p);
            return Array.from(d.value).map(Number);
        };
        const lengths = need('entry/data/unit_cell_lengths');
        const angles = need('entry/data/unit_cell_angles');
        const cells = need('entry/data/unit_cells');
        if (cells.some(x => !(x > 0))) throw new Error('unified structure: invalid unit_cells');
        return {
            lengths: [lengths[0] / cells[0], lengths[1] / cells[1], lengths[2] / cells[2]],
            angles,
            supercell: cells,
        };
    }

    // ----------------------------------------------------------------- readers

    function detectH5Kind(f) {
        if (f.get('scattering/data/data')) return 'unified';
        if (f.get('entry/data/data_values')) return 'unified';
        if (f.get('data') && f.get('lower_limits') && f.get('unit_cell')) return 'yell';
        if (f.get('entry/data/atom_position') || f.get('entry/data/unit_cells')) return 'structure';
        return 'unknown';
    }

    function readUnifiedData(f) {
        if (f.get('scattering/data/data')) return readScatteringGroup(f);
        if (f.get('entry/data/data_values')) return readEntryGroup(f);
        throw new Error('no unified diffuse data group found');
    }

    function readScatteringGroup(f) {
        const g = 'scattering/data/';
        const ds = f.get(g + 'data');
        const shape = ds.shape;                     // [nl, nk, nh]
        if (shape.length !== 3) throw new Error('scattering data must be rank 3');
        const [nl, nk, nh] = shape;
        const values = Float64Array.from(ds.value); // already h fastest
        const corner = Array.from(f.get(g + 'lower_limits').value).map(Number);
        const sv = Array.from(f.get(g + 'step_vectors').value).map(Number);
        const vectors = [0, 1, 2].map(axis => [0, 1, 2].map(comp => sv[axis * 3 + comp]));
        const lengths = Array.from(f.get(g + 'unit_cell_lengths').value).map(Number);
        const angles = Array.from(f.get(g + 'unit_cell_angles').value).map(Number);
        let axes = [1, 2, 3];
        const axesDs = f.get(g + 'data_axes');
        if (axesDs) axes = Array.from(axesDs.value).map(Number);
        let radiation = 'unknown';
        const attrs = f.get('scattering/data').attrs;
        if (attrs && attrs.radiation) {
            const v = attrs.radiation.value;
            radiation = String(Array.isArray(v) ? v[0] : v);
        }
        return {
            dims: [nh, nk, nl], corner, vectors, values,
            cellLengths: lengths, cellAngles: angles, radiation, axes,
        };
    }

    function readEntryGroup(f) {
        const g = 'entry/data/';
        const dims = Array.from(f.get(g + 'data_dimension').value).map(Number); // [nh,nk,nl]
        const [nh, nk, nl] = dims;
        const ds = f.get(g + 'data_values');       // C dims [nh,nk,nl], l fastest
        const flat = ds.value;
        const values = new Float64Array(nh * nk * nl);
        for (let il = 0; il < nl; il++)
            for (let ik = 0; ik < nk; ik++)
                for (let ih = 0; ih < nh; ih++)
                    values[(il * nk + ik) * nh + ih] = flat[(ih * nk + ik) * nl + il];
        const corner = Array.from(f.get(g + 'data_corner').value).map(Number);
        const iv = Array.from(f.get(g + 'data_increment_vector').value).map(Number);
        const vectors = [0, 1, 2].map(axis => [0, 1, 2].map(comp => iv[comp * 3 + axis]));
        const lengths = Array.from(f.get(g + 'unit_cell_lengths').value).map(Number);
        const angles = Array.from(f.get(g + 'unit_cell_angles').value).map(Number);
        let axes = [1, 2, 3];
        const axesDs = f.get(g + 'data_axes');
        if (axesDs) axes = Array.from(axesDs.value).map(Number);
        let radiation = 'unknown';
        const radDs = f.get(g + 'data_radiation');
        if (radDs) {
            const v = radDs.value;
            radiation = String(Array.isArray(v) ? v[0] : v).trim() || 'unknown';
        }
        return {
            dims, corner, vectors, values,
            cellLengths: lengths, cellAngles: angles, radiation, axes,
        };
    }

    function readYell(f) {
        const need = p => {
            const d = f.get(p);
            if (!d) throw new Error('Yell file: missing ' + p);
            return d;
        };
        const dirDs = f.get('is_direct');
        if (dirDs) {
            const v = Number(Array.isArray(dirDs.value) ? dirDs.value[0] : dirDs.value);
            if (v !== 0) throw new Error('Yell file holds direct-space data; only reciprocal space is supported');
        }
        const ds = need('data');
        const shape = ds.shape;                    // [nh, nk, nl], l fastest
        if (shape.length !== 3) throw new Error('Yell data must be rank 3');
        const [nh, nk, nl] = shape;
        const flat = ds.value;
        const values = new Float64Array(nh * nk * nl);
        for (let il = 0; il < nl; il++)
            for (let ik = 0; ik < nk; ik++)
                for (let ih = 0; ih < nh; ih++)
                    values[(il * nk + ik) * nh + ih] = flat[(ih * nk + ik) * nl + il];
        const corner = Array.from(need('lower_limits').value).map(Number);
        let vectors;
        if (f.get('step_sizes_abs') && f.get('step_sizes_ord') && f.get('step_sizes_top')) {
            vectors = ['step_sizes_abs', 'step_sizes_ord', 'step_sizes_top']
                .map(p => Array.from(f.get(p).value).map(Number));
        } else {
            const s = Array.from(need('step_sizes').value).map(Number);
            vectors = [[s[0], 0, 0], [0, s[1], 0], [0, 0, s[2]]];
        }
        const cell = Array.from(need('unit_cell').value).map(Number);
        return {
            dims: [nh, nk, nl], corner, vectors, values,
            cellLengths: cell.slice(0, 3), cellAngles: cell.slice(3, 6),
            radiation: 'unknown', axes: [1, 2, 3],
        };
    }

    // ------------------------------------------------------------- old text .dat

    function parseOldDat(text, parentCell) {
        const A = cellToLattice(parentCell.lengths, parentCell.angles);
        const tokens = text.trim().split(/\s+/);
        let p = 0;
        const npoints = parseInt(tokens[p++], 10);
        const nsec = parseInt(tokens[p++], 10);
        if (!(npoints > 0) || !(nsec >= 1)) throw new Error('bad npoints/nsec header');
        const perRow = 3 + 3 * nsec + 1;
        if (tokens.length < 2 + npoints * perRow) throw new Error('old-format data file is truncated');

        const pix = new Int32Array(3 * npoints);
        const hklAll = new Float64Array(3 * npoints);
        const vals = new Float64Array(npoints);
        let dims = [0, 0, 0];
        for (let n = 0; n < npoints; n++) {
            const i = parseInt(tokens[p], 10), j = parseInt(tokens[p + 1], 10), k = parseInt(tokens[p + 2], 10);
            const q = [Number(tokens[p + 3]), Number(tokens[p + 4]), Number(tokens[p + 5])];
            vals[n] = Number(tokens[p + perRow - 1]);
            p += perRow;
            if (!(i >= 1 && j >= 1 && k >= 1)) throw new Error('pixel coordinates must be positive');
            pix[3 * n] = i; pix[3 * n + 1] = j; pix[3 * n + 2] = k;
            const hkl = qToHkl(A, q);
            hklAll[3 * n] = hkl[0]; hklAll[3 * n + 1] = hkl[1]; hklAll[3 * n + 2] = hkl[2];
            if (i > dims[0]) dims[0] = i;
            if (j > dims[1]) dims[1] = j;
            if (k > dims[2]) dims[2] = k;
        }
        if (dims[0] * dims[1] * dims[2] !== npoints) {
            throw new Error('old-format data does not cover a full pixel grid');
        }
        const [nh, nk, nl] = dims;
        const values = new Float64Array(npoints);
        let corner = null;
        const stepPix = [null, null, null];
        for (let n = 0; n < npoints; n++) {
            const i = pix[3 * n], j = pix[3 * n + 1], k = pix[3 * n + 2];
            values[((k - 1) * nk + (j - 1)) * nh + (i - 1)] = vals[n];
            const hkl = [hklAll[3 * n], hklAll[3 * n + 1], hklAll[3 * n + 2]];
            if (i === 1 && j === 1 && k === 1) corner = hkl;
            if (i === 2 && j === 1 && k === 1) stepPix[0] = hkl;
            if (i === 1 && j === 2 && k === 1) stepPix[1] = hkl;
            if (i === 1 && j === 1 && k === 2) stepPix[2] = hkl;
        }
        if (!corner) throw new Error('pixel (1,1,1) is missing');
        const vectors = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
        for (let axis = 0; axis < 3; axis++) {
            if (dims[axis] > 1 && stepPix[axis]) {
                vectors[axis] = stepPix[axis].map((x, comp) => x - corner[comp]);
            }
        }
        return {
            dims, corner, vectors, values,
            cellLengths: parentCell.lengths.slice(), cellAngles: parentCell.angles.slice(),
            radiation: 'unknown',
            axes: pickAxes(vectors, dims),
            nsecOriginal: nsec,
        };
    }

    // Mirrors set_axes_from_increment_vectors.
    function pickAxes(vectors, dims) {
        const axes = [1, 2, 3];
        for (let axis = 0; axis < 3; axis++) {
            if (dims[axis] > 1) {
                let best = -1, comp = axis;
                for (let cIdx = 0; cIdx < 3; cIdx++) {
                    if (Math.abs(vectors[axis][cIdx]) > best) {
                        best = Math.abs(vectors[axis][cIdx]);
                        comp = cIdx;
                    }
                }
                axes[axis] = comp + 1;
            }
        }
        for (let axis = 0; axis < 3; axis++)
            for (let other = axis + 1; other < 3; other++)
                if (dims[axis] > 1 && dims[other] > 1 && axes[axis] === axes[other]) return [1, 2, 3];
        return axes;
    }

    function writeOldDat(model, cell) {
        const A = cellToLattice(cell.lengths, cell.angles);
        const B = reciprocalBasis(A);
        const [nh, nk, nl] = model.dims;
        const out = [`${nh * nk * nl} 1`];
        for (let k = 0; k < nl; k++) {
            for (let j = 0; j < nk; j++) {
                for (let i = 0; i < nh; i++) {
                    const hkl = [0, 1, 2].map(c =>
                        model.corner[c] + i * model.vectors[0][c] + j * model.vectors[1][c] + k * model.vectors[2][c]);
                    const q = hklToQ(B, hkl);
                    const v = model.values[(k * nk + j) * nh + i];
                    out.push(`${i + 1} ${j + 1} ${k + 1} ${q[0].toExponential(16)} ${q[1].toExponential(16)} ${q[2].toExponential(16)} ${v.toExponential(16)}`);
                }
            }
        }
        return out.join('\n') + '\n';
    }

    // ------------------------------------------------------------- Scatty VTK

    function isVtk(text) {
        return /^#\s*vtk/i.test(text.trimStart());
    }

    function parseVtk(text, parentCell) {
        const A = cellToLattice(parentCell.lengths, parentCell.angles);
        const lines = text.split(/\r?\n/);
        if (!/^#\s*vtk/i.test((lines[0] || '').trim())) throw new Error('not a VTK file');
        let dims = null, origin = null, spacing = null, dataStart = -1;
        // line 1 is the version comment, line 2 is a free-text title
        for (let i = 2; i < lines.length; i++) {
            const tok = lines[i].trim().split(/\s+/);
            const key = (tok[0] || '').toUpperCase();
            if (key === 'BINARY') throw new Error('binary VTK is not supported (Scatty writes ASCII)');
            if (key === 'DATASET' && (tok[1] || '').toUpperCase() !== 'STRUCTURED_POINTS') {
                throw new Error('only STRUCTURED_POINTS VTK is supported');
            }
            if (key === 'DIMENSIONS') dims = tok.slice(1, 4).map(Number);
            if (key === 'ORIGIN') origin = tok.slice(1, 4).map(Number);
            if (key === 'SPACING' || key === 'ASPECT_RATIO') spacing = tok.slice(1, 4).map(Number);
            if (key === 'LOOKUP_TABLE') { dataStart = i + 1; break; }
        }
        if (!dims || !origin || !spacing || dataStart < 0) {
            throw new Error('VTK header is missing DIMENSIONS/ORIGIN/SPACING/LOOKUP_TABLE');
        }
        if (dims.some(d => !(d >= 1))) throw new Error('bad VTK dimensions');
        const [nh, nk, nl] = dims;
        const npoints = nh * nk * nl;
        const values = new Float64Array(npoints);
        let n = 0;
        for (let r = dataStart; r < lines.length && n < npoints; r++) {
            const t = lines[r].trim();
            if (t === '') continue;
            for (const s of t.split(/\s+/)) {
                if (n >= npoints) break;
                values[n++] = Number(s);
            }
        }
        if (n < npoints) throw new Error(`VTK data is truncated (${n} of ${npoints} values)`);
        // ORIGIN/SPACING are cartesian Q (2*pi/Angstrom); VTK x-fastest point
        // order equals the internal h-fastest layout, so values copy through.
        const corner = qToHkl(A, origin);
        const vectors = [0, 1, 2].map(axis => {
            if (dims[axis] <= 1 || spacing[axis] === 0) return [0, 0, 0];
            const q = [0, 0, 0];
            q[axis] = spacing[axis];
            return qToHkl(A, q);
        });
        return {
            dims, corner, vectors, values,
            cellLengths: parentCell.lengths.slice(), cellAngles: parentCell.angles.slice(),
            radiation: 'unknown', axes: pickAxes(vectors, dims),
        };
    }

    function writeVtk(model, cell) {
        const A = cellToLattice(cell.lengths, cell.angles);
        const B = reciprocalBasis(A);
        const [nh, nk, nl] = model.dims;
        const origin = hklToQ(B, model.corner);
        const spacing = [0, 0, 0];
        for (let axis = 0; axis < 3; axis++) {
            if (model.dims[axis] <= 1) continue;
            const q = hklToQ(B, model.vectors[axis]);
            const len = Math.hypot(q[0], q[1], q[2]);
            for (let comp = 0; comp < 3; comp++) {
                if (comp !== axis && Math.abs(q[comp]) > 1e-6 * Math.max(1e-12, len)) {
                    throw new Error('grid axes are not aligned with cartesian x/y/z in Q; ' +
                        'VTK STRUCTURED_POINTS cannot represent this grid');
                }
            }
            if (q[axis] <= 0) {
                throw new Error('grid step along axis ' + (axis + 1) + ' is not positive in Q; ' +
                    'VTK STRUCTURED_POINTS needs an ascending grid');
            }
            spacing[axis] = q[axis];
        }
        const out = [
            '# vtk DataFile Version 2.0',
            'TITLE diffuse scattering',
            'ASCII',
            'DATASET STRUCTURED_POINTS',
            `DIMENSIONS ${nh} ${nk} ${nl}`,
            'ORIGIN ' + origin.map(x => x.toExponential(16)).join(' '),
            'SPACING ' + spacing.map(x => x.toExponential(16)).join(' '),
            `POINT_DATA ${nh * nk * nl}`,
            'SCALARS diffuse_scattering float',
            'LOOKUP_TABLE default',
        ];
        for (let n = 0; n < model.values.length; n++) {
            out.push(model.values[n].toExponential(16));
        }
        return out.join('\n') + '\n';
    }

    // ----------------------------------------------------------------- writers

    function fixedStr(s) {
        const t = String(s);
        return { data: [t], shape: [1], dtype: 'S' + Math.max(1, t.length) };
    }

    function axisNames(model) {
        const basis = ['h', 'k', 'l'];
        return [0, 1, 2].map(a => basis[(model.axes[a] || a + 1) - 1]);
    }

    function writeUnifiedData(f, model, cell, meta) {
        meta = meta || {};
        const [nh, nk, nl] = model.dims;
        const names = axisNames(model);
        const today = new Date().toISOString().slice(0, 10);
        const method = meta.creationMethod || 'RMCProfile web format converter';
        const author = meta.authorName || 'RMCProfile';
        const radiation = model.radiation && model.radiation !== '' ? model.radiation : 'unknown';

        // ---- /scattering/data (current unified layout) ----
        f.create_attribute('audit_conform_dict_name', 'Disorder scattering');
        f.create_attribute('audit_conform_dict_version', '0.0.0');
        f.create_attribute('audit_creation_date', today);
        f.create_attribute('audit_creation_method', method);
        f.create_attribute('audit_author_name', author);
        f.create_attribute('default', 'scattering');

        const scat = f.create_group('scattering');
        scat.create_attribute('NX_class', 'NXentry');
        scat.create_attribute('default', 'data');
        const sd = scat.create_group('data');
        sd.create_attribute('NX_class', 'NXdata');
        sd.create_attribute('signal', 'data');
        sd.create_attribute('axes', names, [3], 'S1');
        sd.create_attribute('indices_abs', 0, [], '<i');
        sd.create_attribute('indices_ord', 1, [], '<i');
        sd.create_attribute('indices_top', 2, [], '<i');
        sd.create_attribute('radiation', radiation);
        sd.create_attribute('space', 'reciprocal');
        sd.create_attribute('content', 'intensity');
        sd.create_attribute('dimension', Math.max(1, model.dims.filter(d => d > 1).length), [], '<i');
        sd.create_attribute('data_type_experiment', meta.experiment || 'unknown');
        sd.create_attribute('data_type_style', 'single_diffraction');
        sd.create_attribute('data_type_with_bragg', 'unknown');
        sd.create_attribute('data_type_symmetrized', 'none');
        sd.create_attribute('data_rad_symbol', 'unknown');

        sd.create_dataset({ name: 'lower_limits', data: model.corner, shape: [3], dtype: '<d' });
        const sv = new Float64Array(9);
        for (let axis = 0; axis < 3; axis++)
            for (let comp = 0; comp < 3; comp++) sv[axis * 3 + comp] = model.vectors[axis][comp];
        sd.create_dataset({ name: 'step_vectors', data: sv, shape: [3, 3], dtype: '<d' });
        sd.create_dataset({ name: 'data_axes', data: Int32Array.from(model.axes), shape: [3], dtype: '<i' });
        for (let axis = 0; axis < 3; axis++) {
            const comp = model.axes[axis] - 1;
            const vals = new Float64Array(model.dims[axis]);
            for (let i = 0; i < vals.length; i++) vals[i] = model.corner[comp] + i * model.vectors[axis][comp];
            sd.create_dataset({ name: names[axis], data: vals, shape: [vals.length], dtype: '<d' });
        }
        sd.create_dataset({ name: 'data', data: model.values, shape: [nl, nk, nh], dtype: '<d' });
        sd.create_dataset({ name: 'unit_cell_lengths', data: cell.lengths, shape: [3], dtype: '<d' });
        sd.create_dataset({ name: 'unit_cell_angles', data: cell.angles, shape: [3], dtype: '<d' });
        sd.create_dataset({ name: 'data_rad_length', data: [0, 0, 0], shape: [3], dtype: '<d' });

        // ---- /entry/data (older unified layout, kept for compatibility) ----
        const entry = f.create_group('entry');
        const ed = entry.create_group('data');
        const eds = (name, spec) => ed.create_dataset(Object.assign({ name }, spec));
        eds('unit_cell_lengths', { data: cell.lengths, shape: [3], dtype: '<d' });
        eds('unit_cell_angles', { data: cell.angles, shape: [3], dtype: '<d' });
        eds('symmetry_space_group_name_H-M', fixedStr('P 1'));
        eds('space_group_origin', { data: [1], shape: [1], dtype: '<i' });
        eds('symmetry_space_group_abc', fixedStr('abc'));
        eds('space_group_symop_number', { data: [1], shape: [1], dtype: '<i' });
        eds('space_group_symop_operation_mat', {
            data: [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0], shape: [3, 4, 1], dtype: '<d',
        });
        eds('data_type_experiment', fixedStr(meta.experiment || 'unknown'));
        eds('data_type_style', fixedStr('single_diffraction'));
        eds('data_type_content', fixedStr('intensity'));
        eds('data_type_reciprocal', fixedStr('reciprocal'));
        eds('data_type_with_bragg', fixedStr('unknown'));
        eds('data_type_symmetrized', fixedStr('none'));
        eds('data_radiation', fixedStr(radiation));
        eds('data_rad_symbol', fixedStr('unknown'));
        eds('data_rad_length', { data: [0, 0, 0], shape: [3], dtype: '<d' });
        eds('data_dimension', { data: Int32Array.from(model.dims), shape: [3], dtype: '<i' });
        eds('data_axes', { data: Int32Array.from(model.axes), shape: [3], dtype: '<i' });
        eds('data_corner', { data: model.corner, shape: [3], dtype: '<d' });
        const iv = new Float64Array(9);
        for (let axis = 0; axis < 3; axis++)
            for (let comp = 0; comp < 3; comp++) iv[comp * 3 + axis] = model.vectors[axis][comp];
        eds('data_increment_vector', { data: iv, shape: [3, 3], dtype: '<d' });
        const lv = new Float64Array(nh * nk * nl);   // l fastest
        for (let il = 0; il < nl; il++)
            for (let ik = 0; ik < nk; ik++)
                for (let ih = 0; ih < nh; ih++)
                    lv[(ih * nk + ik) * nl + il] = model.values[(il * nk + ik) * nh + ih];
        eds('data_values', { data: lv, shape: [nh, nk, nl], dtype: '<d' });
        eds('audit_conform_dict_name', fixedStr('Disorder scattering'));
        eds('audit_conform_dict_version', fixedStr('0.0.0'));
        eds('audit_creation_date', fixedStr(today));
        eds('audit_creation_method', fixedStr(method));
        eds('audit_author_name', fixedStr(author));
    }

    function writeYell(f, model, cell) {
        const [nh, nk, nl] = model.dims;
        f.create_dataset({ name: 'format', data: ['Yell 1.0'], shape: [1], dtype: 'S8' });
        f.create_dataset({ name: 'is_direct', data: [0], shape: [1], dtype: '<b' });
        f.create_dataset({ name: 'lower_limits', data: model.corner, shape: [3], dtype: '<d' });
        f.create_dataset({
            name: 'step_sizes',
            data: [model.vectors[0][0], model.vectors[1][1], model.vectors[2][2]],
            shape: [3], dtype: '<d',
        });
        f.create_dataset({ name: 'step_sizes_abs', data: model.vectors[0], shape: [3], dtype: '<d' });
        f.create_dataset({ name: 'step_sizes_ord', data: model.vectors[1], shape: [3], dtype: '<d' });
        f.create_dataset({ name: 'step_sizes_top', data: model.vectors[2], shape: [3], dtype: '<d' });
        f.create_dataset({
            name: 'unit_cell',
            data: cell.lengths.concat(cell.angles), shape: [6], dtype: '<d',
        });
        const lv = new Float64Array(nh * nk * nl);   // l fastest, h slowest
        for (let il = 0; il < nl; il++)
            for (let ik = 0; ik < nk; ik++)
                for (let ih = 0; ih < nh; ih++)
                    lv[(ih * nk + ik) * nl + il] = model.values[(il * nk + ik) * nh + ih];
        f.create_dataset({ name: 'data', data: lv, shape: [nh, nk, nl], dtype: '<d' });
    }

    return {
        cellToLattice, latticeToCell, reciprocalBasis, hklToQ, qToHkl, isUnitMetric,
        parseRmc6f, readUnifiedStructure,
        detectH5Kind, readUnifiedData, readYell,
        parseOldDat, writeOldDat,
        isVtk, parseVtk, writeVtk,
        writeUnifiedData, writeYell,
    };
}));
