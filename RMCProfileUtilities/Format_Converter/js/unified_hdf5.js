/*
 * Shared browser helpers for the current RMCProfile/DISCUS unified HDF5
 * structure and data contract. The caller supplies an initialized h5wasm
 * module so this works with both the vendored bundle and an ESM CDN module.
 */
(function (root, factory) {
    root.UnifiedH5 = factory();
}(typeof self !== 'undefined' ? self : this, function () {
    'use strict';

    let sequence = 0;
    const DATA_DICTIONARY = 'Disorder unified data';
    const STRUCTURE_DICTIONARY = 'Disorder structure';
    const LEGACY_DATA_DICTIONARY = 'Disorder scattering';
    const AXES_TYPES = new Set([
        'hkl', 'Q', '2theta', 'dstar', 'sin(theta)/lambda', 'theta',
        'xyz', 'uvw', 'r',
    ]);

    function fixedStr(value) {
        const text = String(value == null ? '' : value);
        return { data: [text], shape: [1], dtype: 'S' + Math.max(1, text.length) };
    }

    function h5Value(obj) {
        if (!obj) return null;
        if (typeof obj.value === 'function') return obj.value();
        return obj.value;
    }

    function textValue(value) {
        if (value == null) return '';
        if (ArrayBuffer.isView(value)) {
            if (value instanceof Uint8Array || value instanceof Int8Array) {
                const bytes = new Uint8Array(value.buffer, value.byteOffset, value.byteLength);
                const end = bytes.indexOf(0);
                return new TextDecoder().decode(end < 0 ? bytes : bytes.subarray(0, end)).trim();
            }
            value = Array.from(value);
        }
        if (Array.isArray(value)) {
            if (!value.length) return '';
            return String(value[0]).replace(/\0+$/g, '').trim();
        }
        return String(value).replace(/\0+$/g, '').trim();
    }

    function get(file, path) {
        try {
            return file.get(path) || file.get('/' + String(path).replace(/^\/+/, '')) || null;
        } catch (_) {
            return null;
        }
    }

    function textDataset(file, path) {
        return textValue(h5Value(get(file, path)));
    }

    function numericDataset(file, path, required) {
        const ds = get(file, path);
        if (!ds) {
            if (required) throw new Error('unified HDF5: missing ' + path);
            return [];
        }
        const value = h5Value(ds);
        const flat = ArrayBuffer.isView(value) ? Array.from(value) :
            Array.isArray(value) ? value.flat(Infinity) : [value];
        return flat.map(Number);
    }

    function datasetShape(file, path) {
        const ds = get(file, path);
        return ds && (Array.isArray(ds.shape) || ArrayBuffer.isView(ds.shape))
            ? Array.from(ds.shape).map(Number) : [];
    }

    function validateH5(h5) {
        if (!h5 || !h5.File || !h5.FS || typeof h5.FS.writeFile !== 'function') {
            throw new Error('h5wasm file API is unavailable');
        }
    }

    function tempPath(tag) {
        return '/tmp/unified_' + tag + '_' + Date.now() + '_' + (sequence++) + '.h5';
    }

    function openBytes(h5, bytes, tag) {
        validateH5(h5);
        const path = tempPath(tag);
        h5.FS.writeFile(path, bytes instanceof Uint8Array ? bytes : new Uint8Array(bytes));
        return { path, file: new h5.File(path, 'r') };
    }

    function closeAndUnlink(h5, opened) {
        try {
            if (opened.file && typeof opened.file.close === 'function') opened.file.close();
        } finally {
            try { h5.FS.unlink(opened.path); } catch (_) {}
        }
    }

    function matrixRows(file, path, rows, columns, required) {
        const flat = numericDataset(file, path, required);
        if (!flat.length && !required) return null;
        const shape = datasetShape(file, path);
        const out = new Array(rows);
        for (let i = 0; i < rows; i++) out[i] = new Array(columns);
        // Current DISCUS/RMC structure files store atom arrays as (N,3).
        // Keep accepting the former RMC disk orientation (3,N).
        if (shape.length === 2 && shape[0] === rows && shape[1] === columns) {
            for (let i = 0; i < rows; i++)
                for (let j = 0; j < columns; j++) out[i][j] = flat[i * columns + j];
        } else if (shape.length === 2 && shape[0] === columns && shape[1] === rows) {
            for (let i = 0; i < rows; i++)
                for (let j = 0; j < columns; j++) out[i][j] = flat[j * rows + i];
        } else if (!shape.length && flat.length === rows * columns) {
            for (let i = 0; i < rows; i++)
                for (let j = 0; j < columns; j++) out[i][j] = flat[i * columns + j];
        } else {
            throw new Error('unified HDF5: ' + path + ' has an incompatible shape');
        }
        return out;
    }

    function readStructure(h5, bytes, name) {
        const opened = openBytes(h5, bytes, 'structure');
        try {
            const file = opened.file;
            const dictionary = textDataset(file, 'entry/data/audit_conform_dict_name');
            if (dictionary && dictionary !== STRUCTURE_DICTIONARY) {
                throw new Error('unified structure: wrong dictionary ' + dictionary);
            }
            const atomCountValue = numericDataset(file, 'entry/data/number_of_atoms', true);
            const atomCount = Math.floor(atomCountValue[0]);
            if (!(atomCount > 0)) throw new Error('unified structure: invalid number_of_atoms');
            const cellLengths = numericDataset(file, 'entry/data/unit_cell_lengths', true).slice(0, 3);
            const cellAngles = numericDataset(file, 'entry/data/unit_cell_angles', true).slice(0, 3);
            const unitCells = numericDataset(file, 'entry/data/unit_cells', true).slice(0, 3).map(Math.round);
            if (cellLengths.length !== 3 || cellAngles.length !== 3 ||
                unitCells.length !== 3 || unitCells.some(v => !(v > 0))) {
                throw new Error('unified structure: invalid cell or unit_cells');
            }
            const atomPosition = matrixRows(file, 'entry/data/atom_position', atomCount, 3, true);
            let atomUnitCell = matrixRows(file, 'entry/data/atom_unit_cell', atomCount, 3, false);
            if (!atomUnitCell) atomUnitCell = Array.from({ length: atomCount }, () => [1, 1, 1]);
            const atomType = numericDataset(file, 'entry/data/atom_type', true).slice(0, atomCount).map(Math.round);
            const typeNames = textDataset(file, 'entry/data/types_names').split(/\s+/).filter(Boolean);
            if (atomType.length !== atomCount || !typeNames.length) {
                throw new Error('unified structure: invalid atom types');
            }
            const elements = atomType.map(index => {
                const element = typeNames[index - 1];
                if (!element) throw new Error('unified structure: atom_type is outside types_names');
                return element;
            });
            const coordinateFlagValues = numericDataset(
                file, 'entry/data/coordinates_are_supercell_fractional', false);
            const coordinatesAreSupercellFractional = coordinateFlagValues.length
                ? Math.round(coordinateFlagValues[0]) : 0;
            if (coordinatesAreSupercellFractional !== 0 &&
                coordinatesAreSupercellFractional !== 1) {
                throw new Error('unified structure: coordinates_are_supercell_fractional must be 0 or 1');
            }
            return {
                file: String(name || ''),
                dictionary: dictionary || STRUCTURE_DICTIONARY,
                legacyContract: !dictionary,
                atomCount,
                cellLengths,
                cellAngles,
                unitCells,
                atomPosition,
                atomUnitCell,
                atomType,
                typeNames,
                elements,
                coordinatesAreSupercellFractional,
            };
        } finally {
            closeAndUnlink(h5, opened);
        }
    }

    function reciprocalKind(axesType) {
        if (axesType === 'xyz') return 'direct';
        if (axesType === 'uvw') return 'patterson';
        return 'reciprocal';
    }

    function writeData(h5, spec) {
        validateH5(h5);
        const dims = Array.from(spec.dims || []).slice(0, 3).map(v => Math.floor(Number(v)));
        if (dims.length !== 3 || dims.some(v => !(v > 0))) {
            throw new Error('unified data: dims must contain three positive values');
        }
        const total = dims[0] * dims[1] * dims[2];
        if (!Number.isSafeInteger(total) || !spec.values || spec.values.length !== total) {
            throw new Error('unified data: values length does not match dims');
        }
        const axesType = String(spec.axesType || 'hkl');
        if (!AXES_TYPES.has(axesType)) {
            throw new Error('unified data: unsupported data_type_axes ' + axesType);
        }
        const corner = Array.from(spec.corner || [0, 0, 0]).slice(0, 3).map(Number);
        const vectors = spec.vectors;
        if (corner.length !== 3 || corner.some(v => !Number.isFinite(v)) ||
            !Array.isArray(vectors) || vectors.length !== 3 ||
            vectors.some(row => !Array.isArray(row) || row.length !== 3 ||
                row.some(v => !Number.isFinite(Number(v))))) {
            throw new Error('unified data: invalid corner or increment vectors');
        }
        const cell = Array.from(spec.cell || [1, 1, 1, 90, 90, 90]).slice(0, 6).map(Number);
        if (cell.length !== 6 || cell.some(v => !Number.isFinite(v)) ||
            cell[0] <= 0 || cell[1] <= 0 || cell[2] <= 0) {
            throw new Error('unified data: invalid unit cell');
        }
        const path = tempPath('data');
        let file = null;
        try {
            file = new h5.File(path, 'w');
            const entry = file.create_group('entry');
            const data = entry.create_group('data');
            const ds = (name, definition) =>
                data.create_dataset(Object.assign({ name }, definition));
            const stringDs = (name, value) => ds(name, fixedStr(value));
            const today = new Date().toISOString().slice(0, 10);
            const iv = new Float64Array(9);
            for (let axis = 0; axis < 3; axis++)
                for (let component = 0; component < 3; component++)
                    iv[component * 3 + axis] = Number(vectors[axis][component]);

            ds('unit_cell_lengths', { data: cell.slice(0, 3), shape: [3], dtype: '<d' });
            ds('unit_cell_angles', { data: cell.slice(3, 6), shape: [3], dtype: '<d' });
            stringDs('symmetry_space_group_name_H-M', spec.symmetryName || 'P 1');
            ds('space_group_origin', { data: [1], shape: [1], dtype: '<i' });
            stringDs('symmetry_space_group_abc', 'abc');
            ds('space_group_symop_number', { data: [1], shape: [1], dtype: '<i' });
            // This is the unchanged unified diffuse-data layout. Structure
            // writers use the separate external (nsym,4,3) contract.
            ds('space_group_symop_operation_mat', {
                data: [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0],
                shape: [3, 4, 1], dtype: '<d',
            });
            stringDs('data_type_experiment', spec.experiment || 'calculated');
            stringDs('data_type_style', spec.style || 'single_diffraction');
            stringDs('data_type_axes', axesType);
            stringDs('data_type_content', spec.content || 'intensity');
            stringDs('data_type_reciprocal', spec.reciprocal || reciprocalKind(axesType));
            stringDs('data_type_with_bragg', spec.withBragg || 'unknown');
            stringDs('data_type_symmetrized', spec.symmetrized || 'none');
            stringDs('data_type_number', 'real');
            stringDs('data_radiation', spec.radiation || 'unknown');
            stringDs('data_rad_symbol', spec.radiationSymbol || 'unknown');
            ds('data_rad_length', { data: [0, 0, 0], shape: [3], dtype: '<d' });
            ds('data_dimension', { data: Int32Array.from(dims), shape: [3], dtype: '<i' });
            ds('data_axes', {
                data: Int32Array.from(spec.axes || [1, 2, 3]), shape: [3], dtype: '<i',
            });
            ds('data_corner', { data: corner, shape: [3], dtype: '<d' });
            ds('data_increment_vector', { data: iv, shape: [3, 3], dtype: '<d' });
            ds('data_values', {
                data: spec.values instanceof Float64Array
                    ? spec.values : Float64Array.from(spec.values),
                shape: dims, dtype: '<d',
            });
            stringDs('audit_conform_dict_name', DATA_DICTIONARY);
            stringDs('audit_conform_dict_version', spec.dictionaryVersion || '0.0.0');
            stringDs('audit_creation_date', today);
            stringDs('audit_creation_method', spec.creationMethod || 'RMCProfile browser utility');
            stringDs('audit_author_name', spec.authorName || 'RMCProfile');
            file.flush();
            file.close();
            file = null;
            const bytes = h5.FS.readFile(path);
            return Uint8Array.from(bytes);
        } finally {
            try { if (file) file.close(); } catch (_) {}
            try { h5.FS.unlink(path); } catch (_) {}
        }
    }

    return Object.freeze({
        DATA_DICTIONARY,
        STRUCTURE_DICTIONARY,
        LEGACY_DATA_DICTIONARY,
        readStructure,
        writeData,
    });
}));
