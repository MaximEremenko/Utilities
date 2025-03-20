    // --- Helper Functions ---
    function weightedMean(data, weights) {
      const d = data.length;
      const n = data[0].length;
      const mean = new Array(d);
      for (let i = 0; i < d; i++) {
        let sum = 0;
        const row = data[i];
        for (let j = 0; j < n; j++) {
          sum += row[j] * weights[j];
        }
        mean[i] = sum;
      }
      return mean;
    }
    function weightedCovariance(data, weights) {
      const d = data.length;
      const n = data[0].length;
      const mean = weightedMean(data, weights);
      let sumW2 = 0;
      for (let j = 0; j < n; j++) {
        sumW2 += weights[j] * weights[j];
      }
      const cov = new Array(d);
      for (let i = 0; i < d; i++) {
        cov[i] = new Array(d).fill(0);
        const rowi = data[i];
        for (let j = 0; j < d; j++) {
          const rowj = data[j];
          let s = 0;
          for (let k = 0; k < n; k++) {
            s += weights[k] * (rowi[k] - mean[i]) * (rowj[k] - mean[j]);
          }
          cov[i][j] = s / (1 - sumW2);
        }
      }
      return cov;
    }
    function choleskyDecomposition(matrix) {
      const n = matrix.length;
      const L = new Array(n);
      for (let i = 0; i < n; i++) {
        L[i] = new Array(n).fill(0);
        for (let j = 0; j <= i; j++) {
          let sum = 0;
          for (let k = 0; k < j; k++) {
            sum += L[i][k] * L[j][k];
          }
          if (i === j) {
            const diag = matrix[i][i] - sum;
            if (diag <= 0) {
              throw new Error("Matrix is not positive definite");
            }
            L[i][j] = Math.sqrt(diag);
          } else {
            L[i][j] = (matrix[i][j] - sum) / L[j][j];
          }
        }
      }
      return L;
    }
    function solveLower(L, b) {
      const n = L.length;
      const y = new Array(n);
      for (let i = 0; i < n; i++) {
        let sum = 0;
        for (let j = 0; j < i; j++) {
          sum += L[i][j] * y[j];
        }
        y[i] = (b[i] - sum) / L[i][i];
      }
      return y;
    }
    function determinantFromCholesky(L) {
      let prod = 1;
      const n = L.length;
      for (let i = 0; i < n; i++) {
        prod *= L[i][i];
      }
      return prod * prod;
    }
    function logSumExp(arr) {
      let maxVal = -Infinity;
      const len = arr.length;
      for (let i = 0; i < len; i++) {
        if (arr[i] > maxVal) maxVal = arr[i];
      }
      let sumExp = 0;
      for (let i = 0; i < len; i++) {
        sumExp += Math.exp(arr[i] - maxVal);
      }
      return maxVal + Math.log(sumExp);
    }
    function erf(x) {
      const sign = x >= 0 ? 1 : -1;
      x = Math.abs(x);
      const a1 = 0.254829592,
            a2 = -0.284496736,
            a3 = 1.421413741,
            a4 = -1.453152027,
            a5 = 1.061405429,
            p = 0.3275911;
      const t = 1 / (1 + p * x);
      const y = 1 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);
      return sign * y;
    }
    function normalCDF(x) {
      return 0.5 * (1 + erf(x / Math.sqrt(2)));
    }
    function gaussianRandom() {
      let u = 0, v = 0;
      while (u === 0) u = Math.random();
      while (v === 0) v = Math.random();
      return Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * v);
    }
    function sampleIndex(cumWeights, r) {
      let lo = 0, hi = cumWeights.length - 1;
      while (lo < hi) {
        const mid = (lo + hi) >> 1;
        if (r > cumWeights[mid]) {
          lo = mid + 1;
        } else {
          hi = mid;
        }
      }
      return lo;
    }

    // --- GPU.js Setup ---
    // Try to initialize GPU.js using either GPU or GPU.GPU.
    let gpuInstance = null;
    try {
      if (typeof GPU === "function") {
        gpuInstance = new GPU();
      } else if (typeof GPU === "object" && typeof GPU.GPU === "function") {
        gpuInstance = new GPU.GPU();
      } else {
        console.warn("GPU.js is not available as a constructor; falling back to CPU.");
        gpuInstance = null;
      }
    } catch (e) {
      console.warn("GPU.js could not be initialized, falling back to CPU.", e);
      gpuInstance = null;
    }

    // --- GaussianKDE Class ---
    class GaussianKDE {
      constructor(dataset, bw_method = "scott", weights = null) {
        if (!Array.isArray(dataset) || dataset.length === 0 || !Array.isArray(dataset[0])) {
          throw new Error("Dataset must be a 2D array with shape [d][n].");
        }
        this.dataset = dataset;
        this.d = dataset.length;
        this.n = dataset[0].length;
        if (this.n < 2) {
          throw new Error("Dataset must contain multiple samples.");
        }
        if (weights) {
          if (weights.length !== this.n) {
            throw new Error("Weights array length must equal the number of samples.");
          }
          let sumWeights = 0;
          for (let i = 0; i < weights.length; i++) {
            sumWeights += weights[i];
          }
          this._weights = weights.map(w => w / sumWeights);
        } else {
          this._weights = new Array(this.n).fill(1 / this.n);
        }
        let sumW2 = 0;
        for (let i = 0; i < this._weights.length; i++) {
          sumW2 += this._weights[i] * this._weights[i];
        }
        this._neff = 1 / sumW2;
        if (this.d > this.n) {
          throw new Error("Number of dimensions is greater than number of samples. Consider transposing the input data.");
        }
        this.gpu = gpuInstance;
        this.setBandwidth(bw_method);
      }
      scottsFactor() {
        return Math.pow(this._neff, -1 / (this.d + 4));
      }
      silvermanFactor() {
        return Math.pow((this._neff * (this.d + 2)) / 4, -1 / (this.d + 4));
      }
      setBandwidth(bw_method) {
        if (bw_method === null || bw_method === undefined) {
          // keep current covarianceFactor
        } else if (bw_method === "scott") {
          this.covarianceFactor = () => this.scottsFactor();
        } else if (bw_method === "silverman") {
          this.covarianceFactor = () => this.silvermanFactor();
        } else if (typeof bw_method === "number") {
          this.covarianceFactor = () => bw_method;
        } else if (typeof bw_method === "function") {
          this.covarianceFactor = () => bw_method(this);
        } else {
          throw new Error("bw_method should be 'scott', 'silverman', a number, or a function.");
        }
        this._computeCovariance();
      }
      _computeCovariance() {
        this.factor = this.covarianceFactor();
        this._dataCovariance = weightedCovariance(this.dataset, this._weights);
        this.covariance = this._dataCovariance.map(row =>
          row.map(val => val * this.factor * this.factor)
        );
        this.choCov = choleskyDecomposition(this.covariance);
        const detCov = determinantFromCholesky(this.choCov);
        this._normConst = Math.pow(2 * Math.PI, this.d / 2) * Math.sqrt(detCov);
        this.whitenedData = new Array(this.n);
        for (let i = 0; i < this.n; i++) {
          const sample = new Array(this.d);
          for (let j = 0; j < this.d; j++) {
            sample[j] = this.dataset[j][i];
          }
          this.whitenedData[i] = solveLower(this.choCov, sample);
        }
      }
      evaluate(points) {
        let pts = Array.isArray(points[0]) ? points : points.map(x => [x]);
        const m = pts[0].length;
        if (pts.length !== this.d) {
          throw new Error(`Points have dimension ${pts.length}, but data has dimension ${this.d}`);
        }
        const pointsArr = new Array(m);
        for (let j = 0; j < m; j++) {
          const point = new Array(this.d);
          for (let i = 0; i < this.d; i++) {
            point[i] = pts[i][j];
          }
          pointsArr[j] = point;
        }
        const L = this.choCov;
        const qWhiteArr = pointsArr.map(point => solveLower(L, point));
        if (this.gpu) {
          const evaluateKernel = this.gpu.createKernel(function(qWhite, whitenedData, weights, d, n, normConst) {
            let sumDensity = 0;
            for (let i = 0; i < n; i++) {
              let quad = 0;
              for (let k = 0; k < d; k++) {
                const diff = qWhite[this.thread.x][k] - whitenedData[i][k];
                quad += diff * diff;
              }
              sumDensity += weights[i] * Math.exp(-0.5 * quad);
            }
            return sumDensity / normConst;
          }).setOutput([m]);
          const result = evaluateKernel(qWhiteArr, this.whitenedData, this._weights, this.d, this.n, this._normConst);
          return Array.from(result);
        } else {
          const results = new Array(m);
          for (let j = 0; j < m; j++) {
            const qWhite = qWhiteArr[j];
            let sumDensity = 0;
            for (let i = 0; i < this.n; i++) {
              let quad = 0;
              const sampleWhite = this.whitenedData[i];
              for (let k = 0; k < this.d; k++) {
                const diff = qWhite[k] - sampleWhite[k];
                quad += diff * diff;
              }
              sumDensity += this._weights[i] * Math.exp(-0.5 * quad);
            }
            results[j] = sumDensity / this._normConst;
          }
          return results;
        }
      }
      logpdf(points) {
        let pts = Array.isArray(points[0]) ? points : points.map(x => [x]);
        const m = pts[0].length;
        if (pts.length !== this.d) {
          throw new Error(`Points have dimension ${pts.length}, but data has dimension ${this.d}`);
        }
        const pointsArr = new Array(m);
        for (let j = 0; j < m; j++) {
          const point = new Array(this.d);
          for (let i = 0; i < this.d; i++) {
            point[i] = pts[i][j];
          }
          pointsArr[j] = point;
        }
        const L = this.choCov;
        const qWhiteArr = pointsArr.map(point => solveLower(L, point));
        let logResults = new Array(m);
        const normLog = Math.log(this._normConst);
        if (this.gpu) {
          const logpdfKernel = this.gpu.createKernel(function(qWhite, whitenedData, weights, d, n) {
            let quad = 0;
            for (let k = 0; k < d; k++) {
              const diff = qWhite[this.thread.x][k] - whitenedData[this.thread.y][k];
              quad += diff * diff;
            }
            return Math.log(weights[this.thread.y]) - 0.5 * quad;
          }).setOutput([m, this.n]);
          const logKernels = logpdfKernel(qWhiteArr, this.whitenedData, this._weights, this.d, this.n);
          for (let j = 0; j < m; j++) {
            const row = [];
            for (let i = 0; i < this.n; i++) {
              row.push(logKernels[j][i]);
            }
            logResults[j] = logSumExp(row) - normLog;
          }
        } else {
          for (let j = 0; j < m; j++) {
            const qWhite = qWhiteArr[j];
            const logKernels = new Array(this.n);
            for (let i = 0; i < this.n; i++) {
              let quad = 0;
              const sampleWhite = this.whitenedData[i];
              for (let k = 0; k < this.d; k++) {
                const diff = qWhite[k] - sampleWhite[k];
                quad += diff * diff;
              }
              logKernels[i] = Math.log(this._weights[i]) - 0.5 * quad;
            }
            logResults[j] = logSumExp(logKernels) - normLog;
          }
        }
        return logResults;
      }
      resample(size = Math.round(this._neff)) {
        const d = this.d, n = this.n;
        const samples = new Array(d);
        for (let i = 0; i < d; i++) {
          samples[i] = new Array(size);
        }
        const cumWeights = new Array(n);
        let acc = 0;
        for (let i = 0; i < n; i++) {
          acc += this._weights[i];
          cumWeights[i] = acc;
        }
        const L = this.choCov;
        for (let i = 0; i < size; i++) {
          const r = Math.random();
          const idx = sampleIndex(cumWeights, r);
          const point = new Array(d);
          const z = new Array(d);
          for (let j = 0; j < d; j++) {
            z[j] = gaussianRandom();
          }
          for (let j = 0; j < d; j++) {
            let noise = 0;
            for (let k = 0; k <= j; k++) {
              noise += L[j][k] * z[k];
            }
            point[j] = this.dataset[j][idx] + noise;
          }
          for (let j = 0; j < d; j++) {
            samples[j][i] = point[j];
          }
        }
        return samples;
      }
      integrateGaussian(mean, cov) {
        if (mean.length !== this.d || cov.length !== this.d || cov[0].length !== this.d) {
          throw new Error("Mean and covariance dimensions must match the KDE.");
        }
        const covSum = new Array(this.d);
        for (let i = 0; i < this.d; i++) {
          covSum[i] = new Array(this.d);
          for (let j = 0; j < this.d; j++) {
            covSum[i][j] = cov[i][j] + this.covariance[i][j];
          }
        }
        const Lsum = choleskyDecomposition(covSum);
        const detCovSum = determinantFromCholesky(Lsum);
        const normConst = Math.pow(2 * Math.PI, this.d / 2) * Math.sqrt(detCovSum);
        let result = 0;
        for (let i = 0; i < this.n; i++) {
          const diff = new Array(this.d);
          for (let k = 0; k < this.d; k++) {
            diff[k] = mean[k] - this.dataset[k][i];
          }
          const y = solveLower(Lsum, diff);
          let quad = 0;
          for (let k = 0; k < this.d; k++) {
            quad += y[k] * y[k];
          }
          result += this._weights[i] * Math.exp(-0.5 * quad);
        }
        return result / normConst;
      }
      integrateKDE(other) {
        if (other.d !== this.d) {
          throw new Error("KDEs must have the same number of dimensions.");
        }
        let small, large;
        if (other.n < this.n) {
          small = other;
          large = this;
        } else {
          small = this;
          large = other;
        }
        const covSum = new Array(this.d);
        for (let i = 0; i < this.d; i++) {
          covSum[i] = new Array(this.d);
          for (let j = 0; j < this.d; j++) {
            covSum[i][j] = small.covariance[i][j] + large.covariance[i][j];
          }
        }
        const Lsum = choleskyDecomposition(covSum);
        const detCovSum = determinantFromCholesky(Lsum);
        const normConst = Math.pow(2 * Math.PI, this.d / 2) * Math.sqrt(detCovSum);
        let result = 0;
        for (let i = 0; i < small.n; i++) {
          const sample_i = new Array(this.d);
          for (let k = 0; k < this.d; k++) {
            sample_i[k] = small.dataset[k][i];
          }
          let innerSum = 0;
          for (let j = 0; j < large.n; j++) {
            const diff = new Array(this.d);
            for (let k = 0; k < this.d; k++) {
              diff[k] = large.dataset[k][j] - sample_i[k];
            }
            const y = solveLower(Lsum, diff);
            let quad = 0;
            for (let k = 0; k < this.d; k++) {
              quad += y[k] * y[k];
            }
            innerSum += large._weights[j] * Math.exp(-0.5 * quad);
          }
          result += small._weights[i] * innerSum;
        }
        return result / normConst;
      }
      marginal(dimensions) {
        let dims = typeof dimensions === "number" ? [dimensions] : dimensions.slice();
        dims = dims.map(d => (d < 0 ? this.d + d : d));
        if (new Set(dims).size !== dims.length) {
          throw new Error("Dimensions must be unique.");
        }
        dims.forEach(d => {
          if (d < 0 || d >= this.d) {
            throw new Error(`Dimension ${d} is invalid for a ${this.d}-dimensional dataset.`);
          }
        });
        const newData = dims.map(d => this.dataset[d].slice());
        return new GaussianKDE(newData, this.covarianceFactor(), this._weights.slice());
      }
      get weights() {
        return this._weights;
      }
      get neff() {
        return this._neff;
      }
    }