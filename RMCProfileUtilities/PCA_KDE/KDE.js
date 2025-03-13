// Helper functions

/**
 * Compute weighted mean for a d x n data array.
 * @param {number[][]} data - Array of arrays with shape [d][n].
 * @param {number[]} weights - Array of length n.
 * @returns {number[]} mean vector of length d.
 */
function weightedMean(data, weights) {
  const d = data.length;
  const n = data[0].length;
  const mean = new Array(d).fill(0);
  for (let i = 0; i < d; i++) {
    let sum = 0;
    for (let j = 0; j < n; j++) {
      sum += data[i][j] * weights[j];
    }
    mean[i] = sum;
  }
  return mean;
}

/**
 * Compute weighted covariance matrix for a d x n data array.
 * Uses the unbiased formula: Cov = Σ w_i (x_i - mean)(x_i - mean)^T / (1 - Σw_i²)
 * @param {number[][]} data - Array of arrays with shape [d][n].
 * @param {number[]} weights - Array of length n.
 * @returns {number[][]} d x d covariance matrix.
 */
function weightedCovariance(data, weights) {
  const d = data.length;
  const n = data[0].length;
  const mean = weightedMean(data, weights);
  const sumW2 = weights.reduce((acc, w) => acc + w * w, 0);
  const cov = Array(d)
    .fill(0)
    .map(() => Array(d).fill(0));
  for (let i = 0; i < d; i++) {
    for (let j = 0; j < d; j++) {
      let s = 0;
      for (let k = 0; k < n; k++) {
        s += weights[k] * (data[i][k] - mean[i]) * (data[j][k] - mean[j]);
      }
      cov[i][j] = s / (1 - sumW2);
    }
  }
  return cov;
}

/**
 * Performs a Cholesky decomposition on a symmetric positive definite matrix.
 * Returns the lower triangular matrix L such that matrix = L * L^T.
 * @param {number[][]} matrix - d x d symmetric positive definite matrix.
 * @returns {number[][]} Lower triangular matrix L.
 */
function choleskyDecomposition(matrix) {
  const n = matrix.length;
  const L = Array(n)
    .fill(0)
    .map(() => Array(n).fill(0));
  for (let i = 0; i < n; i++) {
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

/**
 * Given lower-triangular L and vector b, solves L*y = b.
 * @param {number[][]} L - Lower triangular matrix.
 * @param {number[]} b - Right-hand side vector.
 * @returns {number[]} Solution vector y.
 */
function solveLower(L, b) {
  const n = L.length;
  const y = new Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    let sum = 0;
    for (let j = 0; j < i; j++) {
      sum += L[i][j] * y[j];
    }
    y[i] = (b[i] - sum) / L[i][i];
  }
  return y;
}

/**
 * Compute the determinant of a matrix from its Cholesky factor.
 * For L lower triangular with matrix = L * L^T, det(matrix) = (product(diag(L)))^2.
 * @param {number[][]} L - Lower triangular matrix.
 * @returns {number} Determinant.
 */
function determinantFromCholesky(L) {
  const prod = L.reduce((acc, row, i) => acc * row[i], 1);
  return prod * prod;
}

/**
 * Computes the log-sum-exp of an array in a numerically stable way.
 * @param {number[]} arr
 * @returns {number}
 */
function logSumExp(arr) {
  const maxVal = Math.max(...arr);
  const sumExp = arr.reduce((acc, x) => acc + Math.exp(x - maxVal), 0);
  return maxVal + Math.log(sumExp);
}

/**
 * Approximation of the error function.
 * Source: Numerical Recipes (using a common approximation)
 * @param {number} x
 * @returns {number}
 */
function erf(x) {
  const sign = x >= 0 ? 1 : -1;
  x = Math.abs(x);
  const a1 = 0.254829592;
  const a2 = -0.284496736;
  const a3 = 1.421413741;
  const a4 = -1.453152027;
  const a5 = 1.061405429;
  const p = 0.3275911;
  const t = 1 / (1 + p * x);
  const y =
    1 -
    ((((a5 * t + a4) * t + a3) * t + a2) * t + a1) *
      t *
      Math.exp(-x * x);
  return sign * y;
}

/**
 * Standard normal cumulative distribution function.
 * @param {number} x
 * @returns {number}
 */
function normalCDF(x) {
  return 0.5 * (1 + erf(x / Math.sqrt(2)));
}

/**
 * Generate a standard normally distributed random number
 * using the Box-Muller transform.
 * @returns {number}
 */
function gaussianRandom() {
  let u = 0,
    v = 0;
  while (u === 0) u = Math.random(); // Convert [0,1) to (0,1)
  while (v === 0) v = Math.random();
  return Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * v);
}

// Main class for Gaussian Kernel Density Estimation

class GaussianKDE {
  /**
   * Constructor for GaussianKDE.
   *
   * @param {number[][]} dataset - Data array of shape [d][n] (d dimensions, n samples).
   * @param {('scott'|'silverman'|number|Function)} [bw_method='scott'] - Bandwidth method.
   * @param {number[]} [weights=null] - Optional weights (length n). If not provided, uniform weights are used.
   */
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
    // Set weights (normalize if provided)
    if (weights) {
      if (weights.length !== this.n) {
        throw new Error("Weights array length must equal the number of samples.");
      }
      const sumWeights = weights.reduce((a, b) => a + b, 0);
      this._weights = weights.map((w) => w / sumWeights);
    } else {
      this._weights = Array(this.n).fill(1 / this.n);
    }
    // Effective number of samples
    this._neff = 1 / this._weights.reduce((acc, w) => acc + w * w, 0);
    if (this.d > this.n) {
      throw new Error(
        "Number of dimensions is greater than number of samples. Consider transposing the input data."
      );
    }
    // Set bandwidth method and compute covariance
    this.setBandwidth(bw_method);
  }

  /**
   * Scott's factor: neff^(-1/(d+4))
   * @returns {number}
   */
  scottsFactor() {
    return Math.pow(this.neff, -1 / (this.d + 4));
  }

  /**
   * Silverman's factor: (neff*(d+2)/4)^(-1/(d+4))
   * @returns {number}
   */
  silvermanFactor() {
    return Math.pow((this.neff * (this.d + 2)) / 4, -1 / (this.d + 4));
  }

  /**
   * Set the bandwidth method.
   * @param {('scott'|'silverman'|number|Function)} bw_method
   */
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

  /**
   * Compute covariance matrix (scaled by the bandwidth factor), its Cholesky factor,
   * and precompute the whitened dataset.
   */
  _computeCovariance() {
    this.factor = this.covarianceFactor();
    this._dataCovariance = weightedCovariance(this.dataset, this._weights);
    // Scale covariance by factor^2
    this.covariance = this._dataCovariance.map((row) =>
      row.map((val) => val * this.factor * this.factor)
    );
    // Cholesky decomposition of the covariance matrix
    this.choCov = choleskyDecomposition(this.covariance);
    const detCov = determinantFromCholesky(this.choCov);
    this._normConst = Math.pow(2 * Math.PI, this.d / 2) * Math.sqrt(detCov);

    // Precompute the whitened dataset: for each sample x, compute y = L^{-1}x.
    this.whitenedData = new Array(this.n);
    for (let i = 0; i < this.n; i++) {
      const sample = new Array(this.d);
      for (let j = 0; j < this.d; j++) {
        sample[j] = this.dataset[j][i];
      }
      this.whitenedData[i] = solveLower(this.choCov, sample);
    }
  }

  /**
   * Evaluate the estimated density on a set of points.
   * @param {number[][]} points - Array of points with shape [d][m]. For a single point, provide a [d] array.
   * @returns {number[]} Density estimates (array of length m).
   */
  evaluate(points) {
    let pts;
    if (!Array.isArray(points[0])) {
      pts = points.map((x) => [x]);
    } else {
      pts = points;
    }
    const m = pts[0].length;
    if (pts.length !== this.d) {
      throw new Error(`Points have dimension ${pts.length}, but data has dimension ${this.d}`);
    }
    const results = new Array(m).fill(0);
    // Cache local variables for speed.
    const cho = this.choCov;
    const d = this.d, n = this.n;
    const weights = this._weights;
    const normConst = this._normConst;
    const whitenedData = this.whitenedData;

    // For each evaluation point:
    const query = new Array(d);
    for (let j = 0; j < m; j++) {
      for (let k = 0; k < d; k++) {
        query[k] = pts[k][j];
      }
      // Compute the whitened query point: L^{-1}(query)
      const qWhite = solveLower(cho, query);
      let sum = 0;
      // Compute squared distance in whitened space to each sample.
      for (let i = 0; i < n; i++) {
        let quad = 0;
        const sampleWhite = whitenedData[i];
        for (let k = 0; k < d; k++) {
          const diff = qWhite[k] - sampleWhite[k];
          quad += diff * diff;
        }
        sum += weights[i] * Math.exp(-0.5 * quad);
      }
      results[j] = sum / normConst;
    }
    return results;
  }

  /**
   * Evaluate the log of the estimated density on a set of points.
   * Uses a numerically stable log-sum-exp over the kernels.
   * @param {number[][]} points - Array of points with shape [d][m] or a single [d] array.
   * @returns {number[]} Log-density estimates.
   */
  logpdf(points) {
    let pts;
    if (!Array.isArray(points[0])) {
      pts = points.map((x) => [x]);
    } else {
      pts = points;
    }
    const m = pts[0].length;
    if (pts.length !== this.d) {
      throw new Error(`Points have dimension ${pts.length}, but data has dimension ${this.d}`);
    }
    const logResults = new Array(m).fill(0);
    const cho = this.choCov;
    const d = this.d, n = this.n;
    const weights = this._weights;
    const normLog = Math.log(this._normConst);
    const whitenedData = this.whitenedData;

    const query = new Array(d);
    for (let j = 0; j < m; j++) {
      for (let k = 0; k < d; k++) {
        query[k] = pts[k][j];
      }
      const qWhite = solveLower(cho, query);
      const logKernels = new Array(n);
      for (let i = 0; i < n; i++) {
        let quad = 0;
        const sampleWhite = whitenedData[i];
        for (let k = 0; k < d; k++) {
          const diff = qWhite[k] - sampleWhite[k];
          quad += diff * diff;
        }
        logKernels[i] = Math.log(weights[i]) - 0.5 * quad;
      }
      logResults[j] = logSumExp(logKernels) - normLog;
    }
    return logResults;
  }

  /**
   * Resample from the estimated density.
   * @param {number} [size] - Number of samples to draw (default is effective sample size).
   * @returns {number[][]} New dataset samples (shape [d][size]).
   */
  resample(size = Math.round(this.neff)) {
    const samples = Array(this.d)
      .fill(0)
      .map(() => Array(size).fill(0));
    // Precompute cumulative weights for sampling.
    const cumWeights = [];
    this._weights.reduce((acc, w) => {
      acc += w;
      cumWeights.push(acc);
      return acc;
    }, 0);
    const sampleOne = () => {
      // Choose an index based on weights.
      const r = Math.random();
      let idx = 0;
      while (r > cumWeights[idx] && idx < cumWeights.length - 1) {
        idx++;
      }
      // Add noise to the chosen sample.
      const point = new Array(this.d).fill(0);
      const z = new Array(this.d).fill(0).map(() => gaussianRandom());
      for (let i = 0; i < this.d; i++) {
        let noise = 0;
        for (let j = 0; j <= i; j++) {
          noise += this.choCov[i][j] * z[j];
        }
        point[i] = this.dataset[i][idx] + noise;
      }
      return point;
    };

    for (let i = 0; i < size; i++) {
      const pt = sampleOne();
      for (let d = 0; d < this.d; d++) {
        samples[d][i] = pt[d];
      }
    }
    return samples;
  }

  /**
   * Integrate the product of the KDE with a Gaussian (given mean and covariance).
   * That is, compute ∫ KDE(x) * N(x; mean, cov) dx.
   * @param {number[]} mean - Array of length d.
   * @param {number[][]} cov - d x d covariance matrix.
   * @returns {number} The value of the integral.
   */
  integrateGaussian(mean, cov) {
    if (mean.length !== this.d || cov.length !== this.d || cov[0].length !== this.d) {
      throw new Error("Mean and covariance dimensions must match the KDE.");
    }
    // Compute cov_sum = cov + this.covariance.
    const covSum = [];
    for (let i = 0; i < this.d; i++) {
      covSum[i] = [];
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
      const quad = y.reduce((acc, val) => acc + val * val, 0);
      result += this._weights[i] * Math.exp(-0.5 * quad);
    }
    return result / normConst;
  }

  /**
   * Integrate the product of this KDE with another KDE.
   * @param {GaussianKDE} other
   * @returns {number} The value of the integral.
   */
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
    // Compute cov_sum = small.covariance + large.covariance.
    const covSum = [];
    for (let i = 0; i < this.d; i++) {
      covSum[i] = [];
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
        const quad = y.reduce((acc, val) => acc + val * val, 0);
        innerSum += large._weights[j] * Math.exp(-0.5 * quad);
      }
      result += small._weights[i] * innerSum;
    }
    return result / normConst;
  }

  /**
   * Return a marginal KDE over the specified dimensions.
   * @param {number|number[]} dimensions - A single index or array of indices to keep.
   * @returns {GaussianKDE} New GaussianKDE instance with the marginal distribution.
   */
  marginal(dimensions) {
    let dims;
    if (typeof dimensions === "number") {
      dims = [dimensions];
    } else {
      dims = dimensions.slice();
    }
    dims = dims.map((d) => (d < 0 ? this.d + d : d));
    if (new Set(dims).size !== dims.length) {
      throw new Error("Dimensions must be unique.");
    }
    dims.forEach((d) => {
      if (d < 0 || d >= this.d) {
        throw new Error(`Dimension ${d} is invalid for a ${this.d}-dimensional dataset.`);
      }
    });
    const newData = dims.map((d) => this.dataset[d].slice());
    return new GaussianKDE(newData, this.covarianceFactor(), this._weights.slice());
  }

  // Getter for weights.
  get weights() {
    return this._weights;
  }

  // Getter for effective sample size.
  get neff() {
    return this._neff;
  }
}

// ----- Example usage -----
// Uncomment the code below to test:

// const nSamples = 200;
// const data1 = [], data2 = [];
// for (let i = 0; i < nSamples; i++) {
//   data1.push(gaussianRandom());
//   data2.push(gaussianRandom() * 0.5);
// }
// const dataset = [data1, data2];
// const kde = new GaussianKDE(dataset, "scott");
// const density = kde.evaluate([0, 0]);
// console.log("Density at [0,0]:", density);
// const newSamples = kde.resample(10);
// console.log("Resampled points:", newSamples);
