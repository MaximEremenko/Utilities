import numpy as np
from scipy.stats import gaussian_kde
from py_mini_racer import MiniRacer

# Load the updated KDE.js code
with open('KDE.js', 'r') as f:
    js_code = f.read()

ctx = MiniRacer()
ctx.eval(js_code)

# Generate 3-dimensional data (shape: [d][n])
np.random.seed(42)
data = np.random.normal(0, 1, size=(2, 100))
points = np.random.normal(0, 1, size=(2, 10))

# Reference SciPy KDE
scipy_kde = gaussian_kde(data, bw_method='scott')
scipy_results = scipy_kde.evaluate(points)

# Convert arrays to lists for JavaScript
data_js = data.tolist()     # list of 3 lists, each with 100 samples
points_js = points.tolist() # list of 3 lists, each with 10 samples

# Define dataset and instantiate KDE in the JS context.
# Also define a helper function evaluateKDE to call kde.evaluate with the proper context.
ctx.eval(f'''
  // Define the dataset (a 2D array with shape [d][n])
  var dataset3D = {data_js};
  // Instantiate GaussianKDE
  var kde = new GaussianKDE(dataset3D, 'scott');
  console.log("GaussianKDE instantiated with d =", kde.d, "and n =", kde.n);
  
  // Define a helper function that preserves the correct context for evaluating the KDE
  function evaluateKDE(points) {{
    return kde.evaluate(points);
  }}
''')

# Call the helper function to evaluate the KDE on the provided points.
js_results = ctx.call("evaluateKDE", points_js)

# Convert JavaScript result to a NumPy array and compare with SciPy's results
js_results = np.array(js_results)
np.testing.assert_allclose(js_results, scipy_results, rtol=1e-4, atol=1e-6)
print("JavaScript GaussianKDE matches SciPy KDE in 3D.")
