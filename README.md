# Multi-Complex

A header-only C++ library for **multicomplex numbers** implemented via recursive template metaprogramming. Enables fast, machine-precision numerical differentiation with minimal overhead.

---

## Overview

Multicomplex numbers are a generalization of complex numbers that introduce multiple independent imaginary units (i₁, i₂, ..., iₙ). This mathematical structure enables **automatic numerical differentiation** to machine precision, avoiding the numerical instabilities inherent in finite difference methods.

### Key Features

- **Header-only**: Just include `MultiComplex.hpp` — no linking required
- **Template metaprogramming**: Compile-time recursion enables aggressive optimization
- **High precision**: Achieves machine-limit precision for derivatives
- **Full math support**: Includes trigonometric, hyperbolic, exponential, logarithmic, and inverse functions
- **Arbitrary order**: Compute N-th order derivatives using level-N multicomplex numbers

### How It Works

For a holomorphic function f, the Taylor series around a real point x₀ can be written as:

```
f(x₀ + h·i₁ + ... + h·iₙ) = Σₖ₌₀^∞ (Σₗ₌₁ⁿ iₗ·h)ᵏ · f⁽ᵏ⁾(x₀)/k!
```

Using the multinomial theorem to expand `(Σₗ₌₁ⁿ iₗ·h)ᵏ`, the n-th derivative is the only term containing `hⁿ·(∏ₗⁿ iₗ)` — this occurs uniquely when k₁ = k₂ = ... = kₙ = 1 in the multinomial expansion.

The function `Im₁...ₙ` retrieves the real component corresponding to x₂ⁿ by successively extracting imaginary parts:

```
Im₁...ₙ(ζₙ) = Im₁(Im₂(...(Imₙ(ζₙ))...))
```

Therefore, the n-th order derivative can be computed as:

```
f⁽ⁿ⁾(x₀) = Im₁...ₙ(f(x₀ + Σₖ₌₁ⁿ h·iₖ)) / hⁿ
```

---

## Usage Example

```cpp
#include "MultiComplex.hpp"
#include <iostream>
#include <cmath>

using namespace MComplex;

// Function to differentiate: f(x) = sin(x)
template<typename T>
T f(const T& x) {
    return sin(x);
}

int main() {
    // Point at which to compute derivative
    double x0 = 1.0;
    
    // Small perturbation (can be very small without numerical issues)
    double h = 1e-100;
    
    // Create multicomplex number: x0 + h*i
    MultiComplex<1, double> z{x0, h};
    
    // Evaluate function
    MultiComplex<1, double> result = f(z);
    
    // Extract derivative: imag(f(z)) / h = f'(x0)
    double derivative = result.imag() / h;
    
    // Compare with analytical derivative: cos(x0)
    double analytical = std::cos(x0);
    
    std::cout << "Computed derivative: " << derivative << std::endl;
    std::cout << "Analytical derivative: " << analytical << std::endl;
    std::cout << "Error: " << std::abs(derivative - analytical) << std::endl;
    
    return 0;
}
```

**Output:**
```
Computed derivative: 0.540302
Analytical derivative: 0.540302
Error: 0
```

---

## References

- Lantoine, G., Russell, R. P., & Dargent, T. (2012). [Using Multicomplex Variables for Automatic Computation of High-Order Derivatives](https://dl.acm.org/doi/abs/10.1145/3378538). *ACM Transactions on Mathematical Software*.

---

## Notes

- Inverse functions (logarithm, inverse trigonometric) use branch cuts that are only bijective for the first complex level (C₁)
- The implementation prioritizes performance through template optimizations
- Some improvements are possible, but overall the implementation is robust and efficient
