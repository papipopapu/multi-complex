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

For a holomorphic function $f$, the Taylor series around a real point $x_0$ can be written as:

$$f(x_0 + h \cdot i_1 + \ldots + h \cdot i_n) = \sum_{k=0}^{\infty} \left( \sum_{l=1}^{n} i_l \cdot h \right)^k \frac{f^{(k)}(x_0)}{k!}$$

Using the multinomial theorem to expand $\left( \sum_{l=1}^{n} i_l \cdot h \right)^k$, the n-th derivative is the only term containing $h^n \cdot \prod_{l=1}^{n} i_l$ — this occurs uniquely when $k_1 = k_2 = \ldots = k_n = 1$ in the multinomial expansion. If we ignore terms $O(h^{n+2})$, the $i_1 \cdot i_2 \cdot \ldots \cdot i_n$ product uniquely appears in the $(i_1 + \ldots + i_n)^n$ term, so the real coefficient of this imaginary direction depends only on $f^{(n)}(x_0)$.

The function $Im_{1 \ldots n}$ retrieves the real component corresponding to $x_{2^n}$ by successively extracting imaginary parts:

$$\text{Im}_{1 \ldots n}(\zeta_n) = \text{Im}_1(\text{Im}_2(\ldots(\text{Im}_n(\zeta_n))\ldots))$$

Noting that $\frac{n!}{1! \cdot \ldots \cdot 1!} = n!$, the n-th order derivative can be computed with approximation error $O(h^2)$:

$$f^{(n)}(x_0) = \frac{\text{Im}_{1 \ldots n}\left( f\left( x_0 + \sum_{k=1}^{n} h \cdot i_k \right) \right)}{h^n} + O(h^2)$$

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

## Notes

- Inverse functions (logarithm, inverse trigonometric) use branch cuts that are only bijective for the first complex level (C₁)
- The implementation prioritizes performance through template optimizations
- Some improvements are possible, but overall the implementation is robust and efficient
