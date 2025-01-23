package com.gnssAug.helper.lambdaNew;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.FastMath;

public class GammaIncompleteInverse {

    /**
     * Computes the inverse of the regularized incomplete gamma function.
     *
     * @param p The probability value (regularized incomplete gamma value, between 0 and 1).
     * @param a The shape parameter of the gamma distribution.
     * @return The value x such that P(a, x) = p, where P(a, x) is the regularized incomplete gamma function.
     * @throws IllegalArgumentException if p is not in the range [0, 1] or if a <= 0.
     */
    public static double gammaincinv(double p, double a) {
        if (p < 0.0 || p > 1.0) {
            throw new IllegalArgumentException("p must be in the range [0, 1]");
        }
        if (a <= 0.0) {
            throw new IllegalArgumentException("Shape parameter a must be positive");
        }

        // Tolerance for convergence
        final double tolerance = 1e-8;

        // Initial guess using an approximation (scaled value)
        double x = initialGuess(p, a);

        // Newton-Raphson iteration to refine the guess
        for (int i = 0; i < 100; i++) { // Limit iterations to prevent infinite loops
            double gammaP = Gamma.regularizedGammaP(a, x);
            double gammaPDeriv = FastMath.exp(-x + (a - 1) * FastMath.log(x) - Gamma.logGamma(a));

            // Update using Newton-Raphson
            double delta = (gammaP - p) / gammaPDeriv;
            x -= delta;

            // Convergence check
            if (Math.abs(delta) < tolerance * x) {
                return x;
            }

            // Prevent x from becoming negative or zero
            if (x <= 0) {
                x = tolerance;
            }
        }

        throw new ArithmeticException("Failed to converge to a solution for gammaincinv");
    }

    /**
     * Provides an initial guess for the inverse of the regularized incomplete gamma function.
     *
     * @param p The probability value (regularized incomplete gamma value).
     * @param a The shape parameter of the gamma distribution.
     * @return Initial guess for the Newton-Raphson method.
     */
    private static double initialGuess(double p, double a) {
        if (p > 0.5) {
            return a * FastMath.pow(-FastMath.log(1 - p), 1 / a); // Approximation for large p
        } else {
            return a * FastMath.pow(-FastMath.log(p), 1 / a); // Approximation for small p
        }
    }

   
}