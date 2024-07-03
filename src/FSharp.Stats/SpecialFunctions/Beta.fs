﻿namespace FSharp.Stats.SpecialFunctions

open System
open FSharp.Stats

/// The beta function B(p,q), or the beta integral (also called the Eulerian integral of the first kind) is defined by
///
/// B(p, q) = (Γ(p) * Γ(q)) / Γ(p+q)
module Beta =

    let private EPS = 3.0e-8    // Precision.DoublePrecision;
    let private FPMIN = 1.0e-30 // 0.0.Increment()/eps

    ///<summary>
    /// Computes an approximation of the real value of the log beta function using approximations for the gamma function using Lanczos Coefficients described in Numerical Recipes (Press et al) 
    ///</summary>
    ///<remarks>
    /// The caller is responsible to handle edge cases such as nan, infinity, and -infinity in the input
    ///</remarks>
    /// <param name="z">The function input for approximating ln(B(z, w))</param>
    /// <param name="w">The function input for approximating ln(B(z, w))</param>
    let _betaLn z w = (Gamma._gammaLn z) + (Gamma._gammaLn w) - (Gamma._gammaLn (z+w))

    ///<summary>
    /// Computes an approximation of the real value of the beta function using approximations for the gamma function using Lanczos Coefficients described in Numerical Recipes (Press et al) 
    ///</summary>
    ///<remarks>
    /// The caller is responsible to handle edge cases such as nan, infinity, and -infinity in the input
    ///</remarks>
    /// <param name="z">The function input for approximating B(z, w)</param>
    /// <param name="w">The function input for approximating B(z, w)</param>
    let _beta z w = exp (_betaLn z w)

    ///<summary>
    /// Computes an approximation of the real value of the log beta function using approximations for the gamma function using Lanczos Coefficients described in Numerical Recipes (Press et al) 
    ///</summary>
    ///<remarks>
    /// Edge cases in the input (nan, infinity, and -infinity) are catched and handled. 
    /// This might be slower than the unchecked version `_betaLn` but does not require input sanitation to get expected results for these cases.
    ///</remarks>    
    /// <param name="z">The function input for approximating ln(B(z, w))</param>
    /// <param name="w">The function input for approximating ln(B(z, w))</param>
    let betaLn z w = (Gamma.gammaLn z) + (Gamma.gammaLn w) - (Gamma.gammaLn (z+w))

    ///<summary>
    /// Computes an approximation of the real value of the beta function using approximations for the gamma function using Lanczos Coefficients described in Numerical Recipes (Press et al) 
    ///</summary>
    ///<remarks>
    /// Edge cases in the input (nan, infinity, and -infinity) are catched and handled. 
    /// This might be slower than the unchecked version `_beta` but does not require input sanitation to get expected results for these cases.
    ///</remarks>
    /// <param name="z">The function input for approximating B(z, w)</param>
    /// <param name="w">The function input for approximating B(z, w)</param>
    let beta z w = exp (betaLn z w)

    // incomplete beta function 
    /// <summary>
    /// Returns the regularized lower incomplete beta function 
    /// </summary>
    /// <param name="a">The first Beta parameter, a positive real number.</param>
    /// <param name="b">The second Beta parameter, a positive real number.</param>
    /// <param name="x">The upper limit of the integral.</param>
    let lowerIncompleteRegularized a b x =       
        if (a < 0.0) then invalidArg "a" "Argument must not be negative"
        if (b < 0.0) then invalidArg "b" "Argument must not be negative"
        if (x < 0.0 || x > 1.0) then invalidArg "x" "Argument XY interval is inclusive"
        let bt = 
            if (x = 0.0 || x = 1.0) then
                0.0
            else
                exp (Gamma._gammaLn (a + b) - Gamma._gammaLn a - Gamma._gammaLn b + (a*Math.Log(x)) + (b*Math.Log(1.0 - x)))

        let isSymmetryTransformation = ( x >= (a + 1.0)/(a + b + 2.0))

        let symmetryTransformation a b x =
            let qab = a + b
            let qap = a + 1.0
            let qam = a - 1.0
            let c = 1.0
            let d = 
                let tmp =  1.0 - (qab * x / qap)
                if (abs tmp < FPMIN) then 1. / FPMIN else 1. / tmp
            let h = d
            let rec loop m mm d h c =                
                let aa = float m * (b - float m)*x/((qam + mm)*(a + mm))
                let d' = 
                    let tmp = 1.0 + (aa*d)
                    if (abs tmp < FPMIN) then 1. / FPMIN else 1. / tmp
                let c' = 
                    let tmp = 1.0 + (aa/c)
                    if (abs tmp < FPMIN) then FPMIN else tmp
                let h' = h * d' * c'
                let aa' = -(a + float m)*(qab + float m)*x/((a + mm)*(qap + mm))
                let d'' = 
                    let tmp = 1.0 + (aa' * d')
                    if (abs tmp < FPMIN) then 1. / FPMIN else 1. / tmp
                let c'' = 
                    let tmp = 1.0 + (aa'/c')
                    if (abs tmp < FPMIN) then FPMIN else tmp
                
                let del = d''*c''
                let h'' = h' * del
                
                if abs (del - 1.0) <= EPS then
                     if isSymmetryTransformation then 1.0 - (bt*h''/a) else bt*h''/a
                else
                    if m < 140 then
                        loop (m+1) (mm+2.) d'' h'' c''
                    else 
                            if isSymmetryTransformation then 1.0 - (bt*h''/a) else bt*h''/a
                
            loop 1 2. d h c             

        if isSymmetryTransformation then
            symmetryTransformation b a (1.0-x)
        else
            symmetryTransformation a b x


    /// <summary>
    /// Returns the lower incomplete (unregularized) beta function
    /// </summary>
    /// <param name="a">The first Beta parameter, a positive real number.</param>
    /// <param name="b">The second Beta parameter, a positive real number.</param>
    /// <param name="x">The upper limit of the integral.</param>
    let lowerIncomplete(a, b, x) =
        (lowerIncompleteRegularized a b x) * (beta a b)
 
 
    /// <summary>
    ///   Power series for incomplete beta integral. Use when b*x
    ///   is small and x not too close to 1.
    /// </summary>
    let powerSeries a b x =
        let ai = 1.0 / a
        let ui = (1.0 - b) * x 
        let t1 = ui / (a + 1.0)
        let z  = Ops.epsilon * ai

        let rec loop u t v s n =
            if (abs v > z) then 
                let u' = (n - b) * x / n
                let t' = t * u'
                let v' = t' / (a + n)
                loop u' t' (t' / (a + n)) (s + v') (n + 1.)
            else
                s + t1 + ai
        
        let s = loop ui ui t1 0. 2.  
        let u = a * log x
        
        if ((a + b) < Gamma.maximum && abs u < Ops.logMax) then
            let t = Gamma.gamma (a + b) / (Gamma.gamma a * Gamma.gamma b)
            s * t * System.Math.Pow(x, a)
        else
            let t = Gamma.gammaLn (a + b) - Gamma.gammaLn a - Gamma.gammaLn b + u + log s
            if (t < Ops.logMin) then
                0.0
            else
                exp t


    //   Inverse of incomplete beta integral.
    //let incompleteInverse aa bb yy0 =
    //    aa


    //TODO: Beta into class to allow [<ParamArray>]

    //   Multinomial Beta function.
    //let multinomial ([<ParamArray>] x:float[]) =
    //    let mutable sum = 0.
    //    let mutable prd = 1.

    //    for i = 0 to x.Length-1 do  
    //        sum <- sum + x[i];
    //        prd <- prd * Gamma.gamma(x[i]);

    //    prd / Gamma.gamma sum
  