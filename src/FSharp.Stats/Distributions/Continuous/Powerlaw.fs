namespace FSharp.Stats.Distributions.Continuous

open System
open FSharp.Stats
open FSharp.Stats.Distributions
open FSharp.Stats.Ops

//x0
//alpha
type Powerlaw =

    static member CheckParam xmin alpha =
        if x0=<0. || alpha=0. then
            failwith "A powerlaw distribution cannot have a negative minimal value or a "

    let kolmogorovSmirnov (xArray: float []) (cdfF:float [] -> float []) (fittedF:float [] -> float []) =
        let fittedData      = xArray        |> fittedF
        let cdfObservations = xArray        |> cdfF
        let cdfFittedData   = fittedData    |> cdfF
        Array.map2 (fun s p -> s-p|>abs) cdfObservations cdfFittedData
        |> Array.max


    static member MaximumLikelihoodEstimator (xmin:float) (xArray:float []) =
        CheckParam xmin 1.
        let sum = 
            xArray
            |> Array.sumBy(fun xi -> 
                let division = xi/xmin
                System.Math.Log(division)
            )

        let n = xArray|>Array.length|> float

        1. + (n*(sum**(-1.)))

    static member Sample xmin alpha =
        CheckParam xmin alpha
        let r = new System.Random()
        let randomN = r.NextDouble()

        xmin * (1. - randomN) ** (-1. / (alpha - 1.))


    static member CDF xmin alpha x =
        CheckParam xmin alpha
        let base = x/xmin

        let exponent = ((-1.*alpha)+1.)

        base ** exponent

    let CDFOfXMinAndData (xmin:float) (xArray:float []) =
        let alpha = MaximumLikelihoodEstimator xmin xArray
        CDF xmin alpha xArray


    let binarySearchMod (arr : 'a []) xArray = 
        if Array.isEmpty arr then failwith "Array cannot be empty." 
        else
            let rec loop lower upper (lastV) count= 
                if lower > upper then ~~~ lower 
                elif count=0 then 
                    let middel = lower + ((upper - lower) / 2)
                    loop (middle + 1) upper (CDFOfXMinAndData arr.[middle] xArray) 1
                else 
                    let middle = lower + ((upper - lower) / 2)
                    let resultATM =  (CDFOfXMinAndData arr.[middle] xArray)
                    let comparisonResult = lastV - resultATM
                    if comparisonResult = 0. then
                        middle
                    elif comparisonResult < 0. then
                        loop lower (middle - 1) comparisonResult count+1
                    else
                        loop (middle + 1) upper comparisonResult count+1
            loop 0 (arr.Length - 1) (0.0) 0

    
    static member Init xmin alpha =

        { new ContinuousDistribution<float,float> with
            member d.Sample x                           =   Powerlaw.Sample xmin alpha
            member d.CDF x                              =   Powerlaw.CDF xmin alpha x         
            member d.MaximumLikelihoodEstimator xArray  =   Powerlaw.MaximumLikelihoodEstimator xmin alpha xArray
            //member d.PDF x          =   Beta.PDF alpha beta x

        }


