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

    let kolmogorovSmirnovStatistic (xArray: float []) (cdfF:float [] -> float []) (fittedF:float [] -> float []) =
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

    let binarySearchXmin (arr : 'a []) (xArray:float []) (yArray:float []) = 
        let xMax = arr |> Array.Last 
        let cdf xmin = CDFOfXMinAndData xmin xArray
        let getKSStatistic cdfFittedData = 
            Array.map2 (fun s p -> s-p|>abs) yArray cdfFittedData
            |> Array.max
        if Array.isEmpty arr then failwith "Array cannot be empty." 
        else
            let rec loop lower upper (lastV) count= 
                if lower > upper then ~~~ lower 
                elif count=0 then 
                    loop (0) upper (getKSStatistic (cdf(arr.[0]))) 1
                else 
                    let middle = lower + ((upper - lower) / 2)
                    let resultATM =  getKSStatistic (cdf(arr.[middle]))

                    let isLower = lastV - resultATM
                    if comparisonResult = 0. then
                        middle
                    elif comparisonResult < 0. then
                        loop lower (middle - 1) comparisonResult count+1
                    else
                        loop (middle + 1) upper comparisonResult count+1
            loop 0 (arr.Length - 1) (0.0) 0
        
    let kolmogorovSmirnovProb z =
        //Nearest integer function
        let Nint (v : float) =    
            int <| Math.Round(v, MidpointRounding.AwayFromZero)

        // This function returns the confidence level for the null hypothesis, where:
        //   z = dn*sqrt(n), and
        //   dn  is the maximum deviation between a hypothetical distribution
        //       function and an experimental distribution with
        //   n    events
        let kolmogorovProb z =
            let fj = [|-2.;-8.;-18.;-32.;|]
            let r  = [|0.;0.;0.;0.;|]
            let w  = 2.50662827
            // c1 - -pi**2/8, c2 = 9*c1, c3 = 25*c1
            let c1 = -1.2337005501361697
            let c2 = -11.103304951225528
            let c3 = -30.842513753404244

            match abs z with
            | u when (u < 0.2)    -> 1.
            | u when (u < 0.755)  -> let v = 1./(u*u)
                                    1. - w*(exp (c1*v) + exp(c2*v) + exp(c3*v))/u
            | u when (u < 6.8116) -> let v = u*u
                                    //Int_t maxj = TMath::Max(1,TMath::Nint(3./u));
                                    let maxj = max 1 (Nint (3./ u)) - 1
                                    for j=0 to maxj do
                                        r.[j] <- exp(fj.[j]*v)
            
                                    2.*(r.[0] - r.[1] + r.[2] - r.[3]);
            | _ -> 0. 
        kolmogorovProb z

    let inline calcZ dn n = dn * sqrt(n)  
    let kolmogorovSmirnovProbCalc dn n =
        calcZ dn n
        |> kolmogorovSmirnovProb

    static member Init xmin alpha =

        { new ContinuousDistribution<float,float> with
            member d.Sample                             =   Powerlaw.Sample xmin alpha
            member d.CDF x                              =   Powerlaw.CDF xmin alpha x         
            member d.MaximumLikelihoodEstimator xArray  =   Powerlaw.MaximumLikelihoodEstimator xmin alpha xArray
            //member d.PDF x          =   Beta.PDF alpha beta x

        }


