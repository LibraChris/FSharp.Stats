﻿namespace FSharp.Stats

module ConfidenceInterval =

    open FSharp.Stats.Distributions

    /// The deviation from the sample mean is calculated by using a t distribution. Common ciLevel = 0.95
    let ciDeviation ciLevel (sample : seq<float>)=
        let n = float (Seq.length sample)
        let stDev = Seq.stDev sample
        let t = Testing.TTest.getCriticalTValue (float n - 1.) (1. - ciLevel) Testing.TTest.TwoTailed
        let delta = t * stDev / sqrt (float n)
        delta

    /// The confidence interval is calculated by using a t distribution. Common ciLevel = 0.95
    let ci ciLevel (sample : seq<float>)=
        let mean  = Seq.mean  sample
        let delta = ciDeviation ciLevel sample
        Intervals.create (mean - delta) (mean + delta)
