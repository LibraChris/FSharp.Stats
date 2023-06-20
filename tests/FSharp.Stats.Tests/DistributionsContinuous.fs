﻿module DistributionsContinuousTests

open Expecto
open System
open FSharp.Stats
open FSharp.Stats.Distributions
open FSharp.Stats.Distributions.Continuous

// Defining an accuracy appropriate for testing random sampling and inference
let fittingAccuracy : Accuracy = {absolute= 0.1 ;relative= 0.1}


[<Tests>]
let GammaDistributionTests =

    testList "Distributions.Continuous.Gamma" [
        let alpha = 0.4 
        let beta  = 4.2
    
        let d     = Gamma.Init alpha beta

        let mean  = d.Mean     
        let var   = d.Variance 
        let cdfs  = [| 0.; 0.251017; 0.328997; 0.38435; 0.428371; 0.465289;
                       0.497226; 0.525426; 0.55069; 0.573571 |] 

        let pdfs = [| 0.987113653; 0.635929273; 0.486870787; 0.400046182; 0.341683319;
                      0.299071263; 0.266235685; 0.239955525; 0.218322701; 0.200126249;
                      0.184555971; 0.171046668; 0.159190450; 0.148684554; 0.139298865;
                      0.130854902; 0.123211796; 0.116256647; 0.109897748; 0.104059710;
                      0.098679897; 0.093705765; 0.089092854; 0.084803247; 0.080804376;
                      0.077068078; 0.073569861; 0.070288299; 0.067204554; 0.064301989;
                      0.061565838; 0.058982949; 0.056541557; 0.054231102; 0.052042076;
                      0.049965886; 0.047994748; 0.046121587; 0.044339960; 0.042643979;
                      0.041028256; 0.039487846; 0.038018205; 0.036615142; 0.035274793;
                      0.033993583; 0.032768200; 0.031595571; 0.030472842; 0.029397355;
                      0.028366635; 0.027378369; 0.026430398; 0.025520703; 0.024647389;
                      0.023808683; 0.023002918; 0.022228528; 0.021484040; 0.020768066;
                      0.020079300; 0.019416507; 0.018778524; 0.018164249; 0.017572643;
                      0.017002719; 0.016453546; 0.015924240; 0.015413961; 0.014921914;
                      0.014447344; 0.013989532; 0.013547795; 0.013121484; 0.012709981;
                      0.012312696; 0.011929068; 0.011558563; 0.011200670; 0.010854903;
                      0.010520795; 0.010197904; 0.009885805; 0.009584092; 0.009292377;
                      0.009010290; 0.008737475; 0.008473592; 0.008218316; 0.007971333;
                      0.007732346; 0.007501068; 0.007277223; 0.007060548; 0.006850789;
                      0.006647704; 0.006451059; 0.006260630; 0.006076203; 0.005897569; |]

        
        //testCase "Mean" <| fun () ->
        //    Expect.floatClose Accuracy.high mean 0.21105527638190955 "Mean should be equal"

        //testCase "Variance" <| fun () ->
        //    Expect.floatClose Accuracy.high var 0.055689279830523512 "Variance should be equal"
                
        testCase "Cdfs" <| fun () ->
            cdfs 
            |> Array.iteri (fun i v ->
                let cdf = d.CDF (float i / 10.0)
                Expect.floatClose Accuracy.low cdf cdfs[i] "Cdf should be equal"
                )
                 
        testCase "Pdfs" <| fun () ->
            cdfs 
            |> Array.iteri (fun i v ->
                let pdf = d.PDF ((float i + 1.) / 10.0)
                Expect.floatClose Accuracy.low pdf pdfs[i] "Cdf should be equal"
                )          
           
        //testCase "Pdf" <| fun () ->
        //    Expect.floatClose Accuracy.high pdf 0.987114 "Pdf should be equal"
        
        testCase "FitTest" <| fun () ->
            let observations = Array.init 999999 (fun _ -> float (Continuous.Gamma.Sample alpha beta))
            let alpha',beta' = Continuous.Gamma.Fit observations
            
            Expect.floatClose fittingAccuracy alpha alpha' 
                "alpha" 
            Expect.floatClose fittingAccuracy beta beta' 
                "beta"
    
        testCase "FitTest_from_observations" <| fun () ->
            let observations = [| 1275.56; 1239.44; 1237.92; 1237.22; 1237.1; 1238.41; 1238.62; 1237.05;
                1237.19; 1236.51; 1264.6; 1238.19; 1237.39; 1235.79; 1236.53; 1236.8; 1238.06; 
                1236.5; 1235.32; 1236.44; 1236.58; 1236.3; 1237.91; 1238.6; 1238.49; 1239.21; 
                1238.57; 1244.63; 1236.06; 1236.4; 1237.88; 1237.56; 1236.66; 1236.59; 1236.53; 
                1236.32; 1238.29; 1237.79; 1237.86; 1236.42; 1236.23; 1236.37; 1237.18; 1237.63; 
                1245.8; 1238.04; 1238.55; 1238.39; 1236.75; 1237.07; 1250.78; 1238.6; 1238.36; 
                1236.58; 1236.82; 1238.4; 1257.68; 1237.78; 1236.52; 1234.9; 1237.9; 1238.58; 
                1238.12; 1237.89; 1236.54; 1236.55; 1238.37; 1237.29; 1237.64; 1236.8; 1237.73; 
                1236.71; 1238.23; 1237.84; 1236.26; 1237.58; 1238.31; 1238.4; 1237.08; 1236.61; 
                1235.92; 1236.41; 1237.89; 1237.98; 1246.75; 1237.92; 1237.1; 1237.97; 1238.69; 
                1237.05; 1236.96; 1239.44; 1238.49; 1237.88; 1236.01; 1236.57; 1236.44; 1235.76; 
                1237.62; 1238; 1263.14; 1237.66; 1237; 1236; 1261.96; 1238.58; 1237.77; 1237.06; 
                1236.31; 1238.63; 1237.23; 1236.85; 1236.23; 1236.46; 1236.9; 1237.85; 1238; 
                1237.02; 1236.19; 1236.05; 1235.73; 1258.3; 1235.98; 1237.76; 1246.93; 1239.1; 
                1237.72; 1237.67; 1236.79; 1237.61; 1238.41; 1238.29; 1238.11; 1237; 1236.52; 
                1236.6; 1236.31; 1237.77; 1238.58; 1237.88; 1247.35; 1236.14; 1236.83; 1236.15; 
                1237.93; 1238.16; 1237.34; 1236.78; 1238.66; 1237.76; 1237.19; 1236.7; 1236.04; 
                1236.66; 1237.86; 1238.54; 1238.05; 1238.41; 1236.94; 1240.95; 1261.01; 1237.72; 
                1237.91; 1238.2; 1235.68; 1236.89; 1235.12; 1271.31; 1236.97; 1270.76; 1238.52; 
                1238.19; 1238.6; 1237.16; 1236.72; 1236.71; 1237.14; 1238.48; 1237.95; 1237.42; 
                1235.86; 1236.39; 1236.13; 1236.58; 1237.95; 1237.76; 1237.39; 1238.16; 1236.31; 
                1236.41; 1236.12; 1238.7; 1236.48; 1237.84; 1236.38; 1237.95; 1238.48; 1236.51; 
                1236.56 |]
            let alpha, beta = Continuous.Gamma.Fit observations
            //let mean = 1238.8734170854279
            let alpha' = 41566.439533445438
            let beta'  = 0.029804655654680219
            
            Expect.floatClose fittingAccuracy alpha alpha'
                "Gamma Distribution Fit" 
            Expect.floatClose fittingAccuracy beta beta'
                "Gamma Distribution Fit" 
   
    ]


[<Tests>]
let BetaDistributionTests =
    testList "Distributions.Continuous.Beta" [
        
        let pdf_expect1 = 0.550369534108
        let pdf_expect2 = 5.58793544769e-08
        let pdf_expect3 = 30.
        let pdf_expect4 = 0.
        let pdf_expect5 = 0.
        let pdf_expect6 = 600
        let pdf_expect7 = 0
        let pdf_expect8 = 0
        let pdf_expect9 = 2.76522710171e-199
        let pdf_expect10 = 0.000725971756359

        let cdf_expect1 = 0.011909896429
        let cdf_expect2 = 0.999999999069
        let cdf_expect3 = 0.
        let cdf_expect4 = 1.
        let cdf_expect5 = 0.
        let cdf_expect6 = 1.
        let cdf_expect7 = 0.
        let cdf_expect8 = 1.
        let cdf_expect9 = 0.544007501411
        let cdf_expect10 = 1.
            
        // tested against R dbeta
        testCase "PDF" <| fun () ->

            let pdf_actual1 = (Beta.Init 50. 30.).PDF 0.5
            let pdf_actual2 = (Beta.Init 1. 30.).PDF 0.5
            let pdf_actual3 = (Beta.Init 1. 30.).PDF 0.
            let pdf_actual4 = (Beta.Init 1. 3.).PDF 1.
            let pdf_actual5 = (Beta.Init 600. 1.).PDF 0.
            let pdf_actual6 = (Beta.Init 600. 1.).PDF 1.
            let pdf_actual7 = (Beta.Init 600. 800.).PDF 0.
            let pdf_actual8 = (Beta.Init 600. 800.).PDF 1.
            let pdf_actual9 = (Beta.Init 600. 800.).PDF 0.11
            let pdf_actual10 = (Beta.Init 600. 800.).PDF 0.49
            Expect.floatClose Accuracy.high pdf_actual1 pdf_expect1 "Beta PDF was not determined correctly."
            Expect.floatClose Accuracy.high pdf_actual2 pdf_expect2 "Beta PDF was not determined correctly."
            Expect.floatClose Accuracy.high pdf_actual3 pdf_expect3 "Beta PDF was not determined correctly."
            Expect.floatClose Accuracy.high pdf_actual4 pdf_expect4 "Beta PDF was not determined correctly."
            Expect.floatClose Accuracy.high pdf_actual5 pdf_expect5 "Beta PDF was not determined correctly."
            Expect.floatClose Accuracy.high pdf_actual6 pdf_expect6 "Beta PDF was not determined correctly."
            Expect.floatClose Accuracy.high pdf_actual7 pdf_expect7 "Beta PDF was not determined correctly."
            Expect.floatClose Accuracy.high pdf_actual8 pdf_expect8 "Beta PDF was not determined correctly."
            Expect.floatClose Accuracy.high pdf_actual9 pdf_expect9 "Beta PDF was not determined correctly."
            Expect.floatClose Accuracy.high pdf_actual10 pdf_expect10 "Beta PDF was not determined correctly."
        
        // tested against R dbeta
        testCase "PDFLn" <| fun () ->
            
            let pdf_actual1 = (Beta.PDFLn 50. 30. 0.5)  |> exp
            let pdf_actual2 = (Beta.PDFLn 1. 30. 0.5)   |> exp
            let pdf_actual3 = (Beta.PDFLn 1. 30. 0.)    |> exp
            let pdf_actual4 = (Beta.PDFLn 1. 3. 1.)     |> exp
            //higher alpha and beta values are called already when PDF is used
            Expect.floatClose Accuracy.high pdf_actual1 pdf_expect1 "Beta PDFLn was not determined correctly."
            Expect.floatClose Accuracy.high pdf_actual2 pdf_expect2 "Beta PDFLn was not determined correctly."
            Expect.floatClose Accuracy.high pdf_actual3 pdf_expect3 "Beta PDFLn was not determined correctly."
            Expect.floatClose Accuracy.high pdf_actual4 pdf_expect4 "Beta PDFLn was not determined correctly."
            
        // tested against R pbeta
        testCase "CDF" <| fun () ->
            let cdf_actual1 = (Beta.Init 50. 30.).CDF 0.5
            let cdf_actual2 = (Beta.Init 1. 30.).CDF 0.5
            let cdf_actual3 = (Beta.Init 1. 30.).CDF 0.
            let cdf_actual4 = (Beta.Init 1. 3.).CDF 1.
            let cdf_actual5 = (Beta.Init 600. 1.).CDF 0.
            let cdf_actual6 = (Beta.Init 600. 1.).CDF 1.
            let cdf_actual7 = (Beta.Init 600. 800.).CDF 0.
            let cdf_actual8 = (Beta.Init 600. 800.).CDF 1.
            let cdf_actual9 = (Beta.Init 600. 800.).CDF 0.43
            let cdf_actual10 = (Beta.Init 600. 800.).CDF 1.49

            Expect.floatClose Accuracy.high cdf_actual1 cdf_expect1 "Beta CDF was not determined correctly."
            Expect.floatClose Accuracy.high cdf_actual2 cdf_expect2 "Beta CDF was not determined correctly."
            Expect.floatClose Accuracy.high cdf_actual3 cdf_expect3 "Beta CDF was not determined correctly."
            Expect.floatClose Accuracy.high cdf_actual4 cdf_expect4 "Beta CDF was not determined correctly."
            Expect.floatClose Accuracy.high cdf_actual5 cdf_expect5 "Beta CDF was not determined correctly."
            Expect.floatClose Accuracy.high cdf_actual6 cdf_expect6 "Beta CDF was not determined correctly."
            Expect.floatClose Accuracy.high cdf_actual7 cdf_expect7 "Beta CDF was not determined correctly."
            Expect.floatClose Accuracy.high cdf_actual8 cdf_expect8 "Beta CDF was not determined correctly."
            Expect.floatClose Accuracy.high cdf_actual9 cdf_expect9 "Beta CDF was not determined correctly."
            Expect.floatClose Accuracy.high cdf_actual10 cdf_expect10 "Beta CDF was not determined correctly."
    ]
