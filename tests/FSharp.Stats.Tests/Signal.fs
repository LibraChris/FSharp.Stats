module SignalTests 


open Expecto
open FSharp.Stats
open FSharp.Stats.Signal
open Signal.Outliers
open TestExtensions

[<Tests>]
let outlierTests =
    let ls = [-1.4; -1.4; -1.3; -7.9; 9.4; -1.5; 5.0; 7.0; 1.1; 1.6]
    let m = mean ls //1.06
    
    let dataRow =
        [
            [20.;  11.];
            [22.;  29.];
            [12.;  27.];
            [13.;  15.];
            [19.;  23.];
            [28.;  18.];
            [16.;  30.];
            [25.;  24.];
            [14.;  21.];
            [17.;  26.]
        ]
        |> matrix

    let dataColumn = 
        [
            [20.;22.;12.;13.;19.;28.;16.;25.;14.;17.];
            [11.;29.;27.;15.;23.;18.;30.;24.;21.;26.]
        ]
        |> matrix


    let compareIntervals a b (str:string) =
        Expect.floatClose Accuracy.high (Intervals.getStart a) (Intervals.getStart b) str
        Expect.floatClose Accuracy.high (Intervals.getEnd a) (Intervals.getEnd b) str
    
    testList "OutlierTests" [
        testList "Z-Score" [
            testCase "Z-Score in a population" <| fun() ->
                let s = stDevPopulation(ls) //4.745144887
                Expect.floatClose Accuracy.high (zScore -1.4 m s) -0.5184246337 "Z-Score in a population was calculated incorrectly"

            testCase "Z-Score in a sample" <| fun()->
                let sSample = stDev(ls)
                Expect.floatClose Accuracy.high (zScore -1.4 m sSample) -0.4918207913 "Z-Score in a sample was calculated incorrectly"
                
            testCase "Z-Scores of a population" <| fun()->
                let zLs = [-0.5184246337; -0.5184246337; -0.4973504616; -1.88824582; 1.757585953; -0.5394988058; 0.8303223808; 1.251805823; 0.008429668841; 0.1138005294]
                TestExtensions.sequenceEqual Accuracy.high (zScoresOfPopulation ls) zLs "Z-Score of a population was calculated incorrectly"
                
            testCase "Z-Scores of a sample" <| fun()->
                let zLsSample = [-0.4918207913; -0.4918207913; -0.4718280762; -1.791347272; 1.667392439; -0.5118135064; 0.7877129747; 1.187567277; 0.007997086037; 0.1079606615]
                TestExtensions.sequenceEqual Accuracy.high (zScoresOfSample ls) zLsSample "Z-Score of a sample was calculated incorrectly"
                
            testCase "Population interval by Z-Score" <| fun()->
                let populationInterval = Intervals.create -0.3635434661 3.432572444
                compareIntervals (populationIntervalByZScore -0.3 0.5 ls) populationInterval "Z-Score interval in a population was calculated incorrectly"

            testCase "Sample interval by Z-Score" <| fun()->
                let sampleInterval = Intervals.create -0.4405465671 3.560910945
                compareIntervals (sampleIntervalByZscore -0.3 0.5 ls) sampleInterval "Z-Score interval in a sample was calculated incorrectly"
        ]

        testList "Mahalanobi's Distance" [
            testCase "Mahalanobi's Distance for an observation in a matrix"<| fun() ->
                let obs = Vector.ofList [20.; 11.] 
                Expect.floatClose Accuracy.high (mahalanobisDistanceOfEntry dataRow    Matrix.Sample     Matrix.RowWise obs) 1.843936618 "Mahalanobi's Distance for an observation(Sample, RowWise) calculated incorrectly"
                Expect.floatClose Accuracy.high (mahalanobisDistanceOfEntry dataColumn Matrix.Sample     Matrix.ColWise obs) 1.843936618 "Mahalanobi's Distance for an observation calculated(Sample, ColWise) incorrectly"
                Expect.floatClose Accuracy.high (mahalanobisDistanceOfEntry dataRow    Matrix.Population Matrix.RowWise obs) 1.943679857 "Mahalanobi's Distance for an observation calculated(Population, RowWise) incorrectly"
                Expect.floatClose Accuracy.high (mahalanobisDistanceOfEntry dataColumn Matrix.Population Matrix.ColWise obs) 1.943679857 "Mahalanobi's Distance for an observation calculated(Population, ColWise) incorrectly"

            testCase "Mahalanobi's Distance for every observation in a matrix"<| fun() ->
                let mahalDistancesSample = [1.843936618; 1.315823162; 1.395764847; 1.698572419; 0.1305760401; 1.862248734; 1.280527036; 1.28097611; 0.934074348; 0.6301069471]
                let mahalDistancesPopulation = [1.943679857; 1.386999396; 1.471265332; 1.790452538; 0.1376392315; 1.962982523; 1.349794013; 1.350267379; 0.9846008145; 0.6641910408]
                TestExtensions.sequenceEqual Accuracy.high (mahalanobisDistances Matrix.Sample       Matrix.RowWise dataRow   ) mahalDistancesSample  "Mahalanobi's Distance for every observation in a matrix(Sample, RowWise) was calculated incorrectly"
                TestExtensions.sequenceEqual Accuracy.high (mahalanobisDistances Matrix.Population   Matrix.RowWise dataRow   ) mahalDistancesPopulation  "Mahalanobi's Distance for every observation in a matrix(Population, RowWise) was calculated incorrectly"
                TestExtensions.sequenceEqual Accuracy.high (mahalanobisDistances Matrix.Sample       Matrix.ColWise dataColumn) mahalDistancesSample "Mahalanobi's Distance for every observation in a matrix(Sample, ColWise) was calculated incorrectly"
                TestExtensions.sequenceEqual Accuracy.high (mahalanobisDistances Matrix.Population   Matrix.ColWise dataColumn) mahalDistancesPopulation  "Mahalanobi's Distance for every observation in a matrix(Population, ColWise) was calculated incorrectly"

        ]
    ]


[<Tests>]
let normalizationTests =
    
    let table = 
        [
            [4.5;4.2;3.6]
            [4.3;4.2;0.5]
            [2.5;4.1;0.6]
        ]
        |> matrix

    let tableWithNan = 
        [
            [4.5;nan;3.6]
            [4.3;4.2;nan]
            [2.5;4.1;0.6]
        ]
        |> matrix

    testList "NormalizationTests" [
        testCase "MedianOfRatios" <| fun() ->

            let expectedNormalizedTable = 
                [
                    [3.29784;2.08239;10.99283]
                    [3.15127;2.08239;1.52678]
                    [1.83213;2.03281;1.83213]
                    
                ]
                |> matrix

            let result = Normalization.medianOfRatios table

            TestExtensions.sequenceEqual 4 result expectedNormalizedTable "Value was not normalized correctly"

        testCase "MedianOfRatiosIgnoreNans" <| fun() ->
           
            let result = Normalization.medianOfRatiosBy (fun x -> if System.Double.IsNaN x then 0.1 else x) tableWithNan

            Expect.hasCountOf  result 2u System.Double.IsNaN "Only initial nan values should be nans afterwards"

        testCase "MedianOfRatioWides" <| fun() ->
        
            let result = Normalization.medianOfRatiosWide table
            let expected = 
                table
                |> Matrix.transpose
                |> Normalization.medianOfRatios
                |> Matrix.transpose
            TestExtensions.sequenceEqual 4 result expected "Wide method should return the same result as the non wide method on a transposed matrix"
    ]

