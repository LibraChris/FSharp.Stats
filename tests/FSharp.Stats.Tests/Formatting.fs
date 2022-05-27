﻿module FormattingTests

open Expecto

open FSharp.Stats
open FSharp.Stats.Algebra
open TestExtensions

[<Tests>]
let formatValueTests =

    testList "Formatting.formatValue" [
        testCase "Format small positive float value" (fun _ -> Expect.equal (Formatting.formatValue 10.5342) "10.534" "Incorrect format for this value")
        testCase "Format large positive float value" (fun _ -> Expect.equal (Formatting.formatValue 10000000.234) "1e+07" "Incorrect format for this value")
        testCase "Format small negative float value" (fun _ -> Expect.equal (Formatting.formatValue -122.2424456) "-122.24" "Incorrect format for this value")
        testCase "Format large negative float value" (fun _ -> Expect.equal (Formatting.formatValue -1000000.345) "-1e+06" "Incorrect format for this value")
        testCase "Format small positive int value" (fun _ -> Expect.equal (Formatting.formatValue 10) "10" "Incorrect format for this value")
        testCase "Format large positive int value" (fun _ -> Expect.equal (Formatting.formatValue 10000000) "10000000" "Incorrect format for this value")
        testCase "Format small negative int value" (fun _ -> Expect.equal (Formatting.formatValue -122) "-122" "Incorrect format for this value")
        testCase "Format large negative int value" (fun _ -> Expect.equal (Formatting.formatValue -1000000) "-1000000" "Incorrect format for this value")
    ]

[<Tests>]
let formatTableTests =    
    
    let m1 = 
        [|
            [|0.1;-1003547672323.2|]
            [|-0.1;1003547672323.2|]
        |] 
        |> JaggedArray.map Formatting.formatValue
        |> array2D
        |> Formatting.formatTable

    let expected =""" 0.100 -1e+12
 -0.10  1e+12
"""

    testList "Formatting.formatTable" [
        testCase "string values formatted as table" (fun _ -> Expect.equal m1 expected "Incorrect format for this value")
    ]


[<Tests>]
let matrixFormattingtests =    
    testList "Formatting.MatrixFormatting" [

        let rnd = new System.Random(69)

        let mDense1 = Matrix.init 10 10 (fun i j -> float i * float j * rnd.NextDouble())
        let mDense2 = Matrix.init 10 100 (fun i j -> float i * float j * rnd.NextDouble())
        let mDense3 = Matrix.init 100 10 (fun i j -> float i * float j * rnd.NextDouble())
        let mDense4 = Matrix.init 100 100 (fun i j -> float i * float j * rnd.NextDouble())

        let mSparse1 = Matrix.initSparse 10 10 [ 1,1,13.37; 2,2,6942013.37 ]

        testCase "dense float matrix full display no info" (fun _ -> 
            let expected = """          0     1      2      3      4      5      6      7      8      9
                                                                         
 0 -> 0.000 0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000
 1 -> 0.000 0.337  0.134  1.101  1.121  0.612  5.203  5.427  4.360  0.286
 2 -> 0.000 1.553  2.613  0.825  3.699  9.059  3.094  6.168  9.688 11.717
 3 -> 0.000 1.885  0.262  6.466  8.832 11.299  3.357  8.493  3.829  3.668
 4 -> 0.000 3.543  4.777  0.981  4.582 11.178 19.498 15.068 28.406 12.276
 5 -> 0.000 4.085  7.661  8.177  0.215 24.933 12.299 30.125 32.935 20.369
 6 -> 0.000 2.288  5.435 15.535  1.371 14.937 10.544 14.631 30.293 28.454
 7 -> 0.000 3.239  8.621  8.593 15.556  3.471  7.527 28.003 47.662 56.915
 8 -> 0.000 5.757  8.846  6.348  7.183 33.733 41.390 27.922 45.290 53.188
 9 -> 0.000 2.932 10.459 19.086 19.290  4.265 23.378 46.242  1.773 63.866
"""
            Expect.equal (mDense1.Format(false)) expected "Incorrect format for this value"
        )        
        
        testCase "dense float matrix full display with info" (fun _ -> 
            let expected = """          0     1      2      3      4      5      6      7      8      9
                                                                         
 0 -> 0.000 0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000
 1 -> 0.000 0.337  0.134  1.101  1.121  0.612  5.203  5.427  4.360  0.286
 2 -> 0.000 1.553  2.613  0.825  3.699  9.059  3.094  6.168  9.688 11.717
 3 -> 0.000 1.885  0.262  6.466  8.832 11.299  3.357  8.493  3.829  3.668
 4 -> 0.000 3.543  4.777  0.981  4.582 11.178 19.498 15.068 28.406 12.276
 5 -> 0.000 4.085  7.661  8.177  0.215 24.933 12.299 30.125 32.935 20.369
 6 -> 0.000 2.288  5.435 15.535  1.371 14.937 10.544 14.631 30.293 28.454
 7 -> 0.000 3.239  8.621  8.593 15.556  3.471  7.527 28.003 47.662 56.915
 8 -> 0.000 5.757  8.846  6.348  7.183 33.733 41.390 27.922 45.290 53.188
 9 -> 0.000 2.932 10.459 19.086 19.290  4.265 23.378 46.242  1.773 63.866

Matrix of 10 rows x 10 columns"""

            Expect.equal (mDense1.Format(true)) expected "Incorrect format for this value"
        )

        testCase "dense float matrix omitted cols no info" (fun _ -> 
            let expected = """          0     1      2      3      4      5      6      7      8      9     10     11     12     13     14 ...      85      86      87      88      89      90      91      92      93      94      95      96      97      98      99
                                                                                                                                                                                                                                        
 0 -> 0.000 0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000 ...   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000
 1 -> 0.000 0.369  0.757  1.854  2.351  3.150  5.493  5.023  1.912  0.203  1.865 10.766  0.515  6.479  7.019 ...  29.955  56.661  49.328  79.149  71.904  53.322  52.839  22.898  29.057  52.277  81.282  10.936  46.563  95.361   9.654
 2 -> 0.000 0.272  1.935  5.916  2.384  8.304  1.447 12.270  6.552  2.409 11.980  4.007  2.968 15.547 22.483 ... 135.912  88.033 133.204  29.141 165.335 173.135  65.986 131.448 134.744  90.628  82.510  22.938  40.723 178.095  37.310
 3 -> 0.000 1.121  0.040  0.767  4.443 10.121  7.645  8.962  1.955 11.729  1.939 15.277  8.734 16.680 36.867 ...  77.372 240.773  50.632  39.029   0.957  41.046 161.119  29.161 206.496 203.983  79.302 235.075 144.997 181.575 163.336
 4 -> 0.000 2.820  0.223  6.922  8.921  3.892  2.546  8.523  5.391  9.224 24.389 11.851 21.129 34.277 36.440 ... 333.564 103.036 265.008  46.566   8.926 236.542 156.746 164.032 332.957  77.500  99.738  43.725 101.465  95.155 333.231
 5 -> 0.000 0.055  6.697  8.003  1.518 23.830 28.729  9.385  7.787 18.752 32.954  3.424 54.686 14.554 59.145 ...  20.397 300.204 242.666 189.720 360.157 344.718  52.648 420.230  60.955  46.569 398.739  22.461 350.424 459.232 318.616
 6 -> 0.000 3.766  0.609 17.559 13.382  3.141 35.365  0.151 44.668 19.828 52.289 47.069 38.330 25.797  1.029 ... 353.873 428.473  76.781  64.695 184.216 527.838 165.337  44.741 276.257  69.219 319.278 476.268 291.911 512.444 425.218
 7 -> 0.000 5.551 12.471  5.276 22.355 23.325 32.585 36.842 30.356 18.261 38.454 59.348 22.269 37.424 74.244 ... 512.009 265.467  92.593 236.426  50.217 434.788 238.223 480.226 285.111 476.742 234.154 182.183 113.785  31.278 338.933
 8 -> 0.000 3.845 13.130 19.085  7.366 14.882 16.434 48.249 50.706  0.798 69.946 25.523 95.253 70.482 14.877 ... 578.999  12.041 466.678 308.962 609.064 378.687  64.786  98.704 654.655 182.558  28.639 716.930 617.667 316.483 611.034
 9 -> 0.000 4.067  2.030  0.151 25.881 37.902 32.760 62.774 51.165 27.323 11.974 50.173 52.893 80.853 48.513 ... 382.937 569.681  86.735 345.749 707.106 663.277  82.880 264.259 811.848 480.188 292.505 179.339 541.266 103.111 469.922
"""
            Expect.equal (mDense2.Format(false)) expected "Incorrect format for this value"
        )        
        
        testCase "dense float matrix omitted cols with info" (fun _ -> 
            let expected = """          0     1      2      3      4      5      6      7      8      9     10     11     12     13     14 ...      85      86      87      88      89      90      91      92      93      94      95      96      97      98      99
                                                                                                                                                                                                                                        
 0 -> 0.000 0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000 ...   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000
 1 -> 0.000 0.369  0.757  1.854  2.351  3.150  5.493  5.023  1.912  0.203  1.865 10.766  0.515  6.479  7.019 ...  29.955  56.661  49.328  79.149  71.904  53.322  52.839  22.898  29.057  52.277  81.282  10.936  46.563  95.361   9.654
 2 -> 0.000 0.272  1.935  5.916  2.384  8.304  1.447 12.270  6.552  2.409 11.980  4.007  2.968 15.547 22.483 ... 135.912  88.033 133.204  29.141 165.335 173.135  65.986 131.448 134.744  90.628  82.510  22.938  40.723 178.095  37.310
 3 -> 0.000 1.121  0.040  0.767  4.443 10.121  7.645  8.962  1.955 11.729  1.939 15.277  8.734 16.680 36.867 ...  77.372 240.773  50.632  39.029   0.957  41.046 161.119  29.161 206.496 203.983  79.302 235.075 144.997 181.575 163.336
 4 -> 0.000 2.820  0.223  6.922  8.921  3.892  2.546  8.523  5.391  9.224 24.389 11.851 21.129 34.277 36.440 ... 333.564 103.036 265.008  46.566   8.926 236.542 156.746 164.032 332.957  77.500  99.738  43.725 101.465  95.155 333.231
 5 -> 0.000 0.055  6.697  8.003  1.518 23.830 28.729  9.385  7.787 18.752 32.954  3.424 54.686 14.554 59.145 ...  20.397 300.204 242.666 189.720 360.157 344.718  52.648 420.230  60.955  46.569 398.739  22.461 350.424 459.232 318.616
 6 -> 0.000 3.766  0.609 17.559 13.382  3.141 35.365  0.151 44.668 19.828 52.289 47.069 38.330 25.797  1.029 ... 353.873 428.473  76.781  64.695 184.216 527.838 165.337  44.741 276.257  69.219 319.278 476.268 291.911 512.444 425.218
 7 -> 0.000 5.551 12.471  5.276 22.355 23.325 32.585 36.842 30.356 18.261 38.454 59.348 22.269 37.424 74.244 ... 512.009 265.467  92.593 236.426  50.217 434.788 238.223 480.226 285.111 476.742 234.154 182.183 113.785  31.278 338.933
 8 -> 0.000 3.845 13.130 19.085  7.366 14.882 16.434 48.249 50.706  0.798 69.946 25.523 95.253 70.482 14.877 ... 578.999  12.041 466.678 308.962 609.064 378.687  64.786  98.704 654.655 182.558  28.639 716.930 617.667 316.483 611.034
 9 -> 0.000 4.067  2.030  0.151 25.881 37.902 32.760 62.774 51.165 27.323 11.974 50.173 52.893 80.853 48.513 ... 382.937 569.681  86.735 345.749 707.106 663.277  82.880 264.259 811.848 480.188 292.505 179.339 541.266 103.111 469.922

Matrix of 10 rows x 100 columns"""

            Expect.equal (mDense2.Format(true)) expected "Incorrect format for this value"
        )

        testCase "dense float matrix omitted rows no info" (fun _ -> 
            let expected = """           0      1       2       3       4       5       6       7       8       9
                                                                                   
  0 -> 0.000  0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000
  1 -> 0.000  0.571   1.271   1.433   2.526   3.604   4.162   0.163   0.665   4.357
  2 -> 0.000  1.998   1.050   5.135   6.291   6.347   0.640   3.722  14.379   9.552
  3 -> 0.000  0.750   1.635   7.650  10.142   5.156   2.425   6.708   3.206  13.910
  4 -> 0.000  0.001   5.580   5.742  10.931   4.947  11.174  23.563  10.828  30.482
  5 -> 0.000  0.922   5.368   9.493   1.369  22.752  19.137  29.700   4.131   0.987
  6 -> 0.000  2.103  10.045  15.888  18.467  28.485  10.763  32.996  30.356  15.542
  7 -> 0.000  2.619  12.456  11.928  16.794  16.055  12.672  38.412   9.667  33.958
  8 -> 0.000  1.694  14.840   1.247  28.229   0.172  34.235  35.823  50.018  67.330
  9 -> 0.000  2.575   3.891   3.014  29.770  24.023   8.809  51.333  51.051  59.432
 10 -> 0.000  6.800   4.182   1.517  10.424  13.719  17.677  67.819   2.644  54.878
 11 -> 0.000  9.394  14.871  18.621  37.522  26.204  27.914  60.457   0.044  75.750
 12 -> 0.000  7.911   0.098  25.102  27.933  17.999  11.375  65.355  71.140   7.121
 13 -> 0.000  6.330   1.947  17.117  30.712  31.796   0.208  79.800  82.302  71.110
 14 -> 0.000 10.422  16.986   7.215  18.300  59.633  51.258  63.886  28.993  44.166
  :      ...    ...     ...     ...     ...     ...     ...     ...     ...     ...
 85 -> 0.000 73.624 118.441  22.103 160.636 323.442 202.228 493.349 211.294 282.734
 86 -> 0.000 30.789 139.132 184.208 186.532 183.586 441.424 106.656 531.883  34.917
 87 -> 0.000 13.214 100.411 114.588  48.062 195.899 394.989 162.252 525.717 553.289
 88 -> 0.000 86.459 168.129  67.189  66.690  70.043 278.246 236.987 693.405 418.388
 89 -> 0.000  0.889 142.293 144.545  20.259 201.941 277.408 113.514 322.253 680.284
 90 -> 0.000 34.090  11.532 262.448 163.846 314.696 156.108 162.335 682.970  17.914
 91 -> 0.000 11.839  13.435 164.919  97.661 234.044 219.892 353.157 381.852 313.723
 92 -> 0.000 43.309  35.368  67.570  99.834 389.181 194.571  22.996 281.246 566.609
 93 -> 0.000 53.427 151.371 252.787  38.585 335.966 512.481 640.097 594.736 410.245
 94 -> 0.000 25.315  81.996 271.519 196.991 302.550 528.002 128.506 205.347 458.976
 95 -> 0.000 91.683 124.960  19.821 360.746  62.771 106.340 544.800 532.564 520.765
 96 -> 0.000 24.325 168.105  76.106  35.045 206.341 182.177 112.741 384.289 471.565
 97 -> 0.000 40.690 146.476  10.006 199.027 305.868  20.139 155.753 559.423 547.863
 98 -> 0.000 15.404 149.478 246.843 248.232 473.225 539.075 510.573 750.066 856.927
 99 -> 0.000  9.818  55.937  56.517  56.745 305.276 233.893 119.121 690.283  83.786
"""
            Expect.equal (mDense3.Format(false)) expected "Incorrect format for this value"
        )        
        
        testCase "dense float matrix omitted rows with info" (fun _ -> 
            let expected = """           0      1       2       3       4       5       6       7       8       9
                                                                                   
  0 -> 0.000  0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000
  1 -> 0.000  0.571   1.271   1.433   2.526   3.604   4.162   0.163   0.665   4.357
  2 -> 0.000  1.998   1.050   5.135   6.291   6.347   0.640   3.722  14.379   9.552
  3 -> 0.000  0.750   1.635   7.650  10.142   5.156   2.425   6.708   3.206  13.910
  4 -> 0.000  0.001   5.580   5.742  10.931   4.947  11.174  23.563  10.828  30.482
  5 -> 0.000  0.922   5.368   9.493   1.369  22.752  19.137  29.700   4.131   0.987
  6 -> 0.000  2.103  10.045  15.888  18.467  28.485  10.763  32.996  30.356  15.542
  7 -> 0.000  2.619  12.456  11.928  16.794  16.055  12.672  38.412   9.667  33.958
  8 -> 0.000  1.694  14.840   1.247  28.229   0.172  34.235  35.823  50.018  67.330
  9 -> 0.000  2.575   3.891   3.014  29.770  24.023   8.809  51.333  51.051  59.432
 10 -> 0.000  6.800   4.182   1.517  10.424  13.719  17.677  67.819   2.644  54.878
 11 -> 0.000  9.394  14.871  18.621  37.522  26.204  27.914  60.457   0.044  75.750
 12 -> 0.000  7.911   0.098  25.102  27.933  17.999  11.375  65.355  71.140   7.121
 13 -> 0.000  6.330   1.947  17.117  30.712  31.796   0.208  79.800  82.302  71.110
 14 -> 0.000 10.422  16.986   7.215  18.300  59.633  51.258  63.886  28.993  44.166
  :      ...    ...     ...     ...     ...     ...     ...     ...     ...     ...
 85 -> 0.000 73.624 118.441  22.103 160.636 323.442 202.228 493.349 211.294 282.734
 86 -> 0.000 30.789 139.132 184.208 186.532 183.586 441.424 106.656 531.883  34.917
 87 -> 0.000 13.214 100.411 114.588  48.062 195.899 394.989 162.252 525.717 553.289
 88 -> 0.000 86.459 168.129  67.189  66.690  70.043 278.246 236.987 693.405 418.388
 89 -> 0.000  0.889 142.293 144.545  20.259 201.941 277.408 113.514 322.253 680.284
 90 -> 0.000 34.090  11.532 262.448 163.846 314.696 156.108 162.335 682.970  17.914
 91 -> 0.000 11.839  13.435 164.919  97.661 234.044 219.892 353.157 381.852 313.723
 92 -> 0.000 43.309  35.368  67.570  99.834 389.181 194.571  22.996 281.246 566.609
 93 -> 0.000 53.427 151.371 252.787  38.585 335.966 512.481 640.097 594.736 410.245
 94 -> 0.000 25.315  81.996 271.519 196.991 302.550 528.002 128.506 205.347 458.976
 95 -> 0.000 91.683 124.960  19.821 360.746  62.771 106.340 544.800 532.564 520.765
 96 -> 0.000 24.325 168.105  76.106  35.045 206.341 182.177 112.741 384.289 471.565
 97 -> 0.000 40.690 146.476  10.006 199.027 305.868  20.139 155.753 559.423 547.863
 98 -> 0.000 15.404 149.478 246.843 248.232 473.225 539.075 510.573 750.066 856.927
 99 -> 0.000  9.818  55.937  56.517  56.745 305.276 233.893 119.121 690.283  83.786

Matrix of 100 rows x 10 columns"""

            Expect.equal (mDense3.Format(true)) expected "Incorrect format for this value"
        )

        testCase "dense float matrix omitted rows and cols no info" (fun _ -> 
            let expected = """             0       1       2       3       4       5       6       7       8       9       10       11       12       13       14 ...       85       86       87       88       89       90       91       92       93       94       95       96       97       98       99
                                                                                                                                                                                                                                                                              
  0 ->   0.000   0.529   0.873   2.800   0.905   3.353   1.289   6.381   3.010   3.252    3.711   10.828    6.136    3.723   13.354 ...    0.000    0.000    0.000    0.000    0.000    0.000    0.000    0.000    0.000    0.000    0.000    0.000    0.000    0.000    0.000
  1 ->   0.000   0.578   0.788   3.013   6.916   1.074   8.196  11.855  13.041   4.023    4.793    7.982    3.674   15.233    4.160 ...   75.801   18.898   57.698   25.467   40.929   77.650   25.796   16.615   85.206   53.250   14.148   59.625   11.502   96.185   65.284
  2 ->   0.000   1.118   4.503   8.788   7.450   3.705  13.049  19.729  13.625  11.977   26.147   31.094   24.031   29.374   17.244 ...  147.471  156.730  102.547    5.338   73.668   49.444   15.806   60.798   45.654  142.259  131.712  144.015   82.009  141.904  121.608
  3 ->   0.000   1.169   4.415   8.561   4.177  17.695   3.308  26.933  29.367   7.066   35.622   22.205   20.743   14.349    9.917 ...  244.644  126.556   43.554  163.956  225.609  173.822  159.796  265.869   76.627   89.714  131.496   16.856  105.691   92.460  209.170
  4 ->   0.000   4.524   7.888   9.607  17.457   0.306  11.439   5.504  20.414  24.204   16.898   38.762    8.731   58.314   69.632 ...   16.713  203.790  286.932  291.102  330.963  336.217  278.777   95.589  244.692  109.670  183.635  318.584   78.349  350.069  274.409
  5 ->   0.000   0.047   4.121   5.761   7.626   8.161   9.619  31.054  28.594  42.659    4.011   23.175   31.881   64.223   60.231 ...  264.545  374.850  338.361  224.019  120.705   75.757  197.089  434.310  408.493  269.776  409.849  415.516  393.740  362.849   79.977
  6 ->   0.000   5.342  12.301   4.290  19.705  23.540  26.782   0.911  13.634  11.160   29.231   40.063    4.352   71.709   44.453 ...  337.878  132.782  301.409  323.353  324.233  461.890  226.968  247.177  112.362  139.772  176.255   66.638  313.296  174.878  514.112
  7 ->   0.000   7.921   5.941  23.229  19.724   7.676  31.045  29.147  51.442  67.763   34.229   85.258   76.647   34.910    4.464 ...   72.548   86.027  334.107  550.781   14.765  330.970  246.846  484.915  127.455  320.049  168.836  285.202  239.102  389.955  547.607
  8 ->   0.000   4.012  17.044   7.450   1.595   4.542  23.906   0.825  39.877  72.870   78.892   97.406   47.274    1.792  123.518 ...  489.975  266.069  548.115   85.488  604.199  467.518  519.075  322.168  719.914  701.634  132.426  745.218  327.117  337.292  684.178
  9 ->   0.000   6.426  18.127   2.540   9.615  43.313  39.628  21.694   8.293  49.229   29.199   68.634   81.025   62.407   17.244 ...  254.456  574.030  366.758  226.845  697.452  529.468  564.104  719.287   18.804  626.013  104.894  625.169  676.855  869.577  725.237
 10 ->   0.000   1.916  15.309  31.563  29.945  34.555   2.164   6.854  74.936  40.873   33.498   83.429   25.666   88.146  146.706 ...  757.597  182.356  636.383  385.724  381.282  795.580  430.808  349.462  886.413  590.337  626.901   46.875  642.853  431.668  839.495
 11 ->   0.000   3.638   3.514   0.386  32.409  55.537  25.510   2.715  46.432 100.824  113.174   52.487   41.905    5.175   97.681 ...  641.387  280.888  127.801  860.904  839.837  187.785  813.098  905.628  421.338  141.643  205.391  703.456  639.768  591.237  229.743
 12 ->   0.000   2.077   0.562  33.381  38.564  63.239  68.992  58.000  30.998  28.500   70.465   12.305   84.680   34.242   20.384 ...  936.192  108.240   72.741  353.286  973.266  421.449 1051.022  164.455  886.421 1014.618  525.604  181.642  603.914  884.115  677.696
 13 ->   0.000   2.111  14.847  18.120  18.011  57.661  15.543  61.431  23.831  70.455  138.801   30.324   89.580  180.465    7.954 ...  532.632  430.912 1032.861 1078.717   76.172  852.085  709.718  628.723  892.311  348.149  347.938  982.502  595.019   52.348   61.909
 14 ->   0.000   0.675   6.770  29.105  53.102  52.425  16.541  32.803  24.284 126.257  106.427  145.659  162.738  178.945  147.643 ...  219.823  398.129  513.556  491.670 1096.574  207.844  161.146  142.099  510.745    1.904 1197.729   74.498  623.045 1110.252  471.717
  :        ...     ...     ...     ...     ...     ...     ...     ...     ...     ...      ...      ...      ...      ...      ... ...      ...      ...      ...      ...      ...      ...      ...      ...      ...      ...      ...      ...      ...      ...      ...
 85 -> 121.806   4.449  16.793 327.578 100.791  57.726 565.264 535.246 122.442 272.689  651.825  356.182  667.890  515.488 1225.488 ... 5416.512 3133.189 2571.566 4014.137 2372.295 3628.895 7118.603 1829.817 4973.430 4516.299 7937.209 1536.300 1698.671  944.049 3283.056
 86 -> 105.891  33.137  64.971 150.283 457.666 369.620 136.517  49.680 639.810  11.223  197.939  216.000  707.567  466.536  439.745 ... 2223.210 5983.848   86.574 6295.561 2568.867 5372.146 2951.421 4205.028 5816.305 1118.186 5890.141 6521.126 5114.603 7118.847 7099.432
 87 ->  94.343  40.123  34.012 305.471 372.831 150.012  68.217 131.583 657.489 575.565  742.622 1117.466  439.937  682.290  609.841 ... 5643.359 6841.265  293.441 3264.109 2945.036 5677.997 3732.846 4721.548 7434.687  196.697 4681.785 1093.296 1442.105 3492.534 7013.545
 88 -> 133.381  82.910 332.077 173.668 422.138 514.408 554.212 430.205 810.682 572.664  840.751  982.747  377.115  510.018  983.716 ... 2505.831 1735.069 5085.546 6239.398   15.014 3472.400 6212.248  190.067 4196.619 4808.894 4986.510 3893.319 5391.448 4941.816 2333.632
 89 -> 164.086  21.792 135.157  80.036 491.708 316.695 278.505 153.091 264.112 262.547  434.116  918.655  271.062  706.297  779.770 ... 6730.465 6188.537 6948.444 2916.550 2282.292 3612.986 4299.206 6133.302 2499.048 5055.278 8098.833 5328.177 7017.662 6406.702 4653.487
 90 ->  90.706 257.853  60.762 107.656 313.884  69.582  26.630  52.496 638.013 810.969   13.808   77.701  988.353  108.548  413.416 ... 5602.842  183.958 1005.078 3279.931   25.988 1512.070 7080.854 1228.568 3129.116 3926.970 2268.851 1590.904 3722.962 7587.112 1147.959
 91 ->  25.725  55.973 136.325 272.049 183.699 190.744  45.324 400.715 461.149 237.234  227.457  956.547 1250.416 1092.977   32.565 ... 5750.600 5195.454  338.658  217.351 5709.195 5744.160 5244.375 6165.506 5148.552 5653.679 2025.197 7827.562 3711.028 8221.018 5210.631
 92 ->  18.988  45.937  36.557  16.693 551.776 487.803 608.113 827.936 105.717 965.658  620.135 1080.083  594.109  122.437  558.322 ... 7279.823 3033.580 4129.001 1125.713 3118.908 6007.869 7618.297 8244.613 4618.580 4940.244 7683.478 5787.929 4239.315 7515.542 7348.586
 93 ->  38.262  38.215 358.976  42.608 314.170 399.376 109.054 817.181  97.519 307.903  702.925 1045.147  402.679 1296.133 1304.952 ... 2798.868 2809.466  162.900 5956.140  789.608 7303.674 3996.364 8541.878 6066.736 7306.379 5101.265 8213.958 7751.730 8010.192 4472.806
 94 ->   9.671 227.695  54.216 458.690 558.747 388.078 508.382 657.524 770.807 224.193 1114.785 1212.846 1119.382  878.040  373.753 ... 4140.978  386.735 7922.451 3623.586 1008.386 5196.351 4615.186 4130.571 7507.452 6277.388  169.601 6402.340 8985.614 8135.705 4531.890
 95 ->  30.141 180.208  83.045 387.100 250.980 509.888 290.869 315.783 645.062 210.198  457.477  447.338  766.976 1319.663  834.371 ... 7960.489 6085.569 4954.773 3795.198 4322.706 5524.015 3729.493 1181.202 3027.043 3526.940 4784.044  190.829 6917.789 6213.412 1154.812
 96 ->  55.889  76.809 277.112  40.456 207.985 281.524 662.676 247.922 271.158 969.988  169.134 1044.289  600.528 1235.160 1440.033 ... 7797.409    6.056  646.073 5182.352 4856.829 1632.947 5694.009 3941.509 8126.620 1438.446 4944.268 6529.732 7587.947 2997.643 8931.945
 97 ->  98.432 241.058 119.992 484.366 550.355 154.677 357.822 262.762 529.444 651.399 1040.328  322.286 1009.843 1282.437 1086.260 ... 1516.961 7299.605 7050.671 4207.465 7573.702 1025.514 1412.649 5393.216 2002.291 4547.475 6371.506 1402.909 3799.099 5236.658 5699.003
 98 ->  85.982 138.055 329.857 341.346 562.946  12.474 425.633 120.145 715.950 938.459   99.861  828.154 1255.346 1314.460 1345.354 ... 4238.059 7536.498 4556.223 4890.874 3096.930  503.061 8510.020  909.925 8394.635 5107.284 4145.574 5300.112 8946.610 5619.715 8810.497
 99 ->  74.229  59.814  90.515 478.340  28.036 347.503 492.858 304.780 489.473 928.537  993.152  717.310  503.962  844.119 1444.956 ... 6259.917 2828.528 8545.715 2231.671 3463.213 8191.683  553.990 4689.324  779.949 1782.211 2840.371 4116.932 1877.515  163.961 6333.645
"""
            Expect.equal (mDense4.Format(false)) expected "Incorrect format for this value"
        )        
        
        testCase "dense float matrix omitted rows and cols with info" (fun _ -> 
            let expected = """             0       1       2       3       4       5       6       7       8       9       10       11       12       13       14 ...       85       86       87       88       89       90       91       92       93       94       95       96       97       98       99
                                                                                                                                                                                                                                                                              
  0 ->   0.000   0.529   0.873   2.800   0.905   3.353   1.289   6.381   3.010   3.252    3.711   10.828    6.136    3.723   13.354 ...    0.000    0.000    0.000    0.000    0.000    0.000    0.000    0.000    0.000    0.000    0.000    0.000    0.000    0.000    0.000
  1 ->   0.000   0.578   0.788   3.013   6.916   1.074   8.196  11.855  13.041   4.023    4.793    7.982    3.674   15.233    4.160 ...   75.801   18.898   57.698   25.467   40.929   77.650   25.796   16.615   85.206   53.250   14.148   59.625   11.502   96.185   65.284
  2 ->   0.000   1.118   4.503   8.788   7.450   3.705  13.049  19.729  13.625  11.977   26.147   31.094   24.031   29.374   17.244 ...  147.471  156.730  102.547    5.338   73.668   49.444   15.806   60.798   45.654  142.259  131.712  144.015   82.009  141.904  121.608
  3 ->   0.000   1.169   4.415   8.561   4.177  17.695   3.308  26.933  29.367   7.066   35.622   22.205   20.743   14.349    9.917 ...  244.644  126.556   43.554  163.956  225.609  173.822  159.796  265.869   76.627   89.714  131.496   16.856  105.691   92.460  209.170
  4 ->   0.000   4.524   7.888   9.607  17.457   0.306  11.439   5.504  20.414  24.204   16.898   38.762    8.731   58.314   69.632 ...   16.713  203.790  286.932  291.102  330.963  336.217  278.777   95.589  244.692  109.670  183.635  318.584   78.349  350.069  274.409
  5 ->   0.000   0.047   4.121   5.761   7.626   8.161   9.619  31.054  28.594  42.659    4.011   23.175   31.881   64.223   60.231 ...  264.545  374.850  338.361  224.019  120.705   75.757  197.089  434.310  408.493  269.776  409.849  415.516  393.740  362.849   79.977
  6 ->   0.000   5.342  12.301   4.290  19.705  23.540  26.782   0.911  13.634  11.160   29.231   40.063    4.352   71.709   44.453 ...  337.878  132.782  301.409  323.353  324.233  461.890  226.968  247.177  112.362  139.772  176.255   66.638  313.296  174.878  514.112
  7 ->   0.000   7.921   5.941  23.229  19.724   7.676  31.045  29.147  51.442  67.763   34.229   85.258   76.647   34.910    4.464 ...   72.548   86.027  334.107  550.781   14.765  330.970  246.846  484.915  127.455  320.049  168.836  285.202  239.102  389.955  547.607
  8 ->   0.000   4.012  17.044   7.450   1.595   4.542  23.906   0.825  39.877  72.870   78.892   97.406   47.274    1.792  123.518 ...  489.975  266.069  548.115   85.488  604.199  467.518  519.075  322.168  719.914  701.634  132.426  745.218  327.117  337.292  684.178
  9 ->   0.000   6.426  18.127   2.540   9.615  43.313  39.628  21.694   8.293  49.229   29.199   68.634   81.025   62.407   17.244 ...  254.456  574.030  366.758  226.845  697.452  529.468  564.104  719.287   18.804  626.013  104.894  625.169  676.855  869.577  725.237
 10 ->   0.000   1.916  15.309  31.563  29.945  34.555   2.164   6.854  74.936  40.873   33.498   83.429   25.666   88.146  146.706 ...  757.597  182.356  636.383  385.724  381.282  795.580  430.808  349.462  886.413  590.337  626.901   46.875  642.853  431.668  839.495
 11 ->   0.000   3.638   3.514   0.386  32.409  55.537  25.510   2.715  46.432 100.824  113.174   52.487   41.905    5.175   97.681 ...  641.387  280.888  127.801  860.904  839.837  187.785  813.098  905.628  421.338  141.643  205.391  703.456  639.768  591.237  229.743
 12 ->   0.000   2.077   0.562  33.381  38.564  63.239  68.992  58.000  30.998  28.500   70.465   12.305   84.680   34.242   20.384 ...  936.192  108.240   72.741  353.286  973.266  421.449 1051.022  164.455  886.421 1014.618  525.604  181.642  603.914  884.115  677.696
 13 ->   0.000   2.111  14.847  18.120  18.011  57.661  15.543  61.431  23.831  70.455  138.801   30.324   89.580  180.465    7.954 ...  532.632  430.912 1032.861 1078.717   76.172  852.085  709.718  628.723  892.311  348.149  347.938  982.502  595.019   52.348   61.909
 14 ->   0.000   0.675   6.770  29.105  53.102  52.425  16.541  32.803  24.284 126.257  106.427  145.659  162.738  178.945  147.643 ...  219.823  398.129  513.556  491.670 1096.574  207.844  161.146  142.099  510.745    1.904 1197.729   74.498  623.045 1110.252  471.717
  :        ...     ...     ...     ...     ...     ...     ...     ...     ...     ...      ...      ...      ...      ...      ... ...      ...      ...      ...      ...      ...      ...      ...      ...      ...      ...      ...      ...      ...      ...      ...
 85 -> 121.806   4.449  16.793 327.578 100.791  57.726 565.264 535.246 122.442 272.689  651.825  356.182  667.890  515.488 1225.488 ... 5416.512 3133.189 2571.566 4014.137 2372.295 3628.895 7118.603 1829.817 4973.430 4516.299 7937.209 1536.300 1698.671  944.049 3283.056
 86 -> 105.891  33.137  64.971 150.283 457.666 369.620 136.517  49.680 639.810  11.223  197.939  216.000  707.567  466.536  439.745 ... 2223.210 5983.848   86.574 6295.561 2568.867 5372.146 2951.421 4205.028 5816.305 1118.186 5890.141 6521.126 5114.603 7118.847 7099.432
 87 ->  94.343  40.123  34.012 305.471 372.831 150.012  68.217 131.583 657.489 575.565  742.622 1117.466  439.937  682.290  609.841 ... 5643.359 6841.265  293.441 3264.109 2945.036 5677.997 3732.846 4721.548 7434.687  196.697 4681.785 1093.296 1442.105 3492.534 7013.545
 88 -> 133.381  82.910 332.077 173.668 422.138 514.408 554.212 430.205 810.682 572.664  840.751  982.747  377.115  510.018  983.716 ... 2505.831 1735.069 5085.546 6239.398   15.014 3472.400 6212.248  190.067 4196.619 4808.894 4986.510 3893.319 5391.448 4941.816 2333.632
 89 -> 164.086  21.792 135.157  80.036 491.708 316.695 278.505 153.091 264.112 262.547  434.116  918.655  271.062  706.297  779.770 ... 6730.465 6188.537 6948.444 2916.550 2282.292 3612.986 4299.206 6133.302 2499.048 5055.278 8098.833 5328.177 7017.662 6406.702 4653.487
 90 ->  90.706 257.853  60.762 107.656 313.884  69.582  26.630  52.496 638.013 810.969   13.808   77.701  988.353  108.548  413.416 ... 5602.842  183.958 1005.078 3279.931   25.988 1512.070 7080.854 1228.568 3129.116 3926.970 2268.851 1590.904 3722.962 7587.112 1147.959
 91 ->  25.725  55.973 136.325 272.049 183.699 190.744  45.324 400.715 461.149 237.234  227.457  956.547 1250.416 1092.977   32.565 ... 5750.600 5195.454  338.658  217.351 5709.195 5744.160 5244.375 6165.506 5148.552 5653.679 2025.197 7827.562 3711.028 8221.018 5210.631
 92 ->  18.988  45.937  36.557  16.693 551.776 487.803 608.113 827.936 105.717 965.658  620.135 1080.083  594.109  122.437  558.322 ... 7279.823 3033.580 4129.001 1125.713 3118.908 6007.869 7618.297 8244.613 4618.580 4940.244 7683.478 5787.929 4239.315 7515.542 7348.586
 93 ->  38.262  38.215 358.976  42.608 314.170 399.376 109.054 817.181  97.519 307.903  702.925 1045.147  402.679 1296.133 1304.952 ... 2798.868 2809.466  162.900 5956.140  789.608 7303.674 3996.364 8541.878 6066.736 7306.379 5101.265 8213.958 7751.730 8010.192 4472.806
 94 ->   9.671 227.695  54.216 458.690 558.747 388.078 508.382 657.524 770.807 224.193 1114.785 1212.846 1119.382  878.040  373.753 ... 4140.978  386.735 7922.451 3623.586 1008.386 5196.351 4615.186 4130.571 7507.452 6277.388  169.601 6402.340 8985.614 8135.705 4531.890
 95 ->  30.141 180.208  83.045 387.100 250.980 509.888 290.869 315.783 645.062 210.198  457.477  447.338  766.976 1319.663  834.371 ... 7960.489 6085.569 4954.773 3795.198 4322.706 5524.015 3729.493 1181.202 3027.043 3526.940 4784.044  190.829 6917.789 6213.412 1154.812
 96 ->  55.889  76.809 277.112  40.456 207.985 281.524 662.676 247.922 271.158 969.988  169.134 1044.289  600.528 1235.160 1440.033 ... 7797.409    6.056  646.073 5182.352 4856.829 1632.947 5694.009 3941.509 8126.620 1438.446 4944.268 6529.732 7587.947 2997.643 8931.945
 97 ->  98.432 241.058 119.992 484.366 550.355 154.677 357.822 262.762 529.444 651.399 1040.328  322.286 1009.843 1282.437 1086.260 ... 1516.961 7299.605 7050.671 4207.465 7573.702 1025.514 1412.649 5393.216 2002.291 4547.475 6371.506 1402.909 3799.099 5236.658 5699.003
 98 ->  85.982 138.055 329.857 341.346 562.946  12.474 425.633 120.145 715.950 938.459   99.861  828.154 1255.346 1314.460 1345.354 ... 4238.059 7536.498 4556.223 4890.874 3096.930  503.061 8510.020  909.925 8394.635 5107.284 4145.574 5300.112 8946.610 5619.715 8810.497
 99 ->  74.229  59.814  90.515 478.340  28.036 347.503 492.858 304.780 489.473 928.537  993.152  717.310  503.962  844.119 1444.956 ... 6259.917 2828.528 8545.715 2231.671 3463.213 8191.683  553.990 4689.324  779.949 1782.211 2840.371 4116.932 1877.515  163.961 6333.645

Matrix of 100 rows x 100 columns"""

            Expect.equal (mDense4.Format(true)) expected "Incorrect format for this value"
        )

        testCase "sparse float matrix full display no info" (fun _ -> 
            let expected = """          0      1       2     3     4     5     6     7     8     9
                                                                    
 0 -> 0.000  0.000   0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
 1 -> 0.000 13.370   0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
 2 -> 0.000  0.000 6.9e+06 0.000 0.000 0.000 0.000 0.000 0.000 0.000
 3 -> 0.000  0.000   0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
 4 -> 0.000  0.000   0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
 5 -> 0.000  0.000   0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
 6 -> 0.000  0.000   0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
 7 -> 0.000  0.000   0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
 8 -> 0.000  0.000   0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
 9 -> 0.000  0.000   0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
"""
            Expect.equal (mSparse1.Format(false)) expected "Incorrect format for this value"
        )

        testCase "sparse float matrix full display with info" (fun _ -> 
            let expected = """          0      1       2     3     4     5     6     7     8     9
                                                                    
 0 -> 0.000  0.000   0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
 1 -> 0.000 13.370   0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
 2 -> 0.000  0.000 6.9e+06 0.000 0.000 0.000 0.000 0.000 0.000 0.000
 3 -> 0.000  0.000   0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
 4 -> 0.000  0.000   0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
 5 -> 0.000  0.000   0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
 6 -> 0.000  0.000   0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
 7 -> 0.000  0.000   0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
 8 -> 0.000  0.000   0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
 9 -> 0.000  0.000   0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000

Matrix of 10 rows x 10 columns"""
            Expect.equal (mSparse1.Format(true)) expected "Incorrect format for this value"
        )
    ]