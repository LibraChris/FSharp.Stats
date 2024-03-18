﻿namespace FSharp.Stats

[<AutoOpen>]
module JaggedArray =

    /// <summary>Creates an jagged array with the given dimensions and a generator function to compute the elements</summary>
    /// <remarks></remarks>
    /// <param name="rowN"></param>
    /// <param name="colN"></param>
    /// <param name="f"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>
    let init rowN colN f =
        Array.init rowN (fun rowi ->
            Array.init colN (fun coli -> f rowi coli)) 

    /// <summary>Creates an jagged array where the entries are initially the default value Unchecked.defaultof</summary>
    /// <remarks></remarks>
    /// <param name="rowN"></param>
    /// <param name="colN"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>
    let zeroCreate rowN colN = 
        Array.init rowN (fun _ -> Array.zeroCreate colN)

    /// <summary>Copies the jagged array</summary>
    /// <remarks></remarks>
    /// <param name="arr"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>
    let copy (arr : _[][]) = 
        Array.init arr.Length (fun i ->
            Array.copy arr.[i]
            )    

    /// <summary>Transposes a jagged array (unchecked)</summary>
    /// <remarks>The resulting row count is determined by first collection length!</remarks>
    /// <param name="arr"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>
    let transpose_ (arr: 'T [][]) =
        if arr.Length > 0 then 
            let colSize = arr.[0].Length
            Array.init colSize (fun rowI -> Array.init (arr.Length) (fun colI -> (arr.[colI].[rowI])))
        else
            arr

    /// <summary>Transposes a jagged array</summary>
    /// <remarks>inner collections (rows) have to be of equal length</remarks>
    /// <param name="arr"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>
    let transpose (arr: 'T [][]) =
        let colSize = 
            let sizes = arr |> Array.distinctBy Array.length
            if sizes.Length <> 1 then failwithf "All inner collections (rows) have to be of equal length!" else sizes.[0].Length
        transpose_ arr
    
    /// <summary>Converts from an Array2D into an jagged array</summary>
    /// <remarks></remarks>
    /// <param name="arr"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>
    let toArray2D (arr: 'T [][]) =
        let n = arr.Length
        if n > 0 then 
            let m = arr.[0].Length
            Array2D.init n m (fun i j -> arr.[i].[j])
        else
            Array2D.zeroCreate 0 0    

    /// <summary>Converts a jagged array into an Array2D</summary>
    /// <remarks></remarks>
    /// <param name="arr"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>
    let ofArray2D (arr:'T[,]) =
        let n,m = Array2D.length1 arr,Array2D.length2 arr
        Array.init n (fun i ->
                Array.init m (fun j -> arr.[i,j])) 

    /// <summary>Converts a jagged list into a jagged array</summary>
    /// <remarks></remarks>
    /// <param name="data"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>
    let ofJaggedList (data: 'T list list) =
        data
        |> List.map Array.ofList
        |> Array.ofList

    /// <summary>Converts a jagged array into a jagged list</summary>
    /// <remarks></remarks>
    /// <param name="arr"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>
    let toJaggedList (arr: 'T [][]) =
        arr
        |> Array.map List.ofArray
        |> List.ofArray

    /// <summary>Converts a jagged Seq into a jagged array</summary>
    /// <remarks></remarks>
    /// <param name="data"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>
    let ofJaggedSeq (data: seq<#seq<'T>>) =
        data
        |> Seq.map Array.ofSeq
        |> Array.ofSeq

    /// <summary>Converts a jagged array into a jagged seq</summary>
    /// <remarks></remarks>
    /// <param name="arr"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>
    let toJaggedSeq (arr: 'T [][]) =
        arr
        |> Seq.map Array.toSeq

    /// <summary>Builds a new jagged array whose inner arrays are the results of applying the given function to each of their elements.</summary>
    /// <remarks></remarks>
    /// <param name="mapping"></param>
    /// <param name="jArray"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>
    let map (mapping: 'T -> 'U) (jArray : 'T[][]) =
        jArray
        |> Array.map (Array.map mapping)

    /// <summary>Builds a new jagged array whose inner arrays are the results of applying the given function to the corresponding elements of the inner arrays of the two jagged arrays pairwise. <br />All corresponding inner arrays must be of the same length, otherwise ArgumentException is raised.</summary>
    /// <remarks></remarks>
    /// <param name="mapping"></param>
    /// <param name="jArray1"></param>
    /// <param name="jArray2"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>
    let map2 (mapping: 'T1 -> 'T2 -> 'U) (jArray1 : 'T1[][]) (jArray2 : 'T2[][]) = 
        jArray1
        |> Array.mapi (fun index x -> (Array.map2 mapping x jArray2.[index]))

    /// <summary>Builds a new jagged array whose inner arrays are the results of applying the given function to the corresponding elements of the inner arrays of the tree jagged arrays triplewise. <br />All corresponding inner arrays must be of the same length, otherwise ArgumentException is raised.</summary>
    /// <remarks></remarks>
    /// <param name="mapping"></param>
    /// <param name="jArray1"></param>
    /// <param name="jArray2"></param>
    /// <param name="jArray3"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>
    let map3 (mapping : 'T1 -> 'T2 -> 'T3 -> 'U) (jArray1 : 'T1[][]) (jArray2 : 'T2[][]) (jArray3 : 'T3[][]) =
        jArray1 |> Array.mapi (fun index x -> (Array.map3 mapping x jArray2.[index] jArray3.[index]))

    /// <summary>Builds a new jagged array whose inner arrays are the results of applying the given function to each of their elements. The integer index passed to the function indicates the index of element in the inner array being transformed.</summary>
    /// <remarks></remarks>
    /// <param name="mapping"></param>
    /// <param name="jArray"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>
    let mapi (mapping: int -> 'T -> 'U) (jArray : 'T[][]) =
        jArray
        |> Array.map (fun x -> x |> Array.mapi mapping)

    /// <summary>Applies a function to each element of the inner arrays of the jagged array, threading an accumulator argument through the computation.</summary>
    /// <remarks></remarks>
    /// <param name="folder"></param>
    /// <param name="state"></param>
    /// <param name="jArray"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>
    let innerFold (folder: 'State -> 'T -> 'State) (state: 'State) (jArray : 'T[][]) =
        jArray
        |> Array.map (fun x -> x |> Array.fold folder state )

    ///Applies a function to each element of the inner arrays of the jagged array, threading an accumulator argument through the computation. 
    /// <summary>A second function is the applied to each result of the predeceding computation, again passing an accumulater through the computation </summary>
    /// <remarks></remarks>
    /// <param name="innerFolder"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>
    let fold (innerFolder : 'State1 -> 'T -> 'State1) (outerFolder : 'State2 -> 'State1 -> 'State2) (innerState : 'State1) (outerState : 'State2) ((jArray : 'T[][])) =
        jArray
        |> innerFold innerFolder innerState
        |> Array.fold outerFolder outerState

    /// <summary>Returns a new jagged array whose inner arrays only contain the elements for which the given predicate returns true</summary>
    /// <remarks></remarks>
    /// <param name="predicate"></param>
    /// <param name="jArray"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>
    let innerFilter (predicate: 'T -> bool) (jArray: 'T[][]) =
        jArray
        |> Array.map (fun x -> x |> Array.filter predicate)

    /// <summary>Applies the given function to each element in the inner arrays of the jagged array. Returns the jagged array whose inner arrays are comprised of the results x for each element where the function returns Some(x)</summary>
    /// <remarks></remarks>
    /// <param name="chooser"></param>
    /// <param name="jArray"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>
    let innerChoose (chooser: 'T -> 'U option) (jArray: 'T[][]) =
        jArray
        |> Array.map (fun x -> x |> Array.choose chooser)


    /// <summary>Shuffles each column of a jagged array separately (method: Fisher-Yates) in place</summary>
    /// <remarks></remarks>
    /// <param name="arr"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>
    let shuffleColumnWiseInPlace (arr: 'T [][]) =
        if arr.Length > 0 then 
            let rowCount    = arr.Length
            let columnCount = arr.[0].Length
            
            for ci = columnCount - 1 downto 0 do 
                for ri = rowCount downto  1 do
                    // Pick random element to swap.
                    let rj = Random.rndgen.NextInt(ri) // 0 <= j <= i-1
                    // Swap.
                    let tmp         =  arr.[rj].[ci]
                    arr.[rj].[ci]     <- arr.[ri - 1].[ci]
                    arr.[ri - 1].[ci] <- tmp
            arr            

        else
            arr


    /// <summary>Shuffles each row of a jagged array separately (method: Fisher-Yates) in place</summary>
    /// <remarks></remarks>
    /// <param name="arr"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>
    let shuffleRowWiseInPlace (arr: 'T [][]) =
        if arr.Length > 0 then 
            let rowCount    = arr.Length
            let columnCount = arr.[0].Length
            
            for ri = rowCount - 1 downto  0 do
                for ci = columnCount downto 1 do 
                    // Pick random element to swap.
                    let cj = Random.rndgen.NextInt(ci) // 0 <= j <= i-1
                    // Swap.
                    let tmp           =  arr.[ri].[cj]
                    arr.[ri].[cj]     <- arr.[ri].[ci - 1]
                    arr.[ri].[ci - 1] <- tmp
            arr            

        else
            arr


    /// <summary>Shuffels a jagged array (method: Fisher-Yates) in place</summary>
    /// <remarks></remarks>
    /// <param name="arr"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>
    let shuffleInPlace (arr: 'T [][]) =
        if arr.Length > 0 then 
            let rowCount    = arr.Length
            let columnCount = arr.[0].Length
            for ri = rowCount downto 1 do
                for ci = columnCount downto 1 do 
                    // Pick random element to swap.
                    let rj = Random.rndgen.NextInt(ri) // 0 <= j <= i-1
                    let cj = Random.rndgen.NextInt(ci)
                    // Swap.
                    let tmp               =  arr.[rj].[cj]
                    arr.[rj].[cj]         <- arr.[ri - 1].[ci - 1]
                    arr.[ri - 1].[ci - 1] <- tmp
            arr            

        else
            arr

    /// <summary>Shuffles each column of a jagged array separately (method: Fisher-Yates)</summary>
    /// <remarks></remarks>
    /// <param name="arr"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>
    let shuffleColumnWise (arr: 'T [][]) =
        let arr' = copy arr
        if arr.Length > 0 then 
            let rowCount    = arr.Length
            let columnCount = arr.[0].Length
            
            for ci = columnCount - 1 downto 0 do 
                for ri = rowCount downto  1 do
                    // Pick random element to swap.
                    let rj = Random.rndgen.NextInt(ri) // 0 <= j <= i-1
                    // Swap.
                    let tmp         =  arr'.[rj].[ci]
                    arr'.[rj].[ci]     <- arr'.[ri - 1].[ci]
                    arr'.[ri - 1].[ci] <- tmp
            arr'            

        else
            arr'


    /// <summary>Shuffles each row of a jagged array separately (method: Fisher-Yates)</summary>
    /// <remarks></remarks>
    /// <param name="arr"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>
    let shuffleRowWise (arr: 'T [][]) =
        let arr' = copy arr
        if arr.Length > 0 then 
            let rowCount    = arr.Length
            let columnCount = arr.[0].Length
            
            for ri = rowCount - 1 downto  0 do
                for ci = columnCount downto 1 do 
                    // Pick random element to swap.
                    let cj = Random.rndgen.NextInt(ci) // 0 <= j <= i-1
                    // Swap.
                    let tmp           =  arr'.[ri].[cj]
                    arr'.[ri].[cj]     <- arr'.[ri].[ci - 1]
                    arr'.[ri].[ci - 1] <- tmp
            arr'            

        else
            arr'


    /// <summary>Shuffels a jagged array (method: Fisher-Yates)</summary>
    /// <remarks></remarks>
    /// <param name="arr"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>
    let shuffle (arr: 'T [][]) =
        let arr' = copy arr
        if arr.Length > 0 then 
            let rowCount    = arr.Length
            let columnCount = arr.[0].Length
            for ri = rowCount downto 1 do
                for ci = columnCount downto 1 do 
                    // Pick random element to swap.
                    let rj = Random.rndgen.NextInt(ri) // 0 <= j <= i-1
                    let cj = Random.rndgen.NextInt(ci)
                    // Swap.
                    let tmp               =  arr'.[rj].[cj]
                    arr'.[rj].[cj]         <- arr'.[ri - 1].[ci - 1]
                    arr'.[ri - 1].[ci - 1] <- tmp
            arr'            

        else
            arr'


[<AutoOpen>]
module JaggedList =
    
    /// <summary>Transposes a jagged list (unchecked)</summary>
    /// <remarks>The resulting row count is determined by first collection length!</remarks>
    /// <param name="data"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>
    let transpose_ (data: 'T list list) =
        let rec transpose = function
            | (_::_)::_ as M -> List.map List.head M :: transpose (List.map List.tail M)
            | _ -> []
        transpose data

    /// <summary>Transposes a jagged list</summary>
    /// <remarks>all inner collections (rows) have to be of equal length</remarks>
    /// <param name="data"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>
    let transpose (data: 'T list list) =
        let colSizes = data |> List.map List.length
        if (colSizes |> List.distinct |> List.length) <> 1 then failwithf "All inner collections (rows) have to be of equal length!"
        transpose_ data

    /// <summary>Converts a jagged list into a jagged array </summary>
    /// <remarks></remarks>
    /// <param name="data"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>   
    let toJaggedList (data: 'T list list) =
        data
        |> List.map (fun l -> l |> Array.ofList)
        |> Array.ofList

    /// <summary>Converts a jagged array into a jagged list</summary>
    /// <remarks></remarks>
    /// <param name="arr"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>   
    let ofJaggedArray (arr: 'T [][]) =
        arr
        |> Array.map (fun a -> a |> List.ofArray)
        |> List.ofArray

    /// <summary>Converts a jagged Seq into a jagged list</summary>
    /// <remarks></remarks>
    /// <param name="data"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>   
    let ofJaggedSeq (data: seq<#seq<'T>>) =
        data
        |> Seq.map (fun s -> s |> List.ofSeq)
        |> List.ofSeq

    /// <summary>Converts a jagged list into a jagged seq</summary>
    /// <remarks></remarks>
    /// <param name="data"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>   
    let toJaggedSeq (data: 'T list list) =
        data
        |> Seq.map (fun s -> s |> List.toSeq) 

    /// <summary>Builds a new jagged list whose inner lists are the results of applying the given function to each of their elements.</summary>
    /// <remarks></remarks>
    /// <param name="mapping"></param>
    /// <param name="jlist"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>
    let map (mapping: 'T -> 'U) (jlist : 'T list list) =
        jlist
        |> List.map (fun x -> x |> List.map mapping)

    /// <summary>Builds a new jagged list whose inner lists are the results of applying the given function to the corresponding elements of the inner lists of the two jagged lists pairwise. <br />All corresponding inner lists must be of the same length, otherwise ArgumentException is raised.</summary>
    /// <remarks></remarks>
    /// <param name="mapping"></param>
    /// <param name="jlist1"></param>
    /// <param name="jlist2"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>
    let map2 (mapping: 'T1 -> 'T2 -> 'U) (jlist1 : 'T1 list list) (jlist2 : 'T2 list list) = 
        jlist1
        |> List.mapi (fun index x -> (List.map2 mapping x jlist2.[index]))

    /// <summary>Builds a new jagged list whose inner lists are the results of applying the given function to the corresponding elements of the inner lists of the tree jagged lists triplewise. <br />All corresponding inner lists must be of the same length, otherwise ArgumentException is raised.</summary>
    /// <remarks></remarks>
    /// <param name="mapping"></param>
    /// <param name="jlist1"></param>
    /// <param name="jlist2"></param>
    /// <param name="jlist3"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>
    let map3 (mapping : 'T1 -> 'T2 -> 'T3 -> 'U) (jlist1 : 'T1 list list) (jlist2 : 'T2 list list) (jlist3 : 'T3 list list) =
        jlist1 |> List.mapi (fun index x -> (List.map3 mapping x jlist2.[index] jlist3.[index]))

    /// <summary>Builds a new jagged list whose inner lists are the results of applying the given function to each of their elements. The integer index passed to the function indicates the index of element in the inner list being transformed.</summary>
    /// <remarks></remarks>
    /// <param name="mapping"></param>
    /// <param name="jlist"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>
    let mapi (mapping: int -> 'T -> 'U) (jlist : 'T list list) =
        jlist
        |> List.map (fun x -> x |> List.mapi mapping)

    /// <summary>Applies a function to each element of the inner lists of the jagged list, threading an accumulator argument through the computation.</summary>
    /// <remarks></remarks>
    /// <param name="folder"></param>
    /// <param name="state"></param>
    /// <param name="jlist"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>
    let innerFold (folder: 'State -> 'T -> 'State) (state: 'State) (jlist : 'T list list) =
        jlist
        |> List.map (fun x -> x |> List.fold folder state )

    ///Applies a function to each element of the inner lists of the jagged list, threading an accumulator argument through the computation. 
    /// <summary>A second function is the applied to each result of the predeceding computation, again passing an accumulater through the computation </summary>
    /// <remarks></remarks>
    /// <param name="innerFolder"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>
    let fold (innerFolder : 'State1 -> 'T -> 'State1) (outerFolder : 'State2 -> 'State1 -> 'State2) (innerState : 'State1) (outerState : 'State2) ((jlist : 'T list list)) =
        jlist
        |> innerFold innerFolder innerState
        |> List.fold outerFolder outerState

    /// <summary>Returns a new jagged list whose inner lists only contain the elements for which the given predicate returns true</summary>
    /// <remarks></remarks>
    /// <param name="predicate"></param>
    /// <param name="jlist"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>
    let innerFilter (predicate: 'T -> bool) (jlist: 'T list list) =
        jlist
        |> List.map (fun x -> x |> List.filter predicate)

    /// <summary>Applies the given function to each element in the inner lists of the jagged List. Returns the jagged list whose inner lists are comprised of the results x for each element where the function returns Some(x)</summary>
    /// <remarks></remarks>
    /// <param name="chooser"></param>
    /// <param name="jlist"></param>
    /// <returns></returns>
    /// <example>
    /// <code>
    /// </code>
    /// </example>
    let innerChoose (chooser: 'T -> 'U option) (jlist: 'T list list) =
        jlist
        |> List.map (fun x -> x |> List.choose chooser)


