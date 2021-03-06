Attribute VB_Name = "main"
Public ct() As Double
Public subseq() As Double
Public Dimension As Integer

Public Const SWAP As Integer = 0
Public Const REINSERTION As Integer = 1
Public Const OR_OPT2 As Integer = 2
Public Const OR_OPT3 As Integer = 3
Public Const TWO_OPT As Integer = 4

Public Const T As Integer = 0
Public Const C As Integer = 1
Public Const W As Integer = 2

Public improv_flag As Boolean

Public Const inf As Double = 10 ^ 16
Public Const EPSILON As Double = 10 ^ (-16)

Public Const SIZE_EMPTY As Long = 0

'Return Value:
'   -1 - Not an Array
'    0 - Empty
Public Function size( _
    ByRef r_values As Variant _
  , Optional ByVal dimensionOneBased As Long = 1 _
) As Long
  Dim Result As Long: Result = SIZE_EMPTY 'Default to Empty

  Dim lowerBound As Long
  Dim upperBound As Long
  
  On Error GoTo NormalExit
  
  lowerBound = LBound(r_values, dimensionOneBased) 'Possibly generates error
  upperBound = UBound(r_values, dimensionOneBased) 'Possibly generates error
  If (lowerBound < upperBound) Then
    Result = upperBound - lowerBound + 1 'Size greater than 1
  Else
    If (lowerBound = upperBound) Then
      Result = 1 'Size equal to 1
    End If
  End If
  
NormalExit:
  size = Result
End Function

Public Function isEmpty( _
    ByRef r_values() As Integer _
  , Optional ByVal dimensionOneBased As Long = 1 _
) As Boolean
  isEmpty = size(r_values, dimensionOneBased) = 0
End Function

Sub printArray(ByRef Arr As Variant, name As String)
    Dim output As String
    output = name & ": "
    For i = 0 To size(Arr) - 1
        output = output & Arr(i) & " "
    Next
    Debug.Print output
End Sub

Sub readData()
    Dim myFile As String, text As String, textline As String, cost As Double
    myFile = "C:\Users\victo\repos\mlp\distance_matrix"
    'myFile = Application.GetOpenFilename()
    Open myFile For Input As #1
    
    Dim i As Integer
    Dim j As Integer
    Dim Index As Integer
    
    Line Input #1, textline
    Dimension = CInt(textline)
    
    ReDim ct(Dimension, Dimension) As Double
    
    i = 0
    Do Until EOF(1)
        Line Input #1, textline
        If StrComp(textline, "EOF", vbBinaryCompare) = 0 Then
            Exit Do
        End If
        
        j = i + 1
        Do While Len(textline) > 0
            'Debug.Print textline
            Index = InStr(textline, " ")
            cost = CDbl(Left(textline, Index - 1))
            ct(i, j) = cost
            ct(j, i) = cost
            
            textline = Right(textline, Len(textline) - Index)
            j = j + 1
        Loop
        
        i = i + 1
    Loop
    
    Close
    
    'For i = 0 To Dimension - 1
        'For j = 0 To Dimension - 1
            'Debug.Print CStr(ct(i, j)) + " ";
        'Next j
        'Debug.Print ""
    'Next i
End Sub

Sub subseq_load(ByRef s() As Integer, ByRef seq() As Double)
    Dim k As Integer
    
    For i = 0 To (Dimension)
        k = 1 - i - IIf(i <> 0, 0, 1)
        
        seq(i, i, T) = 0#
        seq(i, i, C) = 0#
        seq(i, i, W) = IIf(i <> 0, 1, 0)
        
        Dim j_prev As Integer
        For j = (i + 1) To (Dimension)
            j_prev = j - 1
            
            seq(i, j, T) = ct(s(j_prev), s(j)) + seq(i, j_prev, T)
            seq(i, j, C) = seq(i, j, T) + seq(i, j_prev, C)
            seq(i, j, W) = j + k
            
        Next
    Next
    
End Sub

Sub sort_until_by(ByRef Arr() As Integer, Index As Integer, r As Integer)
    Dim tmp As Integer
    
    For i = 0 To Index
        For j = size(Arr) - 1 To i + 1 Step -1
            If ct(r, Arr(j)) < ct(r, Arr(j - 1)) Then
                tmp = Arr(j)
                Arr(j) = Arr(j - 1)
                Arr(j - 1) = tmp
            End If
        Next
    Next
    
End Sub

Function construction(alpha As Double) As Integer()
    'initial solution
    Dim s() As Integer
    ReDim s(0) As Integer
    s(0) = 0
    
    Dim cList() As Integer
    ReDim cList(Dimension - 2) As Integer
    
    'Debug.Print ""
    For j = 1 To Dimension - 1
        cList(j - 1) = j
    Next
    'Debug.Print ""
    'printArray cList, "cList"
    'Debug.Print ""
    Dim count As Integer
    count = 0
    
    Dim r As Integer
    Dim item As Integer
    Dim Index As Long
    Dim cN As Long
    r = 0
    Do While IsArrayEmpty(cList) = False
        item = CInt(size(cList) * alpha) + 1
        sort_until_by cList, item, r
        
        If size(cList) = 1 Then
            item = 0
        End If
        
        ' random between item and 0
        Index = CInt((item * Rnd) + 0)
        'Debug.Print "index ", index
        cN = cList(Index)
        'Debug.Print "cN ", cN
        InsertElementIntoArray s, UBound(s) + 1, cN
        r = cN
        DeleteArrayElement cList, Index, True
        
        'Debug.Print cN
        'Debug.Print "tamanho ", size(cList)
        
        count = count + 1
    Loop
    
    
    InsertElementIntoArray s, UBound(s) + 1, 0
    
    construction = s
    
End Function

Sub swap_f(s() As Integer, i As Integer, j As Integer)
    Dim tmp As Integer
    tmp = s(i)
    s(i) = s(j)
    s(j) = tmp
End Sub

Sub reverse(s() As Integer, i As Integer, j As Integer)
    Dim a As Integer
    Dim b As Integer
    b = j
    For a = i To CInt((j + i) \ 2)
        swap_f s, a, b
        b = b - 1
    Next a
End Sub

Sub reinsert(s() As Integer, i As Integer, j As Integer, pos As Integer)
    If i < pos Then
        For k = i To j
            tmp = s(i)
            InsertElementIntoArray s, CLng(pos), tmp
            'printArray s, "insert"
            DeleteArrayElement s, CLng(i), True
            'printArray s, "delete"
        Next
    Else
        For k = i To j
            tmp = s(j)
            DeleteArrayElement s, CLng(j), True
            InsertElementIntoArray s, CLng(pos), tmp
        Next
    End If
End Sub

Sub search_swap(s() As Integer, seq() As Double)
    Dim cost_best As Double
    Dim cost_new As Double
    
    ' concat costs vars
    Dim cost_concat_1 As Double
    Dim cost_concat_2 As Double
    Dim cost_concat_3 As Double
    Dim cost_concat_4 As Double
    
    ' pos vars
    Dim i_prev As Integer
    Dim i_next As Integer
    Dim j_next As Integer
    Dim j_prev As Integer
    
    Dim I_best As Integer
    Dim J_best As Integer
    
    cost_best = inf
    
    For i = 1 To Dimension - 2
        i_next = i + 1
        i_prev = i - 1
        
        cost_concat_1 = seq(0, i_prev, T) + ct(s(i_prev), s(i_next))
        cost_concat_2 = cost_concat_1 + seq(i, i_next, T) + ct(s(i), s(i_next + 1))
        
        cost_new = seq(0, i_prev, C) _
                + seq(i, i_next, W) * cost_concat_1 + ct(s(i_next), s(i)) _
                + seq(i_next + 1, Dimension, W) * cost_concat_2 + seq(i_next + 1, Dimension, C)
                
        If cost_new < cost_best Then
            cost_best = cost_new - EPSILON
            I_best = i
            J_best = i_next
        End If
        
        For j = i + 2 To Dimension - 1
            j_next = j + 1
            j_prev = j - 1
            
            cost_concat_1 = seq(0, i_prev, T) + ct(s(i_prev), s(j))
            cost_concat_2 = cost_concat_1 + ct(s(j), s(i_next))
            cost_concat_3 = cost_concat_2 + seq(i_next, j_prev, T) + ct(s(j_prev), s(i))
            cost_concat_4 = cost_concat_3 + ct(s(i), s(j_next))
            
            cost_new = seq(0, i_prev, C) _
                    + cost_concat_1 _
                    + seq(i_next, j_prev, W) * cost_concat_2 + seq(i_next, j_prev, C) _
                    + cost_concat_3 _
                    + seq(j_next, Dimension, W) * cost_concat_4 + seq(j_next, Dimension, C)
                    
            If cost_new < cost_best Then
                cost_best = cost_new - EPSILON
                I_best = i
                J_best = j
            End If
                    
        Next
        
    Next
    
    If cost_best < seq(0, Dimension, C) Then
        'Debug.Print "Swap"
        'Debug.Print cost_best, I_best, J_best
        swap_f s, I_best, J_best
        subseq_load s, seq
        'Debug.Print seq(0, Dimension, C)
        improv_flag = True
    End If
        
End Sub

Sub search_two_opt(s() As Integer, seq() As Double)
    Dim cost_best As Double
    Dim cost_new As Double
    
    ' concat costs vars
    Dim cost_concat_1 As Double
    Dim cost_concat_2 As Double
    Dim rev_seq_cost As Double
    
    Dim i_prev As Integer
    Dim j_next As Integer
    
    Dim I_best As Integer
    Dim J_best As Integer
    
    cost_best = inf
    
    For i = 1 To Dimension - 2
        i_prev = i - 1
        rev_seq_cost = seq(i, i + 1, T)
        
        
        For j = i + 2 To Dimension - 1
            j_next = j + 1
            
            rev_seq_cost = rev_seq_cost + ct(s(j - 1), s(j)) * (seq(i, j, W) - 1)
            
            cost_concat_1 = seq(0, i_prev, T) + ct(s(j), s(i_prev))
            
            cost_concat_2 = cost_concat_1 + seq(i, j, T) + ct(s(j_next), s(i))
            
            
            cost_new = seq(0, i_prev, C) _
                    + seq(i, j, W) * cost_concat_1 + rev_seq_cost _
                    + seq(j_next, Dimension, W) * cost_concat_2 + seq(j_next, Dimension, C)
            'Debug.Print cost_new, i, j
            'Debug.Print (seq(0, i_prev, C)), seq(i, j, W) * cost_concat_1, rev_seq_cost, seq(j_next, Dimension, W) * cost_concat_2, seq(j_next, Dimension, C)
            If cost_new < cost_best Then
                cost_best = cost_new - EPSILON
                I_best = i
                J_best = j
            End If
            'End
            
        Next
    Next
    
    If cost_best < seq(0, Dimension, C) Then
        'Debug.Print "two opt"
        'Debug.Print cost_best
        'Debug.Print I_best, J_best
        'printArray s, "antes"
        reverse s, I_best, J_best
        'printArray s, "depois"
        subseq_load s, seq
        'Debug.Print seq(0, Dimension, C)
        improv_flag = True
    End If
    
    
    
End Sub

Sub search_reinsertion(s() As Integer, seq() As Double, opt As Integer)
    Dim cost_best As Double
    Dim cost_new As Double
    cost_best = inf
    
    ' concat costs vars
    Dim cost_concat_1 As Double
    Dim cost_concat_2 As Double
    Dim cost_concat_3 As Double
    
    Dim j As Integer
    Dim i_prev As Integer
    Dim j_next As Integer
    Dim k_next As Integer
    
    Dim I_best As Integer
    Dim J_best As Integer
    Dim POS_best As Integer
    
    For i = 1 To Dimension - opt
        j = opt + i - 1
        i_prev = i - 1
        j_next = j + 1
        
        For k = 0 To i_prev - 1
            k_next = k + 1
            
            cost_concat_1 = seq(0, k, T) + ct(s(k), s(i))
            cost_concat_2 = cost_concat_1 + seq(i, j, T) + ct(s(j), s(k_next))
            cost_concat_3 = cost_concat_2 + seq(k_next, i_prev, T) + ct(s(i_prev), s(j_next))
            
            cost_new = seq(0, k, C) _
                    + seq(i, j, W) * cost_concat_1 + seq(i, j, C) _
                    + seq(k_next, i_prev, W) * cost_concat_2 + seq(k_next, i_prev, C) _
                    + seq(j_next, Dimension, W) * cost_concat_3 + seq(j_next, Dimension, C)
                    
            If cost_new < cost_best Then
                cost_best = cost_new - EPSILON
                I_best = i
                J_best = j
                POS_best = k
            End If
            
        Next
        
        For k = i + opt To Dimension - opt - 2
            k_next = k + 1
            
            cost_concat_1 = seq(0, i_prev, T) + ct(s(i_prev), s(j_next))
            cost_concat_2 = cost_concat_1 + seq(j_next, k, T) + ct(s(k), s(i))
            cost_concat_3 = cost_concat_2 + seq(i, j, T) + ct(s(j), s(k_next))
            
            cost_new = seq(0, i_prev, C) _
                    + seq(j_next, k, W) * cost_concat_1 + seq(j_next, k, C) _
                    + seq(i, j, W) * cost_concat_2 + seq(i, j, C) _
                    + seq(k_next, Dimension, W) * cost_concat_3 + seq(k_next, Dimension, C)
                    
            If cost_new < cost_best Then
                cost_best = cost_new - EPSILON
                I_best = i
                J_best = j
                POS_best = k
            End If
            
        Next
        
    Next
    
    If cost_best < seq(0, Dimension, C) Then
        'Debug.Print "reinsertion"
        'Debug.Print cost_best, I_best, J_best, POS_best
        reinsert s, I_best, J_best, POS_best + 1
        subseq_load s, seq
        'Debug.Print seq(0, Dimension, C)
        improv_flag = True
    End If
    
End Sub

Sub RVND(s() As Integer, subseq() As Double)
    Dim neighbd_list() As Variant
    neighbd_list = Array(SWAP, REINSERTION, OR_OPT2, OR_OPT3, TWO_OPT)
    'neighbd_list = Array(SWAP)
    
    'Dim s() As Integer
    'ReDim s(0) As Integer
    's(0) = 0
    'For j = (Dimension - 1) To 1 Step -1
    '    InsertElementIntoArray s, UBound(s) + 1, j
    'Next
    
    'InsertElementIntoArray s, UBound(s) + 1, 0
    'subseq_load s, subseq
    'printArray s, "s"
    'Debug.Print "s_cost"; subseq(0, Dimension, C)
    
    
    Dim i As Long
    Do While IsArrayEmpty(neighbd_list) = False
        i = CInt((size(neighbd_list) - 1) * Rnd)
        
        improv_flag = False
        'printArray neighbd_list, "NL "
        'Debug.Print "i", i, size(neighbd_list)
        Select Case neighbd_list(i)
            Case SWAP
                'Debug.Print "swap"
                search_swap s, subseq
            Case REINSERTION
                'Debug.Print "reinsert"
                search_reinsertion s, subseq, REINSERTION
            Case OR_OPT2
                'Debug.Print "or_opt2"
                search_reinsertion s, subseq, OR_OPT2
            Case OR_OPT3
                'Debug.Print "or_opt3"
                search_reinsertion s, subseq, OR_OPT3
            Case TWO_OPT
                'Debug.Print "two_opt"
                search_two_opt s, subseq
        End Select
        'End
        
        If improv_flag = True Then
            neighbd_list = Array(SWAP, REINSERTION, OR_OPT2, OR_OPT3, TWO_OPT)
        Else
            DeleteArrayElement neighbd_list, i, True
        End If
            
        'Exit Sub
    Loop
    
End Sub

Function perturb(sl As Variant)
    Dim s() As Integer
    Dim A_start As Integer, A_end As Integer
    Dim B_start As Integer, B_end As Integer
    Dim size_max As Integer, size_min As Integer, max As Integer
    
    A_start = 1
    A_end = 1
    B_start = 1
    B_end = 1
    
    size_max = size(sl) / 10
    size_max = IIf(size_max >= 2, size_max, 2)
    size_min = 2
    
    Do While (A_start <= B_start And B_start <= A_end) Or (B_start <= A_start And A_start <= B_end)
        max = size(sl) - 2 - size_max
        A_start = CInt((max * Rnd) + 1)
        A_end = A_start + CInt(((size_max - size_min) * Rnd) + size_min)
        
        B_start = CInt((max * Rnd) + 1)
        B_end = B_start + CInt(((size_max - size_min) * Rnd) + size_min)
    Loop
    
    s = sl
    
    If A_start < B_start Then
        reinsert s, B_start, B_end - 1, A_end
        reinsert s, A_start, A_end - 1, B_end
    Else
        reinsert s, A_start, A_end - 1, B_end
        reinsert s, B_start, B_end - 1, A_end
    End If
    
    perturb = s
    
End Function

Sub solve()
    readData
    Dim s() As Integer
    Dim sl() As Integer
    'ReDim s(dimension + 1) As Integer
    ReDim subseq(Dimension + 1, Dimension + 1, 3) As Double
    Dim rvnd_cost_best As Double
    Dim rnvd_cost_crnt As Double
    
    
    Dim R_size As Integer
    Dim R_table()
    R_table = Array(0#, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25)
    R_size = 26
    Dim Imax As Integer
    Imax = 10
    Dim Iils As Integer
    Iils = IIf(Dimension < 100, Dimension, 100)
    
    'printArray s, "Sinit"
    
    'Debug.Print subseq(0, Dimension, C)

    TotalTime = 0
    TotalTime = TotalTime - CDbl(MicroTimer)
    Dim s_best() As Integer
    Dim cost_best As Double
    cost_best = inf
    
    tabspace = "        "
    
    For i = 1 To Imax
        Debug.Print "[+] Local Search " & i
        Debug.Print tabspace & "[+] Constructing Inital Solution.."
        Dim alpha As Double
    
        alpha = R_table(CInt(((R_size - 1) * Rnd) + 0))
        s = construction(alpha)
        sl = s
        
        Debug.Print tabspace & "[+] Looking for the best Neighbor.."
        subseq_load s, subseq
        
        rvnd_cost_best = subseq(0, Dimension, C) - EPSILON
        Dim iterILS As Integer
        iterILS = 0
        Do While iterILS < Iils
            RVND s, subseq
            rvnd_cost_crnt = subseq(0, Dimension, C) - EPSILON
            If rvnd_cost_crnt < rvnd_cost_best Then
                'Debug.Print rvnd_cost_crnt
                rvnd_cost_best = rvnd_cost_crnt
                sl = s
                iterILS = 0
            End If
            'Debug.Print iterILS
            s = perturb(sl)
            subseq_load s, subseq
            iterILS = iterILS + 1
        Loop
        
        subseq_load sl, subseq
        sl_cost = subseq(0, Dimension, C) - EPSILON
        
        If sl_cost < cost_best Then
            cost_best = sl_cost
            s_best = sl
        End If
        
        Debug.Print tabspace & "Current best solution cost: " & cost_best
        
    Next
    
    Debug.Print "COST: " & cost_best
    TotalTime = TotalTime + MicroTimer
    Debug.Print "TIME: " & TotalTime
End Sub
