Attribute VB_Name = "Módulo1"
Public ct() As Double
Public subseq() As Double
Public Dimension As Integer

Public Const T As Integer = 0
Public Const c As Integer = 1
Public Const W As Integer = 2


Public Const SIZE_EMPTY As Long = 0

'Return Value:
'   -1 - Not an Array
'    0 - Empty
Public Function size( _
    ByRef r_values() As Integer _
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
        
        seq(i, i, T) = 0
        seq(i, i, W) = 0
        seq(i, i, c) = IIf(i <> 0, 1, 0)
        
        Dim j_prev As Integer
        For j = (i + 1) To (Dimension)
            j_prev = j - 1
            
            seq(i, j, T) = ct(s(j_prev), s(j)) + seq(i, j_prev, T)
            seq(i, j, c) = seq(i, j, T) + seq(i, j_prev, c)
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
    For j = 0 To size(s) - 1
        Debug.Print s(j)
    Next
    Debug.Print ""
    For j = 1 To Dimension - 1
        cList(j - 1) = j
    Next
    Debug.Print ""
    For j = 0 To size(cList) - 1
        Debug.Print cList(j)
    Next
    Debug.Print ""
    Dim count As Integer
    count = 0
    
    Dim r As Integer
    Dim item As Integer
    Dim Index As Long
    Dim cN As Long
    r = 0
    Do While IsArrayEmpty(cList) <> True
        item = CInt(size(cList) * alpha) + 1
        sort_until_by cList, item, r
        
        If size(cList) = 1 Then
            item = 0
        End If
        
        Index = CInt((item * Rnd) + 0)
        'Debug.Print "index ", index
        cN = cList(Index)
        'Debug.Print "cN ", cN
        InsertElementIntoArray s, UBound(s) + 1, cN
        r = cN
        DeleteArrayElement cList, Index, True
        
        'Debug.Print cN
        'Debug.Print "tamanho ", size(cList)
        
        If count = 20 Then
            Exit Function
        End If
        
        count = count + 1
    Loop
    
    
    InsertElementIntoArray s, UBound(s) + 1, 0
    
    construction = s
    
End Function


Sub solve()
    readData
    Dim s() As Integer
    'ReDim s(dimension + 1) As Integer
    ReDim subseq(Dimension + 1, Dimension + 1, 3) As Double
    
    Dim Imax As Integer
    Dim Iils As Integer
    Iils = IIf(Dimension < 100, Dimension, 100)
    
    
    s = construction(0.12)
    subseq_load s, subseq
    Debug.Print "Sinit:"
    For i = 0 To size(s) - 1
        Debug.Print s(i)
    Next
    
    Debug.Print subseq(0, Dimension, c)
End Sub
