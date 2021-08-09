Attribute VB_Name = "Módulo1"
Public ct() As Double
Public subseq() As Double
Public dimension As Integer

Public Const T As Integer = 0
Public Const C As Integer = 1
Public Const W As Integer = 2


Public Const SIZE_EMPTY As Long = 0

'Return Value:
'   -1 - Not an Array
'    0 - Empty
Public Function size( _
    ByRef r_values() As Integer _
  , Optional ByVal dimensionOneBased As Long = 1 _
) As Long
  Dim result As Long: result = SIZE_EMPTY 'Default to Empty

  Dim lowerBound As Long
  Dim upperBound As Long
  
  On Error GoTo NormalExit
  
  lowerBound = LBound(r_values, dimensionOneBased) 'Possibly generates error
  upperBound = UBound(r_values, dimensionOneBased) 'Possibly generates error
  If (lowerBound < upperBound) Then
    result = upperBound - lowerBound + 1 'Size greater than 1
  Else
    If (lowerBound = upperBound) Then
      result = 1 'Size equal to 1
    End If
  End If
  
NormalExit:
  size = result
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
    Dim index As Integer
    
    Line Input #1, textline
    dimension = CInt(textline)
    
    ReDim ct(dimension, dimension) As Double
    
    i = 0
    Do Until EOF(1)
        Line Input #1, textline
        If StrComp(textline, "EOF", vbBinaryCompare) = 0 Then
            Exit Do
        End If
        
        j = i + 1
        Do While Len(textline) > 0
            'Debug.Print textline
            index = InStr(textline, " ")
            cost = CDbl(Left(textline, index - 1))
            ct(i, j) = cost
            ct(j, i) = cost
            
            textline = Right(textline, Len(textline) - index)
            j = j + 1
        Loop
        
        i = i + 1
    Loop
    
    Close
    
    For i = 0 To dimension - 1
        For j = 0 To dimension - 1
            Debug.Print CStr(ct(i, j)) + " ";
        Next j
        Debug.Print ""
    Next i
End Sub

Sub subseq_load(ByRef s() As Integer, ByRef seq() As Double)
    Dim k As Integer
    
    For i = 0 To (dimension)
        k = 1 - i - IIf(i <> 0, 0, 1)
        
        seq(i, i, T) = 0
        seq(i, i, W) = 0
        seq(i, i, C) = IIf(i <> 0, 1, 0)
        
        Dim j_prev As Integer
        For j = (i + 1) To (dimension)
            j_prev = j - 1
            
            seq(i, j, T) = ct(s(j_prev), s(j)) + seq(i, j_prev, T)
            seq(i, j, C) = seq(i, j, T) + seq(i, j_prev, C)
            seq(i, j, W) = j + k
            
        Next
    Next
    
End Sub



Sub solve()
    readData
    Dim s() As Integer
    ReDim s(dimension + 1) As Integer
    ReDim subseq(dimension + 1, dimension + 1, 3) As Double
    
    Dim Imax As Integer
    Dim Iils As Integer
    Iils = IIf(dimension < 100, dimension, 100)
    
    
    For i = 0 To dimension - 1
        s(i) = i
        Debug.Print s(i)
    Next
    s(dimension) = 0
    
    subseq_load s, subseq
    Debug.Print subseq(0, dimension, C)
End Sub
