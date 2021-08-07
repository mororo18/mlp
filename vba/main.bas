Attribute VB_Name = "Módulo1"
Public c() As Double
Public dimension As Integer

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
    
    ReDim c(dimension, dimension) As Double
    
    i = 1
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
            c(i, j) = cost
            c(j, i) = cost
            
            textline = Right(textline, Len(textline) - index)
            j = j + 1
        Loop
        
        i = i + 1
    Loop
    
    Close
    
    For i = 1 To dimension
        For j = 1 To dimension
            Debug.Print CStr(c(i, j)) + " ";
        Next j
        Debug.Print ""
    Next i
End Sub

Sub solve()
    readData
    
End Sub
