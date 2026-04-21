unit LinearAlgebra;

{$mode ObjFPC}{$H+}{$I proprietary.txt}

interface

{ WesStat7 September 24, 2025, by Wesley R. Parsons, wespar@bellouth.net, www.wespar.com.}

uses
  BasicFuncs,
  Crt,
  DataDisplay,
  DataManager,
  Globals,
  Math,
  SysUtils;

var
  WRank: Integer;
  WDet, WCond, WTrace: FloatType;
  WSym, WSquare: Boolean;
  WInput: String;

{ Copy Functions }
function DeepCopyMatrix(const A: WMatrixType): WMatrixType; overload;
function DeepCopyMatrix(const A: CMatrixType): CMatrixType; overload;
function ConvertWMatrix(const A: WMatrixType): CMatrixType;

{ Multiply Matrices }
function TransposeWMatrices(const A: WMatrixType): WMatrixType;
function TransposeCMatrices(const A: CMatrixType): CMatrixType;
function VectorDotProduct(const Col1, Col2: WVectorType): FloatType;
function MultiplyMatrices(const A, B: WMatrixType): WMatrixType;
function MultiplyWMatrixVector(const A: WMatrixType; const B: WVectorType): CVectorType;
function MultiplyCMatrixVector(const A: CMatrixType; const B: CVectorType): CVectorType;
function MultiplyWMatrixCVector(const A: WMatrixType; const B: CVectorType): WVectorType;
function MultiplyWCMatrices(const A: WMatrixType; B: CMatrixType): CMatrixType;
function ProduceWMatrixTypeVector(const A: CMatrixType; const B: CVectorType): CVectorType;
procedure SwapMatrixRows(var A: WMatrixType; Row1, Row2: Integer);

{ Intermediate Products }
function CreateATA(const A: WMatrixType): CMatrixType;
function CreateCATA(const A: CMatrixType): CMatrixType;
function CreateXtYC(const A: WMatrixType; B: WVectorType): CVectorType;

{ Decompositions }
procedure LUDecomposition(const A: WMatrixType; var Lower, Upper: WMatrixType);
procedure QRDecomposition(const A: WMatrixType; var Q, R: WMatrixType);
procedure CholeskyDecomposition(const A: WMatrixType; var L: WMatrixType);
procedure ReportEigenResults;
procedure JacobiEigenDecomposition(const A: CMatrixType; var V: CMatrixType; var D: CVectorType; var SpecRad: Float);
procedure FillSigmaDiagonal(var A: CMatrixType; const V: CVectorType);
procedure ReportPSD(const EigenValues: CVectorType; var MinEig: FloatType);
procedure SVD(const A: WMatrixType; var U, S, VT: CMatrixType);

{ Inverses }
function TestInvertCMatrix(const A: CMatrixType): Boolean;
procedure InvertCMatrixType(const A: CMatrixType; var AInv: CMatrixType; var Okay: Boolean);
procedure InvertWMatrixType(const A: WMatrixType; var AInv: WMatrixType; var Okay: Boolean);
function InvertCMatrix(const A: CMatrixType; var Okay: Boolean): CMatrixType;
function InvertWMatrix(const A: WMatrixType; var Okay: Boolean): WMatrixType;

{ Matrix Characteristics}
function Determinant(const A: WMatrixType): FloatType;
function MatrixRank(const Data: WMatrixType): Integer;
function IsSquare(const A: WMatrixType): Boolean;
function IsSymmetric(const A: WMatrixType): Boolean;
procedure ConditionNumber(const A: WMatrixType; out Cond: FloatType; out Okay: Boolean);
function FrobeniusNormC(const A: CMatrixType): FloatType;
function FrobeniusNormW(const A: WMatrixType): FloatType;
function Trace(const A: WMatrixType): FloatType;

{ Matrix Spaces }
function RankFromR(const R: WMatrixType): Integer;
procedure RowSpace(const A: WMatrixType; var Qcol: CMatrixType);
procedure Spaces;

{ Matrix Options }
procedure MatrixCharacteristics;

implementation

{ Copy Functions }
// Copy W matrix.
function DeepCopyMatrix(const A: WMatrixType): WMatrixType; overload;
var
  i: Integer;
begin
  SetLength(Result, Length(A));
  for i := 0 to High(A) do
    DeepCopyVector(A[i], Result[i]);
end;

// Copy C matrix.
function DeepCopyMatrix(const A: CMatrixType): CMatrixType; overload;
var
  i: Integer;
begin
  SetLength(Result, Length(A));
  for i := 0 to High(A) do
    DeepCopyVector(A[i], Result[i]);
end;

// Convert W matrix to C matrix.
function ConvertWMatrix(const A: WMatrixType): CMatrixType;
var
  i: Integer;
begin
  SetLength(Result, Length(A));
  for i := 0 to High(A) do
    DeepCopyVector(A[i], Result[i]);
end;

{ Multiply Matrices }
// Transpose for regression W matrices.
function TransposeWMatrices(const A: WMatrixType): WMatrixType;
var
  i, j, m, n: Integer;
begin
  // A has n variables/columns.
  // A has m rows.
  n := Length(A);
  m := Length(A[0].Value);
  SetLength(Result, m);
  for j := 0 to m - 1 do
    SetLength(Result[j].Value, n);

  for i := 0 to n - 1 do
    for j := 0 to m - 1 do
      Result[j].Value[i] := A[i].Value[j];
  if DebugOn then
    DisplayMatrixM('Transposed Matrix', Result, T2);
end;

// Transpose for regression C matrices.
function TransposeCMatrices(const A: CMatrixType): CMatrixType;
var
  i, j, m, n: Integer;
begin
  // A has n variables/columns.
  // A has m rows.
  n := Length(A);
  m := Length(A[0]);
  SetLength(Result, m);
  for j := 0 to m - 1 do
    SetLength(Result[j], n);

  for i := 0 to n - 1 do
    for j := 0 to m - 1 do
      Result[j, i] := A[i, j];
  if DebugOn then
    DisplayMatrixM('Transposed Matrix', Result, T2);
end;

// Create dot produce of vectors, for regression, W Matrices.
function VectorDotProduct(const Col1, Col2: WVectorType): FloatType;
var
  j: Integer;
  SumCol: FloatType;
begin
  SumCol := 0.0;
  for j := 0 to High(Col1.Value) do
    SumCol := SumCol + Col1.Value[j] * Col2.Value[j];
  Result := SumCol;
end;

// Create A' * B, vectors, for regression, W Matrices.
function MultiplyWMatrixVector(const A: WMatrixType; const B: WVectorType): CVectorType;
var
  i, j, m, n: Integer;
begin
  n := Length(A);
  m := Length(A[0].Value);
  SetLength(Result, n);
  // Compute A^T * B
  for j := 0 to m - 1 do begin
    Result[j] := 0.0;
    for i := 0 to n - 1 do
      Result[j] := Result[j] + A[i].Value[j] * B.Value[i];
  end;
end;

// Create A' * B, vectors, for regression, C Matrices.
function MultiplyCMatrixVector(const A: CMatrixType; const B: CVectorType): CVectorType;
var
  i, j, m, n: Integer;
begin
  n := Length(A);
  m := Length(A[0]);
  SetLength(Result, m);
  // Compute A' * B.
  for j := 0 to m - 1 do begin
    Result[j] := 0.0;
    for i := 0 to n - 1 do
      Result[j] := Result[j] + A[i, j] * B[i];
  end;
end;

// Create A' * B for regression, W and C Matrices.
function MultiplyWMatrixCVector(const A: WMatrixType; const B: CVectorType): WVectorType;
var
  i, j, m, n: Integer;
begin
  n := Length(A);
  m := Length(A[0].Value);
  SetLength(Result.Value, m);
  // Compute A' * B.
  for j := 0 to m - 1 do begin
    Result.Value[j] := 0.0;
    for i := 0 to n - 1 do
      Result.Value[j] := Result.Value[j] + A[i].Value[j] * B[i];
  end;
end;

// Create A' * B for regression, C Matrices.
function ProduceWMatrixTypeVector(const A: CMatrixType; const B: CVectorType): CVectorType;
var
  i, j, m, n: Integer;
begin
  n := Length(A);
  m := Length(A[0]);
  SetLength(Result, n);
  // Compute A' * B.
  for j := 0 to m - 1 do begin
    Result[j] := 0.0;
    for i := 0 to n - 1 do
      Result[j] := Result[j] + A[i, j] * B[i];
  end;
end;

// Create A * B for regression, W Matrices.
function MultiplyMatrices(const A, B: WMatrixType): WMatrixType;
var
  m, n, p, i, j, k: Integer;
begin
  m := Length(A);           // Rows of A.
  n := Length(A[0].Value);  // Columns of A = Rows of B.
  p := Length(B);           // Columns of B.

  if Length(B) <> n then begin
    AddToErrorStack('Error: matrix dimensions are incompatible.');
    Result := nil;
    Exit;
  end;

  // A is m x n, B is n x p; Result is m x p.
  SetLength(Result, m);
  for i := 0 to m - 1 do
    SetLength(Result[i].Value, p);

  // Multiply.
  for j := 0 to m - 1 do
    for k := 0 to p - 1 do begin
      Result[j].Value[k] := 0.0;
      for i := 0 to n - 1 do
        Result[j].Value[k] := Result[j].Value[k] + A[j].Value[i] * B[i].Value[k];
    end;
end;

// Create A * B for regression, W and C Matrices.
function MultiplyWCMatrices(const A: WMatrixType; B: CMatrixType): CMatrixType;
var
  m, n, p, i, j, k: Integer;
begin
  m := Length(A);           // Rows of A.
  n := Length(A[0].Value);  // Columns of A = Rows of B.
  p := Length(B);           // Columns of B.

  if Length(B) <> n then begin
    AddToErrorStack('Error: matrix dimensions are incompatible.');
    Result := nil;
    Exit;
  end;
  // A is m x n, B is n x p.
  // Result is m x p.
  SetLength(Result, m);
  for i := 0 to m - 1 do
    SetLength(Result[i], p);
  // Multiply.
  for j := 0 to m - 1 do
    for k := 0 to p - 1 do begin
      Result[j, k] := 0.0;
      for i := 0 to n - 1 do
        Result[j, k] := Result[j, k] + A[j].Value[i] * B[i, k];
    end;
  end;

// Swap rowa of W Matrix.
procedure SwapMatrixRows(var A: WMatrixType; Row1, Row2: Integer);
var
  Temp: FloatType;
  i: Integer;
begin
  for i := 0 to High(A[0].Value) do begin
    Temp := A[Row1].Value[i];
    A[Row1].Value[i] := A[Row2].Value[i];
    A[Row2].Value[i] := Temp;
  end;
end;

{ Intermediate Products }
// Create X'X for regression, W Matrices.
function CreateAtA(const A: WMatrixType): CMatrixType;
var
  i, k, j, m, n: Integer;
begin
  n := Length(A);
  m := Length(A[0].Value);
  SetLength(Result, n, n);
  for i := 0 to n - 1 do
    for j := i to n - 1 do begin
      Result[i, j] := 0.0;
      for k := 0 to m - 1 do
        Result[i, j] := Result[i, j] + A[i].Value[k] * A[j].Value[k];
      Result[j, i] := Result[i, j];  // symmetry
    end;
end;

// Create X'X for regression, C Matrices.
function CreateCAtA(const A: CMatrixType): CMatrixType;
var
  i, k, j, m, n: Integer;
begin
  n := Length(A);
  m := Length(A[0]);
  SetLength(Result, n, n);
  for i := 0 to n - 1 do
    for j := i to n - 1 do begin
      Result[i, j] := 0.0;
      for k := 0 to m - 1 do
        Result[i, j] := Result[i, j] + A[i, k] * A[j, k];
      Result[j, i] := Result[i, j];  // symmetry
    end;
end;

// Create X'Y for regression, W Matrices.
function CreateXtYC(const A: WMatrixType; B: WVectorType): CVectorType;
var
  i, j, m, n: Integer;
begin
  n := Length(A);
  m := Length(A[0].Value);
  SetLength(Result, n);
  for i := 0 to n - 1 do begin
    Result[i] := 0.0;
    for j := 0 to m - 1 do
      Result[i] := Result[i] + A[i].Value[j] * B.Value[j];
  end;
end;

{ Decompositions }
// Normalize colums on a C Matrix, used in SVD.
procedure NormalizeCColumns(var A: CMatrixType);
var
  i, j, m, n: Integer;
  MeanVal, StdDev, Sum, SumSq: FloatType;
begin
  n := Length(A);
  m := Length(A[0]);

  for i := 0 to n - 1 do begin
    // Compute mean of column j.
    Sum := 0.0;
    for j := 0 to m - 1 do
      Sum := Sum + A[i, j];
    MeanVal := Sum / m;

    // Compute standard deviation of column j.
    SumSq := 0.0;
    for j := 0 to m - 1 do
      SumSq := SumSq + Sqr(A[i, j] - MeanVal);
    StdDev := Sqrt(SumSq / m - 1);

    // Normalize column j.
    if StdDev > 0 then
      for j := 0 to m - 1 do
        A[i, j] := (A[i, j] - MeanVal) / StdDev
    else
      for j := 0 to m - 1 do
        A[i, j] := 0.0; // Constant column.
  end;
end;

// LU decomposition.
procedure LUDecomposition(const A: WMatrixType; var Lower, Upper: WMatrixType);
  var
    i, j, k, m, n: Integer;
    Sum: FloatType;
    FileName: String;
  begin
    DisplayMatrixM('Before LU Decomp',A, T1);
    m := Length(A);
    n := Length(A[0].Value);
    if m <> n then begin
      AddToErrorStack('Matrix must be square.');
      Lower := nil;
      Upper := nil;
      Exit;
    end;
    SetLength(Lower, n);
    SetLength(Upper, n);

    for i := 0 to n - 1 do begin
      SetLength(Lower[i].Value, n);
      SetLength(Upper[i].Value, n);
    end;

    for i := 0 to n - 1 do begin
      // Upper Triangular.
      for j := i to n - 1 do begin
        Sum := 0;
        for k := 0 to i - 1 do
          Sum := Sum + Lower[i].Value[k] * Upper[k].Value[j];
        Upper[i].Value[j] := A[i].Value[j] - Sum;
      end;

      // Lower Triangular.
      for j := i to n - 1 do begin
        if i = j then
          Lower[i].Value[i] := 1.0
        else begin
          Sum := 0;
          for k := 0 to i - 1 do
            Sum := Sum + Lower[j].Value[k] * Upper[k].Value[i];
          if Upper[i].Value[i] = 0.0 then begin
            AddToErrorStack('Error: Zero pivot encountered.');
            Lower := nil;
            Upper := nil;
            Exit;
          end;
          Lower[j].Value[i] := (A[j].Value[i] - Sum) / Upper[i].Value[i];
        end;
      end;
    end;

    if DebugOn then begin
      DisplayMatrixM('End LU Matrix L', Lower, T1);
      DisplayMatrixM('End LU Matrix U', Upper, T1);
    end;

    if (Lower <> nil) and (Upper <> nil) then begin
      writeln('A = L * U');
      DisplayMatrixM('Matrix L', Lower, T2);
      DisplayMatrixM('Matrix U', Upper, T2);
    end;

    // Save files.
    write('Save L and U?) (y/n)');
    Readln(WInput);
    if UpCase(WInput) = 'Y' then begin
      write('Name of file for L: ');
      readln(FileName);
      SaveFile(Lower, FileName);
      write('Name of file for U: ');
      readln(FileName);
      SaveFile(Upper, FileName);
    end;
 end;

// QR decomposition.
procedure QRDecomposition(const A: WMatrixType; var Q, R: WMatrixType);
var
  i, j, k, n: Integer;
  V: WVectorType;
  QT: WMatrixType;
  FileName: String;
  LocalTolerance: FloatType;
begin
  n := Length(A);
  if not IsSquare(A) then begin
    AddToErrorStack('Error on QR Decomposition: Matrix must be square.');
    Q := nil;
    R := nil;
    Exit;
  end;
  SetLength(Q, n);
  SetLength(R, n);
  for i:= 0 to n - 1 do begin
    SetLength(Q[i].Value, n);
    SetLength(R[i].Value, n);
  end;

  for j := 0 to n - 1 do begin
    SetLength(V.Value, n);
    for i := 0 to n - 1 do
      V.Value[i] := A[i].Value[j];
    for k := 0 to j - 1 do begin
      R[k].Value[j] := VectorDotProduct(Q[k], V);
      for i := 0 to n - 1 do
        V.Value[i] := V.Value[i] - R[k].Value[j] * Q[k].Value[i];
    end;

    R[j].Value[j] := Sqrt(VectorDotProduct(V, V));

    LocalTolerance := IfThen(UserTolerance > 0, UserTolerance, Tolerance);
    if R[j].Value[j] < LocalTolerance then begin
      AddToErrorStack('Error in QR Decomposition: Columns are not linearly independent.');
      Q := nil;
      R := nil;
      Exit;
    end;

    for i := 0 to n - 1 do
      Q[j].Value[i] := V.Value[i] / R[j].Value[j];
  end;

  // Transpose Q to match A = Q * R.
  SetLength(QT, n);
  for i := 0 to n - 1 do begin
    SetLength(QT[i].Value, n);
    for j := 0 to n - 1 do
      QT[i].Value[j] := Q[j].Value[i];
  end;
  Q := DeepCopyMatrix(QT);

  if DebugOn then begin
    DisplayMatrixM('End QR Matrix Q', Q, T1);
    DisplayMatrixM('End QR Matrix R', R, T1);
  end;

  if (Q <> nil) and (R <> nil) then begin
    writeln('A = Q * R');
    DisplayMatrixM('Matrix Q', Q, T2);
    writeln;
    DisplayMatrixM('Matrix R', R, T2);
  end;

  // Save files.
  write('Save Q and R?) (y/n)');
  Readln(WInput);
  if UpCase(WInput) = 'Y' then begin
    write('Name of file for Q: ');
    readln(FileName);
    SaveFile(Q, FileName);
    write('Name of file for R: ');
    readln(FileName);
    SaveFile(R, FileName);
  end;
end;

// Cholesky decomposition. Work on zero value.
procedure CholeskyDecomposition(const A: WMatrixType; var L: WMatrixType);
var
  i, j, k, n: Integer;
  Sum, Val, LocalTolerance: FloatType;
  FileName: String;
begin
  n := Length(A);
  if not IsSymmetric(A) then begin
     writeln('Error on cholesky decomposition: matrix must be symmetric.');
     L := nil;
     Exit;
   end;

  SetLength(L, n);
  for i := 0 to n - 1 do
    SetLength(L[i].Value, n);

  if DebugOn then
    DisplayMatrixM('Start Cholesky Matrix A', A, T1);

  LocalTolerance := IfThen(UserTolerance > 0, UserTolerance, Tolerance);
  for i := 0 to n - 1 do
    for j := 0 to i do begin
      Sum := 0.0;
      for k := 0 to j - 1 do
        Sum := Sum + L[i].Value[k] * L[j].Value[k];
      Val := A[i].Value[j] - Sum;

      if Val < LocalTolerance then begin
        AddToErrorStack('Error on Cholesky decomposition: matrix not positive definite.');
        L := nil;
        Exit;
      end
      else if Val < 0.0 then Val := 0.0;

      if i = j then
        L[i].Value[j] := Sqrt(Val)
      else
        L[i].Value[j] := Val / L[j].Value[j];
  end;

  if DebugOn then
     DisplayMatrixM('End Cholesky Matrix L', L, T1);

  if L <> nil then begin
    writeln('A = L * L-Transpose');
    DisplayMatrixM('Matrix L', A, T2);
    DisplayMatrixM('Matrix L-Transpose', L, T2);
  end;

  // Save files.
  write('Save L and L-Transpose?) (y/n)');
  Readln(WInput);
  if UpCase(WInput) = 'Y' then begin
    write('Name of file for L: ');
    readln(FileName);
    SaveFile(A, FileName);
    write('Name of file for L-Transpose: ');
    readln(FileName);
    SaveFile(L, FileName);
  end;
end;

// Report the eigen decomposition results.
procedure ReportEigenResults;
var i, j, p: Integer;
  TotalVar, Percent, MinEig: FloatType;
  Loadings: CMatrixType;
  SpecRad: FloatType;
  Eigenvectors: CMatrixType;
  Eigenvalues: CVectorType;
begin
  JacobiEigenDecomposition(CreateATA(WData), Eigenvectors, Eigenvalues, SpecRad);
  p := Length(Eigenvalues);  // Covariance matrix = number of variables.
  SetLength(Loadings, p, p);

  // Compute total variance (should be p for correlation matrix).
  TotalVar := 0.0;
  for i := 0 to p - 1 do
    TotalVar := TotalVar + Eigenvalues[i];
  if TotalVar < 0.0 then begin
    writeln('Variance < 0.');
    Exit;
  end;

  WriteLn('Eigenvalues and Variance Explained:');
  for i := 0 to p - 1 do begin
    Percent := 100.0 * Eigenvalues[i] / TotalVar;
    WriteLn(('  PC' + IntToStr(i + 1)) : 5, ': Lambda = ', Eigenvalues[i] : Width : Precision,
            ', Variance = ', Percent : 5 : 2, '%');
  end;
  writeln;

  WriteLn('Eigenvectors (columns = PCs):');
  for i := 0 to p - 1 do begin
    Write(WData[i].Name : 12, ': ');
    for j := 0 to p - 1 do
      Write(Eigenvectors[i, j] : Width : Precision);
    WriteLn;
  end;

  ReportPSD(Eigenvalues, MinEig);
  if MinEig > 0 then begin
    Writeln;
    WriteLn('Loadings (Eigenvectors * Sqrt(Eigenvalues)):');
    for i := 0 to p - 1 do begin
      Write(WData[i].Name : 12, ': ');
      for j := 0 to p - 1 do begin
        Loadings[i, j] := Eigenvectors[i][j] * Sqrt(Eigenvalues[j]);
        Write(Loadings[i, j] : Width : Precision);
      end;
      Writeln;
    end;
    writeln('SpectralRadius = ', SpecRad);
  end;

  SavePartialData('Save eigen decomposition results to a file? (y/n) ', @ReportEigenResults);
  Pause;
end;

// Eigen decomposition routine.
procedure JacobiEigenDecomposition(const A: CMatrixType; var V: CMatrixType; var D: CVectorType; var SpecRad: Float);
var
  n, j, i, p, q, iter, LocalMaxIter: Integer;
  amax, theta, t, c, s, tau, temp, LocalTolerance: FloatType;
begin
  n := Length(A[0]);
  SetLength(V, n, n);
  SetLength(D, n);
  If DebugOn then
    DisplayMatrixM('Start Jacobi Matrix A ', A, T1);

  { Initialize eigenvectors to identity. }
  for i := 0 to n - 1 do
    for j := 0 to n - 1 do
      if i = j then
        V[j, i] := 1.0
      else
        V[j, i] := 0.0;

  { Initialize D to diagonal of A. }
  for i := 0 to n - 1 do
    D[i] := A[i, i];

  iter := 0;
  LocalMaxIter := IfThen(UserMaxIter > 0, UserMaxIter, MediumMaxIter);

  while iter < LocalMaxIter do begin
    { Find largest off-diagonal element. }
    amax := 0.0;
    p := 0; q := 1;
    for i := 0 to n - 1 do
      for j := i + 1 to n - 1 do
        if Abs(A[j, i]) > amax then begin
          amax := Abs(A[j, i]);
          p := i;
          q := j;
        end;

    { Convergence check. }
    LocalTolerance := IfThen(UserTolerance > 0, UserTolerance, Tolerance);
    if amax < LocalTolerance then Break;

    { Compute Jacobi rotation. }
    if A[q, p] <> 0.0 then begin
      theta := 0.5*(A[p, p] - A[q, q]) / A[q, p];
      t := 1.0 / (Abs(theta) + Sqrt(1.0 + theta * theta));
      if theta < 0.0 then t := -t;
      c := 1.0 / Sqrt(1.0 + t * t);
      s := t * c;
      tau := s / (1.0 + c);
    end
    else begin
      c := 1.0;
      s := 0.0;
      tau := 0.0;
    end;

    { Rotate matrix A. }
    temp := A[p, p];
    A[p, p] := temp - t * A[q, p];
    A[q, q] := A[q, q] + t * A[q, p];
    A[q, p] := 0.0;
    A[p, q] := 0.0;

    for i := 0 to n - 1 do
      if (i <> p) and (i <> q) then begin
        temp := A[i, p];
        A[i, p] := temp - s * (A[i, q] + tau * temp);
        A[p, i] := A[i, p];
        temp := A[i, q];
        A[i, q] := temp + s * (temp * A[i, p] + tau * temp);    { Optional, careful indexing. }
        A[q, i] := A[i, q];
      end;

    { Update eigenvectors V. }
    for i := 0 to n - 1 do begin
      temp := V[i, p];
      V[i, p] := temp - s * (V[i, q] + tau * temp);
      V[i, q] := V[i, q] + s * (temp - tau * V[i, q]);
    end;

    Inc(iter);
  end;

  { Extract eigenvalues from diagonal. }
  for i := 0 to n - 1 do
    D[i] := A[i, i];
  if Verbose or DebugOn then begin
    DisplayMatrixM('Jacobi End Matrix V ', V, T2);
    DisplayVectorM('Jacobi End Matrix D ', D, T2);
  end;

  // Compute spectral radius.
  SpecRad := 0;
  for i := 0 to n - 1 do
    if Abs(D[i]) > SpecRad then
      SpecRad := Abs(D[i]);
end;

// For SVD.
procedure FillSigmaDiagonal(var A: CMatrixType; const V: CVectorType);
var
  i, j, n, m: Integer;
begin
  n := Length(A);            // Accepts non-square matrices.
  SetLength(A, n);
  // Initialize matrix to zero.
  for i := 0 to n - 1 do
    m := Length(A[i]);
    for j := 0 to m - 1 do
      A[i, i] := 0.0;
  // Fill diagonal with square roots of V.
  for i := 0 to n - 1 do
    if V[i] >= 0.0 then
      A[i, i] := Sqrt(V[i])
    else
      A[i, i] := 0.0;        // Handle numerical noise.
end;

// For reporting eigenvector and eigenvalue results.
procedure ReportPSD(const EigenValues: CVectorType; var MinEig: FloatType);
begin
  MinEig := MinValue(EigenValues);
  if MInEig < 0 then begin
    writeln('Smallest eigenvalue less than zero: ', MinEig);
    writeln('Matrix not positive semi-definite. PCA not possible.');
  end else
    writeln('Matrix is positive semi-definite. PCA will proceed.');
end;

procedure SVD(const A: WMatrixType; var U, S, VT: CMatrixType);
var
  AtA, EigenVectorsV: CMatrixType;
  EigenValuesV: CVectorType;
  i: Integer;
  WInput: Char;
  FileName: String;
  SpecRad: FloatType;
begin
  // Compute ATA = A^T * A.
  AtA := CreateAtA(A);

  // TransposeProduct(A, ATA).
  if Verbose or DebugOn then
    DisplayMatrixM('SVD AtA', AtA, T1);

  // Compute eigenvalues and eigenvectors of ATA.
  JacobiEigenDecomposition(AtA, EigenVectorsV, EigenValuesV, SpecRad);
  if Verbose or DebugOn then begin
    DisplayMatrixM('SVD Eigenvectors ', EigenVectorsV, T1);
    DisplayVectorM('SVD Eigenvalues ', EigenValuesV, T2);
  end;

  // Compute singular values matrix S.
  SetLength(S, Length(A));
  for i := 0 to High(A) do
    SetLength(S[i], Length(A));
  FillSigmaDiagonal(S, EigenValuesV);  // Squart root of eigenvalues.

  // Transpose V.
  SetLength(VT, Length(EigenVectorsV));
  VT := TransposeCMatrices(EigenVectorsV);

  // Compute U = A * V * S^{-1}.
  SetLength(U, Length(A));
  for i := 0 to High(A) do
    SetLength(U[i], Length(A));
  U := MultiplyWCMatrices(A, EigenVectorsV);
  NormalizeCColumns(U);                // Ensure orthonormality.

  // Final output: U, S, VT.
  writeln;
  writeln('A = Sigma * S * V-Transpose');
  writeln('Matrix Sigma');
  DisplayMatrixM('', U, T2);
  writeln('Matrix S');
  DisplayMatrixM('', S, T2);
  writeln('Matrix V-Transpose');
  DisplayMatrixM('', VT, T2);

  // Save files.
  write('Save U, S, and V-Transpose?) (y/n)');
  Readln(WInput);
  if UpCase(WInput) = 'Y' then begin
    write('Name of file for Sigma: ');
    readln(FileName);
    SaveFile(U, FileName);
    write('Name of file for S: ');
    readln(FileName);
    SaveFile(S, FileName);
    write('Name of file for V-Tranpose: ');
    readln(FileName);
    SaveFile(VT, FileName);
  end;
end;

{ Inverses }
// Do I need both functions and procs? Settle on one.
// Invert C Matrix - Function.
procedure InvertCMatrixType(const A: CMatrixType; var AInv: CMatrixType; var Okay: Boolean);
var
  i, j, k, n, pivRow: Integer;
  Aug: CMatrixType;
  maxVal, tmpVal, pivot, factor, LocalTolerance: FloatType;
begin
  Okay := True;
  n := Length(A);

  // Must be square.
  if (n = 0) or (Length(A[0]) <> n) then begin
    Okay := False;
    writeln('Matrix not square');
    Exit;
  end;

  // Create augmented matrix [A | I].
  SetLength(Aug, n);
  for i := 0 to n - 1 do begin
    SetLength(Aug[i], 2 * n);
    for j := 0 to n - 1 do
      Aug[i][j] := A[i][j];
    for j := 0 to n - 1 do
      Aug[i][j + n] := Ord(i = j);
  end;

  // Gauss–Jordan with partial pivoting.
  for i := 0 to n - 1 do begin

    // Find best pivot row.
    pivRow := i;
    maxVal := Abs(Aug[i][i]);
    for k := i + 1 to n - 1 do begin
      tmpVal := Abs(Aug[k][i]);
      if tmpVal > maxVal then begin
        maxVal := tmpVal;
        pivRow := k;
      end;
    end;

    // Singular?
    LocalTolerance := IfThen(UserTolerance > 0, UserTolerance, HighTolerance);
    if maxVal < LocalTolerance then begin
      Okay := False;
      writeln('Matrix singular or nearly singular.');
      Exit;
    end;

    // Swap rows if needed.
    if pivRow <> i then begin
      for j := 0 to 2 * n - 1 do begin
        tmpVal := Aug[i][j];
        Aug[i][j] := Aug[pivRow][j];
        Aug[pivRow][j] := tmpVal;
      end;
    end;

    // Normalize pivot row.
    pivot := Aug[i][i];
    for j := 0 to 2 * n - 1 do
      Aug[i][j] := Aug[i][j] / pivot;

    // Eliminate column.
    for k := 0 to n - 1 do begin
      if k <> i then begin
        factor := Aug[k][i];
        for j := 0 to 2 * n - 1 do
          Aug[k][j] := Aug[k][j] - factor * Aug[i][j];
      end;
    end;
  end;

  // Extract inverse.
  SetLength(AInv, n);
  for i := 0 to n - 1 do begin
    SetLength(AInv[i], n);
    for j := 0 to n - 1 do
      AInv[i][j] := Aug[i][j + n];
  end;
end;

// Only tests if C Matrix is invertible. Do I need?
function TestInvertCMatrix(const A: CMatrixType): Boolean;
var
  i, j, n: Integer;
  Aug: CMatrixType;
  pivot: FloatType;
begin
  Result := True;
  n := Length(A);
  if Length(A[0]) <> n then //square
    Exit(False);

  // Initialize augmented matrix [A | I]
  SetLength(Aug, n, 2 * n);
  for j := 0 to n - 1 do
    for i := 0 to n - 1 do begin
      Aug[j, i] := A[j, i];
      Aug[j, i + n] := Ord(j = i); // Identity matrix
    end;

  // Forward elimination
  for i := 0 to n - 1 do begin
    pivot := Aug[i][i];
    if Abs(pivot) < 1e-12 then
      Result := False;
  end;
end;

// Invert C Matrix - Function.
function InvertCMatrix(const A: CMatrixType; var Okay: Boolean): CMatrixType;
var
  i, j, k, n: Integer;
  Aug: CMatrixType;
  factor, pivot: FloatType;
begin
  n := Length(A);
  Okay := True;
  if Length(A[0]) <> n then begin  //square
    if DebugOn or Verbose then
      writeln('Matrix is not square.');
    Okay := False;
    Exit;
  end;

  // Initialize augmented matrix [A | I].
  SetLength(Aug, n, 2 * n);
  for j := 0 to n - 1 do
    for i := 0 to n - 1 do begin
      Aug[j, i] := A[j, i];
      Aug[j, i + n] := Ord(j = i); // Identity matrix.
    end;

  // Forward elimination.
  for i := 0 to n - 1 do begin
    pivot := Aug[i, i];
    if Abs(pivot) < 1e-12 then begin
      if DebugOn or Verbose then
        writeln('Matrix is singular or nearly singular.');
        Okay := False;
      Exit;
    end;

    // Normalize pivot row.
    for j := 0 to 2 * n - 1 do
      Aug[i, j] := Aug[i, j] / pivot;

    // Eliminate other rows.
    for k := 0 to n - 1 do begin
      if k <> i then begin
        factor := Aug[k, i];
        for j := 0 to 2 * n - 1 do
          Aug[k, j] := Aug[k, j] - factor * Aug[i, j];
      end;
    end;
  end;

  // Extract inverse from augmented part.
  SetLength(Result, n, n);
  for i := 0 to n - 1 do
    for j := 0 to n - 1 do
      Result[i, j] := Aug[i, j + n];
end;

// Invert W Matrix - Function.
function InvertWMatrix(const A: WMatrixType; var Okay: Boolean): WMatrixType;
var
  i, j, k, n: Integer;
  Aug: CMatrixType;
  factor, pivot: FloatType;
begin
  n := Length(A);
  Okay := True;
  if Length(A[0].Value) <> n then begin
    if DebugOn or Verbose then
      writeln('Matrix is not square.');
    Okay := False;
    Exit;
  end;

  // Initialize augmented matrix [A | I].
  SetLength(Aug, n, 2 * n);
  for j := 0 to n - 1 do
    for i := 0 to n - 1 do begin
      Aug[j, i] := A[j].Value[i];
      Aug[j, i + n] := Ord(j = i); // Identity matrix.
    end;

  // Forward elimination.
  for i := 0 to n - 1 do begin
    pivot := Aug[i, i];
    if Abs(pivot) < 1e-12 then begin
      if DebugOn or Verbose then
        writeln('Matrix is singular or nearly singular.');
      Okay := False;
      Exit;
    end;

    // Normalize pivot row.
    for j := 0 to 2 * n - 1 do
      Aug[i, j] := Aug[i, j] / pivot;

    // Eliminate other rows.
    for k := 0 to n - 1 do begin
      if k <> i then begin
        factor := Aug[k, i];
        for j := 0 to 2 * n - 1 do
          Aug[k, j] := Aug[k, j] - factor * Aug[i, j];
      end;
    end;
  end;

   // Extract inverse from augmented part.
  SetLength(Result, n);
  for i := 0 to n - 1 do begin
    SetLength(Result[i].Value, n);
    for j := 0 to n - 1 do
      Result[i].Value[j] := Aug[i, j + n];
  end;
end;

// Invert W Matrix - Procedure.
procedure InvertWMatrixType(const A: WMatrixType; var AInv: WMatrixType; var Okay: Boolean);
var
  i, j, k, n: Integer;
  Aug: CMatrixType;
  factor, pivot: FloatType;
begin
  Okay := True;
  n := Length(A);
  if Length(A[0].Value) <> n then begin
      Okay := False;
      writeln('Matrix is not square.');
      Exit;
  end;

  // Initialize augmented matrix [A | I].
  SetLength(Aug, n, 2 * n);
  for j := 0 to n - 1 do
    for i := 0 to n - 1 do begin
      Aug[j, i] := A[j].Value[i];
      Aug[j, i + n] := Ord(j = i); // Identity matrix.
    end;

  // Forward elimination.
  for i := 0 to n - 1 do begin
    pivot := Aug[i][i];
    if Abs(pivot) < 1e-12 then
      writeln('Matrix is singular or nearly singular.');

    // Normalize pivot row.
    for j := 0 to 2 * n - 1 do
      Aug[i][j] := Aug[i][j] / pivot;

    // Eliminate other rows.
    for k := 0 to n - 1 do begin
      if k <> i then begin
        factor := Aug[k][i];
        for j := 0 to 2 * n - 1 do
          Aug[k][j] := Aug[k][j] - factor * Aug[i][j];
      end;
    end;
  end;

  // Extract inverse from augmented part.
  SetLength(AInv, n);
  for i := 0 to n - 1 do begin
    SetLength(AInv[i].Value, n);
    for j := 0 to n - 1 do
      AInv[i].Value[j] := Aug[i, j + n];
  end;
end;

{ Matrix Characteristics }
// Determinant.
function Determinant(const A: WMatrixType): FloatType;
var
  X: WMatrixType;
  i, j, k, n: Integer;
  Factor, Det, LocalTolerance: FloatType;
  PivotFound: Boolean;
begin
  if not IsSquare(A) then begin
    AddToErrorStack('Error on determinant of matrix: not square.');
    Result := NaN;
    Exit;
  end;

  X := DeepCopyMatrix(A);
  n := Length(X);
  Det := 1.0;

  LocalTolerance := IfThen(UserTolerance > 0, UserTolerance, HighTolerance);
  for k := 0 to n - 1 do begin
    pivotFound := False;
    for j := k to n - 1 do
      if Abs(X[j].Value[k]) > LocalTolerance then begin
        PivotFound := True;
        Break;
      end;
    if not pivotFound then begin
      Result := 0.0;
      Exit;
    end;

    if j <> k then begin
      for i := 0 to n - 1 do begin
        Factor := X[j].Value[i];
        X[j].Value[i] := X[k].Value[i];
        X[k].Value[i] := Factor;
      end;
      Det := -Det;
    end;

    Det := Det * X[k].Value[k];

    for j := k + 1 to n - 1 do begin
      Factor := X[j].Value[k] / X[k].Value[k];
      for i := k to n - 1 do
        X[j].Value[i] := X[j].Value[i] - Factor * X[k].Value[i];
    end;
  end;

  Result := Det;
end;

// Compute rank of a matrix using Gaussian elimination.
function MatrixRank(const Data: WMatrixType): Integer;
var
  i, j, rows, cols, lead, pivotRow: Integer;
  factor, pivot, LocalTolerance: FloatType;
  A: WMatrixType;
begin
  if Length(Data) = 0 then begin
    AddToErrorStack('Error on Matrix Rank: Matrix has no rows.');
    Exit(0);
  end;

  A := DeepCopyMatrix(Data); // This is opposite of usual use.
  rows := Length(A);
  cols := Length(A[0].Value);

  lead := 0;
  LocalTolerance := IfThen(UserTolerance > 0, UserTolerance, LowTolerance);
  for pivotRow := 0 to rows - 1 do  begin
    if lead >= cols then Break;

    // Find pivot in current column.
    i := pivotRow;
    while (Abs(A[i].Value[lead]) < LocalTolerance) do begin
      Inc(i);
      if i = rows then begin
        i := pivotRow;
        Inc(lead);
        if lead = cols then Break;
      end;
    end;
    if lead >= cols then Break;

    // Swap rows.
    if i <> pivotRow then
      SwapMatrixRows(A, i, pivotRow);

    // Normalize pivot row.
    pivot := A[pivotRow].Value[lead];
    for j := 0 to cols - 1 do
      A[pivotRow].Value[j] := A[pivotRow].Value[j] / pivot;

    // Eliminate column.
    for i := 0 to rows - 1 do
      if i <> pivotRow then begin
        factor := A[i].Value[lead];
        for j := 0 to cols - 1 do
          A[i].Value[j] := A[i].Value[j] - factor * A[pivotRow].Value[j];
      end;
    Inc(lead);
  end;

  // Count non-zero rows.
  Result := 0;
  for i := 0 to rows - 1 do begin
    for j := 0 to cols - 1 do
      if Abs(A[i].Value[j]) > 1e-10 then begin
        Inc(Result);
        Break;
      end;
  end;
end;

// Square.
function IsSquare(const A: WMatrixType): Boolean;
var
  i, n: Integer;
begin
  n := Length(A);
  for i := 0 to n - 1 do
    if n <> Length(A[i].Value) then
      Exit(False);
  Result := True;
end;

// Symmetric.
function IsSymmetric(const A: WMatrixType): Boolean;
var
  i, j, m: Integer;
begin
  m := Length(A);
  if m = 0 then Exit(False);
  if Length(A[0].Value) <> m then begin
    AddToErrorStack('Error: matrix dimensions are incompatible.');
    Result := False;
    Exit;
  end;

  for j := 0 to m - 1 do
    for i := j + 1 to m - 1 do
      if Abs(A[i].Value[j] - A[j].Value[i]) > Tolerance then
        Exit(False);
  Result := True;
end;

// Infinity Norm, use in condition number.
function MatrixInfNorm(const A: WMatrixType): FloatType;
var
  i, j, n: Integer;
  rowsum, maxsum: FloatType;
begin
  n := Length(A);
  maxsum := 0.0;

  for i := 0 to n - 1 do begin
    rowsum := 0.0;
    for j := 0 to High(A[i].Value) do
      rowsum := rowsum + Abs(A[i].Value[j]);

    if rowsum > maxsum then
      maxsum := rowsum;
  end;

  Result := maxsum;
end;

// New Version 12/30/2025. ChatGPT.
procedure ConditionNumber(const A: WMatrixType; out Cond: FloatType; out Okay: Boolean);
var
  AInv: WMatrixType;
  normA, normAInv: FloatType;
begin
  Okay := True;

  normA := MatrixInfNorm(A);
  if normA = 0.0 then begin
    Okay := False;
    Cond := 0.0;
    Exit;
  end;

  InvertWMatrixType(A, AInv, Okay);
  if not Okay then begin
    Cond := 0.0;
    Exit;
  end;

  normAInv := MatrixInfNorm(AInv);

  Cond := normA * normAInv;
end;

// Old Version: Condition number: norm(A) * norm(inv(A)).
function ConditionNumber(const A: WMatrixType): FloatType;
var
  ATemp, Ainv: CMatrixType;
  Okay: Boolean;                     //Issue c and w matrices
  i: Integer;
begin
  ATemp := ConvertWMatrix(A);
  SetLength(AInv, Length(A));

  for i := 0 to Length(ATemp) - 1 do
    SetLength(AInv[i], Length(ATemp[i]));

  InvertCMatrixType(ATemp, AInv, Okay);

  if not Okay then begin
    AddToErrorStack('Error on Condition Number: matrix not invertible.');
    Result := 0.0;
    Exit;
  end
  else Result := FrobeniusNormW(A) * FrobeniusNormC(Ainv);
end;

// Frobenius Norm, C Matrix (used for condition number).
function FrobeniusNormC(const A: CMatrixType): FloatType;
var
  i, j: Integer;
begin
  Result := 0;

  for j := 0 to High(A) do
    for i := 0 to High(A[j]) do
      Result := Result + Sqr(A[j, i]);

  Result := Sqrt(Result);
end;

// Frobenius Norm, W Matrix.
function FrobeniusNormW(const A: WMatrixType): FloatType;
var
  i, j: Integer;
begin
  Result := 0;

  for j := 0 to High(A) do
    for i := 0 to High(A[j].Value) do
      Result := Result + Sqr(A[j].Value[i]);

  Result := Sqrt(Result);
end;

// Trace.
function Trace(const A: WMatrixType): FloatType;
var
  i, n: Integer;
begin
  Result := 0.0;
  n := Min(Length(A), Length(A[0].Value));  // Handles rectangular matrices safely.

  for i := 0 to n - 1 do
    Result := Result + A[i].Value[i];
end;

{ Matrix Spaces }
// Function used in Spaces.
function RankFromR(const R: WMatrixType): Integer;
var
  i, n: Integer;
  LocalTolerance: FloatType;
begin
  Result := 0;
  n := Min(Length(R), Length(R[0].Value));       // Number of diagonals.

  LocalTolerance := IfThen(UserTolerance > 0, UserTolerance, LowTolerance);
  for i := 0 to n - 1 do
    if Abs(R[i].Value[i]) > LowTolerance then
      Inc(Result);
end;

// Column.
procedure ColumnSpace(const A: WMatrixType; var QCol: CMatrixType);
var
  Q, R: WMatrixType;
  Rank, i, j, m: Integer;
begin
  QRDecomposition(A, Q, R);
  Rank := RankFromR(R);      // Count non-zero diagonals.
  m := Length(Q);
  SetLength(QCol, m);

  for i := 0 to m - 1 do begin
    SetLength(QCol[i], Rank);
    for j := 0 to rank - 1 do
      Qcol[i, j] := Q[i].Value[j];
  end;
end;

// Row.
procedure RowSpace(const A: WMatrixType; var QCol: CMatrixType);
begin
  ColumnSpace(TransposeWMatrices(A), Qcol);
end;

// Null.
procedure NullSpace(const A: WMatrixType; var NS: CMatrixType);
var
  AtA, V: CMatrixType;
  Eigvals: CVectorType;
  i, j, n, r: Integer;
  SpecRad: FloatType;
  LocalTolerance: FloatType;
begin
  // Compute Aᵀ * A.
  AtA := CreateAtA(A);

  // Compute eigen decomposition of AtA.
  JacobiEigendecomposition(AtA, V, Eigvals, SpecRad);

  // Select eigenvectors with near-zero eigenvalues.
  LocalTolerance := IfThen(UserTolerance > 0, UserTolerance, LowTolerance);
  n := Length(Eigvals);
  r := 0;
  for i := 0 to n - 1 do
    if Eigvals[i] < LowTolerance then
      Inc(r);
  SetLength(NS, n, r);
  r := 0;

  for i := 0 to n - 1 do
    if Eigvals[i] < LowTolerance then begin
      for j := 0 to n - 1 do
        NS[j, r] := V[j, i];
      Inc(r);
    end;
end;

// Left Null.
procedure LeftNullSpace(const A: WMatrixType; var LNS: CMatrixType);
begin
  NullSpace(TransposeWMatrices(A), LNS);
end;

// Report each space.
procedure ReportSpace(const Name: string; const Space: CMatrixType);
var
  m, n: Integer;
begin
  m := Length(Space);
  if m > 0 then
    n := Length(Space[0])
  else
    n := 0;

  if (m = 0) or (n = 0) then
    writeln(Name, ' is trivial (dimension = 0).')
  else begin
    writeln(Name, ' basis has dimension ', n, '.');
    DisplayMatrixM(Name, Space, T3);
  end;
end;

// Report all the spaces.
procedure Spaces;
var
  NS, LNS, RS, CS: CMatrixType;
begin
  writeln('Bases for spaces of A.');
  ColumnSpace(WData, CS);
  ReportSpace('Column Space', CS);
  RowSpace(WData, RS);
  ReportSpace('Row Space', RS);
  NullSpace(WData, NS);
  ReportSpace('Null Space', NS);
  LeftNullSpace(WData, LNS);
  ReportSpace('Left Null Space', LNS);

  SavePartialData('Save spaces to a file? (y/n) ', @Spaces);
  Pause;
end;

{ Matrix Options }
// Report many matrix characteristics.
procedure MatrixCharacteristics;
var
  m, n : Integer;
  Okay: Boolean;
begin
  m := Length(WData);
  n := Length(WData[0].Value);

  writeln('Matrix A has ', n, ' columns and ', m, ' rows. It is:');
  DisplayMatrixM('Matrix ', WData, T2);
  writeln;

  WRank := MatrixRank(WData);
  writeTabln('The rank is ', WRank);

  WDet := Determinant(WData);
  writeTabln('The determinant is ', WDet);

  ConditionNumber(WData, WCond, Okay);
  if not Okay then writeln('Condition number cannot be computed.')
  else writeTabln('The matrix has condition ', WCond);

  WSquare := IsSquare(WData);
  write('The matrix is ');
  If WSquare = False then write('not ');
  writeln('square.');

  WTrace := Trace(WData);
  if WSquare then writeTabln('The matrix has trace ', WTrace);

  WSym := IsSymmetric(WData);
  write('The matrix is ');
  If WSym = False then write('not ');
  writeln('symmetric.');

  if FirstPass then
    SavePartialData('Save matrix characteristics to a file? (y/n) ', @MatrixCharacteristics);
  Pause;
end;

end.

