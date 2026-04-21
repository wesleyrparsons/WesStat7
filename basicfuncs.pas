unit BasicFuncs;

{$mode ObjFPC}{$H+}{$I proprietary.txt}

interface

{ WesStat7 September 24, 2025, by Wesley R. Parsons, wespar@bellouth.net, www.wespar.com.}

uses
  Crt,
  DataDisplay,
  Globals,
  Math,
  SysUtils;

{Cleaning Routines}
procedure CleanUpZero(var x: FloatType);
procedure CleanUpOne(var x: FloatType);
procedure CleanUpBoth(var x: FloatType);
procedure CleanUpWVector(var Col: WVectorType);
function IsEffectivelyEqual(const a, b: FloatType): Boolean;

{Sorting Routines}
procedure Swap(var A, B: Integer); overload;
procedure Swap(var A, B: FloatType); overload;
procedure QuickSort(var A: WVectorType; Low, High: Integer);
procedure QuickSortIndex(const keys: WVectorType; var Indices: IVectorType);
function Partition(const lo, hi: Integer; const keys: WVectorType; var indices: IVectorType): Integer;
procedure QuickSelect(var A: WVectorType; l, r, k: Integer);
function MedianOfSorted(const Arr: CVectorType): FloatType;

{Matrix Routines}
procedure DeepCopyVector(const Source: WVectorType; var Target: WVectorType); overload;
procedure DeepCopyVector(const Source: WVectorType; var Target: CVectorType); overload;
procedure DeepCopyVector(const Source: CVectorType; var Target: CVectorType); overload;
procedure DeepCopyWMatrix(const Source: WMatrixType; var Target: WMatrixType);

{Random Routines}
function RandomNormal(const mu, sigma: FloatType): FloatType;

{Utility Routines}
procedure Pause;
procedure ConvertNaNinData;
function Factorial(const n: Integer): QWord;
procedure ShowOptions;

implementation

{Cleaning Routines}
procedure CleanUpZero(var x: FloatType);
begin
  if IsNan(x) or IsInfinite(x) then Exit;
  if Abs(x) < HighTolerance then
    x := 0.0;
end;

procedure CleanUpOne(var x: FloatType);
begin
  if IsNan(x) or IsInfinite(x) then Exit;
  if Abs(x - 1.0) < HighTolerance then
    x := 1.0
  else if Abs(x + 1.0) < HighTolerance then
    x := -1.0;
end;

procedure CleanUpBoth(var x: FloatType);
begin
  CleanUpZero(x);
  CleanUpOne(x);
end;

procedure CleanUpWVector(var Col: WVectorType);
var
  j: Integer;
begin
  for j := 0 to High(Col.Value) do
    CleanUpBoth(Col.Value[j]);
  CleanUpBoth(Col.Mean);
  CleanUpBoth(Col.Median);
  CleanUpBoth(Col.StdDev);
  CleanUpBoth(Col.Min);
  CleanUpBoth(Col.Max);
end;

function IsEffectivelyEqual(const a, b: FloatType): Boolean;
begin
 if IsNan(a) or IsInfinite(a) or IsNan(b) or IsInfinite(b) then Exit;
 Result := Abs(a - b) <= HighTolerance * Max(1.0, Max(Abs(a), Abs(b)));
end;

{Sorting Routines}
// Swap for floats.
procedure Swap(var A, B: FloatType); overload;
var
  Temp: FloatType;
begin
  Temp := A;
  A := B;
  B := Temp;
end;

// Swap for integers.
procedure Swap(var A, B: Integer); overload;
var
  Temp: Integer;
begin
  Temp := A;
  A := B;
  B := Temp;
end;

// QuickSort procedure.
procedure QuickSort(var A: WVectorType; Low, High: Integer);
var
  i, j: Integer;
  pivot, temp: FloatType;
begin
  if Low >= High then Exit;
  i := Low;
  j := High;
  pivot := A.Value[(Low + High) div 2];
  repeat
    while A.Value[i] < pivot do Inc(i);
    while A.Value[j] > pivot do Dec(j);
    if i <= j then begin
      temp := A.Value[i];
      A.Value[i] := A.Value[j];
      A.Value[j] := temp;
      Inc(i);
      Dec(j);
    end;
  until i > j;
  if Low < j then QuickSort(A, Low, j);
  if i < High then QuickSort(A, i, High);
end;

// QuickSort with index procedure.
procedure QuickSortIndex(const keys: WVectorType; var Indices: IVectorType);

  procedure Sort(lo, hi: Integer);
  var
    p: Integer;
  begin
    if lo < hi then begin
      p := Partition(lo, hi, keys, indices);
      sort(lo, p - 1);
      sort(p + 1, hi);
    end;
  end;

var
  len: Integer;
begin
  len := Length(keys.Value);
  if len > 1 then
    Sort(0, len - 1);
end;

// Partition procedure for QuickSelect.
function Partition(const lo, hi: Integer; const keys: WVectorType; var indices: IVectorType): Integer;
var
  Pivotkey: FloatType;
  i, j, temp: Integer;
begin
  Pivotkey := keys.Value[indices[hi]];
  i := lo - 1;
  for j := lo to hi - 1 do
    if keys.Value[indices[j]] <= Pivotkey then begin
      Inc(i);
      if i <> j then begin
        temp := indices[i];
        indices[i] := indices[j];
        indices[j] := temp;
      end;
    end;
  temp := indices[i + 1];
  indices[i + 1] := indices[hi];
  indices[hi] := temp;
  Result := i + 1;
end;

// Quick routine to find kth or median value.
procedure QuickSelect(var A: WVectorType; l, r, k: Integer);
var
  i, storeIndex, pivotIndex: Integer;
  pivotValue, tmp: FloatType;
begin
  while l < r do begin
    // Choose pivot randomly.
    pivotIndex := l + Random(r - l + 1);
    pivotValue := A.Value[pivotIndex];

    // Move pivot to end.
    tmp := A.Value[pivotIndex];
    A.Value[pivotIndex] := A.Value[r];
    A.Value[r] := tmp;

    // Partition (Lomuto).
    storeIndex := l;
    for i := l to r - 1 do
      if A.Value[i] < pivotValue then begin
        tmp := A.Value[i];
        A.Value[i] := A.Value[storeIndex];
        A.Value[storeIndex] := tmp;
        Inc(storeIndex);
      end;

    // Move pivot to final place.
    tmp := A.Value[storeIndex];
    A.Value[storeIndex] := A.Value[r];
    A.Value[r] := tmp;

    // Now: storeIndex = pivot final position.
    if k = storeIndex then Exit
    else if k < storeIndex then
      r := storeIndex - 1    // Recurse into left partition.
    else
      l := storeIndex + 1;   // Recurse into right partition.
  end;
end;

// Finds median of already-sorted vector.
function MedianOfSorted(const Arr: CVectorType): FloatType;
var
  len: Integer;
begin
  len := Length(Arr);
  if len = 0 then Exit(NaN);
  if Odd(len) then
    Result := Arr[len div 2]
  else
    Result := 0.5 * (Arr[(len div 2) - 1] + Arr[len div 2]);
end;

{Matrix Routines}
procedure DeepCopyVector(const Source: WVectorType; var Target: WVectorType); overload;
var
  m: Integer;
begin
  m := Length(Source.Value);
  SetLength(Target.Value, m);
  Target.Name := Source.Name;
  Target.Median := Source.Median;
  Target.Mean := Source.Mean;
  Target.StdDev := Source.StdDev;
  Target.Min := Source.Min;
  Target.Max := Source.Max;
  Target.Value := Copy(Source.Value, 0, Length(Source.Value));
  Target.IsNan := Copy(Source.IsNan, 0, Length(Source.IsNan));
{ // Number of bytes in the Value array.
  nBytes := Length(Source.Value) * SizeOf(FloatType);
  // Fast block copy of the entire float array.
  if nBytes > 0 then
     Move(Source.Value[0], Target.Value[0], nBytes);
  // Number of bytes in the IsNan array.
  nBytes := Length(Source.Value) * SizeOf(FloatType);
  // Fast block copy of the entire boolean array.
  if nBytes > 0 then
     Move(Source.IsNan[0], Target.IsNan[0], nBytes);}
end;

// Copies WVector into CVector.
procedure DeepCopyVector(const Source: WVectorType; var Target: CVectorType); overload;
var
  m: Integer;
begin
  m := Length(Source.Value);
  SetLength(Target, m);
  Target := Copy(Source.Value, 0, Length(Source.Value));
end;

procedure DeepCopyVector(const Source: CVectorType; var Target: CVectorType); overload;
var
  m: Integer;
begin
  m := Length(Source);
  SetLength(Target, m);
  Target := Copy(Source, 0, Length(Source));
end;

procedure DeepCopyWMatrix(const Source: WMatrixType; var Target: WMatrixType);
var
  i, n: Integer;
begin
  n := Length(Source);
  SetLength(Target, n);
  for i := 0 to n - 1 do begin
    SetLength(Target[i].Value, Length(Source[i].Value));
    DeepCopyVector(Source[i], Target[i]);
  end;
end;

{ Random Number Routines }
function RandomNormal(const mu, sigma: FloatType): FloatType;
var
  u1, u2: FloatType;
begin
  repeat
    u1 := Random;
  until u1 > 0.0;
  u2 := Random;
  Result := mu + sigma * sqrt(-2.0 * ln(u1)) * cos(TwoPi * u2);
end;

{ Utility Routines }
procedure Pause;
begin
  Write('Hit <CR> to continue.');
  readln;
end;

procedure ConvertNaNinData;
var
  i, j: Integer;
begin
  for i := 0 to nCol - 1 do
    for j := 0 to nRow - 1 do begin
      Case MissingData of
        MDMean: if IsNan(WData[i].Value[j]) then
          WData[i].Value[j] := WData[i].Mean;
        MDMedian: if IsNan(WData[i].Value[j]) then
          WData[i].Value[j] := WData[i].Median;
        MDZero: if IsNan(WData[i].Value[j]) then
          WData[i].Value[j] := 0.0;
      end;

  if IsNan(WData[i].Value[j]) then writeln('hitnan ',           WData[i].Value[j]);
    end;
end;

function Factorial(const n: Integer): QWord;
var
  i: Integer;
  res: QWord;
begin
  if n < 0 then begin
    AddToErrorStack('Error on factorial: negative argument');
    Exit(0);
  end;
  res := 1;
  for i := 2 to n do
    res := res * QWord(i);
  Result := res;
end;

function CheckRangedNumber(const S: String; const x, R1, R2: Integer): Boolean; overload;
begin
  Result := False;
  if (x >= R1) and (x <= R2) then
    Result := True
  else
    writeln(S, ' is not in the required range ', R1 : 3, ' to ', R2 : 3, '.');
end;

function CheckRangedNumber(const S: String; const x, R1, R2: FloatType): Boolean; overload;
begin
  Result := False;
  if (x >= R1) and (x <= R2) then
    Result := True
  else
    writeln(S, ' is not in the required range ', R1 : 3, ' to ', R2 : 3, '.');
end;

procedure ShowOptions;
begin
  writeln('Settings & Options: ');
  writeln('  Precision: ', Precision: 8);
  writeln('  Width: ', Width: 8);
  writeln('  Use intercept for regression: ', UseIntercept);
  writeln('  Matrix inversion epsilon: ', SmartFloat(HighTolerance));
  writeln('  Matrix rank tolerance: ', SmartFloat(LowTolerance));
  writeln('  Regression convergence epsilon: ', SmartFloat(HighTolerance));
  writeln('  Regression ML iterations: ', SlowMaxIter: 8); // Check.
  writeln('  Global tolerance: ', SmartFloat(Tolerance));  // Add other tols?
  writeln('  Global maximum iterations: ', MediumMaxIter); // Add two other iters?
  writeln('  Verbose: ', Verbose);
  writeln('  DebugOn: ', DebugOn);
  writeln('  Large threshold for scientific notation: ', SmartFloat(LargeThreshold));
  writeln('  Small threshold for scientific notation: ', SmartFloat(SmallThreshold));
  writeln('  Tab for dsiplay: ', Tab: 8);
  writeln('  Show significance value: ', ShowSignificance);
  writeln('  Use sample statistic: ', Sample);
  write('  Default ML estimation type: ');
  if RegressionMode = MLEGD then
    writeln('Gradient Descent Method')
  else if RegressionMode = MLENR then
    writeln('Newton-Raphson Method')
  else if RegressionMode = MLENest then
    writeln('Nesterov Method');
end;

end.

