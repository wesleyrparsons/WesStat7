unit Nominal;

{$mode ObjFPC}{$H+}{$I proprietary.txt}

interface

{ WesStat7 September 24, 2025, by Wesley R. Parsons, wespar@bellouth.net, www.wespar.com.}

uses
  BasicFuncs,
  Correlations,
  DataDisplay,
  DataManager,
  Globals,
  Math,
  SysUtils;

procedure ChiSquareTest(const Data: WMatrixType; out Chi2, PValue: FloatType);
procedure NominalStatistics(const WVar: IVectorType);

implementation

function Lambda(const Data: WMatrixType): FloatType;
var
  Count, m, n, i, j, Total, MaxCol, MaxOverall, SumMaxCol: Integer;
begin
  // Initialize.
  m := Length(Data[0].Value);
  n := Length(Data);
  Total := 0;
  MaxOverall := 0;
  SumMaxCol := 0;

  // Count.
  for j := 0 to m - 1 do begin
    MaxCol := 0;
    for i := 0 to n - 1 do begin
      Count := Round(Data[i].Value[j]);
      Total := Total + Count;
      if Count > MaxCol then MaxCol := Count;
      if Count > MaxOverall then MaxOverall := Count;
    end;
    SumMaxCol := SumMaxCol + MaxCol;
  end;

  if Total = 0 then begin
    AddToErrorStack('Error on Lambda: zero total count.');
    Exit(NaN);
  end;

  Result := (SumMaxCol - MaxOverall) / (Total - MaxOverall);
end;

function PhiCoefficient(const Data: WMatrixType): FloatType;
var
  Chi2, PValue: FloatType;
  i, j, m, n, p: Integer;
begin
  // Initialize.
  m := Length(Data[0].Value);
  n := Length(Data);
  if (m <> 2) or (n <> 2) then begin
    AddToErrorStack('Error on Phi: coefficient requires a 2 x 2 table.');
    Exit(NaN);
  end;

  // Count.
  p := 0;
  for j := 0 to 1 do
    for i := 0 to 1 do
      p := p + Round(Data[i].Value[j]);

  ChiSquareTest(Data, Chi2, PValue);
  if IsNan(Chi2) or (p = 0) then begin
    AddToErrorStack('Error on Phi: invalid input or zero total.');
    Exit(NaN);
  end;

  Result := Sqrt(Chi2 / n);
end;

function CramersV(const Data: WMatrixType): FloatType;
var
  Chi2, PValue: FloatType;
  i, j, m, n, p, df: Integer;
begin
  // Initialize.
  m := Length(Data[0].Value);
  n := Length(Data);
  p := 0;

  // Count.
  for j := 0 to m - 1 do
    for i := 0 to n - 1 do
      p := p + Round(Data[i].Value[j]);

  ChiSquareTest(Data, Chi2, PValue);
  if IsNan(Chi2) or (p = 0) then begin
    AddToErrorStack('Error on Cramér’s V: invalid input or zero total.');
    Exit(NaN);
  end;

  df := Min(m - 1, n - 1);
  if df <= 0 then begin
    AddToErrorStack('Error on Cramér''s V: degrees of freedom must be positive.');
    Exit(NaN);
  end;

  Result := Sqrt(Chi2 / (N * df));
end;

{ Chi-Square Test for m x n contingency table }
procedure ChiSquareTest(const Data: WMatrixType; out Chi2, PValue: FloatType);
var
  m, n, j, i, df: Integer;
  RowTotals, ColTotals: CVectorType;
  GrandTotal, Expected, LocalTolerance: FloatType;
begin
  if DebugOn then
    DisplayMatrixM('Chi Square', Data, T1);

  // Get dimensions.
  n := Length(Data);
  m := Length(Data[0].Value);
  if m < 2 then begin
    AddToErrorStack('Error on chi square: too few rows.');
    Chi2 := NaN;
    PValue := NaN;
    Exit;
  end;
  if n < 2 then begin
    AddToErrorStack('Error on chi square: too few columns.');
    Chi2 := NaN;
    PValue := NaN;
    Exit;
  end;

  // Initialize arrays for row and column totals.
  SetLength(RowTotals, m);
  SetLength(ColTotals, n);
  GrandTotal := 0.0;

  // Initialize column totals.
  for i := 0 to n - 1 do begin
    ColTotals[i] := 0.0;
    if Length(Data[i].Value) <> m then begin
      AddToErrorStack('Error on chi square: inconsistent row length.');
      Chi2 := NaN;
      PValue := NaN;
      Exit;
    end;
  end;

  // Compute row and column totals.
  for j := 0 to m - 1 do begin
    RowTotals[j] := 0.0;
    for i := 0 to n - 1 do begin
      if Data[i].Value[j] < 0 then begin
        AddToErrorStack('Error on chi square: negative values.');
        Exit;
      end;
      RowTotals[j] := RowTotals[j] + Data[i].Value[j];
      ColTotals[i] := ColTotals[i] + Data[i].Value[j];
    end;
    GrandTotal := GrandTotal + RowTotals[j];
  end;

  // Check for zero totals.
  if GrandTotal <= 0 then begin
    AddToErrorStack('Error on chi square: zero grand total.');
    Chi2 := NaN;
    PValue := NaN;
    Exit;
  end;

  // Compute chi-square statistic.
  Chi2 := 0.0;
  LocalTolerance := IfThen(UserTolerance > 0.0, UserTolerance, Tolerance);
  for j := 0 to m - 1 do
    for i := 0 to n - 1 do begin
      Expected := (RowTotals[j] * ColTotals[i]) / GrandTotal;
      if Expected < LocalTolerance then begin
        AddToErrorStack('Error on chi square: expected frequency too small.');
        Inc(iES);
        Exit;
      end;
      Chi2 := Chi2 + Sqr(Data[i].Value[j] - Expected) / Expected;
    end;

  // Compute degrees of freedom.
  df := (m - 1) * (n - 1);
  if df <= 0 then begin
    AddToErrorStack('Error on chi square: invalid degrees of freedom.');
    Chi2 := NaN;
    PValue := NaN;
    Exit;
  end;

  PValue := 1 - IncompleteGamma(df / 2.0, chi2 / 2.0);
  CleanUpZero(Chi2);
  CleanUpOne(PValue);
end;

procedure ChiSquareAnalysis(const Data: WMatrixType);
var
  m, n: Integer;
  Chi2, PValue: FloatType;
begin
  // Initialize.
  m := Length(Data[0].Value);
  n := Length(Data);
  if (m < 2) or (n < 2) then begin
    AddToErrorStack('Error on chi-square: matrix must have at least 2 rows and 2 columns.');
    Exit;
  end;

  // Perform chi-square test.
  ChiSquareTest(Data, Chi2, PValue);
  if IsNan(Chi2) then begin
    AddToErrorStack('Error on chi-square: invalid input (negative values, zero totals, or small expected frequencies).');
    Exit;
  end
  else begin
    Write('Chi-square statistic: ', SmartFloat(Chi2));
    Write('  P-value: ', SmartFloat(PValue));
    WriteLn('  Degrees of freedom: ', (m - 1 ) * (n - 1));
  end;
end;

procedure NominalStatistics(const WVar: IVectorType);
var
  VData: WMatrixType;
  i, m, n: Integer;
  GKL, CramerV, Phi: FloatType;
begin
  writeln('Nominal Statistics');
  // Create VData from WData, using only the variables in WVar.
  SetLength(VData, Length(WVar));
  for i := 0 to High(WVar) do begin
    SetLength(VData[i].Value, Length(WData[WVar[i]].Value));
    DeepCopyVector(WData[WVar[i]], VData[i]);
  end;
  m := Length(VData[0].Value);
  n := Length(VData);

  // Analyses to perform.
  ChiSquareAnalysis(VData);
  GKL := Lambda(VData);
  writeln('Goodman & Kruskal''s Lambda: ', SmartFloat(GKL));
  CramerV := CramersV(VData);
  writeln('Cramer''s V: ', SmartFloat(CramerV));
  if (m = 2) and (n = 2) then begin
    Phi := PhiCoefficient(VData);
    writeln;
    writeln('Phi Coefficient: ', SmartFloat(Phi));
  end;

  DisplayErrorStack;
  if FirstPass then
    SavePartialData('Save nominal statistics to a file? (y/n) ', WVar, @NominalStatistics);
  Pause;
end;

end.

