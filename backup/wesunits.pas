unit WesUnits;

{$mode ObjFPC}{$H+}

interface

uses XAlglib, Sysutils;

procedure QuickSort(var A: TVector; L, R: Integer);
function HarmonicMean(const Data.Value: TVector): Double;
function GeometricMean(const Data.Value: TVector): Double;
function CoefficientOfVariation(const WMean, WSD: Double): Double;
function MedianAbsoluteDeviationFromMean(const Data.Value: TVector): Double;
function MeanAbsoluteDeviation(const Data.Value: TVector): Double;
function GiniCoefficient(const Data.Value: TVector): Double;
function HerfindahlIndex(const Shares: TVector): Double;

implementation

procedure QuickSort(var A: TVector; L, R: Integer);
var
  I, J: Integer;
  P, T: Double;
begin
  I := L;
  J := R;
  P := A[(L+R) div 2];
  repeat
    while A[I] < P do Inc(I);
    while A[J] > P do Dec(J);
    if I <= J then
    begin
      T := A[I];
      A[I] := A[J];
      A[J] := T;
      Inc(I);
      Dec(J);
    end;
  until I > J;
  if L < J then QuickSort(A, L, J);
  if I < R then QuickSort(A, I, R);
end;

function HarmonicMean(const Data.Value: TVector): Double;
var
  i: Integer;
  SumRecip: Double;
begin
  if Length(Data.Value) = 0 then
    raise Exception.Create('Empty Data.Value array.');

  SumRecip := 0;
  for i := 0 to High(Data.Value) do
  begin
    if Data.Value[i] <= 0 then
      raise Exception.Create('Term is undefined for non-positive values.');
    SumRecip := SumRecip + (1.0 / Data.Value[i]);
  end;

  Result := Length(Data.Value) / SumRecip;
end;

function GeometricMean(const Data.Value: TVector): Double;
var
  i: Integer;
  SumLog: Double;
begin
  if Length(Data.Value) = 0 then
    raise Exception.Create('Empty Data.Value array.');

  SumLog := 0;
  for i := 0 to High(Data.Value) do
  begin
    if Data.Value[i] <= 0 then
      raise Exception.Create('Geometric mean undefined for non-positive values.');
    SumLog := SumLog + Ln(Data.Value[i]);
  end;

  Result := Exp(SumLog / Length(Data.Value));
end;

function CoefficientOfVariation(const WMean, WSD: Double): Double;
begin
  if (WMean = 0) then writeln('Because mean is zero, term is undefined.')
  else Result := WSD / WMean;
end;

function MedianAbsoluteDeviationFromMean(const Data.Value: TVector): Double;
var
  i: Integer;
  MeanVal, MedianVal: Double;
  Deviations: TVector;
begin
  if Length(Data.Value) = 0 then
    raise Exception.Create('Empty Data.Value array.');

  MeanVal := SampleMean(Data.Value);
  SetLength(Deviations, Length(Data.Value));

  for i := 0 to High(Data.Value) do
    Deviations[i] := Abs(Data.Value[i] - MeanVal);
  samplemedian(Deviations, MedianVal);
  Result := MedianVal;
end;

function MeanAbsoluteDeviation(const Data.Value: TVector): Double;
var
  i: Integer;
  MeanVal, SumAbsDev: Double;
begin
  if Length(Data.Value) = 0 then
    raise Exception.Create('Empty Data.Value array.');

  MeanVal := 0;
  for i := 0 to High(Data.Value) do
    MeanVal := MeanVal + Data.Value[i];
  MeanVal := MeanVal / Length(Data.Value);

  SumAbsDev := 0;
  for i := 0 to High(Data.Value) do
    SumAbsDev := SumAbsDev + Abs(Data.Value[i] - MeanVal);

  Result := SumAbsDev / Length(Data.Value);
end;

function GiniCoefficient(const Data.Value: TVector): Double;
var
  Sorted: TVector;
  i, n: Integer;
  Sum, WeightedSum: Double;
begin
  n := Length(Data.Value);
  if n = 0 then
    raise Exception.Create('Empty Data.Value array.');
  // Copy and sort the Data.Value
  Sorted := Copy(Data.Value, 0, n);
  QuickSort(Sorted, 0, High(Sorted));
  Sum := 0;
  WeightedSum := 0;
  for i := 0 to n - 1 do
  begin
    Sum := Sum + Sorted[i];
    WeightedSum := WeightedSum + (i + 1) * Sorted[i];
  end;
  if Sum = 0 then
    raise Exception.Create('Gini undefined for zero total.');
  Result := (2 * WeightedSum) / (n * Sum) - (n + 1) / n;
end;

end.

