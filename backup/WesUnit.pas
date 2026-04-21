unit WesUnit;

{$mode ObjFPC}{$H+}{$I proprietary.txt}

interface

{ WesStat7 September 24, 2025, by Wesley R. Parsons, wespar@bellouth.net, www.wespar.com.}

uses
  BasicFuncs,
  Globals,
  Math;

{PARAMETRIC}
{Measures of Central Tendency}
function ArithmeticMean(const Col: WVectorType): FloatType;
function GeometricMean(const Col: WVectorType): FloatType;
function HarmonicMean(const Col: WVectorType): FloatType;
function MeanPairwiseDistance(const Col: WVectorType): FloatType;
function RootMeanSquare(const Col: WVectorType): FloatType;
function CubicMean(const Col: WVectorType): Real;

{Measures of Scale}
function Variance(const Col: WVectorType; const Sample: Boolean): FloatType;
function StandardDeviation(const Col: WVectorType; const Sample: Boolean): FloatType;
function Skewness(const Col: WVectorType; const Sample: Boolean): FloatType;
function YuleSkewness(const Col: WVectorType; const WQ1, WQ3: FloatType): FloatType;  //FixThis
function MedCouple(const Col: WVectorType): FloatType;
function Kurtosis(const Col: WVectorType; const Sample: Boolean): FloatType;
function ExcessKurtosis(const Col: WVectorType; const Sample: Boolean): FloatType;
function BimodalityCoefficient(const Col: WVectorType): FloatType;
function CoefficientOfVariation(const Col: WVectorType; const Sample: Boolean): FloatType;
function MeanAbsoluteDeviation(const Col: WVectorType): FloatType;
function StandardErrorOfMean(const Col: WVectorType; const Sample: Boolean): FloatType;
function SymmetryIndex(const Col: WVectorType): FloatType;
function HillEstimator(const Col: WVectorType; const k: Integer): FloatType;

{Measures of Equality}
function GiniCoefficient(const Col: WVectorType): FloatType;
function LorenzMuenzner(Col: WVectorType): FloatType;
function TheilIndex(const Col: WVectorType): FloatType;
function ShannonEntropy(const Col: WVectorType; BinType: EntropyType): FloatType;

{NON-PARAMETRIC STATISTICS}
{Measures of Central Tendency}
function Median(const Col: WVectorType): FloatType;
function MedianPairwiseDistance(const Col: WVectorType): FloatType;
function HighMedian(const Col: WVectorType): FloatType;
function LowMedian(const Col: WVectorType): FloatType;
function Trimedian(const Col: WVectorType): FloatType;
function Mode(const Col: WVectorType): FloatType;
function HuberEstimator(const Col: WVectorType): FloatType;
function HodgesLehmann(const Col: WVectorType): FloatType;
function TukeyBiweight(const Col: WVectorType): FloatType;
function SnEstimator(const Col: WVectorType): FloatType;
function QnEstimator(const Col: WVectorType): FloatType;

{Measures of Range}
function MinmValue(const Col: WVectorType): FloatType;
function MaxmValue(const Col: WVectorType): FloatType;
function Range(const Col: WVectorType): FloatType;
function Midrange(const Col: WVectorType): FloatType;
function Percentile(const D: WVectorType; P: FloatType): FloatType;
procedure Quartiles(const Col: WVectorType; out Q1, Q2, Q3: FloatType);
//procedure Quartiles(const A: WVectorType; out Q1, Q2, Q3: FloatType);
function InterQuartileRange(const Col: WVectorType): FloatType;
function SemiInterquartileRange(const Q1, Q3: FloatType): FloatType;
function Midhinge(const Q1, Q3: FloatType): FloatType;
function Q3Q1Ratio(const Q1, Q3: FloatType): FloatType;
function P90P10Ratio(const X: WVectorType): FloatType;
function P95P5Ratio(const X: WVectorType): FloatType;
function P75P25Spread(const X: WVectorType): FloatType;
function P95P5Spread(const X: WVectorType): FloatType;
function RangeRatio(const X: WVectorType): FloatType;
function QuartileCoefficientOfDispersion(const Q1, Q3: FloatType): FloatType;

{Measures of Scale}
function MeanAbsoluteDeviationFromMean(const Col: WVectorType): FloatType;
function MedianAbsoluteDeviationFromMean(const Col: WVectorType): FloatType;
function MeanAbsoluteDeviationFromMedian(const Col: WVectorType): FloatType;
function MedianAbsoluteDeviationFromMedian(const Col: WVectorType): FloatType;
function StandardErrorOfMedian(const Col: WVectorType): FloatType;
function SpearmanRankSkewness(const Col: WVectorType): FloatType;
function BowleySkewness(const Q1, Q2, Q3: FloatType): FloatType;
function KellySkewness(const Col: WVectorType): FloatType;
procedure LMoments(const Col: WVectorType; out L1, L2, L3, L4: FloatType);

{ Summing Routines }
function SumOfSquares(const Col: WVectorType): FloatType; overload;
function SumOfSquares(const Col: CVectorType): FloatType; overload;

{ADDITIONAL}
function GeneralizedMean(const Col: WVectorType; const p: FloatType): FloatType;
function HarrellDavisQuantile(const Col: WVectorType; p: FloatType): FloatType;
function Quantile(const Col: WVectorType; q: FloatType): FloatType;
function PearsonKthMoment(const Col: WVectorType; const k: Integer): FloatType;
function DurbinWatson(const R: WVectorType): FloatType;
procedure PercentOutliers(const Col: WVectorType; const Cut: FloatType; out Bot, Top: FloatType);
procedure PercentOutliersMAD(const Col: WVectorType; const Threshold: FloatType; out Bot, Top: FloatType);

implementation

{PARAMETRIC}
{Measures of Central Tendency}

function ArithmeticMean(const Col: WVectorType): FloatType;
var
  j, m: Integer;
  Sum: FloatType;
begin
  m := Length(Col.Value);
  if m = 0 then begin
    AddToErrorStack('Error on Mean: empty vector.');
    Exit(NaN);
  end;
  Sum := 0;
  for j := 0 to m do
    Sum := Sum + Col.Value[j];
  Result := Sum / m;
  CleanUpZero(Result);
end;

function GeometricMean(const Col: WVectorType): FloatType;
var
  j, m: Integer;
  SumLog: FloatType;
begin
  m := Length(Col.Value);
  SumLog := 0.0;
  for j := 0 to m - 1 do begin
    if Col.Value[j] <= 0 then begin
      AddToErrorStack('Error on Geometric Mean: undefined for non-positive values.');
      Exit(NaN);
    end;
    SumLog := SumLog + Ln(Col.Value[j]);
  end;
  Result := Exp(SumLog / m);
  CleanUpZero(Result);
end;

function HarmonicMean(const Col: WVectorType): FloatType;
var
  j, m: Integer;
  SumRecip: FloatType;
begin
  m := Length(Col.Value);
  SumRecip := 0.0;
  for j := 0 to m - 1 do begin
    if Col.Value[j] <= 0 then begin
      AddToErrorStack('Error on Harmonic Mean: undefined for non-positive values.');
      Exit(NaN);
    end;
    SumRecip := SumRecip + (1.0 / Col.Value[j]);
  end;
  Result := m / SumRecip;
  CleanUpZero(Result);
end;

function MeanPairwiseDistance(const Col: WVectorType): FloatType;
var
  i, j, m: Integer;
  PairSum: FloatType;
begin
  m := length(Col.Value);
  if m < 2 then begin
    AddToErrorStack('Error on Mean Pairwise Distance: length of vector < 2.');
    Exit(NaN);
  end;
  PairSum := 0.0;
  for i := 0 to m - 2 do
    for j := i + 1 to m - 1 do
      PairSum := PairSum + Abs(Col.Value[i] - Col.Value[j]);
  Result := (2 * PairSum) / (m * (m - 1));  // Normalize by number of pairs
end;

function RootMeanSquare(const Col: WVectorType): FloatType;
var
  j, m: Integer;
  Sum: FloatType;
begin
  m := Length(Col.Value);
  if m = 0 then begin
    AddToErrorStack('Error on Root Mean Square: empty vector.');
    Exit(NaN);
  end;

  Sum := 0;
  for j := 0 to m - 1 do
    Sum := Sum + Col.Value[j] * Col.Value[j];
  if m > 0 then
    Result := Sqrt(Sum / m)
  else
    Result := 0;
  CleanUpZero(Result);
end;

function CubicMean(const Col: WVectorType): Real;
var
  j, m: Integer;
  SumCube: FloatType;
begin
  m := Length(Col.Value);
  if m = 0 then begin
    AddToErrorStack('Error on Cubic Mean: empty vector.');
    Exit(NaN);
  end;

  SumCube := 0;
  for j := 0 to m - 1 do
    SumCube := SumCube + Power(Col.Value[j], 3);
  if m > 0 then
    Result := Sign(SumCube) * Power(Abs(SumCube / m), 1.0 / 3.0)
  else
    Result := 0.0;
  CleanUpZero(Result);
end;

{Measures of Scale}
function Variance(const Col: WVectorType; const Sample: Boolean): FloatType;
var
  j, m: Integer;
  mu, sumSq: FloatType;
begin
  m := Length(Col.Value);
  if m < 2 then begin
    AddToErrorStack('Error on Variance: length of vector < 2.');
    Exit(NaN);
  end;

  mu := ArithmeticMean(Col);
  sumSq := 0.0;

  for j := 0 to m - 1 do
    sumSq := sumSq + Sqr(Col.Value[j] - mu);
  if Sample then
    Result := sumSq / (m - 1)  //sample
  else
    Result := sumSq / m;  //population
  CleanUpBoth(Result);
end;

function StandardDeviation(const Col: WVectorType; const Sample: Boolean): FloatType;
var
  v: FloatType;
  m: Integer;
begin
  m := Length(Col.Value);
  if m < 2 then begin
    AddToErrorStack('Error on Standard Deviation: length of vector < 2.');
    Exit(NaN);
  end;

  v := Variance(Col, Sample);
  if v < 0.0 then begin
    AddToErrorStack('Error on Standard Deviation: negative variance.');
    Exit(NaN);
  end;
  Result := Sqrt(v);
  CleanUpBoth(Result);
end;

function Skewness(const Col: WVectorType; const Sample: Boolean): FloatType;
var
  mu, sigma, SumSkew, Correction: FloatType;
  j, m: Integer;
begin
  m := Length(Col.Value);
  if m < 3 then begin
    AddToErrorStack('Error on Skewness: length of vector < 3.');
    Exit(NaN);
  end;

  mu := ArithmeticMean(Col);
  sigma := StandardDeviation(Col, Sample);
  if sigma = 0.0 then begin
    AddToErrorStack('Error on Skewness: standard deviation is zero.');
    Exit(NaN);
  end;

  SumSkew := 0.0;
  for j := 0 to m - 1 do
    SumSkew := SumSkew + Power(Col.Value[j] - mu, 3);
  SumSkew := SumSkew / m / Power(sigma, 3);

  if Sample then
    Correction := Sqrt((m * (m - 1)) / (m - 2))
  else
    Correction := 1.0;
  Result := SumSkew * Correction;
  CleanUpZero(Result);
end;

function YuleSkewness(const Col: WVectorType; const WQ1, WQ3: FloatType): FloatType;
var
  med: FloatType;
  SortCol: WVectorType;
begin
  if Length(Col.Value) = 0 then begin
    AddToErrorStack('Error on Yule Skewness: empty vector.');
    Exit(NaN);
  end;

  SetLength(SortCol.Value, Length(Col.Value));
  DeepCopyVector(Col, SortCol);
  med := Median(SortCol);
  if (WQ3 = WQ1) then begin
    AddToErrorStack('Error on Yule Skewness: Q1 = Q3.');
    Exit(NaN);
  end
  else
    Result := (WQ3 + WQ1 - 2 * med) / (WQ3 - WQ1);
end;

function MaxAbs(const Col: WVectorType): FloatType;
var
  x: FloatType;
  j: Integer;
begin
  if Length(Col.Value) = 0 then begin
    AddToErrorStack('Error on Maximum: empty vector.');
    Exit(NaN);
  end;
Result := Abs(Col.Value[0]);
  for j := 1 to High(Col.Value) do begin
    x := Abs(Col.Value[j]);
    if x > Result then Result := x;
  end;
end;

{ From Copilot 11/28/2025, new improved MedianofKernel routine }
function MedianOfKernel(const Dplus, Dminus: CVectorType): FloatType;
var
  lowp, highp, mid: FloatType;
  i, j, count, total: Integer;
begin
  // Bounds of kernel values
  lowp := -1.0;
  highp := 1.0;
  total := Length(Dplus) * Length(Dminus);

  // Binary search for median
  while (highp - lowp > HighTolerance) do begin
    mid := 0.5 * (lowp + highp);
    count := 0;

    // Count kernel values <= mid without building them
    for i := 0 to High(Dplus) do
      for j := 0 to High(Dminus) do
        if (Dplus[i] - Dminus[j]) / (Dplus[i] + Dminus[j]) <= mid then
          Inc(count);

    if count >= (total div 2) then
      highp := mid
    else
      lowp := mid;
  end;

  Result := 0.5 * (lowp + highp);
end;

{ This routine came from Copilot 11/28/2025, new mediankernel routine }
function Medcouple(const Col: WVectorType): FloatType;
var
  Data: WVectorType;
  Dplus, Dminus: CVectorType;
  i, n, k: Integer;
  Xm: FloatType;

begin
  n := Length(Col.Value);
  if n < 3 then begin
    AddToErrorStack('Error on Medcouple: length of vector < 3.');
    Exit(NaN);
  end;

  // Sort data
  DeepCopyVector(Col, Data);
  QuickSort(Data, 0, High(Data.Value));

  // Median
  Xm := MedianOfSorted(Data.Value);

  // Split into Dminus and Dplus
  k := 0;
  while (k < n) and (Data.Value[k] < Xm) do
    Inc(k);

  SetLength(Dminus, k);
  SetLength(Dplus, n - k);

  for i := 0 to k - 1 do
    Dminus[k - 1 - i] := Xm - Data.Value[i];   // descending
  for i := 0 to (n - k) - 1 do
    Dplus[i] := Data.Value[k + i] - Xm;      // ascending

  // Optimized median search
  Result := MedianOfKernel(Dplus, Dminus);
end;

{ This version of Medcouple comes from ChatGPT. 11/10/2025. }
{function MedCouple(const Col: WVectorType): FloatType;
var
  Data, Dplus, Dminus, H: WVectorType;
  n, i, j, k, idxH: Integer;
  Xm, A, B: FloatType;

  function MedianOfSorted(const Arr: WVectorType): FloatType;
  var
    len: Integer;
  begin
    len := Length(Arr.Value);
    if len = 0 then Exit(NaN);
    if Odd(len) then
      Result := Arr.Value[len div 2]
    else
      Result := 0.5 * (Arr.Value[(len div 2) - 1] + Arr.Value[len div 2]);
  end;

begin
  n := Length(Col.Value);
  if n < 3 then Exit(NaN);

  // Copy and sort input
  DeepCopyVector(Col, Data);
  QuickSort(Data, 0, High(Data.Value));

  // Median
  Xm := MedianOfSorted(Data);

  // Split data into below and above median
  k := 0;
  while (k < n) and (Data.Value[k] < Xm) do Inc(k);

  SetLength(Dminus.Value, k);
  SetLength(Dplus.Value, n - k);

  // Fill centered arrays
  for i := 0 to k - 1 do
    Dminus.Value[k - 1 - i] := Xm - Data.Value[i];   // descending order
  for i := 0 to (n - k) - 1 do
    Dplus.Value[i] := Data.Value[k + i] - Xm;        // ascending order

  // Compute kernel matrix
  SetLength(H.Value, Length(Dplus.Value) * Length(Dminus.Value));
  idxH := 0;

  for i := 0 to High(Dplus.Value) do
    for j := 0 to High(Dminus.Value) do begin
      A := Dplus.Value[i];
      B := Dminus.Value[j];

      if (Abs(A) < Epsilon15) and (Abs(B) < Epsilon15) then
        H.Value[idxH] := 0.0
      else
        H.Value[idxH] := (A - B) / (A + B);

      Inc(idxH);
    end;

  // Sort H and take median
  QuickSort(H, 0, High(H.Value));
  Result := MedianOfSorted(H);
end;  }

function Kurtosis(const Col: WVectorType; const Sample: Boolean): FloatType;
var
  mu, sigma, SumKurt, Correction: FloatType;
  j, m: Integer;
begin
  m := Length(Col.Value);
  if m < 4 then begin
    AddToErrorStack('Error on Kurtosis: vector length < 4.');
    Exit(NaN);
  end;
  mu := ArithmeticMean(Col);
  sigma := StandardDeviation(Col, Sample);
  if sigma = 0.0 then begin
    AddToErrorStack('Error on Kurtosis: standard deviation is zero.');
    Exit(0.0);
  end;
  SumKurt := 0.0;
  for j := 0 to m - 1 do
    SumKurt := SumKurt + Power((Col.Value[j] - mu), 4);
  SumKurt := SumKurt / (m * Power(sigma, 4));
  if Sample then begin
    Correction := (m * (m + 1)) / ((m - 1) * (m - 2) * (m - 3));
    Result := Correction * SumKurt;
  end else
    Result := SumKurt;
  CleanUpZero(Result);
end;

// Bias Corrected version.
function ExcessKurtosis(const Col: WVectorType; const Sample: Boolean): FloatType;
var
  m: Integer;
begin
  m := Length(Col.Value);
  if m < 4 then begin
    AddToErrorStack('Error on Excess Kurtosis: Vector length < 4.');
    Exit(NaN);
  end;

  Result := Kurtosis(Col, Sample) - 3.0;
  CleanUpZero(Result);
end;

function BimodalityCoefficient(const Col: WVectorType): FloatType;
var
  g1, g2, m: FloatType;
begin
  m := Length(Col.Value);
  if m <= 3 then begin
    AddToErrorStack('Error on Bimodality Coefficient: vector length < 3.');
    Exit(NaN);
  end;

  g1 := Skewness(Col, True);
  g2 := ExcessKurtosis(Col, True);
  Result := (Sqr(g1) + 1.0) / (g2 + 3.0 + (3.0 * Sqr(m - 1)) / ((m - 2) * (m - 3)));
end;

function CoefficientOfVariation(const Col: WVectorType; const Sample: Boolean): FloatType;
var
  mu: FloatType;
begin
  if Length(Col.Value) = 0 then begin
    AddToErrorStack('Error on Mean Absolute Deviation: empty vector.');
    Exit(NaN);
  end;
  mu := ArithmeticMean(Col);

  if (mu = 0) then begin
    AddToErrorStack('Error on Coefficient of Variation: mean is zero.');
    Exit(NaN);
  end
  else
    Result := StandardDeviation(Col, Sample) / mu;
  CleanUpBoth(Result);
end;

function MeanAbsoluteDeviation(const Col: WVectorType): FloatType;
var
  j, m: Integer;
  mu, AbsSum: FloatType;
begin
  m := Length(Col.Value);
  if m = 0 then begin
    AddToErrorStack('Error on Mean Absolute Deviation: empty vector.');
    Exit(NaN);
  end;

  mu := ArithmeticMean(Col);
  AbsSum := 0.0;
  for j := 0 to m - 1 do
    AbsSum := AbsSum + Abs(Col.Value[j] - mu);
  Result := AbsSum / m;
  CleanUpZero(Result);
end;

function StandardErrorOfMean(const Col: WVectorType; const Sample: Boolean): FloatType;
var
  m: Integer;
  sd: FloatType;
begin
  m := Length(Col.Value);
  if m = 0 then begin
    AddToErrorStack('Error on Standard Error of Mean: empty vector.');
    Exit(NaN);
  end;
  if Sample and (m < 2) then begin
    AddToErrorStack('Error on Standard Error of Mean: sample vector length < 2.');
    Exit(NaN);
  end;

  sd := StandardDeviation(Col, Sample);
  Result := sd / Sqrt(m);
  CleanUpBoth(Result);
end;

function SymmetryIndex(const Col: WVectorType): FloatType;
var
  MeanVal, MedianVal, ModeVal, Denom: FloatType;
  SortCol: WVectorType;
begin
  if length(Col.Value) < 2 then begin
    AddToErrorStack('SymmetryIndex error: insufficient data.');
    Exit(NaN);
  end;

  SetLength(SortCol.Value, Length(Col.Value));
  DeepCopyVector(Col, SortCol);

  MeanVal := ArithmeticMean(Col);
  MedianVal := Median(SortCol);
  ModeVal := Mode(Col);

  Denom := MeanVal - MedianVal;
  if Abs(Denom) < 1e-12 then begin
    AddToErrorStack('SymmetryIndex error: mean equals median, division by near-zero.');
    Exit(NaN);
  end;

  Result := (MeanVal - ModeVal) / Denom;
end;

function HillEstimator(const Col: WVectorType; const k: Integer): FloatType;
var
  i, j, m: Integer;
  Sorted: WVectorType;
  LogSum: FloatType;
begin
  m := length(Col.Value);
  if (m < 2) or (k <= 1) or (k >= m) then begin
    AddToErrorStack('Hill Estimator error: invalid sample size or k.');
    Exit(NaN);
  end;

  SetLength(Sorted.Value, m);
  for i := 0 to m - 1 do
    Sorted.Value[i] := Col.Value[i];

  QuickSort(Sorted, 0, m);  // Assumes descending sort routine
  //So reverse sort
  j := High(Sorted.Value);
  for i := 0 to (j div 2) do
    Swap(Sorted.Value[i], Sorted.Value[j - i]);

  // Check for negative values
  LogSum := 0;
  for i := 0 to k - 1 do begin
    if Sorted.Value[i] <= 0 then begin
      AddToErrorStack('Hill Estimator error: non-positive value encountered.');
      Exit(NaN);
    end;
    LogSum := LogSum + Ln(Sorted.Value[i]);
  end;

  if Sorted.Value[k] <= 0 then begin
    AddToErrorStack('Hill Estimator error: kth value not positive.');
    Exit(NaN);
  end;

  LogSum := 0.0;
  for i := 0 to k - 1 do
    LogSum := LogSum + Ln(Sorted.Value[i]);

  Result := (LogSum / k) - Ln(Sorted.Value[k]);
end;

{Measures of Equality}
function GiniCoefficient(const Col: WVectorType): FloatType;
var
  CopiedCol: WVectorType;
  j, m: Integer;
  SumCol, WeightedSum: FloatType;
begin
  DeepCopyVector(Col, CopiedCol);
  m := Length(CopiedCol.Value);
  if m < 2 then begin
    AddToErrorStack('Error on Gini Coefficient: vector length < 2.');
    Exit(0.0);
  end;
  QuickSort(CopiedCol, 0, m - 1);
  SumCol := 0.0;
  WeightedSum := 0.0;
  for j := 0 to m - 1 do begin
    SumCol := SumCol + CopiedCol.Value[j];
    WeightedSum := WeightedSum + (j + 1) * CopiedCol.Value[j];  // Correct rank weight
  end;
  if SumCol = 0.0 then begin
    AddToErrorStack('Error for Gini Coefficient: undefined for zero total.');
    Exit(NaN);
  end;
  Result := (2 * WeightedSum) / (m * SumCol) - (m + 1) / m;
  CleanUpZero(Result);
end;

function LorenzMuenzner(Col: WVectorType): FloatType;
var
  m, j: Integer;
  Total, CumuSum, CumuShare, SumY: FloatType;
  xVector, yVector: CVectorType;
begin
  m := Length(Col.Value);
  if m = 0 then begin
    AddToErrorStack('Error on Lorenz-Muenzner: empty vector.');
    Exit(NaN);
  end;

  SetLength(xVector, m);
  SetLength(yVector, m);
  QuickSort(Col, 0, High(Col.Value));

  Total := 0;
  for j := 0 to m do Total := Total + Col.Value[j];
  if Total = 0.0 then begin
    AddToErrorStack('Error Lorenz-Muenzner: total sum is zero.');
    Exit(NaN);
  end;

  CumuSum := 0.0;
  SumY := 0.0;
  for j := 0 to m - 1 do begin
    CumuSum := CumuSum + Col.Value[j];
    CumuShare := CumuSum / Total;
    xVector[j] := j / m;
    yVector[j] := CumuShare;
    SumY := SumY + CumuShare;
  end;
  Result := (Succ(m) - 2 * SumY) / Pred(m);
end;

function TheilIndex(const Col: WVectorType): FloatType;
var
  j, m: Integer;
  mu, SumCol: FloatType;
begin
  m := Length(Col.Value);
  mu := ArithmeticMean(Col);

  if (m = 0) then begin
    AddToErrorStack('Error on Theil Index: empty vector.');
    Exit(NaN);
  end;

  // Compute Theil index
  SumCol := 0;
  for j := 0 to m - 1 do begin
    if Col.Value[j] <= 0 then begin
      AddToErrorStack('Error on Theil Index: all data must be positive.');
      Exit(NaN);
    end
    else
      SumCol := SumCol + (Col.Value[j] / mu) * Ln(Col.Value[j] / mu);
  end;

  Result := SumCol / m;
  CleanUpBoth(Result);
end;

function ShannonEntropy(const Col: WVectorType; EntropyAlso: EntropyType): FloatType;
var
  MinVal, MaxVal, BinWidth: FloatType;
  j, k, BinIndex, m, BinCount: Integer;
  Bins: IVectorType;
  p, Entropy, LocalTolerance: FloatType;
begin
  m := Length(Col.Value);
  if (m = 0) or (BinCount <= 0) then begin
    AddToErrorStack('Error on Shannon Entropy: empty vector or bin count.');
    Exit(NaN)
  end
  else if (m < 3) then begin
    AddToErrorStack('Error on Shannon Entropy: length of vector < 3.');
    Exit(NaN);
  end;
  Case EntropyAlgo of
    SqRoot: BinCount := Ceil(Sqrt(m));
    Scott: begin
      k := Ceil((MaxmValue(Col) - MinmValue(Col)) / 2.0);
      BinCount :=  Ceil(3.5 * StandardDeviation(Col, Sample) / Power(m, 1.0 / 3.0));
    end;
    Sturges: begin
      BinCount := Ceil(Log2(m) + 1.0);
      if BinCount < 1 then BinCount := 1;
    end;
    Rice: begin
      BinCount := Ceil(2.0 * power(m, 1.0 / 3.0));
    end;
    FD: begin
      BinWidth := Ceil(2.0 * InterQuartileRange(Col) / Power(m, 1.0 / 3.0));
    end;
  end;
  // Step 1: Find min and max.
  MinVal := Col.Value[0];
  MaxVal := Col.Value[0];
  for j := 1 to m - 1 do begin
    if Col.Value[j] < MinVal then MinVal := Col.Value[j];
    if Col.Value[j] > MaxVal then MaxVal := Col.Value[j];
  end;

  // Step 2: Initialize bins.
  if EntropyAlgo = FD then
    BinCount := Ceil((MaxVal - MinVal) / BinWidth)
  else
    BinWidth := (MaxVal - MinVal) / BinCount;
  if BinWidth = 0.0 then begin
    AddToErrorStack('Error on Shannon Entropy: all values are identical.');
    Exit(NaN);
  end;
  if BinCount = 0 then begin
    AddToErrorStack('Error on Shannon Entropy:: bin count is zero.');
    Exit(NaN);
  end;
  SetLength(Bins, BinCount);

  // Step 3: Count values per bin.
  for j := 0 to m - 1 do begin
    BinIndex := Trunc((Col.Value[j] - MinVal) / BinWidth);
    if BinIndex >= BinCount then
      BinIndex := binCount - 1; // Edge case.
    Inc(Bins[BinIndex]);
  end;

  // Step 4: Compute entropy.
  Entropy := 0;
  LocalTolerance := IfThen(UserTolerance > 0, UserTolerance, HighTolerance);
  for j := 0 to BinCount - 1 do
    if Bins[j] > 0 then begin
      p := Bins[j] / m;
      if p = 0.0 then
        p := LocalTolerance;
      Entropy := Entropy - p * Ln(p);  // Log base 2.
  end;

  Result := Entropy;
  CleanUpBoth(Result);
end;

{NON-PARAMETRIC STATISTICS}

{Measures of Central Tendency}

// This Median function finds the median by using QuickSelect.
// More efficient.
function Median(const Col: WVectorType): FloatType;
var
  m, mid: Integer;
  QSCol: WVectorType;
begin
  m := Length(Col.Value);
  if m = 0 then begin
    AddToErrorStack('Error on Median: empty vector.');
    Exit(NaN);
  end;

  // Create vector to send to QuickSelect.
  SetLength(QSCol.Value, m);
  DeepCopyVector(Col, QSCol);
  mid := m div 2;
  if (m mod 2) = 1 then begin
    // Odd length: find middle element
    QuickSelect(QSCol, 0, m - 1, mid);
    Result := QSCol.Value[mid];
  end
  else begin
    // Even length: find two middle elements
    QuickSelect(QSCol, 0, m - 1, mid - 1);
    QuickSelect(QSCol, 0, m - 1, mid);
    Result := (QSCol.Value[mid - 1] + QSCol.Value[mid]) / 2.0;
  end;

end;

// This Median function finds the median by sorting the vector.
// Not efficient.
function MedianBySorting(const Col: WVectorType): FloatType;
var
  CopiedCol: WVectorType;
  m: Integer;
begin
  m := Length(Col.Value);
  if m = 0 then begin
    AddToErrorStack('Error on Median: empty vector.');
    Exit(NaN);
  end
  else begin
    DeepCopyVector(Col, CopiedCol);
    QuickSort(CopiedCol, 0, m - 1);
    if (m mod 2) = 1 then
      Result := CopiedCol.Value[(m div 2)]
    else
      Result := (CopiedCol.Value[m div 2 - 1] + CopiedCol.Value[m div 2]) / 2;
  end;
  CleanUpZero(Result);
end;

function MedianPairwiseDistance(const Col: WVectorType): FloatType;
var
  i, j, m, k: Integer;
  Distances: WVectorType;
begin
  m := length(Col.Value);
  if m < 2 then begin
    AddToErrorStack('Error on Median Pairwise Distance: length of vector < 2.');
    Exit(NaN);
  end;

  SetLength(Distances.Value, (m * (m - 1)) div 2);
  k := 0;
  for i := 0 to m - 2 do
    for j := i + 1 to m - 1 do begin
      Distances.Value[k] := Abs(Col.Value[i] - Col.Value[j]);
      Inc(k);
    end;

  QuickSort(Distances, 0, High(Distances.Value));
  k := Length(Distances.Value);
  if (k mod 2) = 1 then
    Result := Distances.Value[k div 2]
  else
    Result := 0.5 * (Distances.Value[(k div 2) - 1] + Distances.Value[k div 2]);
end;

function HighMedian(const Col: WVectorType): FloatType;
var
  m, mid: Integer;
  CopiedCol: WVectorType;
begin
  DeepCopyVector(Col, CopiedCol);
  m := Length(CopiedCol.Value);
  if m = 0 then begin
    AddToErrorStack('Error on High Median: empty vector.');
    Exit(NaN);
  end;

  mid := m div 2;

  if Odd(m) then begin
    // Odd length: standard median
    QuickSelect(CopiedCol, 0, m - 1, mid);
    Result := CopiedCol.Value[mid];
  end
  else begin
    // Even length: take the higher of the two middle values
    QuickSelect(CopiedCol, 0, m - 1, mid);
    Result := CopiedCol.Value[mid];
  end;
end;

{function HighMedian(const Col: WVectorType): FloatType;
var
  CopiedCol: WVectorType;
  m, mid: Integer;
begin
  DeepCopyVector(Col, CopiedCol);
  m := Length(CopiedCol.Value);
  if m = 0 then begin
    AddToErrorStack('Error on High Median: empty vector.');
    Exit(NaN);
  end;
  QuickSort(CopiedCol, 0, m - 1);
  mid := m div 2;
  Result := CopiedCol.Value[mid]  // High median for even N
end;}

function LowMedian(const Col: WVectorType): FloatType;
var
  m, mid: Integer;
  CopiedCol: WVectorType;
begin
  DeepCopyVector(Col, CopiedCol);
  m := Length(CopiedCol.Value);
  if m = 0 then begin
    AddToErrorStack('Error on High Median: empty vector.');
    Exit(NaN);
  end;

  mid := m div 2;

  if Odd(m) then begin
    QuickSelect(CopiedCol, 0, m - 1, mid);
    Result := CopiedCol.Value[mid];
  end
  else begin
    // Even length: take the lower of the two middle values
    QuickSelect(CopiedCol, 0, m - 1, mid-1);
    Result := CopiedCol.Value[mid-1]; // lower middle element
  end;
end;

{function LowMedian(const Col: WVectorType): FloatType;
var
  CopiedCol: WVectorType;
  m, mid: Integer;
begin
  DeepCopyVector(Col, CopiedCol);
  m := Length(CopiedCol.Value);
  if m = 0 then begin
    AddToErrorStack('Error on Low Median: empty vector.');
    Exit(NaN);
  end;
  QuickSort(CopiedCol, 0, m - 1);
  mid := m div 2;
  if (m mod 2 = 0) then
    Result := CopiedCol.Value[mid - 1]  // lower median for even
  else
    Result := CopiedCol.Value[mid];     // middle value for odd
end;}

function Trimedian(const Col: WVectorType): FloatType;
var
  Q1, Q2, Q3: FloatType;
begin
  if Length(Col.Value) = 0 then begin
    AddToErrorStack('Error on Trimedian: empty vector.');
    Exit(NaN);
  end;
  Quartiles(Col, Q1, Q2, Q3);
  Result := (Q1 + 2 * Q2 + Q3) / 4;
  CleanUpZero(Result);
end;

function Mode(const Col: WVectorType): FloatType;
var
  i, j, Count, MaxCount, m: Integer;
begin
  m := Length(Col.Value);
  if m = 0 then begin
    AddToErrorStack('Error on Mode: empty vector.');
    Exit(NaN);
  end;

  MaxCount := 0;
  Result := Col.Value[0];
  for j := 0 to m - 1 do begin
    Count := 0;
    for i := 0 to m - 1 do
      if Col.Value[i] = Col.Value[j] then Inc(Count);
    if Count > MaxCount then begin
      MaxCount := Count;
      Result := Col.Value[j];
    end;
  end;
end;

function HuberEstimator(const Col: WVectorType): FloatType;
const
  C = 1.345;
var
  Mu, OldMu, Wi, SumWeights, SumWeightedValues: FloatType;
  j, Iter: Integer;
  SortCol: WVectorType;
begin
  if Length(Col.Value) = 0 then begin
    AddToErrorStack('Error on Huber Estimator: empty vector.');
    Exit(NaN);
  end;

  SetLength(SortCol.Value, Length(Col.Value));
  DeepCopyVector(Col, SortCol);
  Mu := Median(SortCol);
  Iter := 0;

  repeat
    OldMu := Mu;
    SumWeights := 0.0;
    SumWeightedValues := 0.0;
    for j := 0 to High(Col.Value) do begin
      if Abs(Col.Value[j] - Mu) <= C then
        Wi := 1.0
      else
        Wi := C / Abs(Col.Value[j] - Mu);
      SumWeights := SumWeights + Wi;
      SumWeightedValues := SumWeightedValues + Wi * Col.Value[j];
    end;
    if SumWeights > 0.0 then
      Mu := SumWeightedValues / SumWeights
    else begin
      Result := OldMu;
      Exit;
    end;
    Inc(Iter);
  until (Abs(Mu - OldMu) < Tolerance) or (Iter >= SlowMaxIter);

  Result := Mu;
  CleanUpBoth(Result);
end;

function HodgesLehmann(const Col: WVectorType): FloatType;
var
  i, j, m, k, mWalsh: Integer;
  Walsh: WVectorType;
begin
  m := Length(Col.Value);
  if m = 0 then begin
    AddToErrorStack('Error on HodgesLehmann: empty vector.');
    Exit(NaN);
  end;

  SetLength(Walsh.Value, (m * (m + 1)) div 2);
  k := 0;
  for i := 0 to m - 1 do
    for j := i to m - 1 do begin
      Walsh.Value[k] := (Col.Value[i] + Col.Value[j]) / 2;
      Inc(k);
    end;

  QuickSort(Walsh, 0, High(Walsh.Value));
  mWalsh := Length(Walsh.Value);
  if mWalsh mod 2 = 1 then
    Result := Walsh.Value[mWalsh div 2]
  else
    Result := (Walsh.Value[mWalsh div 2 - 1] + Walsh.Value[mWalsh div 2]) / 2;
  CleanUpZero(Result);
end;

function TukeyBiweight(const Col: WVectorType): FloatType;
const
  c = 4.685;  // Tuning constant
var
  SortCol: WVectorType;
  m, j: Integer;
  num, denom, u, med, w: FloatType;
begin
  m := Length(Col.Value);
  if m = 0 then begin
    AddToErrorStack('Error on Tukey Biweight: empty vector.');
    Exit(NaN);
  end;

  SetLength(SortCol.Value, m);
  DeepCopyVector(Col, SortCol);

  med := Median(SortCol);
  num := 0.0;
  denom := 0.0;

  for j := 0 to m - 1 do begin
    u := Col.Value[j] - med;
    if Abs(u) < c then begin
      w := Sqr(1.0 - Sqr(u / c));
      num := num + Col.Value[j] * w;
      denom := denom + w;
    end;
  end;

  if denom = 0.0 then
    Result := med  // fallback if all weights are zero (all values are outliers)
  else
    Result := num / denom;
  CleanUpBoth(Result);
end;

{Measures of Range}
function MinmValue(const Col: WVectorType): FloatType;
var
  j: Integer;
  MinVal: FloatType;
begin
  if Length(Col.Value) = 0 then begin
    AddToErrorStack('Error on Minimum: empty vector.');
    Exit(NaN);
  end;

  MinVal := Col.Value[0];
  for j := 1 to High(Col.Value) do
    if Col.Value[j] < MinVal then MinVal := Col.Value[j];
  Result := MinVal;
  CleanUpZero(Result);
end;

function MaxmValue(const Col: WVectorType): FloatType;
var
  j: Integer;
  MaxVal: FloatType;
begin
  if Length(Col.Value) = 0 then begin
    AddToErrorStack('Error on Maximum: empty vector.');
    Exit(NaN);
  end;

  MaxVal := Col.Value[0];
  for j := 1 to High(Col.Value) do
    if Col.Value[j] > MaxVal then MaxVal := Col.Value[j];
  Result := MaxVal;
  CleanUpZero(Result);
end;

function Range(const Col: WVectorType): FloatType;
var
  i, m: Integer;
  MinVal, MaxVal: FloatType;
begin
  m := Length(Col.Value);
  if Length(Col.Value) = 0 then begin
    AddToErrorStack('Error on Range: empty vector.');
    Exit(0.0);
  end;

  MinVal := Col.Value[0];
  MaxVal := Col.Value[0];

  for i := 1 to m do begin
    if Col.Value[i] < MinVal then MinVal := Col.Value[i];
    if Col.Value[i] > MaxVal then MaxVal := Col.Value[i];
  end;

  Result := MaxVal - MinVal;
end;

function Midrange(const Col: WVectorType): FloatType;
begin
  Result := Abs(Range(Col)) / 2.0;
end;

function Percentile(const D: WVectorType; P: Float): FloatType;
var
  x: FloatType;
  k0, k1: Integer;
begin
  if Length(D.Value) = 0 then begin
    AddToErrorStack('Error on Percentile: undefined for non-positive values.');
    Exit(NaN);
  end;

  x := P * (Length(D.Value) - 1);
  k0 := Floor(x);
  k1 := Ceil(x);
  if k0 = k1 then
    Result := D.Value[k0]
  else
    Result := D.Value[k0] + (x - k0) * (D.Value[k1] - D.Value[k0]);
end;

procedure Quartiles(const Col: WVectorType; out Q1, Q2, Q3: FloatType);
var
  m, mid: Integer;
  a: WVectorType;
begin
  m := Length(Col.Value);
  if m = 0 then begin
    AddToErrorStack('Error on Quartiles: empty vector.');
    Q1 := NaN;
    Q2 := NaN;
    Q3 := NaN;
    Exit;
  end;

  // Create vector to send to QuickSelect.
  SetLength(A.Value, Length(Col.Value));
  DeepCopyVector(Col, A);

  // Median (Q2).
  mid := m div 2;
  QuickSelect(A, 0, m - 1, mid);
  if Odd(m) then
    Q2 := A.Value[mid]
  else begin
    QuickSelect(A, 0, m - 1, mid-1);
    Q2 := (A.Value[mid-1] + A.Value[mid]) / 2.0;
  end;

  // First quartile (Q1).
  QuickSelect(A, 0, m - 1, m div 4);
  Q1 := A.Value[m div 4];

  // Third quartile (Q3).
  QuickSelect(A, 0, m - 1, (3 * m) div 4);
  Q3 := A.Value[(3 * m) div 4];
end;

{procedure Quartiles(const A: WVectorType; out Q1, Q2, Q3: Float);
var
  Sorted: WVectorType;
  j, m: Integer;
begin
  m := Length(A.Value);
  if m = 0 then begin
    AddToErrorStack('Error on Quartiles: empty vector.');
    Q1 := NaN;
    Q2 := NaN;
    Q3 := NaN;
    Exit;
  end;

  SetLength(Sorted.Value, m);
  for j := 0 to m - 1 do
    Sorted.Value[j] := A.Value[j];
  QuickSort(Sorted, 0, High(Sorted.Value));
  Q1 := Percentile(Sorted, 0.25);
  Q2 := Percentile(Sorted, 0.50);
  Q3 := Percentile(Sorted, 0.75);
end;}

function InterquartileRange(const Col: WVectorType): FloatType;
var
  q1, q2, q3: FloatType;
begin
  if Length(Col.Value) = 0 then begin
    AddToErrorStack('Error on Interquartile Range: empty vector.');
    Exit(NaN);
  end;

  Quartiles(Col, q1, q2, q3);
  Result := q3 - q1;
  CleanUpZero(Result);
end;

function SemiInterquartileRange(const Q1, Q3: FloatType): FloatType;
begin
  Result := (Q3 - Q1) / 2;
  CleanUpZero(Result);
end;

function P90P10Ratio(const X: WVectorType): FloatType;
var
  p10, p90: Float;
begin
  if Length(X.Value) = 0 then begin
    AddToErrorStack('Error on Ratio: empty vector.');
    Exit(NaN);
  end;

  p10 := Percentile(X, 0.10);
  p90 := Percentile(X, 0.90);
  if IsEffectivelyEqual(p10, 0.0) or IsEffectivelyEqual(p90, 0.0) then
    Exit(0.0)
  else
    Result := Percentile(X, 0.90) / Percentile(X, 0.10);
end;

function P95P5Ratio(const X: WVectorType): FloatType;
var
  p05, p95: Float;
begin
  if Length(X.Value) = 0 then begin
    AddToErrorStack('Error on Ratio: empty vector.');
    Exit(NaN);
  end;

  p05 := Percentile(X, 0.05);
  p95 := Percentile(X, 0.95);
  if IsEffectivelyEqual(p05, 0.0) or IsEffectivelyEqual(p95, 0.0) then
    Exit(0.0)
  else
    Result := Percentile(X, 0.95) / Percentile(X, 0.05);
end;

function P75P25Spread(const X: WVectorType): FloatType;
begin
  Result := Percentile(X, 0.75) - Percentile(X, 0.25);
end;

function P95P5Spread(const X: WVectorType): FloatType;
begin
  Result := Percentile(X, 0.95) - Percentile(X, 0.05);
end;

function RangeRatio(const X: WVectorType): FloatType;
var
  p: FloatType;
begin
  if Length(X.Value) = 0 then begin
    AddToErrorStack('Error on Range Ratio: empty vector.');
    Exit(NaN);
  end;

  p := MinmValue(X);
  if IsEffectivelyEqual(p, 0.0) then
    Exit(0.0)
  else
    Result := MaxmValue(X) / p;
end;

function QuartileCoefficientOfDispersion(const Q1, Q3: FloatType): FloatType;
begin
  if (Q3 + Q1) <> 0.0 then
    Result := (Q3 - Q1) / (Q3 + Q1)
  else begin
    AddToErrorStack('Error on Quartile Coefficient of Dispersion: division by zero.');
    Exit(NaN);
  end;
  CleanUpBoth(Result);
end;

function Midhinge(const Q1, Q3: FloatType): FloatType;
begin
  Result := (Q1 + Q3) / 2;
  CleanUpZero(Result);
end;

function Q3Q1Ratio(const Q1, Q3: FloatType): FloatType;
begin
  if Q1 <> 0 then
    Result := Q3 / Q1
  else begin
    AddToErrorStack('Error on Q3/Q1 ratio: division by zero.');
    Exit(NaN);
  end;
  CleanUpBoth(Result);
end;

function SnEstimator(const Col: WVectorType): FloatType;
const
  C = 1.1926; // Consistency constant for normal distribution
var
  m, i, j, DiffCount: Integer;
  Sorted, Differences: WVectorType;
begin
  m := Length(Col.Value);
  if m < 2 then begin
    AddToErrorStack('Error on SnExtimator: vector length < 2.');
    Exit(NaN);
  end;

  // Sort the data
  SetLength(sorted.Value, m);
  DeepCopyVector(Col, Sorted);
  {for j := 0 to m - 1 do
    sorted.Value[j] := Col.Value[j];}
  QuickSort(Sorted, 0, m - 1);

  // Compute pairwise absolute differences (i < j)
  DiffCount := (m * (m - 1)) div 2;
  SetLength(Differences.Value, DiffCount);
  DiffCount := 0;
  for i := 0 to m - 2 do
    for j := i + 1 to m - 1 do begin
      Differences.Value[DiffCount] := Abs(Sorted.Value[j] - Sorted.Value[i]);
      Inc(DiffCount);
    end;

  // Compute median of differences
  Result := C * Median(Differences);
  CleanUpZero(Result);
end;

function QnEstimator(const Col: WVectorType): FloatType;
const
  C = 2.2219;  // Consistency constant
var
  m, i, j, k, h, DiffCount: Integer;
  Diffs: WVectorType;
  Temp: FloatType;
begin
  m := Length(Col.Value);
  if m < 2 then begin
    AddToErrorStack('Error on QnExtimator: vector length < 2.');
    Exit(NaN);
  end;

  // --- Step 1: Compute all pairwise absolute differences ---
  DiffCount := (m * (m - 1)) div 2;
  SetLength(Diffs.Value, DiffCount);

  k := 0;
  for i := 0 to m - 2 do
    for j := i + 1 to m - 1 do begin
      Diffs.Value[k] := Abs(Col.Value[i] - Col.Value[j]);
      Inc(k);
    end;

  // --- Step 2: Sort differences ---
  for i := 0 to DiffCount - 2 do
    for j := i + 1 to diffCount - 1 do
      if Diffs.Value[i] > Diffs.Value[j] then begin
        Temp := Diffs.Value[i];
        Diffs.Value[i] := Diffs.Value[j];
        Diffs.Value[j] := Temp;
      end;

  // --- Step 3: Compute h and k ---
  h := (m div 2) + 1;
  k := (h * (h - 1)) div 2;  // Index of the k-th smallest difference

  // --- Step 4: Qn = c * d_(k) ---
  if k >= diffCount then
    Result := C * Diffs.Value[DiffCount - 1]  // fallback
  else
    Result := C * Diffs.Value[k];
end;

{ Measures of Scale }
function MeanAbsoluteDeviationFromMean(const Col: WVectorType): FloatType;
var
  m, j: Integer;
  mu, SumAbsDev: FloatType;
begin
  m := Length(Col.Value);
  if m = 0 then begin
    AddToErrorStack('Error on Mean Absolute Deviation from Mean: empty vector.');
    Exit(NaN);
  end;

  // Compute mean
  mu := ArithmeticMean(Col);

  // Compute mean absolute deviation from mean
  SumAbsDev := 0.0;
  for j := 0 to m - 1 do
    SumAbsDev := SumAbsDev + Abs(Col.Value[j] - mu);

  Result := SumAbsDev / m;
  CleanUpZero(Result);
end;

function MedianAbsoluteDeviationFromMean(const Col: WVectorType): FloatType;
var
  j: Integer;
  Mu, MedianVal: FloatType;
  Deviations: WVectorType;
begin
  if Length(Col.Value) = 0 then begin
    AddToErrorStack('Error on Median Absolute Deviation from Mean: empty vector.');
    Exit(NaN);
  end;

  mu := ArithmeticMean(Col);
  SetLength(Deviations.Value, Length(Col.Value));
  for j := 0 to High(Deviations.Value) do
    Deviations.Value[j] := Abs(Col.Value[j] - mu);
  MedianVal := Median(Deviations);

  Result := MedianVal;
  CleanUpZero(Result);
end;

function MeanAbsoluteDeviationFromMedian(const Col: WVectorType): FloatType;
var
  j, m: Integer;
  mu, sum: FloatType;
  SortCol: WVectorType;
begin
  m := Length(Col.Value);
  if m = 0 then begin
    AddToErrorStack('Error on Mean Absolute Deviation fromMedian: empty vector.');
    Exit(NaN);
  end;

  SetLength(SortCol.Value, m);
  DeepCopyVector(Col, SortCol);
    mu := Median(SortCol);  // mu is the median
  sum := 0.0;
  for j := 0 to m - 1 do
    sum := sum + Abs(Col.Value[j] - mu);

    Result := sum / m;
  CleanUpZero(Result);
end;

function MedianAbsoluteDeviationFromMedian(const Col: WVectorType): FloatType;
var
  Med: FloatType;
  AbsDeviations, SortCol: WVectorType;
  j, m: Integer;
begin
  m := Length(Col.Value);
  if m = 0 then begin
    AddToErrorStack('Error on Median Absolute Deviations from Median: empty vector.');
    Exit(Nan);
  end;

  SetLength(SortCol.Value, m);
  DeepCopyVector(Col, SortCol);
  Med := Median(SortCol);

  SetLength(AbsDeviations.Value, m);
  for j := 0 to m - 1 do
    AbsDeviations.Value[j] := Abs(Col.Value[j] - Med);

  Result := Median(AbsDeviations);
  CleanUpZero(Result);
end;

function StandardErrorOfMedian(const Col: WVectorType): FloatType;
const
  EfficiencyFactor = SqrtPi2; // ≈ sqrt(π/2), for normal distribution
var
  m: Integer;
begin
  m := Length(Col.Value);
  if m < 2 then begin
    AddToErrorStack('Error on standard error of mean: vector length < 2.');
    Exit(NaN);
  end;
  Result := EfficiencyFactor * StandardDeviation(Col, True) / Sqrt(m);
  CleanUpZero(Result);
end;

function FindRank(const Sorted: WVectorType; const x: FloatType): FloatType;
var
  j, m, startIdx, endIdx: Integer;
begin
  m := Length(Sorted.Value);
  if m = 0 then begin
    AddToErrorStack('Error on FindRank: empty vector.');
    Exit(NaN);
  end;

  startIdx := -1;
  endIdx := -1;

  // Locate all positions where Sorted.Value[j] = x
  for j := 0 to m - 1 do begin
    if Sorted.Value[j] = x then begin
      if startIdx = -1 then startIdx := j;
      endIdx := j;
    end;
  end;

  if startIdx = -1 then begin
    AddToErrorStack('Error on FindRank: value not found in sorted vector.');
    Exit(NaN);
  end;

  // Average rank for ties: (startIdx + endIdx) / 2 + 1
  Result := (startIdx + endIdx) / 2 + 1;
end;

function SpearmanRankSkewness(const Col: WVectorType): FloatType;

  // Assign average ranks (1..n) with ties handled by average rank
  procedure AssignAverageRanks(const SortedVals: WVectorType; const OrigVals: WVectorType;
    out Ranks: array of FloatType);
  var
    n, i, j: Integer;
    startIdx, endIdx: Integer;
    sumRanks: Double;
    avgRank: Double;
  begin
    n := Length(SortedVals.Value);
    // map each value in OrigVals to its average rank from SortedVals
    // SortedVals must be sorted ascending

    for i := 0 to High(OrigVals.Value) do begin
      // find the block of equal values in SortedVals
      // binary search to locate first occurrence, then scan forward for ties
      // simple linear search from 0 is okay for moderate n; replace with binary search if desired
      startIdx := -1;
      for j := 0 to n - 1 do
        if SortedVals.Value[j] = OrigVals.Value[i] then begin
          startIdx := j;
          Break;
        end;
      if startIdx = -1 then begin
        // value not found (shouldn't happen with exact equality), set NaN rank
        Ranks[i] := NaN;
        Continue;
      end;
      endIdx := startIdx;

      while (endIdx + 1 < n) and (SortedVals.Value[endIdx + 1] = SortedVals.Value[startIdx]) do
        Inc(endIdx);
      // average rank for positions startIdx..endIdx (ranks are 1-based)
      sumRanks := 0.0;
      for j := startIdx to endIdx do
        sumRanks := sumRanks + (j + 1);
      avgRank := sumRanks / (endIdx - startIdx + 1);
      Ranks[i] := avgRank;
    end;
  end;

var
  m, j: Integer;
  Ranks: array of FloatType;
  Sorted: WVectorType;
  MeanRank, SumSq, SumCube: FloatType;
begin
  m := Length(Col.Value);
  if m < 3 then begin
    AddToErrorStack('Error on Spearman''s Rank Skewness: vector length < 3.');
    Exit(NaN);
  end;

  // Step 1: copy and sort the data (we need sorted values for rank assignment)
  DeepCopyVector(Col, Sorted);
  QuickSort(Sorted, 0, High(Sorted.Value)); // sort ascending

  // allocate rank array (use FloatType to allow fractional (average) ranks)
  SetLength(Ranks, m);

  // assign average ranks (tie-aware)
  AssignAverageRanks(Sorted, Col, Ranks);

  // verify no NaN ranks (defensive)
  for j := 0 to m - 1 do
    if IsNan(Ranks[j]) then begin
      AddToErrorStack('Error on Spearman''s Rank Skewness: rank assignment failed.');
      Exit(NaN);
    end;

  // Step 2: compute skewness from ranks
  MeanRank := (m + 1) / 2.0;
  SumSq := 0.0;
  SumCube := 0.0;
  for j := 0 to m - 1 do begin
    SumSq := SumSq + Sqr(Ranks[j] - MeanRank);
    SumCube := SumCube + Power(Ranks[j] - MeanRank, 3);
  end;

  if SumSq = 0.0 then begin
    AddToErrorStack('Error on Spearman''s Rank Skewness: zero variance in ranks.');
    Exit(NaN);
  end;

  // bias-corrected formula for sample skewness on ranks (Spearman version)
  Result := (m * (m * m - 1)) / ((m - 1) * (m - 2)) * (SumCube / Power(SumSq, 1.5));

  CleanUpZero(Result);
end;

function BowleySkewness(const Q1, Q2, Q3: FloatType): FloatType;
begin
  if (Q3 - Q1) = 0.0 then begin
    AddToErrorStack('Error on Bowley''s Skewness: division by zero.');
    Exit(NaN);
  end;

  Result := (Q3 + Q1 - 2 * Q2) / (Q3 - Q1);
  CleanUpZero(Result);
end;

function KellySkewness(const Col: WVectorType): FloatType;
var
  P10, P50, P90, Denominator: FloatType;
begin
  if Length(Col.Value) = 0 then begin
    AddToErrorStack('Error on Kelly''s Skewness: empty vector.');
    Exit(NaN);
  end;

  P10 := Percentile(Col, 0.10);
  P50 := Percentile(Col, 0.50);
  P90 := Percentile(Col, 0.90);
  Denominator := P90 - P10;

  if Denominator = 0.0 then begin
    AddToErrorStack('Error on Kelly''s Skewness: division by zero.');
    Exit(NaN);
  end;

  Result := (P90 + P10 - 2 * P50) / Denominator;
  CleanUpZero(Result);
end;

{CORRELATION}
{ Procedure to calculate probability-weighted moments (b_r) }
procedure CalculatePWMs(const Col: WVectorType; var b0, b1, b2, b3: FloatType);
var
  j, m: Integer;
  CopiedCol: WVectorType;
begin
  m := Length(Col.Value);
  if m = 0 then begin
    AddToErrorStack('Error on Calculatins PWMs: empty vector.');
    b0 := NaN;
    b1 := NaN;
    b2 := NaN;
    b3 := NaN;
    Exit;
  end;

  DeepCopyVector(Col, CopiedCol);
  QuickSort(CopiedCol, 0, m - 1);

  b0 := 0.0; b1 := 0.0; b2 := 0.0; b3 := 0.0;
  for j := 0 to m - 1 do begin
    b0 := b0 + CopiedCol.Value[j] / m;
    if m > 1 then
      b1 := b1 + CopiedCol.Value[j] * (j) / (m * (m - 1));
    if m > 2 then
      b2 := b2 + CopiedCol.Value[j] * (j) * (j - 1) / (m * (m - 1) * (m - 2));
    if m > 3 then
      b3 := b3 + CopiedCol.Value[j] * (j) * (j - 1) * (j - 2) /
             (m * (m - 1) * (m - 2) * (m - 3));
  end;
end;

procedure LMoments(const Col: WVectorType; out L1, L2, L3, L4: FloatType);
var
  n, i: Integer;
  sorted: WVectorType;
  b0, b1, b2, b3, coef: FloatType;
begin
  n := Length(Col.Value);
  if n < 4 then begin
    AddToErrorStack('Error on  L Moments: length of vector <4.');
    L1 := NaN;
    L2 := NaN;
    L3 := NaN;
    L4 := NaN;
    Exit;
  end;

  DeepCopyVector(Col, sorted);
  QuickSort(sorted, 0, n - 1);

  b0 := 0.0; b1 := 0.0; b2 := 0.0; b3 := 0.0;
  for i := 0 to n - 1 do begin
    coef := 1.0 / n;
    b0 := b0 + sorted.Value[i] * coef;
    b1 := b1 + sorted.Value[i] * i * coef / (n - 1);
    b2 := b2 + sorted.Value[i] * i * (i - 1) * coef / ((n - 1) * (n - 2));
    b3 := b3 + sorted.Value[i] * i * (i - 1) * (i - 2) * coef /
           ((n - 1) * (n - 2) * (n - 3));
  end;

  L1 := b0;
  L2 := 2 * b1 - b0;
  L3 := 6 * b2 - 6 * b1 + b0;
  L4 := 20 * b3 - 30 * b2 + 12 * b1 - b0;

  CleanUpZero(L1);
  CleanUpZero(L2);
  CleanUpZero(L3);
  CleanUpZero(L4);
end;

{Summing Routines}
function SumOfSquares(const Col: WVectorType): FloatType; overload;
var
  j: Integer;
  SumSqCol: FloatType;
begin
  SumSqCol := 0.0;
  for j := 0 to High(Col.Value) do
    SumSqCol := SumSqCol + sqr(Col.Value[j]);
  Result := SumSqCol;
end;

function SumOfSquares(const Col: CVectorType): FloatType; overload;
var
  j: Integer;
  SumSqCol: FloatType;
begin
  SumSqCol := 0.0;
  for j := 0 to High(Col) do
    SumSqCol := SumSqCol + sqr(Col[j]);
  Result := SumSqCol;
end;


{ADDITIONAL}
{If p=1 you get the arithmetic mean.
 If p=2 you get the quadratic mean (root mean square).
 If p=3 you get the harmonic mean.
 If p=0 you get the geometric mean (code handles this with the special case).
 If p=0.45 you get fractional order mean, somewhere between geometric and arithmetic.}
function GeneralizedMean(const Col: WVectorType; const p: FloatType): FloatType;
var
  j, m: Integer;
  Sum: FloatType;
begin
  m := Length(Col.Value);
  if m = 0 then begin
    AddToErrorStack('Error on Generalized Mean: empty vector.');
    Exit(NaN);
  end;

  if p = 0.0 then
    Exit(GeometricMean(Col));

  Sum := 0.0;
  for j := 0 to m - 1 do
    Sum := Sum + Power(Abs(Col.Value[j]), p);

  Result := Power(Sum / m, 1.0 / p);
  CleanUpZero(Result);
end;

function HarrellDavisQuantile(const Col: WVectorType; p: FloatType): FloatType;
var
  CopiedCol: WVectorType;
  m, idx: Integer;
  f, frac: FloatType;
begin
  m := Length(Col.Value);
  if m = 0 then begin
    AddToErrorStack('Error on Harrell Davis Quantile: empty vector.');
    Exit(NaN);
  end;

  DeepCopyVector(Col, CopiedCol);
  QuickSort(CopiedCol, 0, m - 1);

  f := (m + 1.0 / 3.0) * p + 1.0 / 3.0;
  idx := Floor(f) - 1;
  frac := f - Floor(f);

  if idx < 0 then idx := 0;
  if idx >= m - 1 then
    Result := CopiedCol.Value[m - 1]
  else
    Result := (1 - frac) * CopiedCol.Value[idx] + frac * CopiedCol.Value[idx + 1];

  CleanUpZero(Result);
end;

function Quantile(const Col: WVectorType; q: FloatType): FloatType;
var
  CopiedCol: WVectorType;
  i, m: Integer;
  Pos: FloatType;
begin
  m := Length(Col.Value);
  if m = 0 then begin
    AddToErrorStack('Error on Quantile: empty vector.');
    Exit(NaN);
  end;
  if (q < 0.0) or (q > 1.0) then begin
    AddToErrorStack('Error on Quantile: q must be between 0 and 1.');
    Exit(NaN);
  end;

  DeepCopyVector(Col, CopiedCol);
  QuickSort(CopiedCol, 0, m - 1);

  Pos := q * (m - 1);
  i := Trunc(Pos);

  if i = m - 1 then
    Result := CopiedCol.Value[i]
  else
    Result := CopiedCol.Value[i] + (Pos - i) * (CopiedCol.Value[i + 1] - CopiedCol.Value[i]);

  CleanUpZero(Result);
end;

function PearsonKthMoment(const Col: WVectorType; const k: Integer): FloatType;
var
  j, m: Integer;
  sum, mu: FloatType;
begin
  m := Length(Col.Value);
  if m = 0 then begin
    AddToErrorStack('Error on Pearson''s Kth Moment: empty vector.');
    Exit(NaN);
  end;
  if k < 1 then begin
    AddToErrorStack('Error on Pearson''s Kth Moment: k must be ≥ 1.');
    Exit(NaN);
  end;

  sum := 0.0;
  mu := ArithmeticMean(Col);
  for j := 0 to m - 1 do
    sum := sum + Power(Col.Value[j] - mu, k);

  Result := sum / m;
  CleanUpZero(Result);
end;

function DurbinWatson(const R: WVectorType): FloatType;
var
  j, m: Integer;
  numerator, denominator, diff: FloatType;
begin
  m := Length(R.Value);
  if m < 2 then begin
    AddToErrorStack('Error on Durbin Watson Statistic: length of vector < 2.');
    Exit(NaN);
  end;

  numerator := 0.0;
  denominator := 0.0;
  for j := 0 to m - 1 do
    denominator := denominator + Sqr(R.Value[j]);
  for j := 1 to m - 1 do begin
    diff := R.Value[j] - R.Value[j - 1];
    numerator := numerator + Sqr(diff);
  end;

  if denominator = 0.0 then begin
    AddToErrorStack('Error on Durbin Watson Statistic: zero denominator.');
    Exit(NaN);
  end;

  Result := numerator / denominator;
  CleanUpZero(Result);
end;

procedure PercentOutliers(const Col: WVectorType; const Cut: FloatType; out Bot, Top: FloatType);
var
  Sorted: WVectorType;
  m, i: Integer;
  pTopIndex, pBotIndex: Integer;
  pTopValue, pBotValue: FloatType;
  CountOutTop, CountOutBot: Integer;
begin
  m := Length(Col.Value);
  if m = 0 then begin
    Top := 0.0;
    Bot := 0.0;
    Exit;
  end;

  // Copy and sort
  SetLength(Sorted.Value, m);
  DeepCopyVector(Col, Sorted);
  QuickSort(Sorted, 0, High(Sorted.Value));

  // Compute indexes and values
  pBotIndex := Floor(Cut * (m - 1));
  pTopIndex := Floor((1 - Cut) * (m - 1));
  pBotValue := Sorted.Value[pBotIndex];
  pTopValue := Sorted.Value[pTopIndex];

  // Count outliers
  CountOutTop := 0;
  CountOutBot := 0;
  for i := 0 to m - 1 do begin
    if (Col.Value[i] <= pBotValue) then
      Inc(CountOutBot);
    if (Col.Value[i] >= pTopValue) then
    Inc(CountOutTop);
  end;

  // Return percentage
  Top := (CountOutTop / m) * 100.0;
  Bot := (CountOutBot / m) * 100.0;
//  if DebugOn then
  //  writeln('PctOut ', pbotindex, ' ',ptopindex, ' ', countoutbot, ' ', countouttop); pause;
end;

procedure PercentOutliersMAD(const Col: WVectorType; const Threshold: FloatType;
  out Bot, Top: FloatType);
var
  m, i: Integer;
  med, madRaw, madSigma, z: FloatType;
  absDevs: WVectorType;
  CountOutTop, CountOutBot: Integer;
begin
  m := Length(Col.Value);

  if m = 0 then begin
    Top := 0.0;
    Bot := 0.0;
    Exit;
  end;

  med := Median(Col);

  SetLength(absDevs.Value, m);
  for i := 0 to m - 1 do
    absDevs.Value[i] := Abs(Col.Value[i] - med);

  madRaw := Median(absDevs);  // median of raw array

  madSigma := madRaw * 1.482602218505601860547076529360;

  // If MAD is zero, no spread => no outliers
  if madSigma = 0.0 then begin
    Top := 0.0;
    Bot := 0.0;
    Exit;
  end;

  CountOutTop := 0;
  CountOutBot := 0;

  for i := 0 to m - 1 do begin
    z := (Col.Value[i] - med) / madSigma;

    if z > Threshold then
      Inc(CountOutTop)
    else if z < -Threshold then
      Inc(CountOutBot);
  end;

  Top := (CountOutTop / m) * 100.0;
  Bot := (CountOutBot / m) * 100.0;
end;

end.
