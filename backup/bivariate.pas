unit Bivariate;

{$mode ObjFPC}{$H+}{$I proprietary.txt}

interface

{ WesStat7 September 24, 2025, by Wesley R. Parsons, wespar@bellouth.net, www.wespar.com.}

uses
  BasicFuncs,
  Correlations,
  Crt,
  DataDisplay,
  DataManager,
  Globals,
  Math,
  SysUtils,
  WesUnit;

function JensenShannonDistance(const P, Q: WVectorType): FloatType;
function JensenShannonDivergence(const P, Q: WVectorType): FloatType;
function KLDivergence(const P, Q: WVectorType): FloatType;
function HellingerDistance(const P, Q: WVectorType): FloatType;
function TotalVariationDistance(const P, Q: WVectorType): FloatType;
function BhattacharyyaCoefficient(const P, Q: WVectorType): FloatType;
function BhattacharyyaDistance(const P, Q: WVectorType): FloatType;
function ChiSquareDivergence(const P, Q: WVectorType): FloatType;
function CanberraDistance(const P, Q: WVectorType): FloatType;
procedure MannWhitneyU(const Col1, Col2: WVectorType);
procedure SignTest(const Col1, Col2: WVectorType);
function DistanceCorrelation(const X, Y: WVectorType): FloatType;
procedure PearsonCI(const r: FloatType; const n: Integer; const alpha: FloatType; out rLower, rUpper: FloatType);
procedure BivariateCorrelations(const Col1, Col2: WVectorType);

implementation

{ Kullback–Leibler divergence (asymmetric, [0, ∞]) }
function KLDivergence(const P, Q: WVectorType): FloatType;
var
  i: Integer;
  sum: FloatType;
begin
  if Length(P.Value) <> Length(Q.Value) then begin
    AddToErrorStack('Error on Hellinger Distance: vectors are different lengths.');
    Exit(0.0);
  end;

  sum := 0.0;
  for i := 0 to High(P.Value) do begin
    if (P.Value[i] > 0.0) then
      if (Q.Value[i] > 0.0) then
        sum := sum + P.Value[i] * Ln(P.Value[i] / Q.Value[i])
      else begin
        AddToErrorStack('Error on Kullback Leibler Divergence: if P is positive, Q cannot be <= 0.');
        Exit(NaN);
      end;
    Result := sum;
  end;
end;

{ Jensen–Shannon Divergence (symmetrised, bounded in [0,1]) }
{ JS(P∥Q) = ½·KL(P∥M) + ½·KL(Q∥M) where M = ½(P+Q) }
function JensenShannonDivergence(const P, Q: WVectorType): FloatType;
var
  i: Integer;
  M: WVectorType;
begin
  if Length(P.Value) <> Length(Q.Value) then begin
    AddToErrorStack('Error on Hellinger Distance: vectors are different lengths.');
    Exit(0.0);
  end;

  SetLength(M.Value, Length(P.Value));
  for i := 0 to High(P.Value) do
    M.Value[i] := 0.5 * (P.Value[i] + Q.Value[i]);
  Result := 0.5 * KLDivergence(P, M) + 0.5 * KLDivergence(Q, M);
end;

{ Jensen–Shannon Distance (square-root version, a true metric) }
{ JSDist(P∥Q) = √JS(P∥Q)  ∈ [0,1] }
function JensenShannonDistance(const P, Q: WVectorType): FloatType;
begin
  Result := Sqrt(JensenShannonDivergence(P, Q));
end;

{ Symmetric Kullback–Leibler (Jeffreys divergence) [0, ∞] }
function JeffreysDivergence(const P, Q: WVectorType): FloatType;
begin
  Result := KLDivergence(P, Q) + KLDivergence(Q, P);
end;

{ Hellinger distance (metric)  [0, 1] }
function HellingerDistance(const P, Q: WVectorType): FloatType;
var
  i: Integer;
  sum: FloatType;
begin
  if Length(P.Value) <> Length(Q.Value) then begin
    AddToErrorStack('Error on Hellinger Distance: vectors are different lengths.');
    Exit(0.0);
  end;
  for i := 0 to High(P.Value) do
    if (P.Value[i] < 0) or (P.Value[i] > 1) or (Q.Value[i] < 0) or (Q.Value[i] > 1) then begin
    AddToErrorStack('Error on Hellinger Distance: values are not [0, 1].');
    Exit(NaN);
  end;

  sum := 0.0;
  for i := 0 to High(P.Value) do
    sum := sum + Sqr(Sqrt(P.Value[i]) - Sqrt(Q.Value[i]));
  Result := Sqrt(sum) / Sqrt2;
end;

{ Total Variation distance  [0, 1] }
function TotalVariationDistance(const P, Q: WVectorType): FloatType;
var
  i: Integer;
  sum: FloatType;
begin
  if Length(P.Value) <> Length(Q.Value) then begin
    AddToErrorStack('Error on Total Variation Distance: vectors are different lengths.');
    Exit(0.0);
  end;

  sum := 0.0;
  for i := 0 to High(P.Value) do
    sum := sum + Abs(P.Value[i] - Q.Value[i]);
  Result := 0.5 * sum;
end;

{ Bhattacharyya distance & coefficient }
function BhattacharyyaCoefficient(const P, Q: WVectorType): FloatType;
var
  i: Integer;
  sum: FloatType;
begin
  if Length(P.Value) <> Length(Q.Value) then begin
    AddToErrorStack('Error on Bhattacharyya Coefficient: vectors are different lengths.');
    Exit(0.0);
  end;
  for i := 0 to High(P.Value) do
    if (P.Value[i] < 0) or (P.Value[i] > 1) or (Q.Value[i] < 0) or (Q.Value[i] > 1) then begin
    AddToErrorStack('Error on Bhattarchayya Coefficient: values are not [0, 1].');
    Exit(NaN);
  end;

  sum := 0.0;
  for i := 0 to High(P.Value) do
    sum := sum + Sqrt(P.Value[i] * Q.Value[i]);
  Result := sum;             // ∈ [0,1]  (1 = identical)
end;

function BhattacharyyaDistance(const P, Q: WVectorType): FloatType;
begin
  if IsNan(BhattacharyyaCoefficient(P, Q)) then begin
    Result := NaN;
    Exit;
  end;
 if BhattacharyyaCoefficient(P, Q) = 0.0 then
   Result := 0.0             // ∈ [0, ∞)
 else
   Result := -Ln(BhattacharyyaCoefficient(P, Q));
 CleanUpBoth(Result);
end;

{ Chi-square divergence (symmetric version)  [0, ∞] }
function ChiSquareDivergence(const P, Q: WVectorType): FloatType;
var
  i: Integer;
  sum: FloatType;
begin
  if Length(P.Value) <> Length(Q.Value) then begin
    AddToErrorStack('Error on Chi Square Divergence: vectors are different lengths.');
    Exit(0.0);
  end;

  sum := 0.0;
  for i := 0 to High(P.Value) do begin
     if Q.Value[i] > 0 then
      sum := sum + Sqr(P.Value[i] - Q.Value[i]) / Q.Value[i]
    else if P.Value[i] > 0 then
      Exit(NaN); // Divergence is infinite if Q[i] = 0 but P[i] > 0.
  end;

  Result := sum;
end;

{ Canberra distance (weighted L1, good for sparse data)  [0, ∞] }
function CanberraDistance(const P, Q: WVectorType): FloatType;
var
  i: Integer;
  denom, sum: FloatType;
begin
  if Length(P.Value) <> Length(Q.Value) then begin
    AddToErrorStack('Error on Canberra Distance: vectors are different lengths.');
    Exit(0.0);
  end;

  sum := 0.0;
  for i := 0 to High(P.Value) do begin
    denom := Abs(P.Value[i]) + Abs(Q.Value[i]);
    if denom > 0 then
      sum := sum + Abs(P.Value[i] - Q.Value[i]) / denom;
    // If both are zero, contribution is 0.
  end;
  Result := sum;
end;

procedure PearsonCI(const r: FloatType; const n: Integer; const alpha: FloatType; out rLower, rUpper: FloatType);
var
  z, zCrit, SE, zLower, zUpper: FloatType;
begin
  if (n <= 3) or (Abs(r) >= 1.0) then begin
    writeln('Error on Confidence Interval: length of vector < 3 or |r| > 1.');
    rLower := 0.0;
    rUpper := 0.0;
    Exit;
  end;

  // Fisher z-transform.
  z := 0.5 * Ln((1 + r) / (1 - r));

  // Standard error.
  SE := 1.0 / Sqrt(n - 3);

  // Critical value from standard normal (e.g., 1.96 for 95%).
 case Round(alpha * 100) of
   90: zCrit := 1.64485362695147;
   95: zCrit := 1.95996398454005;
   99: zCrit := 2.57582930354890;
  else
    zCrit := Ln((1 + alpha) / (1 - alpha)) / 2.0;
  end;

  // Confidence interval in z-space.
  zLower := z - zCrit * SE;
  zUpper := z + zCrit * SE;

  // Back-transform to r-space.
  rLower := (Exp(2 * zLower) - 1) / (Exp(2 * zLower) + 1);
  rUpper := (Exp(2 * zUpper) - 1) / (Exp(2 * zUpper) + 1);

  CleanUpZero(rLower);
  CleanUpZero(rUpper);
end;

function DistanceCorrelation(const X, Y: WVectorType): FloatType;
var
  m, i, j: Integer;      // Can I use my covariance routines here?
  A, B: CMatrixType;
  AbarRow, AbarCol, BbarRow, BbarCol: array of FloatType;
  AbarAll, BbarAll: FloatType;
  dCov, dVarX, dVarY: FloatType;
begin
  m := Length(X.Value);
  if (m = 0) or (m <> Length(Y.Value)) then begin
    AddToErrorStack('Error in Distance Correlation: length of vectors are zero or are not equal length.');
    Exit(NaN);
  end;

  SetLength(A, m, m);
  SetLength(B, m, m);
  SetLength(AbarRow, m);
  SetLength(AbarCol, m);
  SetLength(BbarRow, m);
  SetLength(BbarCol, m);

  // Step 1: Compute pairwise distance matrices.
  for i := 0 to m - 1 do
    for j := 0 to m - 1 do begin
      A[i, j] := Abs(X.Value[i] - X.Value[j]);
      B[i, j] := Abs(Y.Value[i] - Y.Value[j]);
    end;

  // Step 2: Compute row/col means and grand mean.
  AbarAll := 0.0;
  BbarAll := 0.0;
  for i := 0 to m - 1 do begin
    for j := 0 to m - 1 do begin
      AbarRow[i] := AbarRow[i] + A[i, j];
      AbarCol[j] := AbarCol[j] + A[i, j];
      BbarRow[i] := BbarRow[i] + B[i, j];
      BbarCol[j] := BbarCol[j] + B[i, j];
      AbarAll := AbarAll + A[i, j];
      BbarAll := BbarAll + B[i, j];
    end;
  end;

  for i := 0 to m - 1 do begin
    AbarRow[i] := AbarRow[i] / m;
    AbarCol[i] := AbarCol[i] / m;
    BbarRow[i] := BbarRow[i] / m;
    BbarCol[i] := BbarCol[i] / m;
  end;
  AbarAll := AbarAll / (m * m);
  BbarAll := BbarAll / (m * m);

  // Step 3: Double-center and compute dCov, dVarX, dVarY.
  dCov := 0.0;
  dVarX := 0.0;
  dVarY := 0.0;
  for i := 0 to m - 1 do
    for j := 0 to m - 1 do begin
      A[i, j] := A[i, j] - AbarRow[i] - AbarCol[j] + AbarAll;
      B[i, j] := B[i, j] - BbarRow[i] - BbarCol[j] + BbarAll;
      dCov := dCov + A[i, j] * B[i, j];
      dVarX := dVarX + Sqr(A[i, j]);
      dVarY := dVarY + Sqr(B[i, j]);
    end;

  dCov := dCov / (m * m);
  dVarX := dVarX / (m * m);
  dVarY := dVarY / (m * m);

  if (dVarX = 0.0) or (dVarY = 0.0) then begin
    AddToErrorStack('Error in Distance Correlation: zero distance variance.');
    Exit(NaN);
  end;

  Result := Sqrt(dCov / Sqrt(dVarX * dVarY));
  CleanUpBoth(Result);
end;

procedure MannWhitneyU(const Col1, Col2: WVectorType);
var
  n1, n2, nTotal, i: Integer;
  Ranks, Combined: CVectorType;
  RankSum1, RankSum2, U1, U2, UMean, UVar, Z: FloatType;

  procedure SortWithRanks(var Data, Rank: CVectorType);
  var
    i, j: Integer;
    tmp: FloatType;
  begin
    // Simple bubble sort with rank assignment.
    for i := 0 to High(Data) do
      Rank[i] := i + 1; // Initial ranks.

    for i := 0 to High(Data) - 1 do
      for j := i + 1 to High(Data) do
        if Data[i] > Data[j] then begin
          tmp := Data[i]; Data[i] := Data[j]; Data[j] := tmp;
          tmp := Rank[i]; Rank[i] := Rank[j]; Rank[j] := tmp;
        end;
  end;

begin
  n1 := Length(Col1.Value);
  n2 := Length(Col2.Value);
  nTotal := n1 + n2;

  if (n1 = 0) or (n2 = 0) then begin
    Writeln('Error on Mann-Whitney U: one of the vectors is empty.');
    Exit;
  end;

  // Combine data.
  SetLength(Combined, nTotal);
  SetLength(Ranks, nTotal);

  for i := 0 to n1 - 1 do
    Combined[i] := Col1.Value[i];
  for i := 0 to n2 - 1 do
    Combined[n1 + i] := Col2.Value[i];

  // Sort and assign ranks.
  SortWithRanks(Combined, Ranks);

  // Compute rank sums.
  RankSum1 := 0.0;
  for i := 0 to n1 - 1 do
    RankSum1 := RankSum1 + Ranks[i];
  RankSum2 := (nTotal * (nTotal + 1)) / 2 - RankSum1;

  // Compute U statistics.
  U1 := RankSum1 - n1 * (n1 + 1) / 2;
  U2 := n1 * n2 - U1;

  // Mean and variance of U.
  UMean := n1 * n2 / 2.0;
  UVar := n1 * n2 * (n1 + n2 + 1) / 12.0;

  // Normal approximation (z-score).
  Z := (U1 - UMean) / Sqrt(UVar);

  // Output
  Writeln('Mann-Whitney U test results:');
  Writeln('  RankSum1 = ', SmartFloat(RankSum1));
  Writeln('  RankSum2 = ', SmartFloat(RankSum2));
  Writeln('  U1 = ', U1:0:4, ', U2 = ', SmartFloat(U2));
  Writeln('  Normal approximation Z = ', SmartFloat(Z));
end;

// Helper function: Binomial coefficient.
function BinomialCoefficient(n, k: Integer): FloatType;
var
  i: Integer;
  ResultValue: FloatType;
begin
  if (k < 0) or (k > n) then Exit(0);
  if (k = 0) or (k = n) then Exit(1);

  // Symmetry for numerical stability.
  if k > n div 2 then
    k := n - k;

  ResultValue := 1.0;
  for i := 1 to k do begin
    ResultValue := ResultValue * (n - k + i) / i;
  end;

  Result := ResultValue;
end;

procedure SignTest(const Col1, Col2: WVectorType);
var
  m, j, nTotal, nPos, nNeg, nZero: Integer;
  Diff, pValue: FloatType;
begin
  m := Length(Col1.Value);
  if m <> Length(Col2.Value) then begin
    Writeln('Error on sign test: vectors have different lengths');
    Exit;
  end;

  nTotal := m;
  nPos := 0;
  nNeg := 0;
  nZero := 0;

  // Count signs of differences.
  for j := 0 to m - 1 do begin
    Diff := Col1.Value[j] - Col2.Value[j];
    if Diff > 0 then
      Inc(nPos)
    else if Diff < 0 then
      Inc(nNeg)
    else
      Inc(nZero);
  end;

  // Binomial test: probability of observing this many positives
  // under null hypothesis p = 0.5.
  // Here use is a simple approximation: two-sided p-value.
  pValue := 2 * Min(
    Power(0.5, nTotal - nZero) * BinomialCoefficient(nTotal - nZero, nPos),
    Power(0.5, nTotal - nZero) * BinomialCoefficient(nTotal - nZero, nNeg)
  );

  // Output results
  Writeln('Sign Test results:');
  Writeln('  Total pairs = ', nTotal);
  Writeln('  Positive differences = ', nPos);
  Writeln('  Negative differences = ', nNeg);
  Writeln('  Zero differences = ', nZero);
  Writeln('  Approximate two-sided p-value = ', SmartFloat(pValue));
end;

procedure BivariateCorrelations(const Col1, Col2: WVectorType);
var
  r, adjr, sr, v1, v2, cov, rho, srho, tau, stau, hoeff, shoeff, MutInfo, NormMutInfo, dc,
    JSDist, JSDiv, KLDiv, HellDist, TVDist, BhattCoef, BhattDist, ChiDiv, CanDist,
    rLower90, rUpper90, rLower95, rUpper95, rLower99, rUpper99: FloatType;
begin
  if (Length(Col1.Value) = 0) or (Length(Col1.Value) = 0) then begin
    writeln('Error in Univariate Statistics: one selected variable is empty.');
    Exit;
  end;

  // Compute Pearson statistics.
  r := PearsonsR(Col1, Col2);
  CleanUpBoth(r);
  sr := PearsonsRSignificance(r, Length(Col1.Value));
  CleanUpBoth(sr);
  if (nRow > 2) then
    adjr := 1 - ((1 - r * r) * (nRow - 1) / (nRow - 2))
  else
    adjr := 0.0;
  CleanUpBoth(adjr);

  // Display Pearson statistics.
  writeTab('Pearson''s R = ', r);
  writeTabln('Pearson''s R Significance = ', sr);
  writeTabln('Pearson''s R-Squared = ', (r * r));
  writeTab('Pearson''s Adjusted R = ', adjr);
  writeTab('Pearson''s Adjusted R-Squared = ', (adjr * adjr));

  // Compute variance statistics.
  v1 := Variance(Col1, True);
  v2 := Variance(Col2, True);
  cov := Covariance(Col1, Col2);

  // Display Pearson statistics.
  writeTab('Variance of first variable: ', v1);
  writeTabln('Variance of second variable: ', v2);
  writeTabln('Covariance of variables: ', cov);

  // Compute confidence intervals.
  PearsonCI(r, nRow, 0.90, rLower90, rUpper90);
  PearsonCI(r, nRow, 0.95, rLower95, rUpper95);
  PearsonCI(r, nRow, 0.99, rLower99, rUpper99);

  // Display confidence intervals.
  writeln('Confidence intervals at 90%.     Upper: ', SmartFloat(rUpper90), '   Lower: ', SmartFloat(rLower90));
  writeln('Confidence intervals at 95%.     Upper: ', SmartFloat(rUpper95), '   Lower: ', SmartFloat(rLower95));
  writeln('Confidence intervals at 99%.     Upper: ', SmartFloat(rUpper99), '   Lower: ', SmartFloat(rLower99));

  // Compute non-parametric statistics.
  rho := SpearmansRho(Col1, Col2);
  CleanUpBoth(rho);
  srho := SpearmansRhoSignificance(Col1, Col2, rho);
  CleanUpBoth(srho);
  tau := KendallsTau(Col1, Col2);
  CleanUpZero(tau);
  stau := KendallsTauSignificance(Col1, Col2);
  CleanUpBoth(stau);
  hoeff := HoeffdingD(Col1, Col2, SampleCutoff, SampleSize);
  CleanUpZero(hoeff);
  shoeff := HoeffdingsDSignificance(Col1, Col2, SampleCutoff, SampleSize, Permutations, hoeff);
  CleanUpBoth(shoeff);
  MutInfo := MutualInformation(Col1, Col2, Bins);
  CleanUpZero(MutInfo);
  NormMutInfo := NormalizedMutualInformation(Col1, Col2, Bins);
  CleanUpBoth(NormMutInfo);
  dc := DistanceCorrelation(Col1, Col2);
  CleanUpBoth(dc);

  // Display non-parametric statistics.
  writeTab('Spearman''s Rho = ', rho);
  writeTabln('Spearman''s Rho Significance = ', srho);
  writeTab('Kendall''s Tau: ', tau);
  writeTabln('Kendall''s Tau Significance = ', stau);
  writeTab('Hoeffding''s D: ', hoeff);
  writeTabln('Hoeffding''s D Significance = ', shoeff);
  writeTabln('Mutual Information = ', MutInfo);
  writeTab('Normalized Mutual Information: ', NormMutInfo);
  writeTabln('Distance Correlation: ', dc);

  // Distance and Divergence statistics.
  JSDist := JensenShannonDistance(Col1, Col2);
  JSDiv := JensenShannonDivergence(Col1, Col2);
  writeTab('Jensen Shannon Distance = ', JSDist);
  writeTabln('Jensen Shannon Divergence = ', JSDiv);
  KLDiv := KLDivergence(Col1, Col2);
  writeTabln('Kullback Leibler Divergence = ', KLDiv);
  HellDist := HellingerDistance(Col1, Col2);
  writeTab('Hellinger Distance = ', HellDist);
  TVDist := TotalVariationDistance(Col1, Col2);
  writeTabln('Total Variation Distance = ', TVDist);
  BhattCoef := BhattacharyyaCoefficient(Col1, Col2);
  writeTab('Bhattarchayya Coefficient = ', BhattCoef);
  BhattDist := BhattacharyyaDistance(Col1, Col2);
  writeTabln('Bhattarchayya Distance = ', BhattDist);
  ChiDiv := ChiSquareDivergence(Col1, Col2);
  writeTab('Chi Square Divergence = ', ChiDiv);
  CanDist := CanberraDistance(Col1, Col2);
  writeTabln('Canberra Distance = ', CanDist);

  // Other statistics.
  MannWhitneyU(Col1, Col2);
  SignTest(Col1, Col2);

  // Final touches.
  Writeln(TwoTailed);
  DisplayErrorStack;
  if FirstPass then
    SavePartialData('Save bivariate statistics to a file? (y/n) ', Col1, Col2, @BivariateCorrelations);
  Pause;
end;

end.
