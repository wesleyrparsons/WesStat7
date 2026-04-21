unit Correlations;

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
  SysUtils,
  WesUnit;

{ Statistics Functions }
function IncompleteGamma(const a, x: FloatType): FloatType;
function Gamma(const x: FloatType): FloatType;
function NewLnGamma(const x: FloatType): FloatType;                       // Not used.
function nLnGamma(z: FloatType): FloatType;
function IncompleteBetaSimpson(const x, a, b: FloatType): FloatType;      // Not used.
function IncompleteBetaCF(const x, a, b: Double): Double;
function Beta(const a, b: FloatType): FloatType;
function Erf(const z: FloatType): FloatType;
function Erfc(const a: FloatType): FloatType;
function CDFNormal(const z: FloatType): FloatType;
function TDistPValue(const t: FloatType; const df: Integer): FloatType;

{ Mutual Information }
//function Entropy(const X: WVectorType; const NumBins: Integer): FloatType;
function MutualInformation(const X, Y: WVectorType; const NumBins: Integer): FloatType;
function NormalizedMutualInformation(const X, Y: WVectorType; const NumBins: Integer): FloatType;

{ Hoeffding }
function HoeffdingD(const Col1, Col2: WVectorType;
  Cutoff: Integer = 0; SampleSize: Integer = 0): FloatType;
function HoeffdingsDSignificance(const Col1, Col2: WVectorType; Cutoff, SampleSize, NumPermutations: Integer;
  out OutObservedD: FloatType): FloatType;

{ Kendall }
function KendallsTau(const Col1, Col2: WVectorType): FloatType;
function KendallsTauSignificance(const X, Y: WVectorType): FloatType;

{ Spearman }
function SpearmansRho(const x, y: WVectorType): FloatType;
function SpearmansRhoSignificance(const X, Y: WVectorType; const rho: Float): FloatType;
function AsymptoticSpearmanP(const m: Integer; const rho: FloatType): FloatType;

{ Covariance & Pearson Correlation }
function PearsonsR(const x, y: WVectorType): FloatType;
function PearsonsRSignificance(const r: FloatType; const n: Integer): FloatType;
function Covariance(const x, y: WVectorType): FloatType;
procedure CreateCovarianceMatrix(const Data: WMatrixType; out CovMatrix: CMatrixType);
procedure CreateCorrelationMatrix(const Data: WMatrixType; out CorrMatrix: CMatrixType);
procedure CreateCorrelationFromCovarianceMatrix(const CovMatrix: CMatrixType; out CorrMatrix: CMatrixType);
procedure MultivariateCorrelations(const WVar: IVectorType);

implementation

type
  TValueIndex = record
    Value: FloatType;
    Index: Integer;
  end;


{ Lower regularized incomplete gamma function P(a, x) = γ(a,x)/Γ(a) }
function IncompleteGamma(const a, x: FloatType): FloatType;
const
  FpMin = 1.0e-30;          { Avoid division by zero }      //makre this minrealnumber
var
  sum, del, ap, gln, b, c, d, h, an: FloatType;
  i, LocalMaxIter: Integer;
begin
  if (x < 0.0) or (a <= 0.0) then begin
    AddToErrorstack('Error on gamma function: argument is zero.');
    Result := NaN;
    Exit;
  end;
  gln := nLnGamma(a);
  if x < (a + 1.0) then begin
    { --- Series representation --- }
    ap := a;
    sum := 1.0 / a;
    del := sum;
    while True do begin
      ap := ap + 1.0;
      del := del * x / ap;
      sum := sum + del;
      if Abs(del) < Abs(sum) * LowTolerance then Break;
    end;
    Result := sum * Exp(-x + a * Ln(x) - gln);
  end
  else begin
    { --- Continued fraction for Q(a,x), return P = 1 - Q --- }
    b := x + 1.0 - a;
    c := 1.0 / FPMIN;
    d := 1.0 / b;
    h := d;
    LocalMaxIter := IfThen(UserMaxIter > 0, UserMaxIter, MediumMaxIter);
    for i := 1 to LocalMaxIter do begin
      an := -i * (i - a);
      b := b + 2.0;
      d := an * d + b;
      if Abs(d) < FpMin then d := FpMin;
      c := b + an / c;
      if Abs(c) < FpMin then c := FpMin;
      d := 1.0 / d;
      del := d * c;
      h := h * del;
      if Abs(del - 1.0) < Tolerance then Break;
    end;
    Result := 1.0 - h * Exp(-x + a * Ln(x) - gln);
  end;
end;

// Gamma function.
function Gamma(const x: FloatType): FloatType;
const
  Coeff: array[0..8] of FloatType = (
    0.99999999999980993,
    676.5203681218851,
   -1259.1392167224028,
    771.32342877765313,
   -176.61502916214059,
    12.507343278686905,
   -0.13857109526572012,
    9.9843695780195716e-6,
    1.5056327351493116e-7
  );
  g = 7.0;
  Sqrt2Pi = 2.506628274631000502415765284811;
  MaxGam = 171.6243769563027;
var
  i: Integer;
  sum, t, y: FloatType;
begin
  if x <= 0 then Exit(NaN);
  if x < 0.5 then begin
    // Reflection formula.
    if IsEffectivelyEqual(Sin(Pi * x), 0.0) then
      Result := 0.0;
    Result := Pi / (Sin(Pi * x) * Gamma(1.0 - x));
    Exit;
  end;

  if x > MaxGam then
    Result := 0.0;

  sum := Coeff[0];
  for i := 1 to High(Coeff) do
    sum := sum + Coeff[i] / (x + i - 1);

  t := x + g - 0.5;
  y := Power(t, x - 0.5) * Exp(-t) * sum * Sqrt2Pi;

  if IsEffectivelyEqual(y, 0.0) or IsNaN(y) then
        Result := 0.0;

  Result := y;
end;

// Log Gamma Function.
// This is the Spouge method. With a=15, 14 decimnals. A>15 doesn't get getter.
function nlnGamma(z: FloatType): FloatType;
const
  a = 20;
var
  c: Array[0..20] of FloatType;
  k: Integer;
  sum_s, term, ln_fact, temp: FloatType;
begin
  if z <= 0.0 then Exit(0.0);

  // c[0] = sqrt(2*pi)
  c[0] := Sqrt2pi;

  // Compute c[k] for k = 1 to a-1.
  ln_fact := 0.0;  // ln((k-1)!) to avoid overflow
  for k := 1 to a - 1 do begin
    temp := a - k;
    if k = 1 then
      ln_fact := ln(1)
    else
      ln_fact := ln_fact + Ln(k - 1);  // Build ln((k-1)!).

    // Compute ln|c[k]| = ln(temp^(k-0.5)) + temp - ln_fact
    // = (k-0.5)*Ln(temp) + temp - ln_fact
    term := (k - 0.5) * Ln(temp) + temp - ln_fact;

    // Determine sign: (-1)^{k+1}.
    if Odd(k) then
      c[k] := Exp(term)
    else
      c[k] := -Exp(term);
  end;

  // Sum: c[0] + Σ c[k] / (z + k - 1).
  sum_s := c[0];
  for k := 1 to a - 1 do
    sum_s := sum_s + c[k] / (z + k - 1.0);

  // Final formula.
  Result := (z - 0.5) * Ln(z + a - 1.0) - (z + a - 1.0) + Ln(sum_s);
end;

// Log Gamma Function. This gives about 13 decimals. 11-8-2025.
function NewLnGamma(const x: FloatType): FloatType;
const
  g = 7.0;
  p0: array[0..8] of FloatType = (
    0.99999999999980993,
    676.5203681218851,
   -1259.1392167224028,
    771.32342877765313,
   -176.61502916214059,
    12.507343278686905,
   -0.13857109526572012,
    9.9843695780195716e-6,
    1.5056327351493116e-7
  );
var
  i: Integer;
  y, xx, sum, t, sinpix: FloatType;
  reflect: Boolean;
begin
  // Handle reflection if x <= 0.
  reflect := x < 0.5;
  if reflect then
  begin
    y := 1.0 - x; // use reflection.
    sinpix := Sin(Pi * x);
    if sinpix = 0.0 then
      Exit(NaN); // pole at negative integer.
  end
  else
    y := x;

  // Lanczos approximation for y >= 0.5.
  xx := y - 1.0;
  sum := p0[0];
  for i := 1 to High(p0) do
    sum := sum + p0[i] / (xx + i);
  t := xx + g + 0.5;
  Result := 0.5 * Ln(2.0 * Pi) + (xx + 0.5) * Ln(t) - t + Ln(sum);

  // Apply reflection correction only once.
  if reflect then
    Result := Ln(Pi) - Ln(Abs(sinpix)) - Result;
end;

// Beta function: B(a, b) = Gamma(a) * Gamma(b) / Gamma(a + b).
function Beta(const a, b: FloatType): FloatType;
begin
  if (a <= 0) or (b <= 0) then begin
    AddToErrorstack('Error on beta function: argument is zero.');
    Exit(NaN);
  end;
  Result := Gamma(a) * Gamma(b) / Gamma(a + b);
end;

// Log Gamma Function. Simpson's Method Version. Use the LMath version instead.
function SafeLnGamma(const z: FloatType): FloatType;
const
  g = 7.0;
  p: array[0..8] of FloatType = (
    0.99999999999980993,
    676.5203681218851,
   -1259.1392167224028,
    771.32342877765313,
   -176.61502916214059,
     12.507343278686905,
     -0.13857109526572012,
      9.984369578019571e-6,
      1.5056327351493116e-7
  );
var
  i: Integer;
  y, t, sum: FloatType;
begin
  if z < 0.5 then
    Exit(Ln(Pi) - Ln(Sin(Pi*z)) - SafeLnGamma(1.0 - z));

  y := z - 1.0;
  sum := p[0];
  for i := 1 to High(p) do
    sum := sum + p[i] / (y + i);
  t := y + g + 0.5;
  Result := 0.5 * Ln(2 * Pi) + (y + 0.5) * Ln(t) - t + Ln(sum);
end;

// Beta Function related to above.
function BetaFunc(const a, b: FloatType): FloatType;
begin
  Result := Exp(SafeLnGamma(a) + SafeLnGamma(b) - SafeLnGamma(a + b));
end;

// Integrand for incomplete beta function.
function BetaIntegrand(const t, a, b: FloatType): FloatType;
begin
  if (t <= 0.0) or (t >= 1.0) then Exit(0.0);
  Result := Power(t, a - 1.0) * Power(1.0 - t, b - 1.0);
end;

// Beta Function. Simpson’s rule numerical integration from 0 to x.
function SimpsonBeta(const x, a, b: FloatType; const n: Integer): FloatType;
var
  i: Integer;
  h, sum: FloatType;
begin
  if (x <= 0.0) then Exit(0.0);
  if (x >= 1.0) then Exit(BetaFunc(a, b));

  h := x / n;
  sum := BetaIntegrand(0.0, a, b) + BetaIntegrand(x, a, b);

  for i := 1 to n - 1 do begin
    if Odd(i) then
      sum := sum + 4.0 * BetaIntegrand(i * h, a, b)
    else
      sum := sum + 2.0 * BetaIntegrand(i * h, a, b);
  end;

  Result := sum * h / 3.0;
end;

// Regularized incomplete beta using numerical integration.
function IncompleteBetaSimpson(const x, a, b: FloatType): FloatType;
const
  N = 5000; // number of subintervals (even)
var
  bx, btot: FloatType;
begin
  if (x <= 0.0) then Exit(0.0);
  if (x >= 1.0) then Exit(1.0);

  bx := SimpsonBeta(x, a, b, N);
  btot := BetaFunc(a, b);
  Result := bx / btot;
end;

// LMath Version of IncBeta. Use this version }
const
  Big = 1.0 / MachineEpsilon;
const
  MaxLog  = 709.782712893384;
  MinLog  = -708.3964185322641;
  MaxGam  = 170.6243769563027;

  { Power series for incomplete beta integral. Use when B*X is small }
function PSeries(x, A, B : FloatType) : FloatType;
  var
    S, T, U, V, T1, Z, Ai : FloatType;
    N : Integer;
  begin
    Ai := 1.0 / A;
    U := (1.0 - B) * X;
    V := U / (A + 1.0);
    T1 := V;
    T := U;
    N := 2;
    S := 0.0;
    Z := MachineEpsilon * Ai;
    while Abs(V) > Z do
      begin
        U := (N - B) * X / N;
        T := T * U;
        V := T / (A + N);
        S := S + V;
        N := N + 1;
      end;
    S := S + T1;
    S := S + Ai;

    U := A * Ln(X);
    if (A + B < MaxGam) and (Abs(U) < MaxLog) then
      begin
        T := Gamma(A + B) / (Gamma(A) * Gamma(B));
        S := S * T * Power(X, A);
      end
    else
      begin
        T := nLnGamma(A + B) - nLnGamma(A) - nLnGamma(B) + U + Ln(S);
        if T < MinLog then
          S := 0.0
        else
          S := Exp(T);
      end;
    PSeries := S;
  end;

function IsZero(const x: Double): Boolean;
const
  MachEp = 2.220446049250313e-16;
begin
  Result := Abs(x) < MachEp;
end;

  function CFrac1(x, A, B : FloatType) : FloatType;
  { Continued fraction expansion #1 for incomplete beta integral }
  var
    Xk, Pk, Pkm1, Pkm2, Qk, Qkm1, Qkm2,
    K1, K2, K3, K4, K5, K6, K7, K8,
    R, T, Ans, Thresh: FloatType;
    N, LocalMaxIter: Integer;
  label
    CDone;
  begin
    K1 := A;
    K2 := A + B;
    K3 := A;
    K4 := A + 1.0;
    K5 := 1.0;
    K6 := B - 1.0;
    K7 := K4;
    K8 := A + 2.0;

    Pkm2 := 0.0;
    Qkm2 := 1.0;
    Pkm1 := 1.0;
    Qkm1 := 1.0;
    Ans := 1.0;
    R := 1.0;
    N := 0;
    Thresh := 3.0 * MachineEpsilon;
    LocalMaxIter := IfThen(UserMaxIter > 0, UserMaxIter, MediumMaxIter);

    repeat
      Xk := - (X * K1 * K2) / (K3 * K4);
      Pk := Pkm1 + Pkm2 * Xk;
      Qk := Qkm1 + Qkm2 * Xk;
      Pkm2 := Pkm1;
      Pkm1 := Pk;
      Qkm2 := Qkm1;
      Qkm1 := Qk;

      Xk := (X * K5 * K6) / (K7 * K8);
      Pk := Pkm1 + Pkm2 * Xk;
      Qk := Qkm1 + Qkm2 * Xk;
      Pkm2 := Pkm1;
      Pkm1 := Pk;
      Qkm2 := Qkm1;
      Qkm1 := Qk;

      if not IsZero(Qk) then R := Pk / Qk;

      if not IsZero(R) then begin
        T := Abs((Ans - R) / R);
        Ans := R;
      end
      else
        T := 1.0;

      if T < Thresh then goto CDone;

      K1 := K1 + 1.0;
      K2 := K2 + 1.0;
      K3 := K3 + 2.0;
      K4 := K4 + 2.0;
      K5 := K5 + 1.0;
      K6 := K6 - 1.0;
      K7 := K7 + 2.0;
      K8 := K8 + 2.0;

      if Abs(Qk) + Abs(Pk) > Big then begin
        Pkm2 := Pkm2 * MachineEpsilon;
        Pkm1 := Pkm1 * MachineEpsilon;
        Qkm2 := Qkm2 * MachineEpsilon;
        Qkm1 := Qkm1 * MachineEpsilon;
      end;

      if (Abs(Qk) < MachineEpsilon) or (Abs(Pk) < MachineEpsilon) then begin
        Pkm2 := Pkm2 * Big;
        Pkm1 := Pkm1 * Big;
        Qkm2 := Qkm2 * Big;
        Qkm1 := Qkm1 * Big;
      end;
      N := N + 1;
    until N > LocalMaxIter;
    writeln('Error in Incomplete Beta: Floating point loss.');

  CDone:
    CFrac1 := Ans;
  end;

  function CFrac2(x, A, B : FloatType) : FloatType;
  { Continued fraction expansion #2 for incomplete beta integral }
  var
    Xk, Pk, Pkm1, Pkm2, Qk, Qkm1, Qkm2,
    K1, K2, K3, K4, K5, K6, K7, K8,
    R, T, Z, Ans, Thresh: FloatType;
    N, LocalMaxIter: Integer;
  label
    CDone;
  begin
    K1 := A;
    K2 := B - 1.0;
    K3 := A;
    K4 := A + 1.0;
    K5 := 1.0;
    K6 := A + B;
    K7 := A + 1.0;
    K8 := A + 2.0;

    Pkm2 := 0.0;
    Qkm2 := 1.0;
    Pkm1 := 1.0;
    Qkm1 := 1.0;
    Z := X / (1.0 - X);
    Ans := 1.0;
    R := 1.0;
    N := 0;
    Thresh := 3.0 * MachineEpsilon;
    LocalMaxIter := IfThen(UserMaxIter > 0, UserMaxIter, MediumMaxIter);

    repeat
      Xk := - (Z * K1 * K2) / (K3 * K4);
      Pk := Pkm1 + Pkm2 * Xk;
      Qk := Qkm1 + Qkm2 * Xk;
      Pkm2 := Pkm1;
      Pkm1 := Pk;
      Qkm2 := Qkm1;
      Qkm1 := Qk;

      Xk := (Z * K5 * K6) / (K7 * K8);
      Pk := Pkm1 + Pkm2 * Xk;
      Qk := Qkm1 + Qkm2 * Xk;
      Pkm2 := Pkm1;
      Pkm1 := Pk;
      Qkm2 := Qkm1;
      Qkm1 := Qk;

      if Qk <> 0.0 then R := Pk / Qk;

      if R <> 0.0 then begin
        T := Abs((Ans - R) / R);
        Ans := R;
      end
      else
        T := 1.0;

      if T < Thresh then goto CDone;

      K1 := K1 + 1.0;
      K2 := K2 - 1.0;
      K3 := K3 + 2.0;
      K4 := K4 + 2.0;
      K5 := K5 + 1.0;
      K6 := K6 + 1.0;
      K7 := K7 + 2.0;
      K8 := K8 + 2.0;

      if Abs(Qk) + Abs(Pk) > Big then  begin
        Pkm2 := Pkm2 * MachineEpsilon;
        Pkm1 := Pkm1 * MachineEpsilon;
        Qkm2 := Qkm2 * MachineEpsilon;
        Qkm1 := Qkm1 * MachineEpsilon;
      end;

      if (Abs(Qk) < MachineEpsilon) or (Abs(Pk) < MachineEpsilon) then begin
        Pkm2 := Pkm2 * Big;
        Pkm1 := Pkm1 * Big;
        Qkm2 := Qkm2 * Big;
        Qkm1 := Qkm1 * Big;
      end;
      N := N + 1;
    until N > LocalMaxIter;
    writeln('Error on Incomplete Beta: Floating point loss.');

  CDone:
    CFrac2 := Ans;
  end;

{ Regularized incomplete beta I_x(a,b) with Lanczos LnGamma + Lentz CF }
// IncBetaCF: signature (x, a, b) -> I_x(a,b).
// This version also works, about as well as LMath.
function IncompleteBetaCF(const x, a, b: FloatType): FloatType;
var
  m, m2, LocalMaxIter: Integer;
  aa, c, d, del, h, qab, qap, qam, bt, threshold: FloatType;
begin
  if x <= 0.0 then Exit(0.0);
  if x >= 1.0 then Exit(1.0);
  if (a <= 0.0) or (b <= 0.0) then Exit(NaN);

  bt := Exp(nLnGamma(a + b) - nLnGamma(a) - nLnGamma(b) + a * Ln(x) + b * Ln(1.0 - x));

  // Choose branch; use <= to avoid mutual recursion at exact boundary.
  threshold := (a + 1.0) / (a + b + 2.0);
  if x <= threshold then begin
    qab := a + b;
    qap := a + 1.0;
    qam := a - 1.0;
    c := 1.0;
    d := 1.0 - qab * x / qap;
    if Abs(d) < VeryTinyNumber then d := VeryTinyNumber;
    d := 1.0 / d;
    h := d;

    LocalMaxIter := IfThen(UserMaxIter > 0, UserMaxIter, MediumMaxIter);
    for m := 1 to LocalMaxIter do begin
      m2 := 2 * m;
      // Even.
      aa := m * (b - m) * x / ((qam + m2) * (a + m2));
      d := 1.0 + aa * d; if Abs(d) < VeryTinyNumber then d := VeryTinyNumber;
      c := 1.0 + aa / c; if Abs(c) < VeryTinyNumber then c := VeryTinyNumber;
      d := 1.0 / d;
      h := h * d * c;
      // Odd.
      aa := -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
      d := 1.0 + aa * d; if Abs(d) < VeryTinyNumber then d := VeryTinyNumber;
      c := 1.0 + aa / c; if Abs(c) < VeryTinyNumber then c := VeryTinyNumber;
      d := 1.0 / d;
      del := d * c;
      h := h * del;
      if Abs(del - 1.0) < HighTolerance then Break;
    end;
    Result := bt * h / a;
    Exit;
  end;

  // Use symmetry (x > threshold).
  Result := 1.0 - IncompleteBetaCF(1.0 - x, b, a);
end;

// Error Function. This is from EBuxbaun, Github. Gives 13 decimals.
function Erf(const z: FloatType): FloatType;
const
  A: array[1..14] of FloatType =
    (1.1283791670955,      0.34197505591854,     0.86290601455206E-1,
     0.12382023274723E-1,  0.11986242418302E-2,  0.76537302607825E-4,
     0.25365482058342E-5, -0.99999707603738,    -1.4731794832805,
    -1.0573449601594,     -0.44078839213875,    -0.10684197950781,
    -0.12636031836273E-1, -0.1149393366616E-8);

  B: array[1..12] of FloatType =
    (-0.36359916427762,     0.52205830591727E-1, -0.30613035688519E-2,
     -0.46856639020338E-4,  0.15601995561434E-4, -0.62143556409287E-6,
      2.6015349994799,      2.9929556755308,      1.9684584582884,
      0.79250795276064,     0.18937020051337,     0.22396882835053E-1);
var
  U, X, S: FloatType;
begin
  X := ABS(Z);
  if Z >= 0.0 then
    S := 1.0
  else S := -1.0;
  if (abs(Z) < 0.0) then
    Erf := 0.0
  else if (X >= 5.5) then
    Erf := S
  else begin
    U := X * X;
    if (X <= 1.5) then
      Erf :=
      (X*EXP(-U)*(A[1] + U*(A[2] + U*(A[3] + U*(A[4] +
       U*(A[5] + U*(A[6] + U*A[7])))))) / (1.0 + U * (B[1] +
       U*(B[2] + U*(B[3] + U*(B[4] + U*(B[5] + U*B[6]))))))) * S
    else
      Erf :=
      (EXP(-U)*(A[8] + X*(A[9] + X*(A[10] + X*(A[11] + X*(A[12] +
       X*(A[13] + X*A[14])))))) / (1.0 + X*(B[7] + X*(B[8] +
       X*(B[9] + X*(B[10] + X*(B[11] + X*B[12])))))) + 1.0) * S;
    end;
end;

// Error Function.
// Abramowitz & Stegun 7.1.26 approximation for erf(x).
// Good, fast approximation (typical absolute error ~1e-7).
// This function has about 7 decimals precision. 11/8/2025.
function ASxErf(const x: FloatType): FloatType;
const
  a1: FloatType = 0.254829592;
  a2: FloatType = -0.284496736;
  a3: FloatType = 1.421413741;
  a4: FloatType = -1.453152027;
  a5: FloatType = 1.061405429;
  p : FloatType = 0.3275911;
var
  sign, z, t, y: FloatType;
begin
  // Save sign and work with absolute value
  if x < 0.0 then sign := -1.0 else sign := 1.0;
  z := Abs(x);

  // Abramowitz & Stegun approximation
  t := 1.0 / (1.0 + p * z);
  y := 1.0 - (((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t) * Exp(-z * z);

  Result := sign * y;
end;

// Complementary error function (erfc),
function Erfc(const a: FloatType): FloatType;
begin
  Result := 1.0 - Erf(a);
end;

// Standard normal CDF using Erf.
function CDFNormal(const z: FloatType): FloatType;
begin
  Result := 0.5 * (1.0 + Erf(z / Sqrt2));
end;

function TDistPValue(const t: FloatType; const df: Integer): FloatType;
var
  tv, x, ib: FloatType;
begin
  if df < 1 then Exit(NaN);

  // Absolute value of t.
  tv := Abs(t);

  // t = 0 -> p = 1 (two-tailed).
  if tv = 0.0 then Exit(1.0);

  // Very large df: normal approximation.
  if df > 200 then begin
    Result := 2.0 * (1.0 - CDFNormal(tv));
    if Result < 0.0 then Result := 0.0;
    if Result > 1.0 then Result := 1.0;
    Exit;
  end;

  // Transform t for incomplete beta.
  x := df / (df + tv * tv);

  // Compute lower-tail probability.
  ib := IncompleteBetaCF(x, df / 2.0, 0.5);

  // Two-tailed: use symmetry.
  if ib < 0.5 then
    Result := ib
  else
    Result := 1.0 - ib;

  // Clamp just in case.
  if Result < 0.0 then Result := 0.0;
  if Result > 1.0 then Result := 1.0;
end;

// Mutual Information.
function MutualInformation(const X, Y: WVectorType; const NumBins: Integer): FloatType;
var
  j, xj, yj, m: Integer;
  MI, MinX, MaxX, MinY, MaxY, BinWidthX, BinWidthY: FloatType;
  CountX, CountY, CountXY: IVectorType;
  PX, PY, PXY: WVectorType;
begin
  m := Length(X.Value);
  // Basic validation.
  if (m <> Length(Y.Value)) or (m < 2) then begin
    AddToErrorStack('Error on Mutual Information: length of vectors < 2.');
    Exit(NaN);
  end;
  if NumBins <= 0 then begin
    AddToErrorStack('Error on Mutual Information: nbins < 1.');
    Exit(NaN);
  end;

  for j := 0 to m-1 do begin
    if IsNan(X.Value[j]) or IsInfinite(X.Value[j]) or
       IsNan(Y.Value[j]) or IsInfinite(Y.Value[j]) then begin
      AddToErrorStack('Error on Mutual Information: NaN or Infinite values.');
      Exit(NaN);
    end;
  end;

  // Find min and max for binning.
  MinX := X.Value[0];
  MaxX := X.Value[0];
  MinY := Y.Value[0];
  MaxY := Y.Value[0];
  for j := 1 to m - 1 do begin
    if X.Value[j] < MinX then MinX := X.Value[j];
    if X.Value[j] > MaxX then MaxX := X.Value[j];
    if Y.Value[j] < MinY then MinY := Y.Value[j];
    if Y.Value[j] > MaxY then MaxY := Y.Value[j];
  end;
  BinWidthX := (MaxX - MinX) / NumBins;
  BinWidthY := (MaxY - MinY) / NumBins;
  SetLength(CountX, NumBins);
  SetLength(CountY, NumBins);
  SetLength(CountXY, NumBins * NumBins);

  // Bin the data.
  for j := 0 to m - 1 do begin
    xj := Trunc((X.Value[j] - MinX) / BinWidthX);
    yj := Trunc((Y.Value[j] - MinY) / BinWidthY);
    if xj >= NumBins then xj := NumBins - 1;
    if yj >= NumBins then yj := NumBins - 1;
    CountX[xj] := CountX[xj] + 1;
    CountY[yj] := CountY[yj] + 1;
    CountXY[xj * NumBins + yj] := CountXY[xj * NumBins + yj] + 1;
  end;
  SetLength(PX.Value, NumBins);
  SetLength(PY.Value, NumBins);
  SetLength(PXY.Value, NumBins * NumBins);

  // Convert counts to probabilities.
  for xj := 0 to NumBins - 1 do
    PX.Value[xj] := CountX[xj] / m;
  for yj := 0 to NumBins - 1 do
    PY.Value[yj] := CountY[yj] / m;
  for j := 0 to NumBins * NumBins - 1 do
    PXY.Value[j] := CountXY[j] / m;

  // Compute mutual information.
  MI := 0.0;
  for xj := 0 to NumBins - 1 do
    for yj := 0 to NumBins - 1 do begin
      if PXY.Value[xj * NumBins + yj] > 0 then
        MI := MI + PXY.Value[xj * NumBins + yj] * Ln(PXY.Value[xj * NumBins + yj] /
          (PX.Value[xj] * PY.Value[yj]));
    end;
  Result := MI; // In nats.
  CleanUpZero(Result);
end;

{function NumBinsFromSturges(const X: WVectorType): Integer;
var
  n: Integer;
begin
  n := Length(X.Value);
  if n <= 1 then Exit(1);

  Result := Ceil(Log2(n) + 1.0);
  if Result < 1 then Result := 1;
end;

// Entropy function, need bins.
function Entropy(const X: WVectorType; const NumBins: Integer): FloatType;
var
  m, i, bin: Integer;
  minX, maxX, width, p, H: FloatType;
  Counts: IVectorType;
begin
  m := Length(X.Value);
  if m = 0 then Exit(0.0);

  // Find range.
  minX := X.Value[0];
  maxX := X.Value[0];
  for i := 1 to m - 1 do begin
    if X.Value[i] < minX then minX := X.Value[i];
    if X.Value[i] > maxX then maxX := X.Value[i];
  end;

  // Allocate bins.
  SetLength(Counts, NumBins);
  for i := 0 to NumBins - 1 do
    Counts[i] := 0;

  width := (maxX - minX) / NumBins;
  if width = 0 then Exit(0.0);

  // Count values.
  for i := 0 to m - 1 do begin
    bin := Trunc((X.Value[i] - minX) / width);
    if bin >= NumBins then bin := NumBins - 1;   // Edge case.
    Inc(Counts[bin]);
  end;

  // Compute entropy.
  H := 0.0;
  for i := 0 to NumBins - 1 do begin
    if Counts[i] > 0 then begin
      p := Counts[i] / m;
      H := H - p * Ln(p);
    end;
  end;

  Result := H;
end;}

// Normalized Mutual Information.
function NormalizedMutualInformation(const X, Y: WVectorType; const NumBins: Integer): FloatType;
var
  I, HX, HY: FloatType;
begin
  I := MutualInformation(X, Y, NumBins);    // Existing MI function.
  HX := ShannonEntropy(X, EntropyAlgo);
  HY := ShannonEntropy(Y, EntropyAlgo);

  // NMI using average marginal entropies.
  if (HX + HY) = 0 then
    Result := 0.0
  else
    Result := 2.0 * I / (HX + HY);
  CleanUpBoth(Result);
end;

{ Hoeffding }
{ Shuffle: Fisher–Yates full permutation of 0..n-1 in place }
procedure ShuffleIndices(var idx: IVectorType);
var
  i, r, tmp: Integer;
begin
  Randomize;
  for i := High(idx) downto 1 do begin
    r := Random(i + 1); // 0..i
    tmp := idx[i];
    idx[i] := idx[r];
    idx[r] := tmp;
  end;
end;

{ Permutation test for Hoeffding D
  - Col1, Col2: input WVectorType (same length)
  - Cutoff, SampleSize: forwarded to HoeffdingD for speed/sample mode
  - NumPermutations: number of Monte Carlo permutations (e.g. 1000)
  - Seed: if >=0, sets RandSeed for reproducible results (optional) -- Not used.
  Returns: two-sided p-value. Also returns the observed D in OutObservedD (if non-nil).}
function HoeffdingsDSignificance(const Col1, Col2: WVectorType;
  Cutoff, SampleSize, NumPermutations: Integer; out OutObservedD: FloatType): FloatType;
var
  n, i, k,   countExtreme: Integer;
  obsD, dperm: FloatType;
  idx: IVectorType;
  Col2Perm: WVectorType;
begin
  Result := NaN;
  OutObservedD := NaN;

  // Basic checks.
  n := Length(Col1.Value);
  if (n <> Length(Col2.Value)) or (n < 5) then Exit;
  if NumPermutations <= 0 then Exit;

  // Compute observed D.
  obsD := HoeffdingD(Col1, Col2, Cutoff, SampleSize);
  OutObservedD := obsD;

  // Prepare index array and permuted Col2 record.
  SetLength(idx, n);
  for i := 0 to n - 1 do
    idx[i] := i;
  Col2Perm := Col2;     // Ensure Name and Value size exists for building permuted version.

  countExtreme := 0;

  for k := 1 to NumPermutations do begin
    // Shuffle indices (gives a random permutation.
    ShuffleIndices(idx);

    // Build permuted Col2Perm.Value using shuffled indices.
    SetLength(Col2Perm.Value, n);
    for i := 0 to n - 1 do
      Col2Perm.Value[i] := Col2.Value[idx[i]];

    // Compute D for permuted data (forward Cutoff/SampleSize along, but usually use exact here).
    dperm := HoeffdingD(Col1, Col2Perm, Cutoff, SampleSize);

    // Two-sided counting on absolute values.
    if Abs(dperm) >= Abs(obsD) then
      Inc(countExtreme);
  end;

  // Monte Carlo p-value with +1 correction.
  Result := (countExtreme + 1) / (NumPermutations + 1);
end;

// Quicksort for TValueIndex (ascending by Value). For Hoeffdings D.
procedure QuickSortValueIndex(var A: array of TValueIndex; L, R: Integer);
var
  i, j: Integer;
  pivot, tmp: TValueIndex;
begin
  if L >= R then Exit;
  i := L; j := R;
  pivot := A[(L + R) shr 1];
  while i <= j do begin
    while A[i].Value < pivot.Value do Inc(i);
    while A[j].Value > pivot.Value do Dec(j);
    if i <= j then begin
      tmp := A[i]; A[i] := A[j]; A[j] := tmp;
      Inc(i); Dec(j);
    end;
  end;
  if L < j then QuickSortValueIndex(A, L, j);
  if i < R then QuickSortValueIndex(A, i, R);
end;

// Midrank tie-corrected ranks for an array of doubles.
// Returns R as 1-based mid-ranks (type Double). For Hoeffdings D.
procedure ComputeMidRanks(const Data: CVectorType; out R:  CVectorType);
var
  n, i, j, startIdx, endIdx: Integer;
  A: array of TValueIndex;
  midrank: FloatType;
begin
  n := Length(Data);
  SetLength(R, n);
  if n = 0 then Exit;
  SetLength(A, n);

  for i := 0 to n - 1 do begin
    A[i].Value := Data[i];
    A[i].Index := i;
  end;

  QuickSortValueIndex(A, 0, n - 1);
  i := 0;
  while i < n do begin
    startIdx := i;
    endIdx := i;
    while (endIdx + 1 < n) and (Abs(A[endIdx + 1].Value - A[startIdx].Value) <= HighTolerance) do
      Inc(endIdx);

    // Midrank (ranks are 1..n).
    midrank := (startIdx + endIdx + 2) / 2.0;

    for j := startIdx to endIdx do
      R[A[j].Index] := midrank;

    i := endIdx + 1;
  end;
end;

// Joint ranks Q: Q[i] = 1 + count{ j : X[j] < X[i] and Y[j] < Y[i] }.
// Strict < used for both coordinates, as in common implementations. For Hoeffdings D.
procedure ComputeJointRanks_Q(const Xvals, Yvals: CVectorType; out Q: IVectorType);
var
  n, i, j: Integer;
begin
  n := Length(Xvals);
  SetLength(Q, n);
  for i := 0 to n - 1 do begin
    Q[i] := 1;     // Starts at 1.
    for j := 0 to n - 1 do
      if (Xvals[j] + HighTolerance < Xvals[i]) and (Yvals[j] + HighTolerance < Yvals[i]) then
        Inc(Q[i]);
  end;
end;

// Sample indices without replacement (Fisher-Yates partial shuffle). For Hoeffdings D.
procedure SampleIndicesWithoutReplacement(nTotal, SampleSize: Integer; out idx:  IVectorType);
var
  i, r, tmp: Integer;
  pool: IVectorType;
begin
  if SampleSize >= nTotal then begin
    SetLength(idx, nTotal);
    for i := 0 to nTotal - 1 do
      idx[i] := i;
    Exit;
  end;

  SetLength(pool, nTotal);
  for i := 0 to nTotal - 1 do
    pool[i] := i;

  // Shuffle first SampleSize entries.
  // Seed uses global Random state; caller may call RandSeed/Randomize if reproducibility needed
  for i := 0 to SampleSize - 1 do begin
    r := Random(nTotal - i) + i;
    tmp := pool[i];
    pool[i] := pool[r];
    pool[r] := tmp;
  end;

  SetLength(idx, SampleSize);
  for i := 0 to SampleSize - 1 do
    idx[i] := pool[i];
end;

// Hoeffdings D Function.
function HoeffdingD(const Col1, Col2: WVectorType; Cutoff: Integer; SampleSize: Integer): FloatType;
var
  nOrig, n, i: Integer;                          // Cutoff indicates when to switch to sampling.
  UseSampling: Boolean;                          // SampleSize is the size of the sample.
  Xvals, Yvals, Rx, Ry: CVectorType;
  Q, idx: IVectorType;
  D1, D2, D3, denom: FloatType;
begin
  Result := NaN;

  // Basic checks.
  nOrig := Length(Col1.Value);
  if (nOrig <> Length(Col2.Value)) or (nOrig < 5) then begin
    AddToErrorStack('Error on Hoeffding''s D: lengths of rows do not match or are insufficient');
    Exit(NaN);
  end;

  UseSampling := (Cutoff > 0) and (nOrig > Cutoff);

  if useSampling then begin
    if SampleSize < 5 then Exit;  // Require at least 5 in a sample.
    if SampleSize > nOrig then
      SampleSize := nOrig;
    SampleIndicesWithoutReplacement(nOrig, SampleSize, idx);
    n := Length(idx);
    SetLength(Xvals, n);
    SetLength(Yvals, n);
    for i := 0 to n - 1 do begin
      Xvals[i] := Col1.Value[idx[i]];
      Yvals[i] := Col2.Value[idx[i]];
    end;
  end
  else begin
    n := nOrig;
    SetLength(Xvals, n);
    SetLength(Yvals, n);
    for i := 0 to n - 1 do begin
      Xvals[i] := Col1.Value[i];
      Yvals[i] := Col2.Value[i];
    end;
  end;

  // Compute mid-ranks for margins.
  ComputeMidRanks(Xvals, Rx);
  ComputeMidRanks(Yvals, Ry);

  // Compute joint ranks Q.
  ComputeJointRanks_Q(Xvals, Yvals, Q);

  // Compute D1, D2, D3 per the documented formula (SAS / Hmisc-compatible).
  D1 := 0.0;
  D2 := 0.0;
  D3 := 0.0;

  for i := 0 to n - 1 do begin
    // Q is integer >= 1.
    D1 := D1 + (Q[i] - 1) * (Q[i] - 2);                                        // Sum (Q-1)(Q-2).
    D2 := D2 + (Rx[i] - 1.0) * (Rx[i] - 2.0) * (Ry[i] - 1.0) * (Ry[i] - 2.0);  // sum (R-1)(R-2)(S-1)(S-2)
    D3 := D3 + (Rx[i] - 2.0) * (Ry[i] - 2.0) * (Q[i] - 1);                     // Sum (R-2)(S-2)(Q-1).
  end;

  // Denom = n*(n-1)*(n-2)*(n-3)*(n-4).
  denom := n * (n - 1) * (n - 2) * (n - 3) * (n - 4);
  if denom = 0 then Exit;

  // Fnal D per SAS/Hmisc (note factor 30).
  // D = 30 * ( (n-2)*(n-3)*D1 + D2 - 2*(n-2)*D3 ).
  Result := 30.0 * ( (n - 2) * (n - 3) * D1 + D2 - 2.0 * (n - 2) * D3 ) / denom;
  CleanUpBoth(Result);
end;

// For Kendall's Tau Function.
procedure GetTieSums(const A: WVectorType; out n1, vt, t2, t3: FloatType);
var
  i, k, len, t: Integer;
  indices: IVectorType;
begin
  len := Length(A.Value);
  SetLength(indices, len);
  for i := 0 to len - 1 do
    indices[i] := i;
  QuickSortIndex(A, indices);

  n1 := 0; vt := 0; t2 := 0; t3 := 0;
  i := 0;

  while i < len do begin
    k := i;
    while (k < len) and (Abs(A.Value[indices[k]] - A.Value[indices[i]]) < HighTolerance) do
      Inc(k);
    t := k - i;
    if t > 1 then begin
      n1 := n1 + t * (t - 1) / 2.0;
      vt := vt + t * (t - 1) * (2 * t + 5);
      t2 := t2 + t * (t - 1);
      t3 := t3 + t * (t - 1) * (t - 2);
    end;
    i := k;
  end;
end;

// Kendall's Tau Function Significance.
function KendallsTauSignificance(const X, Y: WVectorType): FloatType;
var
  n, i, j, C, D: Integer;
  n1x, n1y, vtx, vty, t2x, t2y, t3x, t3y: FloatType;
  n0, varS, z, p: FloatType;
begin
  Result := NaN;
  n := Length(X.Value);
  if (n <> Length(Y.Value)) or (n < 2) then begin
    AddToErrorStack('Error on Kendall tau p-value: lengths of vectors unequal or n < 2.');
    Exit;
  end;

  // Count concordant / discordant.
  C := 0;
  D := 0;
  for i := 0 to n-2 do
    for j := i+1 to n-1 do begin
      if (Abs(X.Value[i] - X.Value[j]) < 1e-10) or
         (Abs(Y.Value[i] - Y.Value[j]) < 1e-10) then
        Continue;
      if (X.Value[j] - X.Value[i]) * (Y.Value[j] - Y.Value[i]) > 0 then
        Inc(C)
      else
        Inc(D);
    end;

  GetTieSums(X, n1x, vtx, t2x, t3x);
  GetTieSums(Y, n1y, vty, t2y, t3y);

  n0 := n * (n - 1) / 2.0;
  varS := (n0 * (2 * n + 5) - vtx - vty) / 18.0
        + (t2x * t2y) / (2.0 * n * (n - 1))
        + (t3x * t3y) / (9.0 * n * (n - 1) * (n - 2));

  if varS <= 0 then Exit(0.0);

  z := (C - D) / Sqrt(varS);
  if IsInfinite(z) or IsNan(z) then Exit(0.0);

  p := 2.0 * (1.0 - CDFNormal(Abs(z)));
  if p > 1.0 then
    p := 1.0;
  if p < 0.0 then
    p := 0.0;

  Result := p;
  CleanUpBoth(Result);
end;

// Kendall''s Tau b Function, the standard function (with ties).
function KendallsTau(const Col1, Col2: WVectorType): FloatType;
var
  i, j, m: Integer;
  concordant, discordant, tieX, tieY, bothTied: Integer;
  diff1, diff2: FloatType;
begin
  Result := NaN;
  m := Length(Col1.Value);

  // Input validation.
  if (m <> Length(Col2.Value)) or (m < 2) then begin
    AddToErrorStack('Error on Kendall''s Tau: length of vectors unequal or < 2.');
    Exit(NaN);
  end;

  concordant := 0;
  discordant := 0;
  tieX := 0;
  tieY := 0;
  bothTied := 0;

  for i := 0 to m - 2 do
    for j := i + 1 to m - 1 do begin
      diff1 := Col1.Value[i] - Col1.Value[j];
      diff2 := Col2.Value[i] - Col2.Value[j];

      if (diff1 * diff2 > 0) then
        Inc(concordant)
      else if (diff1 * diff2 < 0) then
        Inc(discordant)
      else begin
        // At least one tie
        if Abs(diff1) < 1e-10 then Inc(tieX);
        if Abs(diff2) < 1e-10 then Inc(tieY);
        if (Abs(diff1) < 1e-10) and (Abs(diff2) < 1e-10) then Inc(bothTied);
      end;
    end;

  // Adjust for double-counting in bothTied.
  tieX := tieX - bothTied;
  tieY := tieY - bothTied;

  // Denominator.
  Result := (concordant - discordant) /
    Sqrt((concordant + discordant + tieX) * (concordant + discordant + tieY));

  // Handle edge case: perfect tie in one variable.
  if IsNan(Result) or IsInfinite(Result) then
    Result := 0.0;
  CleanUpBoth(Result);
end;

// Rank a vector (average ranks for ties). Used for Spearman's Rho.
function RankVector(const X: WVectorType): WVectorType;
var
  n, i, j, k: Integer;
  idx: IVectorType;
  rank: FloatType;
  sorted: CVectorType;
begin
  n := Length(X.Value);
  SetLength(Result.Value, n);
  SetLength(idx, n);
  SetLength(sorted, n);

  for i := 0 to n - 1 do begin
    idx[i] := i;
    sorted[i] := X.Value[i];
  end;

  // Sort indices by X value.
  for i := 0 to n - 2 do
    for j := i + 1 to n - 1 do
      if sorted[i] > sorted[j] then begin
        Swap(sorted[i], sorted[j]);
        Swap(idx[i], idx[j]);
      end;

  // Assign average ranks.
  i := 0;
  while i < n do begin
    j := i;
    while (j < n - 1) and (sorted[j + 1] = sorted[i]) do
      Inc(j);
    rank := (i + j + 2) / 2;      // Average rank (1-based).
    for k := i to j do
      Result.Value[idx[k]] := rank;
    i := j + 1;
  end;
end;

// Exact Spearman Significance.
function ExactSpearmanP(const X, Y: WVectorType): FloatType;
var
  j, m, extreme, total: Integer;
  obsRho, rho: FloatType;
  RX, RY, Yperm: WVectorType;
  perm: IVectorType;

  procedure RandomPermuteTest(iterations: Integer = 100);
  var
    i, j, t: Integer;
  begin
    total := 0;
    extreme := 0;
    for i := 0 to iterations - 1 do begin

      // Fisher-Yates
      for j := m - 1 downto 1 do begin
        t := Random(j + 1);
        Swap(perm[j], perm[t]);
      end;
      for j := 0 to m - 1 do
        Yperm.Value[j] := RY.Value[perm[j]];

      rho := SpearmansRho(RX, Yperm);
      if Abs(rho) >= Abs(obsRho) then
        Inc(extreme);
      Inc(total);
    end;
  end;

  procedure Permute(k: Integer);
  var
    i, t: Integer;
  begin
    if k = m - 1 then begin
      for i := 0 to m - 1 do Yperm.Value[i] := RY.Value[perm[i]];
      rho := SpearmansRho(RX, Yperm);
      if Abs(rho) >= Abs(obsRho) then
        Inc(extreme);
      Inc(total);
    end
    else begin
      for i := k to m - 1 do begin
        t := perm[k]; perm[k] := perm[i]; perm[i] := t;
        Permute(k + 1);
        t := perm[k]; perm[k] := perm[i]; perm[i] := t;
      end;
    end;
  end;

begin
  m := Length(X.Value);
  if (m <> Length(Y.Value)) or (m < 3) then begin
    AddToErrorStack('Error on Exact Spearman''s Rhp Significance: length of vectors unequal or length < 3.');
    Exit(NaN);
  end;

  // Rank and compute observed rho.
  RX := RankVector(X);
  RY := RankVector(Y);
  obsRho := SpearmansRho(RX, RY);

  SetLength(perm, m);
  for j := 0 to m - 1 do
    perm[j] := j;
  SetLength(Yperm.Value, m);

  // Choose method.
  if m <= 9 then
    Permute(0)                    // Exact
  else if m <= 20 then
    RandomPermuteTest(10000)      // High precision.
  else
    Result := AsymptoticSpearmanP(m, obsRho);    // Fast

  if total > 0 then
    Result := extreme / total;
  CleanUpBoth(Result);
end;

// Spearman's Rho Significance, Asymptotic p-value for large n.
function AsymptoticSpearmanP(const m: Integer; const rho: FloatType): FloatType;
var
  t, p: FloatType;
begin
  t := rho * Sqrt((m - 2) / (1 - Sqr(rho)));
  p := 2 * (1 - 0.5 * (1 + Erf(Abs(t) / Sqrt2)));
  Result := p;
end;

// Spearman's Rho Significance.
function SpearmansRhoSignificance(const X, Y:WVectorType; const rho: Float): FloatType;
var
  m: Integer;
  p: FloatType;
begin
  m := Length(X.Value);
  // Basic validation.
  if IsEffectivelyEqual(rho, 1.0) then Exit(0.0);
  if IsEffectivelyEqual(rho, 0.0) then Exit(1.0);
  //  Rho := SpearmansRho(X, Y).
  if m <= 30 then
    p := ExactSpearmanP(X, Y)
  else
    p := AsymptoticSpearmanP(m, rho);
  Result := p;
  CleanUpBoth(Result);
end;

// Spearman's Rho.
function SpearmansRho(const x, y: WVectorType): FloatType;
var
  rx, ry: WVectorType;
  m: Integer;
begin
  m := Length(x.Value);
  // Basic validation.
  if m <> Length(y.Value) then begin
    AddToErrorStack('Error in Spearman''s Rho: unequal vectors.');
    Exit(NaN);
  end;
  if m < 2 then begin
    AddToErrorStack('Error in Spearman''s Rho: length of vector < 2.');
    Exit(NaN);
  end;

  // Compute ranks once (RankVector should return average ranks for ties).
  rx := RankVector(x);
  ry := RankVector(y);

  // If either rank vector has zero variance, rho undefined.
  if (Abs(Math.Variance(rx.Value)) <= 0.0) or (Abs(Math.Variance(ry.Value)) <= 0.0) then begin
    AddToErrorStack('Error in Spearman''s Rho: variance undefined.');
    Exit(NaN);
  end;

  // Use Pearson's R on the rank vectors.
  Result := PearsonsR(rx, ry);
  CleanUpBoth(Result);
end;

// Covariance Function.
function Covariance(const x, y: WVectorType): FloatType;
var
  i, n: Integer;
  meanX, meanY, sumXY: FloatType;
begin
  n := Length(x.Value);
  // Validation.
  if (n <> Length(Y.Value)) or (n < 2) then begin
    AddToErrorStack('Error on Covariance: length of vectors unequal or length < 2.');
    Exit(NaN);
  end;

  // NaN / Infinite guard.
  meanX := 0.0; meanY := 0.0;
  for i := 0 to n - 1 do begin
    if IsNan(X.Value[i]) or IsInfinite(X.Value[i]) or
      IsNan(Y.Value[i]) or IsInfinite(Y.Value[i]) then begin
      AddToErrorStack('Error on Covariance: NaN or Infinite values.');
      Exit(NaN);
    end;
    meanX := meanX + X.Value[i];
    meanY := meanY + Y.Value[i];
  end;
  MeanX := 0.0;
  MeanY := 0.0;

  for i := 0 to n - 1 do begin
    MeanX := MeanX + x.Value[i];
    MeanY := MeanY + y.Value[i];
  end;

  MeanX := MeanX / n;
  MeanY := MeanY / n;

  // Compute covariance.
  SumXY := 0.0;
  for i := 0 to n - 1 do
    SumXY := SumXY + (x.Value[i] - MeanX) * (y.Value[i] - MeanY);
  Result := SumXY / (n - 1);
end;

// Compute significance for Pearson's r.
function PearsonsRSignificance(const r: FloatType; const n: Integer): FloatType;
var
  t, p: FloatType;
begin
  if Abs(r) = 1.0 then Exit(0.0);
  if (n < 3) or (Abs(r) > 1.0) then begin
    AddToErrorStack('Error on Pearson''s R Significance: length of vectors < 3 or R > 1.');
    Exit(NaN);
  end;

  // Student's t statistic.
  t := r * Sqrt((n - 2) / (1.0 - r * r));

  // Two-tailed p-value from t-distribution.
  // Use TDistPValue function that calls IncBetaCF, etc.
  p := TDistPValue(t, n - 2);
  Result := p;
  CleanUpBoth(Result);
end;

// Pearson's R Function.
function PearsonsR(const x, y: WVectorType): FloatType;
var
  j, m, n: Integer;
  sumX, sumY, sumXY, sumX2, sumY2: FloatType;
begin
  m := Length(x.Value);
  n := Length(y.Value);

  if (m <= 1) or (n <= 1) then begin
    AddToErrorStack('Error on Pearson''s R: length of vector < 1.');
    Exit(0.0);
  end;
  if (m <> n) then begin
    AddToErrorStack('Error on Pearson''s R: lengths of vectors unequals.');
    Exit(0.0);
  end;

  // Initialize sums.
  sumX := 0.0;
  sumY := 0.0;
  sumXY := 0.0;
  sumX2 := 0.0;
  sumY2 := 0.0;

  // Calculate sums.
  for j := 0 to m do begin
    sumX := sumX + x.Value[j];
    sumY := sumY + y.Value[j];
    sumXY := sumXY + x.Value[j] * y.Value[j];
    sumX2 := sumX2 + x.Value[j] * x.Value[j];
    sumY2 := sumY2 + y.Value[j] * y.Value[j];
  end;

  // R = (n * Σxy - Σx * Σy) / sqrt((n * Σx² - (Σx)²) * (n * Σy² - (Σy)²)).
  if ((m * sumX2) = (SumX * SumX)) or ((m * sumY2) = (SumY * SumY)) then begin
    AddToErrorStack('Error on Pearson''s R: vector has no variance.');
    Exit(0.0);
  end;

  Result := (m * sumXY - sumX * sumY) /
    Sqrt((m * sumX2 - sumX * sumX) * (m * sumY2 - sumY * sumY));
  CleanUpBoth(Result);
end;

// Create covariance matrix.
procedure CreateCovarianceMatrix(const Data: WMatrixType; out CovMatrix: CMatrixType);
var
  m, n, i, j, k: Integer;
  Mean: WVectorType;
  Cov: FloatType;
begin
  // Dimensions: n = number of variables, m = number of observations.
  n := Length(Data);
  if n = 0 then begin
    AddToErrorStack('Error on Covariance: empty data.');
    CovMatrix := nil;
    Exit;
  end;

  m := Length(Data[0].Value);
  if m = 0 then begin
    AddToErrorStack('Error on Covariance: no observations.');
    CovMatrix := nil;
    Exit;
  end;

  // Validate that all vectors have equal length
  for i := 0 to n - 1 do
    if Length(Data[i].Value) <> m then begin
      AddToErrorStack('Error on Covariance: length of vectors unequal.');
      CovMatrix := nil;
      Exit;
    end;

  // Allocate arrays
  SetLength(Mean.Value, n);
  SetLength(CovMatrix, n, n);

  // Step 1: Compute means for each variable.
  for i := 0 to n - 1 do begin
    Mean.Value[i] := 0.0;
    for k := 0 to m - 1 do
      Mean.Value[i] := Mean.Value[i] + Data[i].Value[k];
    Mean.Value[i] := Mean.Value[i] / m;
  end;

  // Step 2: Compute covariance matrix.
  for i := 0 to n - 1 do
    for j := i to n - 1 do begin
      Cov := 0.0;
      for k := 0 to m - 1 do
        Cov := Cov + (Data[i].Value[k] - Mean.Value[i]) *
                     (Data[j].Value[k] - Mean.Value[j]);
      if m > 1 then
        Cov := Cov / (m - 1)
      else
        Cov := 0.0;

      CovMatrix[i, j] := Cov;
      CovMatrix[j, i] := Cov;     // Symmetric.
    end;
end;

// Create correlationb matrix.
procedure CreateCorrelationMatrix(const Data: WMatrixType; out CorrMatrix: CMatrixType);
var
  m, n, i, j, k: Integer;
  Mean, StdDev: WVectorType;
  Cov, tmp: FloatType;
begin
  // Basic validation.
  n := Length(Data);
  m := Length(Data[0].Value);
  if (n < 1) or (m < 1) then begin
    AddToErrorStack('Error on Correlation: no data.');
    CorrMatrix := nil;
    Exit;
  end;
  for i := 0 to n - 1 do
    if Length(Data[i].Value) <> m then begin
      AddToErrorStack('Error on Correlation: length of vectors unequal.');
      CorrMatrix := nil;
      Exit;
    end;

  // Allocate.
  SetLength(Mean.Value, n);
  SetLength(StdDev.Value, n);
  SetLength(CorrMatrix, n, n);

  // Step 1: compute column means (mean of each variable across observations).
  for i := 0 to n - 1 do begin
    Mean.Value[i] := 0.0;
    StdDev.Value[i] := 0.0;
    for j := 0 to m - 1 do begin
      Mean.Value[i] := Mean.Value[i] + Data[i].Value[j];
      CorrMatrix[i, j] := 0.0;
    end;
    Mean.Value[i] := Mean.Value[i] / m;
  end;

  // Step 2: compute (sum of squared deviations) for stddev.
  for i := 0 to n - 1 do begin
    for j := 0 to m - 1 do
      StdDev.Value[i] := StdDev.Value[i] + Sqr(Data[i].Value[j] - Mean.Value[i]);
    if m > 1 then
      StdDev.Value[i] := Sqrt(StdDev.Value[i] / (m - 1))
    else
      StdDev.Value[i] := 0.0;
  end;

  // Step 3: compute covariance and fill correlation matrix.
  for i := 0 to n - 1 do
    for j := i to n - 1 do begin
      Cov := 0.0;
      for k := 0 to m - 1 do
        Cov := Cov + (Data[i].Value[k] - Mean.Value[i]) * (Data[j].Value[k] - Mean.Value[j]);
      if m > 1 then
        Cov := Cov / (m - 1)
      else
        Cov := 0.0;

      // If either stddev is (near) zero, define correlation by convention:
      // Diagonal = 1.0 for constant variable, off-diagonal = 0.0.
      if (StdDev.Value[i] > 0.0) and (StdDev.Value[j] > 0.0) then
        tmp := Cov / (StdDev.Value[i] * StdDev.Value[j])
      else if i = j then
        tmp := 1.0
      else
        tmp := 0.0;

      // Guard against tiny numerical overshoot.
      if tmp > 1.0 then tmp := 1.0;
      if tmp < -1.0 then tmp := -1.0;
      CleanUpBoth(tmp);

      CorrMatrix[i, j] := tmp;
      CorrMatrix[j, i] := tmp;    // Symmetry.
    end;
end;

// Easy creation of correlation matrix if have covariance matrix.
procedure CreateCorrelationFromCovarianceMatrix(const CovMatrix: CMatrixType; out CorrMatrix: CMatrixType);
var
  n, i, j: Integer;
  StdDev: WVectorType;
  v: FloatType;
begin
  // Validate square matrix.
  n := Length(CovMatrix);
  if n = 0 then begin
    AddToErrorStack('Error on Correlation: no data.');
    CorrMatrix := nil;
    Exit;
  end;
  for i := 0 to n - 1 do
    if Length(CovMatrix[i]) <> n then begin
      AddToErrorStack('Error on Correlation: lengths of vectors unequal.');
      CorrMatrix := nil;
      Exit;
    end;

  // Allocate outputs.
  SetLength(CorrMatrix, n, n);
  SetLength(StdDev.Value, n);

  // Compute standard deviations defensively.
  for i := 0 to n - 1 do begin
    v := CovMatrix[i, i];
    // If numerical noise produces a tiny negative, clamp to zero.
    if IsNan(v) or IsInfinite(v) then
      StdDev.Value[i] := NaN
    else if v < 0 then begin
      if Abs(v) <= HighTolerance then
        StdDev.Value[i] := 0.0
      else
        StdDev.Value[i] := NaN;   // Clearly invalid covariance diagonal.
    end
    else
      StdDev.Value[i] := Sqrt(v);
  end;

  // Build correlation matrix.
  for i := 0 to n - 1 do
    for j := 0 to n - 1 do begin
      // Default.
      CorrMatrix[i, j] := 0.0;

      // If either stddev is NaN propagate NaN.
      if IsNan(StdDev.Value[i]) or IsNan(StdDev.Value[j]) then begin
        CorrMatrix[i, j] := NaN;
        Continue;
      end;

      // If both stddev > HighTolerance compute ratio; otherwise leave as 0 (or NaN by policy).
      if (StdDev.Value[i] > HighTolerance) and (StdDev.Value[j] > HighTolerance) then begin
        CorrMatrix[i, j] := CovMatrix[i, j] / (StdDev.Value[i] * StdDev.Value[j]);

        // Fix numerical overshoot into [-1,1].
        if CorrMatrix[i, j] > 1.0 then
          CorrMatrix[i, j] := 1.0
        else if CorrMatrix[i, j] < -1.0 then
          CorrMatrix[i, j] := -1.0;
      end
      else begin
        // If both variables have (near) zero variance they are constant.
        // By convention set correlation = 1 on diagonal, 0 off-diagonal.
        if (i = j) then
          CorrMatrix[i, j] := 1.0
        else
          CorrMatrix[i, j] := 0.0;
      end;
    end;

  // Ensure exact unit diagonal where possible.
  for i := 0 to n - 1 do
    if not IsNan(CorrMatrix[i, i]) then
      CorrMatrix[i, i] := 1.0;
end;

// Various kinds of multivariate correlations. Pearson, Spearman, Kendall, Hoeffding, MI.
procedure MultivariateCorrelations(const WVar: IVectorType);
var
  i, j, nVar: Integer;
  VData: WMatrixType;
  CData, SData: CMatrixType;
  FileName: String;
begin
  nVar := Length(WVar);
  SetLength(CData, nVar, nVar);
  SetLength(SData, nVar, nVar);

  // Create VData from WData, using only the variables in WVar.
  SetLength(VData, nVar);
  for i := 0 to nVar - 1 do begin
    SetLength(VData[i].Value, Length(WData[WVar[i]].Value));
    DeepCopyVector(WData[WVar[i]], VData[i]);
  end;

  case Corr of
    Pearson: begin
      for i := 0 to nVar - 1 do for j := 0 to nVar - 1 do begin
        CData[i, j] := PearsonsR(VData[i], VData[j]);
        SData[i, j] := PearsonsRSignificance(CData[i, j], Length(VData[i].Value));
      end;
      DisplayCorrelationMatrix(Corr, CData, SData, WVar);
    end;
    Spearman: begin
      for i := 0 to nVar - 1 do for j := 0 to nVar - 1 do begin
        CData[i, j] := SpearmansRho(VData[i], VData[j]);
        SData[i, j] := SpearmansRhoSignificance(VData[i], VData[j], CData[i, j]);
      end;
      DisplayCorrelationMatrix(Corr, CData, SData, WVar);
    end;
    Kendall: begin
      for i := 0 to nVar - 1 do for j := 0 to nVar - 1 do begin
        CData[i, j] := KendallsTau(VData[i], VData[j]);
        SData[i, j] := KendallsTauSignificance(VData[i], VData[j]);;
      end;
      DisplayCorrelationMatrix(Corr, CData, SData, WVar);
    end;
    Hoeffding: begin
      for i := 0 to nVar - 1 do for j := 0 to nVar - 1 do begin
        CData[i, j] := HoeffdingD(VData[i], VData[j], 1000, 500);
        SData[i, j] := HoeffdingsDSignificance(VData[i], VData[j], 1000, 500, 1000, CData[i, j]);
      end;
      DisplayCorrelationMatrix(Corr, CData, SData, WVar);
    end;
    MI: begin
      for i := 0 to nVar - 1 do for j := 0 to nVar - 1 do begin
        CData[i, j] := MutualInformation(VData[i], VData[j], Trunc(Sqrt(Length(VData[i].Value))));
      end;
      DisplayCorrelationMatrix(Corr, CData, SData, WVar);
    end;
    NMI: begin
      for i := 0 to nVar - 1 do for j := 0 to nVar - 1 do begin
        CData[i, j] := NormalizedMutualInformation(VData[i], VData[j], Trunc(Sqrt(Length(VData[i].Value))));
      end;
      DisplayCorrelationMatrix(Corr, CData, SData, WVar);
    end;
    Covar: begin
      for i := 0 to nVar - 1 do for j := 0 to nVar - 1 do begin
        CData[i, j] := Covariance(VData[i], VData[j]);
      end;
      DisplayCovarianceMatrix(CData, WVar);
    end;
  end;

  // Save files.
  if Corr <> Covar then begin
    write('Save correlations and significances?) (y/n)');
    Readln(WInput);
    if UpCase(WInput) = 'Y' then begin
      write('Name of file for correlations: ');
      readln(FileName);
      SaveFile(CData, FileName);
      write('Name of file for significances: ');
      readln(FileName);
      SaveFile(SData, FileName);
    end;
  end
  else begin
    write('Save covariances?) (y/n)');
    Readln(WInput);
    if UpCase(WInput) = 'Y' then begin
      write('Name of file for covariances: ');
      readln(FileName);
      SaveFile(CData, FileName);
    end;
  end;
end;
end.
