UNIT Deskript;

{ descriptive statistics on vectors and matrices }

INTERFACE

USES Math, mathfunc, Vector, Matrix;

CONST
  DeskriptError: BOOLEAN = FALSE;

{ ************************** Position ******************************* }

FUNCTION ArithmeticMean(Data.Value: VectorTyp): float;

FUNCTION GeometricMean(Data.Value: VectorTyp): float;

FUNCTION HarmonicMean(Data.Value: VectorTyp): float;

FUNCTION GeneralMean(Data.Value, Gewichte: VectorTyp; Exponent: float): float;

{ ********************************** Scale ****************************** }

FUNCTION Gini(Data.Value: VectorTyp; Mean: float): float;

FUNCTION Covariance(CONST X, Y: VectorTyp): float;

FUNCTION MeanDeviationFromMean(Data.Value: VectorTyp): float;

FUNCTION MedianDeviationFromMean(Data.Value: VectorTyp): float;

FUNCTION Variance(Data.Value: VectorTyp): float;

FUNCTION StandardDeviation(VAR Data.Value: VectorTyp): double;

FUNCTION CoefficientOfVariation(Mean, v: float; n: WORD): float;

{ ***************************** Moments *********************************** }

FUNCTION Mue(CONST Data.Value: VectorTyp; Mean: float; k: WORD): float;

FUNCTION Skewness(CONST Data.Value: VectorTyp; Mean, StaDev: float): float;

FUNCTION ExcessKurtosis(CONST Data.Value: VectorTyp; Mean, StaDev: float): float;

{ ***************************** Grouped Data.Value ************************ }

FUNCTION WeightedMean(Means, Vars, Lengths: VectorTyp): float;

FUNCTION WeightedStandardDeviation(Means, Vars, Lengths: VectorTyp): float;

{ ***************************** Concentration *************************** }

FUNCTION LorenzMuenzner(Data.Value: VectorTyp; VAR xVektor, yVektor: VectorTyp): float;

FUNCTION HerfindahlIndex(Data.Value: VectorTyp): float;

{ *********************** non-parametric position ************************* }

FUNCTION HiMed(CONST SortedData.Value: VectorTyp): float;

FUNCTION LoMed(CONST SortedData.Value: VectorTyp): float;

FUNCTION Median(VAR Data.Value: VectorTyp): float;

FUNCTION WeightedHiMed(CONST SortedData.Value: MatrixTyp): float;

FUNCTION WeightedLoMed(CONST SortedData.Value: MatrixTyp): float;

FUNCTION WeightedMedian(VAR Data.Value: MatrixTyp): float;

FUNCTION Quantile(VAR Data.Value: VectorTyp; q: float): float;

FUNCTION TriMedian(Data.Value: VectorTyp): float;

FUNCTION HodgesLehmann(Data.Value: VectorTyp): float;

FUNCTION NaiveHodgesLehmann(Data.Value: VectorTyp): float;

{ *********************** non-parametric scale ************************* }

FUNCTION InterQuantilDistance(Q1, Q3: float): float;

FUNCTION MAD(Data.Value: VectorTyp): double;

FUNCTION StandardErrorOfMedian(Data.Value: VectorTyp): float;

FUNCTION QuantileDispersionCoefficient(Q1, Q3: float): float;

FUNCTION Sn(VAR Data.Value: VectorTyp): float;

FUNCTION NaiveSn(VAR Data.Value: VectorTyp): float;

FUNCTION Qn(VAR Data.Value: VectorTyp): float;

FUNCTION NaiveQn(VAR Data.Value: VectorTyp): float;

function SnStatistic(const X: VectorTyp): Float;

function QnStatistic(const X: VectorTyp): Float;

{ *********************** non-parametric moments ************************* }

FUNCTION Ell2(CONST SortedData.Value: VectorTyp): float;

FUNCTION Ell3(CONST SortedData.Value: VectorTyp): float;

FUNCTION Ell4(CONST SortedData.Value: VectorTyp): float;

FUNCTION QuartileCoefficientOfSkewness(Q1, Q2, Q3: float): float;

FUNCTION CentilCoeffKurtosis(Data.Value: VectorTyp): float;

{ ********************** Standardise and normalise vector ******************** }

PROCEDURE MeanNormalise(VAR Data.Value: VectorTyp);
{ subtract arithmetic mean from all elements }

PROCEDURE Z_Standardise(VAR Data.Value: VectorTyp);
{ subtract mean and divide by standard deviation }

PROCEDURE RobustStandardise(VAR Data.Value: VectorTyp);
{ subtract median and divide by Qn }

{ ************************ Description of matrices *********************** }

PROCEDURE VarCovarMatrix(CONST Data.Value: MatrixTyp; VAR VarCovar: MatrixTyp);

PROCEDURE VarCov(CONST Data.Value: MatrixTyp; VAR VarCovar: MatrixTyp);

PROCEDURE MeanVector(CONST Data.Value: MatrixTyp; VAR Mean: VectorTyp);

PROCEDURE StaVector(CONST Data.Value: MatrixTyp; VAR Sta: VectorTyp);

{ ***************** Standardise and normalise matrix columns *************** }

PROCEDURE CentreMatrix(VAR A: MatrixTyp);

PROCEDURE StandardiseMatrix(VAR A: MatrixTyp);

PROCEDURE RobustStandardiseMatrix(VAR A: MatrixTyp);

{ *********************** Distances and outliers *********************** }

PROCEDURE MahalanobisDistance(CONST Data.Value: MatrixTyp; VAR Dm: VectorTyp);

PROCEDURE RobustDistance(CONST Data.Value: MatrixTyp; VAR Dr: VectorTyp);


IMPLEMENTATION

VAR
  ch: CHAR;

FUNCTION ArithmeticMean(Data.Value: VectorTyp): float;

VAR
  i: WORD;

BEGIN
  i := VectorLength(Data.Value);
  IF i > 0
    THEN
      Result := NeumaierSum(Data.Value) / ActualElements(Data.Value)
    ELSE
      BEGIN
        ch := WriteErrorMessage('Arithmetic mean from vector of length 0');
        DeskriptError := TRUE;
      END;
END;


FUNCTION GeometricMean(Data.Value: VectorTyp): float;

VAR
  Sum, x: float;
  n, i: WORD;

BEGIN
  n := VectorLength(Data.Value);
  IF (n = 0)
    THEN
      BEGIN
        ch := WriteErrorMessage('Geometric mean from vector of length 0');
        DeskriptError := TRUE;
        EXIT;
      END;
  CASE VectorSignum(Data.Value) OF
    -1: BEGIN
          ch := WriteErrorMessage('Geometric mean from negative Data.Value');
          DeskriptError := TRUE;
          EXIT;
        end;
    0:  Result := 0;  // product IS 0 if any of the factors is 0
    1: BEGIN
         n := 0;
         Sum := 0;
         FOR i := 1 TO VectorLength(Data.Value) DO
           BEGIN
             x := GetVectorElement(Data.Value, i);
             IF IsNaN(x)
               THEN
               ELSE
                 BEGIN
                   INC(n);
                   Sum := Sum + Ln(x);
                 END;
             Result := Exp(Sum / n);
           END;
       END;
  END; { case }
END;


FUNCTION HarmonicMean(Data.Value: VectorTyp): float;

VAR
  Sum, x: float;
  i, j, n: WORD;

BEGIN
  n := VectorLength(Data.Value);
  IF (n = 0)
    THEN
      BEGIN
        CH := WriteErrorMessage('Harmonic mean from vector of length 0');
        DeskriptError := TRUE;
        EXIT;
      END;
  CASE VectorSignum(Data.Value) OF
    -1: BEGIN
          WriteErrorMessage('Harmonic mean from negative Data.Value');
          DeskriptError := TRUE;
          EXIT;
        end;
     0: Result := 0;
     1: BEGIN
          Sum := 0;
          j := 0;
          n := VectorLength(Data.Value);
          FOR i := 1 TO n DO
            BEGIN
              x := GetVectorElement(Data.Value, i);
              IF IsNaN(x)
                THEN // ignore
                ELSE
                  BEGIN
                    Sum := Sum + 1 / x;
                    INC(j);
                  END;
            END;
          Result := j / Sum;
        END;
  END; // CASE
END;



FUNCTION GeneralMean(Data.Value, Gewichte: VectorTyp; Exponent: float): float;

VAR
  i, n: WORD;
  Sum: float;

BEGIN
  Sum := 0;
  n := VectorLength(Data.Value);
  IF (n = 0)
    THEN
      BEGIN
        CH := WriteErrorMessage('General mean from vector of length 0');
        DeskriptError := TRUE;
        EXIT;
      END;
  IF (n <> VectorLength(Gewichte))
    THEN
      BEGIN
        CH := WriteErrorMessage( 'General mean: unequal vector lengths of Data.Value and weights');
        DeskriptError := TRUE;
        EXIT;
      END;
  FOR i := 1 TO n DO
    Sum := Sum + GetVectorElement(Gewichte, i) *
      pot(GetVectorElement(Data.Value, i), Exponent);
  Result := pot(Sum, 1 / Exponent);
END;



{ *************************** Dispersion ***********************************}

FUNCTION Covariance(CONST X, Y: VectorTyp): float;

VAR
  n, Count, i, j: WORD;
  CH: CHAR;
  xi, yi, xAv, yAv, s: float;

BEGIN
  n := VectorLength(X);
  IF (n = 0)
    THEN
      BEGIN
        CH := WriteErrorMessage('Covariance from vector of length 0');
        DeskriptError := TRUE;
        EXIT;
      END;
  IF (n <> VectorLength(Y))
    THEN
      BEGIN
        CH := WriteErrorMessage('Covariance of vectors of unequal length ');
        DeskriptError := TRUE;
        EXIT;
      END;
  j := 0;
  REPEAT                            // find first complete pair OF Data.Value
    INC(j);
    xi := GetVectorElement(X, j);
    yi := GetVectorElement(Y, j);
  UNTIL NOT (IsNaN(xi) OR IsNaN(yi));
  xAv := xi;
  yAv := yi;
  s := 0.0;
  Count := 1;
  FOR i := Succ(j) TO n DO
    BEGIN
      xi := GetVectorElement(X, i);
      yi := GetVectorElement(Y, i);
      IF (IsNan(xi) OR IsNaN(yi))
        THEN
        ELSE
          BEGIN
            INC(Count);
            s := s * ((Count - 2) / (Count - 1)) + (xi - xAv) * (yi - yAv) / Count;
            xAv := xAv + (xi - xAv) / Count;
            yAv := yAv + (yi - yAv) / Count;
          END;
    END;
  Result := s;
END;


FUNCTION Variance(Data.Value: VectorTyp): float;

BEGIN
  Result := Covariance(Data.Value, Data.Value);
END;

FUNCTION StandardDeviation(VAR Data.Value: VectorTyp): float;

BEGIN
  Result := Sqrt(Covariance(Data.Value, Data.Value));
END;


FUNCTION StandardErrorOfMean(v: float; n: WORD): float;

BEGIN
  IF n > 0
    THEN
      StandardErrorOfMean := Sqrt(v) / Sqrt(n)
    ELSE
      BEGIN
        WriteErrorMessage(' Standard error of mean for n = 0');
        DeskriptError := TRUE;
      END;
END;


FUNCTION CoefficientOfVariation(Mean, v: float; n: WORD): float;

BEGIN
  IF Abs(mean) > Zero
    THEN
      Result := Sqrt(v) / Abs(Mean)
    ELSE
      BEGIN
        WriteErrorMessage(' Coefficient of variation with mean 0');
        DeskriptError := TRUE;
      end;
END;


FUNCTION MeanDeviationFromMean(Data.Value: VectorTyp): float;

VAR
  Mean, Sum: float;
  i, j, n: WORD;

BEGIN
  n := VectorLength(Data.Value);
  IF (n = 0)
    THEN
      BEGIN
        CH := WriteErrorMessage('Mean deviation of mean from vector of length 0');
        DeskriptError := TRUE;
        EXIT;
      END;
  Mean := ArithmeticMean(Data.Value);
  Sum := 0;
  FOR i := 1 TO n DO
    Sum := Sum + Abs(Mean - GetVectorElement(Data.Value, i));
  Result := Sum / n;
END;


FUNCTION MedianDeviationFromMean(Data.Value: VectorTyp): float;

VAR
  Mean: float;
  Abweichungen: VectorTyp;
  i: WORD;

BEGIN
  CreateVector(Abweichungen, VectorLength(Data.Value), 0.0);
  Mean := ArithmeticMean(Data.Value);
  CreateVector(Abweichungen, VectorLength(Data.Value), 0.0);
  FOR i := 1 TO VectorLength(Data.Value) DO
    SetVectorElement(Abweichungen, i, Abs(Mean - GetVectorElement(Data.Value, i)));
  ShellSort(Abweichungen);
  MedianDeviationFromMean := Median(Abweichungen);
  DestroyVector(Abweichungen);
END;


{ *************************** Pearson Moments ****************************** }

FUNCTION Mue(CONST Data.Value: VectorTyp; Mean: float; k: WORD): float;

VAR
  Sum: float;
  i, n: WORD;

BEGIN
  Sum := 0;
  n := VectorLength(Data.Value);
  IF (n = 0)
    THEN
      BEGIN
        CH := WriteErrorMessage('Mue from vector of length 0');
        DeskriptError := TRUE;
        EXIT;
      END;
  FOR i := 1 TO n DO
    Sum := Sum + pot(GetVectorElement(Data.Value, i) - Mean, k);
  Result := Sum / n;
END;


FUNCTION Skewness(CONST Data.Value: VectorTyp; Mean, StaDev: float): float;

BEGIN
  Result := Mue(Data.Value, Mean, 3) / pot(StaDev, 3);
END;


FUNCTION ExcessKurtosis(CONST Data.Value: VectorTyp; Mean, StaDev: float): float;

BEGIN
  Result := (Mue(Data.Value, Mean, 4) / pot(StaDev, 4)) - 3.0;
END;

{ ************************************************************************* }

FUNCTION WeightedMean(Means, Vars, Lengths: VectorTyp): float;

VAR
  Sum1, Sum2: float;
  i, n: WORD;

BEGIN
  Sum1 := 0;
  Sum2 := 0;
  n := VectorLength(Means);
  IF (n = 0)
    THEN
      BEGIN
        CH := WriteErrorMessage('Weighted mean from vector of length 0');
        DeskriptError := TRUE;
        EXIT;
      END;
  FOR i := 1 TO n DO
    BEGIN
      Sum1 := Sum1 + (GetVectorElement(Lengths, i) * GetVectorElement(Means, i) /
                      GetVectorElement(Vars, i));
      Sum2 := Sum2 + (GetVectorElement(Lengths, i) / GetVectorElement(Vars, i));
    END;
  Result := Sum1 / Sum2;
END;


FUNCTION WeightedStandardDeviation(Means, Vars, Lengths: VectorTyp): float;

VAR
  Sum1, Sum2: float;
  i, n: WORD;

BEGIN
  Sum1 := 0;
  Sum2 := 0;
  n := VectorLength(Means);
  IF (n = 0)
    THEN
      BEGIN
        CH := WriteErrorMessage('Weighted standard deviation from vector of length 0');
        DeskriptError := TRUE;
        EXIT;
      END;
  FOR i := 1 TO n DO
    BEGIN
      Sum1 := Sum1 + (GetVectorElement(Lengths, i) - 1) * GetVectorElement(Vars, i);
      Sum2 := Sum2 + GetVectorElement(Lengths, i);
    END;
  Result := Sum1 / (Sum2 - n);
END;


FUNCTION Gini(Data.Value: VectorTyp; Mean: float): float;
  { mittlere Abweichung aller Data.Value voneinander }

VAR
  i, j, n: WORD;
  Sum: float;

BEGIN
  Sum := 0;
  n := VectorLength(Data.Value);
  IF (n = 0)
    THEN
      BEGIN
        CH := WriteErrorMessage('Gini coefficient from vector of length 0');
        DeskriptError := TRUE;
        EXIT;
      END;
  FOR i := 1 TO Pred(n) DO
    FOR j := Succ(i) TO n DO
      Sum := Sum + Abs(GetVectorElement(Data.Value, i) - GetVectorElement(Data.Value, j));
  Result := Sum / (2 * n * n * Mean);
END;


FUNCTION LorenzMuenzner(Data.Value: VectorTyp; VAR xVektor, yVektor: VectorTyp): float;

VAR
  n, i: WORD;
  u, v, Sum, Gesamt, SumV: float;

BEGIN
  n := VectorLength(Data.Value);
  IF (n = 0)
    THEN
      BEGIN
        CH := WriteErrorMessage('Lorenz-Muenzner from vector of length 0');
        DeskriptError := TRUE;
        EXIT;
      END;
  CreateVector(xVektor, n, 0.0);
  CreateVector(yVektor, n, 0.0);
  Gesamt := NeumaierSum(Data.Value);
  v := 0;
  Sum := 0;
  SumV := 0;
  FOR i := 1 TO n DO
    BEGIN
      u := i / n;
      Sum := Sum + GetVectorElement(Data.Value, i);
      v := Sum / Gesamt;
      SetVectorElement(xVektor, i, u);
      SetVectorElement(yVektor, i, v);
      SumV := SumV + v;
    END;
  Result := (Succ(n) - 2 * SumV) / (Pred(n));
END;


FUNCTION HerfindahlIndex(Data.Value: VectorTyp): float;

VAR
  n, i: WORD;
  Sum1, Sum2: float;

BEGIN
  n := VectorLength(Data.Value);
  IF (n = 0)
    THEN
      BEGIN
        CH := WriteErrorMessage('Herfindahl index from vector of length 0');
        DeskriptError := TRUE;
        EXIT;
      END;
  Sum1 := NeumaierSum(Data.Value);
  Sum2 := 0;
  FOR i := 1 TO n DO
    Sum2 := Sum2 + Sqr(GetVectorElement(Data.Value, i) / Sum1);
  Result := Sum2;
END;

{ ********************* non-parametric *********************************** }

FUNCTION HiMed(CONST SortedData.Value: VectorTyp): float;

VAR
  n: WORD;

BEGIN
  n := VectorLength(SortedData.Value);
  IF (n = 0)
    THEN
      BEGIN
        CH := WriteErrorMessage('Hi median from vector of length 0');
        DeskriptError := TRUE;
      END
  ELSE
    Result := GetVectorElement(SortedData.Value, Succ(n DIV 2));
END;

FUNCTION LoMed(CONST SortedData.Value: VectorTyp): float;

VAR
  n: WORD;

BEGIN
  n := VectorLength(SortedData.Value);
  IF (n = 0)
    THEN
      BEGIN
        CH := WriteErrorMessage('Low median from vector of length 0');
        DeskriptError := TRUE;
      END
    ELSE
      Result := GetVectorElement(SortedData.Value, Succ(n) DIV 2);
END;

FUNCTION Median(VAR Data.Value: VectorTyp): float;

VAR
  n: WORD;
  Sorted: VectorTyp;

BEGIN
  n := VectorLength(Data.Value);
  IF (n = 0)
    THEN
      BEGIN
        CH := WriteErrorMessage('Median from vector of length 0');
        DeskriptError := TRUE;
        EXIT;
      END;
  ShellSort(Data.Value);
  IF Odd(n)
    THEN Median := LoMed(Data.Value)  // LOW AND Hi median are identical
    ELSE Median := (LoMed(Data.Value) + HiMed(Data.Value)) / 2;
END;


FUNCTION WeightedHiMed(CONST SortedData.Value: MatrixTyp): float;

VAR
  i: WORD;
  Sum: float;

BEGIN
  Sum := 0.0;
  i := 0;
  REPEAT
    INC(i);
    sum := sum + GetMatrixElement(SortedData.Value, i, 2);
  UNTIL (sum >= 0.5);
  WeightedHiMed := GetMatrixElement(SortedData.Value, i, 1);
END;

FUNCTION WeightedLoMed(CONST SortedData.Value: MatrixTyp): float;

VAR
  i: WORD;
  Sum: float;

BEGIN
  Sum := 0.0;
  i := 0;
  REPEAT
    INC(i);
    sum := sum + GetMatrixElement(SortedData.Value, i, 2);
    IF (sum > 0.5) THEN DEC(i);
  UNTIL (sum >= 0.5);
  WeightedLoMed := GetMatrixElement(SortedData.Value, i, 1);
END;

FUNCTION WeightedMedian(VAR Data.Value: MatrixTyp): float;

VAR
  n: WORD;

BEGIN
  n := MatrixRows(Data.Value);
  IF (n = 0)
    THEN
      BEGIN
        CH := WriteErrorMessage('Weighted median from vector of length 0');
        DeskriptError := TRUE;
        EXIT;
      END;
  ShellSortMatrix(Data.Value, 1);
  WeightedMedian := (WeightedLoMed(Data.Value) + WeightedHiMed(Data.Value)) / 2.0;
END;


FUNCTION Quantile(VAR Data.Value: VectorTyp; q: float): float;

VAR
  n: WORD;

BEGIN
  n := VectorLength(Data.Value);
  IF (n = 0)
    THEN
      BEGIN
        CH := WriteErrorMessage('Quantile from vector of length 0');
        DeskriptError := TRUE;
        EXIT;
      END;
  ShellSort(Data.Value);
  Result := GetVectorElement(Data.Value, Round(n * q + 0.5)); // always Round up
END;


FUNCTION TriMedian(Data.Value: VectorTyp): float;

VAR
  n: WORD;

BEGIN
  n := VectorLength(Data.Value);
  IF (n = 0) THEN
  BEGIN
    CH := WriteErrorMessage('Trimedian from vector of length 0');
    DeskriptError := TRUE;
    EXIT;
  END;
  Result := (Quantile(Data.Value, 0.25) + 2 * Quantile(Data.Value, 0.5) +
    Quantile(Data.Value, 0.75)) / 4;
END;


FUNCTION NaiveHodgesLehmann(Data.Value: VectorTyp): float;

VAR
  n, i, j, k: WORD;
  x, y: float;
  Averages: VectorTyp;

BEGIN
  n := VectorLength(Data.Value);
  IF (n = 0)
    THEN
      BEGIN
        CH := WriteErrorMessage('Hodges-Lehmann estimator from vector of length 0');
        DeskriptError := TRUE;
        EXIT;
      END;
  n := Round(BinomialCoef(n, 2));
  CreateVector(Averages, n, 0.0);
  k := 0;
  FOR i := 1 TO VectorLength(Data.Value) DO
    BEGIN
      x := GetVectorElement(Data.Value, i);
      FOR j := 1 TO Pred(i) DO
        BEGIN
          INC(k);
          y := GetVectorElement(Data.Value, j);
          SetVectorElement(Averages, k, (x + y) / 2);
        END;
    END;
  Result := Median(Averages);
  DestroyVector(Averages);
END;


FUNCTION HodgesLehmann(Data.Value: VectorTyp): float;

VAR
  min, max, Step, Sum, x, y: float;
  Counts: VectorTyp;
  i, j, n, NrBins, Bin: WORD;

BEGIN
  n := VectorLength(Data.Value);
  IF (n = 0)
    THEN
      BEGIN
        CH := WriteErrorMessage('Hodges-Lehmann estimator from vector of length 0');
        DeskriptError := TRUE;
        EXIT;
      END;
  NrBins := Round(Sqrt(BinomialCoef(n, 2)));
  min := FindSmallest(Data.Value);
  max := FindLargest(Data.Value);
  IF (min = max)
    THEN
      BEGIN
        CH := WriteErrorMessage('Hodges-Lehmann estimator: all Data.Value are the same');
        DeskriptError := TRUE;
        EXIT;
      END;
  Step := (max - min) / NrBins;
  CreateVector(Counts, Succ(NrBins), 0.0);  // because vector starts at 1, rather than 0
  FOR i := 1 TO n DO                        // count frequency OF grouped averages
    BEGIN
      x := GetVectorElement(Data.Value, i);
      FOR j := 1 TO Pred(i) DO
        BEGIN
          y := GetVectorElement(Data.Value, j);
          Bin := Succ(Round(((x + y) / 2 - min) / Step));
          SetVectorElement(Counts, Bin, GetVectorElement(Counts, Bin) + 1.0);
        END;
    END;
  Sum := 0.0;
  i := 0;
  REPEAT                                    // identify bin OF median
    INC(i);
    Sum := Sum + GetVectorElement(Counts, i);
  UNTIL sum >= (BinomialCoef(n, 2) / 2);
  Result := Pred(i) * Step + min;
  DestroyVector(Counts);
END;


FUNCTION StandardErrorOfMedian(Data.Value: VectorTyp): float;

VAR
  i, j, n: WORD;

BEGIN
  n := VectorLength(Data.Value);
  IF (n = 0)
    THEN
      BEGIN
        CH := WriteErrorMessage('Weighted low median from vector of length 0');
        DeskriptError := TRUE;
        EXIT;
      END;
  i := Round(n / 2 - Sqrt(3 * n) / 2);
  j := Round(n / 2 + Sqrt(3 * n) / 2);
  Result := (GetVectorElement(Data.Value, j) - GetVectorElement(Data.Value, i)) / 3.4641;
END;


FUNCTION InterQuantilDistance(Q1, Q3: float): float;

BEGIN
  Result := 0.5 * (Q3 - Q1);
END;


FUNCTION QuantileDispersionCoefficient(Q1, Q3: float): float;

BEGIN
  QuantileDispersionCoefficient := (Q3 - Q1) / (Q3 + Q1);
END;


FUNCTION MAD(Data.Value: VectorTyp): double;

CONST
  Scale = 1.4826; //   1 / (Sqrt(2) * qnorm(5/8))

VAR
  m, Sum, cn: float;
  n, i: WORD;
  Dev: VectorTyp;

BEGIN
  n := VectorLength(Data.Value);
  IF (n = 0)
    THEN
      BEGIN
        CH := WriteErrorMessage('MAD from vector of length 0');
        DeskriptError := TRUE;
        EXIT;
      END;
  m := Median(Data.Value);
  CreateVector(Dev, n, 0.0);
  FOR i := 1 TO n DO
    SetVectorElement(Dev, i, Abs(GetVectorElement(Data.Value, i) - m));
  CASE n OF
    2: cn := 1.196;
    3: cn := 1.495;
    4: cn := 1.363;
    5: cn := 1.206;
    6: cn := 1.200;
    7: cn := 1.140;
    8: cn := 1.129;
    9: cn := 1.107
    ELSE
      cn := n / (n - 0.8);
  END; { case }
  Result := cn * Scale * Median(Dev);
  DestroyVector(Dev);
END;

procedure QuickSort(var A: array of Float; L, R: Word);
var
  I, J: Word;
  P, T: Float;
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

procedure SortArray(var A: array of Float);
begin
  if Length(A) > 1 then
    QuickSort(A, 0, High(A));
end;

function MedianArray(var A: array of Float): Float;
var
  B: array of Float;
  n: Word;
begin
  n := Length(A);
  SetLength(B, n);
  Move(A[0], B[0], n * SizeOf(Double));
  SortArray(B);
  if Odd(n) then
    Result := B[n div 2]
  else
    Result := 0.5 * (B[n div 2 - 1] + B[n div 2]);
end;

function SnStatistic(const X: VectorTyp): Float;

var
  n, i, j: Word;
  Dij, Mi, temp: array of Float;
  c : Float;
begin
  n := VectorLength(X);
  SetLength(Mi, n);

  // Compute m_i = median_j |x_i - x_j|
  for i := 0 to VectorLength(X) - 1 do
  begin
    SetLength(Dij, n);
    for j := 0 to VectorLength(X) - 1 do
    Dij[j] := Abs(GetVectorElement(X, i) - GetVectorElement(X, j));
    Mi[i] := MedianArray(Dij);
  end;

  // Correction constant (asymptotic)
  c := 1.1926;

  Result := c * MedianArray(Mi);
end;

function QnStatistic(const X: VectorTyp): Float;
var
  h, n, i, j, m, k: Word;
  Diffs: array of Float;
  c: Float;
begin
  n := VectorLength(X);
  m := n * (n - 1) div 2;
  SetLength(Diffs, m);

  // collect all pairwise absolute differences
  k := 0;
  for i := 0 to VectorLength(X) - 2 do
    for j := i+1 to VectorLength(X) - 1 do
    begin
      Diffs[k] := Abs(GetVectorElement(X, i) - GetVectorElement(X, j));
      Inc(k);
    end;

  // sort the differences
  SortArray(Diffs);

  // find index h (approx first quartile of distances)
  h := ( (n div 2) * ((n div 2) + 1) ) div 2;

  // correction constant (asymptotic)
  c := 2.2219;

  Result := c * Diffs[h-1];  // adjust because arrays are 0-based
end;

FUNCTION NaiveSn(VAR Data.Value: VectorTyp): float;

VAR
  i, j, n: WORD;
  a1, a2: VectorTyp;
  xi, cn: float;

BEGIN
  n := VectorLength(Data.Value);
  IF (n = 0)
    THEN
      BEGIN
        CH := WriteErrorMessage('Sn from vector of length 0');
        DeskriptError := TRUE;
        EXIT;
      END;
  CreateVector(a1, n, 0.0);
  CreateVector(a2, n, 0.0);
  FOR i := 1 TO n DO
    BEGIN
      xi := GetVectorElement(Data.Value, i);
      FOR j := 1 TO n DO
        SetVectorElement(a1, j, Abs(xi - GetVectorElement(Data.Value, j)));
      ShellSort(a1);
      SetVectorElement(a2, i, GetVectorElement(a1, Succ(n DIV 2))); // high median
    END;
  DestroyVector(a1);
  ShellSort(a2);
  CASE n OF
    2: cn := 0.743;
    3: cn := 1.851;
    4: cn := 0.954;
    5: cn := 1.351;
    6: cn := 0.993;
    7: cn := 1.198;
    8: cn := 1.005;
    9: cn := 1.131
    ELSE IF Odd(n)
           THEN cn := n / (n - 0.9)
           ELSE cn := 1;
  END; { case }
  Result := cn * 1.1926 * GetVectorElement(a2, Succ(n) DIV 2);  // LOW median
END;

FUNCTION Sn(VAR Data.Value: VectorTyp): float;

VAR
  n, i, j, l, nA, nB, rightA, rightB, leftA, leftB, tryA, tryB, diff,
  Amin, Amax, even, half: WORD;
  cn, medA, medB: float;
  a2: VectorTyp;

BEGIN
  n := VectorLength(Data.Value);
  IF (n = 0)
    THEN
      BEGIN
        CH := WriteErrorMessage('Sn from vector of length 0');
        DeskriptError := TRUE;
        EXIT;
      END;
  ShellSort(Data.Value);
  CreateVector(a2, n, 0.0);
  SetVectorElement(a2, 1, GetVectorElement(Data.Value, Succ(n DIV 2)) -
    GetVectorElement(Data.Value, 1));
  FOR i := 2 TO (Succ(n DIV 2)) DO // a2(i) = lomed_{j<>i} | xi - xj |, i = 2...Succ(n)/2  Wes Fix
    BEGIN
      nA := Pred(i);
      nB := n - i;
      diff := nB - nA;
      leftA := 1;
      leftB := 1;
      rightA := nB;
      rightB := nB;
      Amin := Succ(diff DIV 2);
      Amax := diff DIV 2 + nA;
      WHILE (leftA < rightA) DO
        BEGIN
          l := Succ(rightA - leftA);
          even := 1 - (l MOD 2);
          half := Pred(l) DIV 2;
          tryA := leftA + half;
          tryB := leftB + half;
          IF (tryA < Amin)
            THEN
              BEGIN
                rightB := tryB;
                leftA := tryA + even;
              END
            ELSE IF (tryA > Amax)
                   THEN
                     BEGIN
                       rightA := tryA;
                       leftB := tryB + even;
                     END
                   ELSE
                     BEGIN
                       medA := GetVectorElement(Data.Value, i) - GetVectorElement(Data.Value, (i - tryA + Amin - 1));
                       medB := GetVectorElement(Data.Value, tryB + i) - GetVectorElement(Data.Value, i);
                       IF (medA >= medB)
                         THEN
                           BEGIN
                             rightA := tryA;
                             leftB := tryB + even;
                           END
                         ELSE
                           BEGIN
                             rightB := tryB;
                             leftA := tryA + even;
                           END;
                     END;
        END; { while }
      IF (leftA > Amax)
        THEN
          SetVectorElement(a2, i, GetVectorElement(Data.Value, leftB + i) -
                           GetVectorElement(Data.Value, i))
        ELSE
          BEGIN
            medA := GetVectorElement(Data.Value, i) - GetVectorElement(Data.Value,
              Pred(i - leftA + Amin));
            medB := GetVectorElement(Data.Value, leftB + i) - GetVectorElement(Data.Value, i);
            SetVectorElement(a2, i, min(medA, medB));
          END;
    END; { for }
  FOR i := Succ(Succ(n) DIV 2) TO Pred(n) DO    // same, but i = Succ(Succ(n)/2) ... n-1
    BEGIN
      nA := n - i;
      nB := Pred(i);
      diff := nB - nA;
      leftA := 1;
      leftB := 1;
      rightA := nB;
      rightB := nB;
      Amin := Succ(diff DIV 2);
      Amax := diff DIV 2 + nA;
      WHILE (leftA < rightA) DO
        BEGIN
          l := Succ(rightA - leftA);
          even := 1 - (l MOD 2);  // 0 OR 1
          half := Pred(l) DIV 2;
          tryA := leftA + half;
          tryB := leftB + half;
          IF (tryA < Amin)
            THEN
              BEGIN
                rightB := tryB;
                leftA := tryA + even;
              END
            ELSE IF (tryA > Amax)
                   THEN
                     BEGIN
                       rightA := tryA;
                       leftB := tryB + even;
                     END
                   ELSE
                     BEGIN
                       medA := GetVectorElement(Data.Value, Succ(i + tryA - Amin)) -
                               GetVectorElement(Data.Value, i);
                       medB := GetVectorElement(Data.Value, i) -
                               GetVectorElement(Data.Value, i - tryB);
                      IF (medA >= medB)
                        THEN
                          BEGIN
                            rightA := tryA;
                            leftB := tryB + even;
                          END
                        ELSE
                          BEGIN
                            rightB := tryB;
                            leftA := tryA + even;
                          END;
                    END;
        END; { while }
      IF (leftA > Amax)
        THEN
          SetVectorElement(a2, i, GetVectorElement(Data.Value, i) -
                           GetVectorElement(Data.Value, i - leftB))
        ELSE
          BEGIN
            medA := GetVectorElement(Data.Value, Succ(i + leftA - Amin)) -
                    GetVectorElement(Data.Value, i);
            medB := GetVectorElement(Data.Value, i) - GetVectorElement(Data.Value, i - leftB);
            SetVectorElement(a2, i, min(medA, medB));
          END;
    END; { for }
  SetVectorElement(a2, n, GetVectorElement(Data.Value, n) -
    GetVectorElement(Data.Value, Succ(n) DIV 2));
  CASE n OF           // finite sample correction factor
    2: cn := 0.743;
    3: cn := 1.851;
    4: cn := 0.954;
    5: cn := 1.351;
    6: cn := 0.993;
    7: cn := 1.198;
    8: cn := 1.005;
    9: cn := 1.131
    ELSE
      IF Odd(n)
        THEN cn := n / (n - 0.9)
        ELSE cn := 1;
  END; { case }
  ShellSort(a2);
  Result := cn * 1.1926 * GetVectorElement(a2, Succ(n) DIV 2);
  DestroyVector(a2);
END; { Sn }

FUNCTION NaiveQn(VAR Data.Value: VectorTyp): float;

VAR
  i, j, n, k: WORD;
  a1: VectorTyp;
  xi, dn: float;

BEGIN
  n := VectorLength(Data.Value);
  IF (n = 0)
    THEN
      BEGIN
        CH := WriteErrorMessage('Qn from vector of length 0');
        DeskriptError := TRUE;
        EXIT;
      END;
  CreateVector(a1, (n - 1) * n DIV 2, 0.0);
  k := 0;
  FOR i := 1 TO n DO
    BEGIN
      xi := GetVectorElement(Data.Value, i);
      FOR j := Succ(i) TO n DO
        BEGIN
          INC(k);
          SetVectorElement(a1, k, Abs(xi - GetVectorElement(Data.Value, j)));
        END;
    END;
  ShellSort(a1);
  CASE n OF   // finite sample correction factor
    2: dn := 0.399;
    3: dn := 0.994;
    4: dn := 0.512;
    5: dn := 0.844;
    6: dn := 0.611;
    7: dn := 0.857;
    8: dn := 0.669;
    9: dn := 0.872;
    ELSE
      IF Odd(n)    // n >= 10
        THEN dn := n / (n + 1.4)
        ELSE dn := n / (n + 3.8);
  END; { case }
  xi := GetVectorElement(a1, Round(BinomialCoef(n, 2)) DIV 4);
  Result := dn * 2.2219 * xi;
  DestroyVector(a1);
END;



FUNCTION Qn(VAR Data.Value: VectorTyp): float;

TYPE
  IntVec = ARRAY [1..MaxVector] OF WORD;

VAR
  work, weight: VectorTyp;
  Left, Right, Q, P: IntVec;
  n, h, i, j, jj, k, knew, kcand, jhelp, nL, nR, sumQ, sumP: WORD;
  found: BOOLEAN;
  trial, dn: float;

  FUNCTION WeightedHiMedian(SortedData.Value, weight: VectorTyp): float;

  VAR
    acand, iwcand: VectorTyp;
    wrest, wleft, wmid, wright, i, n, nn: WORD;
    wtotal: float;

  BEGIN
    n := VectorLength(SortedData.Value);
    nn := n;
    wtotal := 0;
    CreateVector(acand, n, 0.0);
    CreateVector(iwcand, n, 0.0);
    FOR i := 1 TO nn DO
      wtotal := wtotal + GetVectorElement(weight, i);
    wrest := 0;
    WHILE TRUE DO    // infinite LOOP EXIT'ed when result known
      BEGIN
        Trial := GetVectorElement(SortedData.Value, succ(nn div 2));
        wleft := 0;
        wmid := 0;
        wright := 0;
        FOR i := 1 TO nn DO
          IF (GetVectorElement(SortedData.Value, i) < trial)
            THEN wleft := wleft + trunc(GetVectorElement(weight, i))
            ELSE IF (GetVectorElement(SortedData.Value, i) > trial)
                   THEN wright := wright + trunc(GetVectorElement(weight, i))
                   ELSE wmid := wmid + trunc(GetVectorElement(weight, i));
        IF ((2 * wrest + 2 * wleft) > wtotal)
          THEN
            BEGIN
              kcand := 0;
              FOR i := 1 TO nn DO
                IF (GetVectorElement(SortedData.Value, i) < trial)
                  THEN
                    BEGIN
                      Inc(kcand);
                      SetVectorElement(acand, kcand, GetVectorElement(SortedData.Value, i));
                      SetVectorElement(iwcand, kcand, GetVectorElement(weight, i));
                    END;
              nn := kcand;
            END
          ELSE
            BEGIN
              IF ((2 * wrest + 2 * wleft + 2 * wmid) > wtotal)  // Endpunkt wird nicht erreicht
                THEN
                  BEGIN
                    WeightedHiMedian := trial;
                    exit;
                  END
                ELSE
                  BEGIN
                    kcand := 0;
                    FOR i := 1 TO nn DO
                      IF (GetVectorElement(SortedData.Value, i) > trial)
                        THEN
                          BEGIN
                            Inc(kcand);
                            SetVectorElement(acand, kcand,
                              GetVectorElement(SortedData.Value, i));
                            SetVectorElement(iwcand, kcand,
                              GetVectorElement(weight, i));
                          END;
                    nn := kcand;
                    wrest := wrest + wleft + wmid;
                  END; { else }
            END; { else }
        FOR i := 1 TO nn DO
          BEGIN
            SetVectorElement(SortedData.Value, i, GetVectorElement(acand, i));
            SetVectorElement(weight, i, GetVectorElement(iwcand, i));
          END;
      END; { while }
    DestroyVector(acand);
    DestroyVector(iwcand);
  END;  { WeightedHiMedian }

BEGIN
  n := VectorLength(Data.Value);
  IF (n = 0)
    THEN
      BEGIN
        ch := WriteErrorMessage('Qn from vector OF Length 0');
        DeskriptError := TRUE;
        EXIT;
      END;
  h := succ(n div 2);
  k := h * pred(h) div 2;
  ShellSort(Data.Value);
  CreateVector(work, n, 0.0);
  CreateVector(weight, n, 0.0);
  FOR i := 1 TO n DO
    BEGIN
      left[i] := n - i + 2;
      right[i] := n;
    END;
  jhelp := n * succ(n) div 2;
  knew := k + jhelp;
  nL := jhelp;
  nR := sqr(n);
  found := False;
  WHILE ((nR - nL > n) and (not found)) DO
    BEGIN
      j := 1;
      FOR i := 2 TO n DO
        IF (left[i] <= right[i])
          THEN
            BEGIN
              SetVectorElement(weight, j, succ(right[i] - left[i]));
              jhelp := left[i] + trunc(GetVectorElement(weight, j)) div 2;
              SetVectorElement(work, j, GetVectorElement(Data.Value, i) -
                GetVectorElement(Data.Value, succ(n) - jhelp));
              Inc(j);
            END;
      trial := WeightedHiMedian(work, weight);
      j := 0;
      FOR i := n DOWNTO 1 DO
        BEGIN
          WHILE ((j < n) AND (GetVectorElement(Data.Value, i) - GetVectorElement(Data.Value, n - j)
               < trial)) DO
            Inc(j);
          P[i] := j;
        END;
      j := succ(n);
      FOR i := 1 TO n DO
        BEGIN
          WHILE (GetVectorElement(Data.Value, i) - GetVectorElement(Data.Value, n - j + 2) > trial) DO
            Dec(j);
          Q[i] := j;
        End;
      sumP := 0;
      sumQ := 0;
      FOR i := 1 TO n DO
        BEGIN
          sumP := sumP + P[i];
          sumQ := sumQ + pred(Q[i]);
        END;
      IF (knew <= sumP)     // problem about here
        THEN
          BEGIN
            FOR i := 1 TO n DO
              right[i] := P[i];
            nR := sumP;
          END
        ELSE IF (knew > sumQ)
               THEN
                 BEGIN
                   FOR i := 1 TO n DO
                     left[i] := Q[i];
                   nL := sumQ;
                 END
               ELSE
                 BEGIN
                   Qn := trial;
                   found := True;
                 END;
    END; { WHILE }
  IF NOT (found)
    THEN
      BEGIN
        j := 1;
        FOR i := 2 TO n DO
          IF (left[i] <= right[i])
            THEN
              FOR jj := left[i] TO right[i] DO
                BEGIN
                  SetVectorElement(work, j, GetVectorElement(Data.Value, i) -
                    GetVectorElement(Data.Value, succ(n - jj)));
                  Inc(j);  // !!!! grows bejond n !!!!
                END;
      END;
  CASE n OF   // finite sample correction factor
    2: dn := 0.399;
    3: dn := 0.994;
    4: dn := 0.512;
    5: dn := 0.844;
    6: dn := 0.611;
    7: dn := 0.857;
    8: dn := 0.669;
    9: dn := 0.872;
    ELSE
      IF odd(n)
        THEN dn := n / (n + 1.4)
        ELSE dn := n / (n + 3.8);
  END; { case }
  ShellSort(work);
  Result := dn * 2.2219 * GetVectorElement(work, knew - nL);
  DestroyVector(work);
  DestroyVector(weight);
END;

{ **************************** L-moments ******************************** }

FUNCTION Ell2(CONST SortedData.Value: VectorTyp): float;

VAR
  sum, factor: float;
  i, n: WORD;

BEGIN
  n := VectorLength(SortedData.Value);
  IF (n = 0)
    THEN
      BEGIN
        CH := WriteErrorMessage('Ell2 from vector OF Length 0');
        DeskriptError := TRUE;
        EXIT;
      END;
  sum := 0;
  FOR i := 1 TO n DO
    BEGIN
      factor := Pred(2 * i - n);
      sum := sum + factor * GetVectorElement(SortedData.Value, i);
    END;
  Factor := 1 / (n * n - n);
  Result := Factor * Sum;
END;

FUNCTION Ell3(CONST SortedData.Value: VectorTyp): float;

VAR
  sum, factor: float;
  i, n: WORD;

BEGIN
  n := VectorLength(SortedData.Value);
  IF (n = 0)
    THEN
      BEGIN
        CH := WriteErrorMessage('Ell3 from vector OF Length 0');
        DeskriptError := TRUE;
        EXIT;
      END;
  sum := 0;
  FOR i := 1 TO n DO
    BEGIN
      factor := 3 * i * i - 3 * i * n - 3 * i + n * n / 2 + 3 * n / 2 + 1;
      sum := sum + factor * GetVectorElement(SortedData.Value, i);
    END;
  Factor := 2 / (n * n * n - 3 * n * n + 2 * n);
  Result := Factor * Sum;
END;


FUNCTION Ell4(CONST SortedData.Value: VectorTyp): float;

VAR
  sum, factor: float;
  i, n: WORD;

BEGIN
  n := VectorLength(SortedData.Value);
  IF (n = 0)
    THEN
      BEGIN
        CH := WriteErrorMessage('Ell4 from vector OF Length 0');
        DeskriptError := TRUE;
        EXIT;
      END;
  sum := 0;
  FOR i := 1 TO n DO
    BEGIN
      factor := 10 * i * i * i / 3 - 5 * i * i * n - 5 * i * i + 2 * i * n * n + 5 * i * n +
        11 * i / 3 - n * n * n / 6 - n * n - 11 * n / 6 - i;
      sum := sum + factor * GetVectorElement(SortedData.Value, i);
    END;
  Factor := 6 / (n * n * n * n - 6 * n * n * n + 11 * n * n - 6 * n);
  Result := Factor * Sum;
END;


FUNCTION QuartileCoefficientOfSkewness(Q1, Q2, Q3: float): float;

BEGIN
  Result := (Q1 - 2 * Q2 + Q3) / (Q3 - Q1);
END;


FUNCTION CentilCoeffKurtosis(Data.Value: VectorTyp): float;

VAR
  Q1, Q3, Q10, Q90: float;
  n: WORD;

BEGIN
  n := VectorLength(Data.Value);
  IF (n = 0)
    THEN
      BEGIN
        CH := WriteErrorMessage('Centile coefficient OF kurtosis from vector OF Length 0');
        DeskriptError := TRUE;
        EXIT;
      END;
  Q1 := Quantile(Data.Value, 0.25);
  Q3 := Quantile(Data.Value, 0.75);
  Q10 := Quantile(Data.Value, 0.10);
  Q90 := Quantile(Data.Value, 0.90);
  Result := 0.5 * (Q3 - Q1) / (Q90 - Q10);
END;

{ *********** Normalisation and standardisation of vectors *********** }

PROCEDURE MeanNormalise(VAR Data.Value: VectorTyp);

VAR
  Mean: double;
  i, n: WORD;

BEGIN
  n := VectorLength(Data.Value);
  IF (n = 0)
    THEN
      BEGIN
        CH := WriteErrorMessage('Normalisation OF vector OF Length 0');
        DeskriptError := TRUE;
        EXIT;
      END;
  Mean := ArithmeticMean(Data.Value);
  FOR i := 1 TO VectorLength(Data.Value) DO
    SetVectorElement(Data.Value, i, (GetVectorElement(Data.Value, i) - Mean));
END;

PROCEDURE Z_Standardise(VAR Data.Value: VectorTyp);

VAR
  Mean, s: double;
  i, n: WORD;

BEGIN
  n := VectorLength(Data.Value);
  IF (n = 0)
    THEN
      BEGIN
        CH := WriteErrorMessage('Standardisation OF vector OF Length 0');
        DeskriptError := TRUE;
        EXIT;
      END;
  Mean := ArithmeticMean(Data.Value);
  s := StandardDeviation(Data.Value);
  FOR i := 1 TO n DO
    SetVectorElement(Data.Value, i, (GetVectorElement(Data.Value, i) - Mean) / s);
END;


PROCEDURE RobustStandardise(VAR Data.Value: VectorTyp);

VAR
  Mean, s: double;
  i, n: WORD;
  Sorted: VectorTyp;

BEGIN
  n := VectorLength(Data.Value);
  IF (n = 0)
    THEN
      BEGIN
        CH := WriteErrorMessage('Standardisation OF a vector OF Length 0');
        DeskriptError := TRUE;
        EXIT;
      END;
  CopyVector(Data.Value, Sorted);
  Mean := Median(Sorted);  // so that order OF elements IS NOT changed
  s := NaiveQn(Sorted);
  FOR i := 1 TO n DO
    SetVectorElement(Data.Value, i, (GetVectorElement(Data.Value, i) - Mean) / s);
  DestroyVector(Sorted);
END;

{ *********************** description of matrices *************************** }

PROCEDURE VarCovarMatrix(CONST Data.Value: MatrixTyp; VAR VarCovar: MatrixTyp);

VAR
  Rows, Columns, i, j: WORD;
  x, y: VectorTyp;

BEGIN
  Rows := MatrixRows(Data.Value);
  Columns := MatrixColumns(Data.Value);
  CreateMatrix(VarCovar, Columns, Columns, 0.0);
  FOR i := 1 TO Columns DO
    BEGIN
      GetColumn(Data.Value, i, x);
      SetMatrixElement(VarCovar, i, i, Covariance(x, x));
      FOR j := Succ(i) TO Columns DO
        BEGIN
          GetColumn(Data.Value, j, y);
          SetMatrixElement(VarCovar, i, j, Covariance(x, y));
          SetMatrixElement(VarCovar, j, i, GetMatrixElement(VarCovar, i, j));
          DestroyVector(y);
        END;
      DestroyVector(x);
    END;
END;

PROCEDURE VarCov(CONST Data.Value: MatrixTyp; VAR VarCovar: MatrixTyp);

VAR
  Columns, Rows: WORD;
  Ones, I1, Dev, DevT: MatrixTyp;

BEGIN
  Rows := MatrixRows(Data.Value);
  Columns := MatrixColumns(Data.Value);
  CreateMatrix(Ones, Rows, Rows, 1.0);
  MatrixInnerProduct(Ones, Data.Value, I1);
  SkalarMultiplikation(I1, 1 / Rows);
  NegativeMatrix(I1);
  MatrixAdd(Data.Value, I1, Dev);          // deviation scores
  DestroyMatrix(I1);
  DestroyMatrix(Ones);
  MatrixTranspose(Dev, DevT);
  MatrixInnerProduct(DevT, Dev, VarCovar); // deviation score sum OF squares
  SkalarMultiplikation(VarCovar, 1 / Pred(Rows));
  DestroyMatrix(Dev);
  DestroyMatrix(DevT);
END;


PROCEDURE MeanVector(CONST Data.Value: MatrixTyp; VAR Mean: VectorTyp);

VAR
  Columns, i, j: WORD;
  x: VectorTyp;

BEGIN
  Columns := MatrixColumns(Data.Value);
  CreateVector(Mean, Columns, 0.0);
  FOR j := 1 TO Columns DO
    BEGIN
      GetColumn(Data.Value, j, x);
      SetVectorElement(Mean, j, NeumaierSum(x) / ActualElements(x)); // arithmetic mean
      DestroyVector(x);
    END;
END;


PROCEDURE StaVector(CONST Data.Value: MatrixTyp; VAR Sta: VectorTyp);

VAR
  Columns, Rows, i, j: WORD;
  x: VectorTyp;

BEGIN
  Rows := MatrixRows(Data.Value);
  Columns := MatrixColumns(Data.Value);
  CreateVector(Sta, Columns, 0.0);
  FOR j := 1 TO Columns DO
    BEGIN
      GetColumn(Data.Value, j, x);
      SetVectorElement(Sta, j, StandardDeviation(x));
      DestroyVector(x);
    END;
END;

{ *********** Normalisation and standardisation of matrix columns *********** }

PROCEDURE CentreMatrix(VAR A: MatrixTyp);

VAR
  Data.Value: VectorTyp;
  j, Columns: WORD;

BEGIN
  Columns := MatrixColumns(A);
  FOR j := 1 TO Columns DO
    BEGIN
      GetColumn(A, j, Data.Value);
      Centre(Data.Value);
      SetColumn(A, Data.Value, j);
      DestroyVector(Data.Value);
    END;
END;

PROCEDURE StandardiseMatrix(VAR A: MatrixTyp);

VAR
  Data.Value: VectorTyp;
  j, Columns: WORD;

BEGIN
  Columns := MatrixColumns(A);
  FOR j := 1 TO Columns DO
    BEGIN
      GetColumn(A, j, Data.Value);
      Z_Standardise(Data.Value);
      SetColumn(A, Data.Value, j);
      DestroyVector(Data.Value);
    END;
END;

PROCEDURE RobustStandardiseMatrix(VAR A: MatrixTyp);

VAR
  Data.Value: VectorTyp;
  j, Columns: WORD;

BEGIN
  Columns := MatrixColumns(A);
  FOR j := 1 TO Columns DO
    BEGIN
      GetColumn(A, j, Data.Value);
      RobustStandardise(Data.Value);
      SetColumn(A, Data.Value, j);
      DestroyVector(Data.Value);
    END;
END;

{ *********************** Distances and outliers *********************** }

PROCEDURE MahalanobisDistance(CONST Data.Value: MatrixTyp; VAR Dm: VectorTyp);

VAR
  S, T, C, Inter, Res: MatrixTyp;
  Mean: VectorTyp;
  i, j, Rows, Columns: WORD;
  Sum: float;

BEGIN
  Rows := MatrixRows(Data.Value);
  Columns := MatrixColumns(Data.Value);
  CreateVector(Dm, Rows, 0.0);
  CreateMatrix(C, Rows, Columns, 0.0);
  MeanVector(Data.Value, Mean);
  VarCovarMatrix(Data.Value, S);
  InverseMatrix(S);
  FOR i := 1 TO Rows DO                  // (x-µ)
    FOR j := 1 TO Columns DO
      SetMatrixElement(C, i, j, GetMatrixElement(Data.Value, i, j) -
        GetVectorElement(Mean, j));
  MatrixTranspose(C, T);                 // (x-µ)^T
  MatrixInnerProduct(C, S, Inter);       // (x-µ) S^-1
  MatrixInnerProduct(Inter, T, Res);     // (x-µ) S^-1 (x-µ)^T
  FOR j := 1 TO Rows DO                  // Sqrt OF diagonal elements
    SetVectorElement(Dm, j, Sqrt(GetMatrixElement(Res, j, j)));
  DestroyMatrix(S);
  DestroyMatrix(T);
  DestroyMatrix(C);
  DestroyMatrix(Inter);
  DestroyMatrix(Res);
  DestroyVector(Mean);
END;


PROCEDURE RobustDistance(CONST Data.Value: MatrixTyp; VAR Dr: VectorTyp);

VAR
  i, j, n, p: WORD;
  Position, Scale, x: VectorTyp;
  dist, Sum: float;

BEGIN
  n := MatrixRows(Data.Value);
  p := MatrixColumns(Data.Value);
  CreateVector(Position, p, 0.0);
  CreateVector(Scale, p, 0.0);
  CreateVector(Dr, n, 0.0);
  FOR j := 1 TO p DO
    BEGIN
      GetColumn(Data.Value, j, x);
      SetVectorElement(Position, j, HodgesLehmann(x));
      SetVectorElement(Scale, j, NaiveQn(x));
      DestroyVector(x);
    END;
  FOR i := 1 TO n DO
    BEGIN
      sum := 0;
      FOR j := 1 TO p DO
        BEGIN
          dist := (GetMatrixElement(Data.Value, i, j) -
            GetVectorElement(Position, j)) / GetVectorElement(Scale, j);
          Sum := Sum + Sqr(dist);
        END;
      SetVectorElement(Dr, i, Sqrt(Sum));
    END;
  DestroyVector(Position);
  DestroyVector(Scale);
END;

END.

