program WesStat2;
{$mode objfpc}{$H+}
//Statistical routines obtained from Dr. Engelbert Buxsbaum,
//github.com/ebuxbaum/statistics, as of September 1, 2025.
//Licensed under Version 3 of the GNU General Public License
//as published by the Free Software Foundation.
//Modified in some cases by Wesley R. Parsons,
//wespar@bellouth.net, wespar.com, Miami, Florida, 2025.

uses
  Crt,
  Math,
  MathFunc,
  Vector,
  Matrix,
  Zufall,
  Deskript,
  Correlations,
  Data.ValueHandling,
  Printing;

var
  WInput: string;
  WData.Value: MatrixTyp;
  nRow, nCol: word;

Procedure MultCorrelations;

var i, j: word;
    CData.Value, SData.Value, CVData.Value: MatrixTyp;
    WVector1, WVector2: VectorTyp;
    CovXY, r: float;
    s: SignificanceType;
    QInput: String;

begin
CreateMatrix(CData.Value, nCol, nCol, 0.0);
CreateMatrix(SData.Value, nCol, nCol, 0.0);
CreateMatrix(CVData.Value, nCol, nCol, 0.0);
   writeln('Correlations');
   writeln('Options:');
   writeln('1 Correlation Matrix (Pearson''s R)');
   writeln('2 Correlation Matrix (Pearson''s R) with Significance');
   writeln('3 Covariance Matrix');
   writeln('4 Rank Correlation Matrix (Spearman''s RhoR)');
   writeln('x Exit');
   readln(WInput);
   Case WInput of
     '1':
       Begin
    for i := 1 to nCol do
    begin
      GetColumn(WData.Value, i, WVector1);
      for j:= 1 to nCol do
        begin
          GetColumn(WData.Value, j, WVector2);
          r := PearsonProductMomentCorrelation(WVector1, WVector2, CovXY, s);
          SetMatrixElement(CData.Value, i, j, r);
        end;
      end;
    PrintCorrelations(nCol, CData.Value, SData.Value, WInput);
    end;
     '2':       Begin
    for i := 1 to nCol do
    begin
      GetColumn(WData.Value, i, WVector1);
      for j:= 1 to nCol do
        begin
          GetColumn(WData.Value, j, WVector2);
          r := PearsonProductMomentCorrelation(WVector1, WVector2, CovXY, s);
          SetMatrixElement(CData.Value, i, j, r);
        end;
      end;
    PrintCorrelations(nCol, CData.Value, SData.Value, WInput);
    end;
     '3':       Begin
    for i := 1 to nCol do
    begin
      GetColumn(WData.Value, i, WVector1);
      for j:= 1 to nCol do
        begin
          GetColumn(WData.Value, j, WVector2);
          r := PearsonProductMomentCorrelation(WVector1, WVector2, CovXY, s);
          SetMatrixElement(CVData.Value, i, j, CovXY);
        end;
      end;
    PrintCovariances(nCol, CVData.Value);
    end;
     '4':;
     'x', 'X': ;
   end;
  readln;
  writeln('Save correlation Data.Value? y/n');
  readln(QInput);
  if (QInput = 'y') or (QInput = 'Y') then
   begin
  //    writeCSVFile('WCorrs', CData.Value);
  //    writeCSVFile('WSCorrs', SData.Value);
  //    writeCSVFile('WCovars', CVData.Value);

   end;
  DestroyMatrix(CData.Value);
DestroyMatrix(SData.Value);
DestroyMatrix(CVData.Value);
end;

Procedure BivariateCorrelations;

var V1, V2: word;
    WVector1, WVector2: VectorTyp;
    CovXY, r: float;
    s: SignificanceType;

begin
  write('Input number of first variable: ');
  readln(V1);
  write('Input number of second variable: ');
  readln(V2);
  GetColumn(WData.Value, V1, WVector1);
  GetColumn(WData.Value, V2, WVector2);
  r := PearsonProductMomentCorrelation(WVector1, WVector2, CovXY, s);
  writeln('Variables ', V1, ' and ', V2);
  writeln('Pearsons R = ', r);
  writeln ('T-Statistic = ', s.TestValue, '  Degrees of Freedom = ', s.Freedom,
    '  Probablity = ', s.P0);
  readln;
end;

  procedure DescriptiveStatsOne;
  const
    nRow = 50;
  var
    i: word;
    q1, q2, q3: Float;
    AM, SD: Float;
    MData.Value: VectorTyp;
  begin
    CreateVector(MData.Value, nRow, 0.0);
    ClrScr;
    writeln('Descriptive Statistics on One Variable');
    writeln;
    for i := 1 to nRow do
    begin
      SetVectorElement(MData.Value, i, RandomNormal(100, 15));
    end;

    {Measures of Central Tendency and Dispersion}
    AM := ArithmeticMean(MData.Value);
    SD := StandardDeviation(MData.Value);
    Writeln('Basic Statistics');
    writeln('Arithmetic Mean: ', AM: 0: 5);
    writeln('Geometric Mean: ', GeometricMean(MData.Value): 0: 5);
    Writeln('Median Deviation from Mean: ', MedianDeviationFromMean(MData.Value): 0: 5);
    Writeln('Coefficient of Variation: ', CoefficientOfVariation(AM, SD*SD, 0): 0: 5);
    writeln('Standard Deviation: ', SD: 0: 5);
    writeln('Variance: ', Variance(MData.Value): 0: 5);
    writeln('Skewness: ', Skewness(MData.Value, SD, AM): 0: 5);
    writeln('Excess Kurtosis: ', ExcessKurtosis(MData.Value, SD, AM): 0: 5);
    writeln('Median: ', Median(MData.Value): 0: 5);
    writeln('Maximum: ', FindLargest(MData.Value): 0: 5);
    writeln('Minimum: ', FindSmallest(MData.Value): 0: 5);

    {Measures of Equality}
    Writeln('Gini: ', Gini(MData.Value, AM): 0: 5);
    Writeln('Herfindahl-Index: ', HerfindahlIndex(MData.Value): 0: 5);

    {Non-Parametric}
    q2 := Median(MData.Value);
    q1 := Quantile(MData.Value, 0.25);
    q3 := Quantile(MData.Value, 0.75);
    Writeln('Non-Parametric Statistics');
    Writeln('Median: ', q2: 0: 5);
    Writeln('Trimedian: ', Trimedian(MData.Value): 0: 5); // Data.Value have been sorted by Median
    Writeln('Inter-quartile distance: ', InterQuantilDistance(q1, q3): 0: 5);
    Writeln('Mean Absolute Deviation: ', MAD(MData.Value): 0: 5);
    Writeln('Mean Deviation from Mean: ', MeanDeviationFromMean(MData.Value): 0: 5);
    Writeln('Standard Error of Median: ', StandardErrorOfMedian(MData.Value): 0: 5);
    Writeln('Naive Hodges-Lehmann Estimator: ', NaiveHodgesLehmann(MData.Value): 0: 5);
    Writeln('Hodges-Lehmann Estimator: ', HodgesLehmann(MData.Value): 0: 5);
    // Data.Value have been sorted by Median
    Writeln('Naive Rousseeuw–Croux Sn Estimator: ', NaiveSn(MData.Value): 0: 5);
    Writeln('Rousseeuw–Croux Sn Estimator: ', Sn(MData.Value): 0: 5);
    Writeln('Naive Rousseeuw–Croux Qn Estimator: ', NaiveQn(MData.Value): 0: 5);
    Writeln('Rousseeuw–Croux Qn Estimator: ', Qn(MData.Value): 0: 5);
    Writeln('High Median: ', HiMed(MData.Value): 0: 5);
    Writeln('Low Median: ', LoMed(MData.Value): 0: 5);
    Writeln('Trimedian: ', TriMedian(MData.Value): 0: 5);
    Writeln('L-scale (2nd moment): ', Ell2(MData.Value): 0: 5);
    Writeln('L-skewness (3rd moment): ', Ell3(MData.Value): 0: 5);
    Writeln('L-scale (4th moment): ', Ell4(MData.Value): 0: 5);

{FUNCTION Quantile(VAR Data.Value: VectorTyp; q: float): float;

FUNCTION InterQuantilDistance(Q1, Q3: float): float;

FUNCTION QuantileDispersionCoefficient(Q1, Q3: float): float;

FUNCTION QuartileCoefficientOfSkewness(Q1, Q2, Q3: float): float;

FUNCTION CentilCoeffKurtosis(Data.Value: VectorTyp): float;
}
    DestroyVector(MData.Value);
    readln;
  end;

  procedure DescriptiveStatsMultiple;
  var
    i, j: word;
    AM, SD: Float;
    WVect: VectorTyp;
  begin
    ClrScr;
    writeln('Descriptive Statistics for Multiple Variables');
    writeln;
    for i := 1 to nRow do
    begin
      for j := 1 to nCol do
      begin
        SetMatrixElement(WData.Value, i, j, RandomNormal(10, 5));
      end;
    end;
    CreateVector(WVect, nRow, 0.0);
    writeln('               Mean         Std Dev      Skew        Kurt          Median       Min          Max');
    for i := 1 to nCol do
    begin
      GetColumn(WData.Value, i, WVect);
      AM := ArithmeticMean(WVect);
      SD := StandardDeviation(WVect);
      writeln('Var ', i, '     ', AM: 12: 5, ' ', SD: 12: 5, ' ',
      Skewness(WVect, SD, AM): 12: 5, ' ',
      ExcessKurtosis(WVect, SD, AM): 12: 5, ' ',
      Median(WVect): 12: 5, ' ',
      FindLargest(WVect): 12: 5, ' ', FindSmallest(WVect): 12: 5);
    end;
    DestroyVector(WVect);
    readln;
  end;

  procedure ObtainData.Value;
  var
    WData.ValueInput: String;
  begin
    writeln('Wes''s Statistics');
    writeln('Obtain Data.Value');
    writeln('September 2025');
    writeln;
    Write('Enter number of columns (variables): ');
    readln(nCol);
    Write('Enter number of rows (observations): ');
    readln(nRow);
    writeln;
    CreateMatrix(WData.Value, nRow, nCol, 0.0);
    writeln('Options: ');
    writeln('  C Create Random Data.Value');
    writeln('  R Read Data.Value');
    writeln('  M Manually Enter Data.Value');
    readln(WData.ValueInput);

    case WData.ValueInput of
      'C', 'c': begin
        CreateMatrix(WData.Value, nRow, nCol, 0.0);
        CreateRandomData.Value(nRow, nCol, WData.Value);
      end;
      'R', 'r': begin
        CreateMatrix(WData.Value, nRow, nCol, 0.0);
        ReadCSVData.Value(nRow, nCol, WData.Value);
      end;
      'M', 'm': begin
        CreateMatrix(WData.Value, nRow, nCol, 0.0);
        ManuallyInputData.Value(nRow, nCol, WData.Value);
        Write('Enter each piece of Data.Value, follow by <Enter>');
      end;
    end;
  end;

Procedure LinearAlgebra;
Var IMatrix: MatrixTyp;
  begin
  CreateMatrix(IMatrix, nCol, nCol, 0.0);
  PrintMatrix(nRow, nCol, WData.Value);
  CopyMatrix(WData.Value, IMatrix);
  writeln('Inverse Matrix');
  PrintMatrix(nRow, nCol, IMatrix);
  write('Determinant: ');
  writeln(Determinante(WData.Value));
  DestroyMatrix(IMatrix);
end;

begin
  while True do
  begin
    ClrScr;
    writeln('Wes''s Statistics');
    writeln('Main Page');
    writeln('September 2025');
    writeln;
    writeln('Options: ');
    writeln('  O Obtain Data.Value');
    writeln('  D Descriptive Statistics on One Variable');
    writeln('  B Bivariate Correlations');
    writeln('  M Descriptive Statistics on Multiple Variable');
    writeln('  C Multivariate Correlations');
    writeln('  L Linear Algebra');
    writeln('  P Print Data.Value');
    writeln('  H Help');
    writeln('  X Exit');
    readln(Winput);
    case WInput of
      'O', 'o': ObtainData.Value;
      'D', 'd': DescriptiveStatsOne;
      'M', 'm': DescriptiveStatsMultiple;
      'C', 'c': MultCorrelations;
      'B', 'b': BivariateCorrelations;
      'L', 'l': LinearAlgebra;
      'P', 'p': PrintMatrix(nRow, nCol, WData.Value);
      'H', 'h': writeln('No help available.');
      'X', 'x': halt;
    end;
  end;
  DestroyMatrix(WData.Value);
end.
