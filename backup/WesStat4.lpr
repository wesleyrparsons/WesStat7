program WesStat4;
{$mode objfpc}{$H+}

//Wesley R. Parsons,
//wespar@bellouth.net, wespar.com, Miami, Florida, 2025.

//Data.Value is contained in WData.Value[0..nRow-1, 0..nCol-1]
//i indexes columns, and j indexes rows
//V1 = nCol-1 and O1 = nRow-1
uses
  Crt,
  Data.ValueHandling,
  Data.ValueInAndOut,
  SysUtils,
  Math,
  WesUnit;

type
  TVector = array of double;
  TMatrix = array of TVector;

const Title: String = 'Wes''s Statistics Version 4 September 2025';

var
  WInput: char;
  i, j: integer;
  nRow, nCol: word;
  WData.Value: TMatrix = ((2, 3, 4, 0, 9), (5, 6, 9, 1, 5), (6, 4, 3, 11, 5),
    (1, 6, 5, 0, 3), (1, 5, 3, 14, 2), (9, 3, 8, 9, 1), (9, 3, 16, 0, 9),
    (1, 4, 9, 1, 5), (10, 5, -3, 14, 2), (9, 13, 18, 9, -1),
    (9, -3, 16, 0, 19), (11, 4, -9, 1, 0));

  function RandomNormal(const m, SD: double): double;
  var
    u1, u2: double;
  begin
    repeat
      u1 := Random;
    until u1 > 0.0;
    u2 := Random;
    Result := m + SD * sqrt(-2.0 * ln(u1)) * cos(2.0 * pi * u2);
  end;

  procedure CorrelationMatrix;
  var
    i, j, V1, O1: integer;
    TwoT, LeftT, RightT: double;
    CData.Value, SData.Value, CVData.Value: TMatrix;
    WVector1, WVector2: TVector;
    CovXY, r: double;
    IsPearson, PrintSig: boolean;
    QInput: char;
  begin
    writeln('Correlations');
    writeln('Options:');
    writeln('1 Correlation Matrix (Pearson''s R)');
    writeln('2 Correlation Matrix (Pearson''s R) with Significance');
    writeln('3 Rank Correlation Matrix (Spearman''s Rho)');
    writeln('4 Rank Correlation Matrix (Spearman''s Rho) with Significance');
    writeln('5 Covariance Matrix');
    writeln('x Exit');
    V1 := nCol - 1;
    O1 := nRow - 1;
    readln(QInput);
    SetLength(CData.Value, nCol, nCol);
    SetLength(SData.Value, nCol, nCol);
    SetLength(CVData.Value, nCol, nCol);
    if (QInput = '2') then begin
      {PearsonRCorrelation(WData.Value, CData.Value);}
      PrintCorrelations(nCol, CData.Value, SData.Value, True, False);
      readln;
    end;
  end;

  procedure MultCorrelations;
  var
    i, j, V1, O1: integer;
    TwoT, LeftT, RightT: double;
    ColArray: TMatrix;
    CData.Value, SData.Value, CVData.Value: TMatrix;
    WVector1, WVector2: TVector;
    CovXY, r: double;
    IsPearson, PrintSig: boolean;
    QInput: string;
  begin
    ClrScr;
    writeln(Title);
    writeln('Main Page');
    writeln('Correlations');
    repeat
    writeln('Options:');
    writeln('1 Correlation Matrix (Pearson''s R)');
    writeln('2 Correlation Matrix (Pearson''s R) with Significance');
    writeln('3 Rank Correlation Matrix (Spearman''s Rho)');
    writeln('4 Rank Correlation Matrix (Spearman''s Rho) with Significance');
    writeln('5 Covariance Matrix');
    writeln('x Exit');
    V1 := nCol - 1;
    O1 := nRow - 1;
    readln(QInput);

    until (QInput='x') or (QInput='X');;
    SetLength(CData.Value, nCol, nCol);
    SetLength(SData.Value, nCol, nCol);
    SetLength(CVData.Value, nCol, nCol);
    SetLength(ColArray, nCol);
    for i := 0 to V1 do SetLength(ColArray[i], nRow);
    for i := 0 to V1 do for j := 0 to O1 do ColArray[i][j] := WData.Value[i, j];
    case QInput of
      '1': begin
        IsPearson := True;
        PrintSig := False;
        for i := 0 to V1 do for j := 0 to V1 do begin
            if IsPearson then begin
              CData.Value[i, j] := PearsonRCorrelation(ColArray[i], ColArray[j]);
              SData.Value[i, j] := PearsonRCorrelationSignificance(CData.Value[i, j], nCol);
              CVData.Value[i, j] := Covariance(ColArray[i], ColArray[j]);
            end
            else
            begin
              CData.Value[i, j] := SpearmanCorrelation(ColArray[i], ColArray[j]);
              SData.Value[i, j] := SpearmanCorrelationSignificance(CData.Value[i, j], nCol);
            end;
          end;
        PrintCorrelations(nCol, CData.Value, SData.Value, IsPearson, PrintSig);
      end;
      '2': begin
        IsPearson := True;
        PrintSig := True;
        for i := 0 to V1 do for j := 0 to V1 do begin
            if IsPearson then begin
              CData.Value[i, j] := PearsonRCorrelation(ColArray[i], ColArray[j]);
              SData.Value[i, j] := PearsonRCorrelationSignificance(CData.Value[i, j], nRow);
              CVData.Value[i, j] := Covariance(ColArray[i], ColArray[j]);
            end
            else
            begin
              CData.Value[i, j] := SpearmanCorrelation(ColArray[i], ColArray[j]);
              SData.Value[i, j] := SpearmanCorrelationSignificance(CData.Value[i, j], nRow);
            end;
          end;
        PrintCorrelations(nCol, CData.Value, SData.Value, IsPearson, PrintSig);
      end;
      '3': begin
        IsPearson := False;
        PrintSig := False;
        for i := 0 to V1 do for j := 0 to V1 do begin
            if IsPearson then begin
              CData.Value[i, j] := PearsonRCorrelation(ColArray[i], ColArray[j]);
              SData.Value[i, j] := PearsonRCorrelationSignificance(CData.Value[i, j], nRow);
              CVData.Value[i, j] := Covariance(ColArray[i], ColArray[j]);
            end
            else
            begin
              CData.Value[i, j] := SpearmanCorrelation(ColArray[i], ColArray[j]);
              SData.Value[i, j] := SpearmanCorrelationSignificance(CData.Value[i, j], nRow);
            end;
          end;
        PrintCorrelations(nCol, CData.Value, SData.Value, IsPearson, PrintSig);
      end;
      '4': begin
        IsPearson := False;
        PrintSig := True;
        if IsPearson then begin
          CData.Value[i, j] := PearsonRCorrelation(ColArray[i], ColArray[j]);
          SData.Value[i, j] := PearsonRCorrelationSignificance(CData.Value[i, j], nRow);
          CVData.Value[i, j] := Covariance(ColArray[i], ColArray[j]);
        end
        else
        begin
          CData.Value[i, j] := SpearmanCorrelation(ColArray[i], ColArray[j]);
          SData.Value[i, j] := SpearmanCorrelationSignificance(CData.Value[i, j], nRow);
        end;

        PrintCorrelations(nCol, CData.Value, SData.Value, IsPearson, PrintSig);
      end;
      '5': begin
        {CovarianceMatrix(WData.Value, nRow, nCol, CVData.Value);}
        PrintCovarianceMatrix(nCol, CVData.Value);
      end;
    end;
    readln;
{    writeln('Save correlation Data.Value? y/n');
    readln(QInput);
    if (QInput = 'y') or (QInput = 'Y') then
    begin
          writeCSVFile('WCorrs', CData.Value);
          writeCSVFile('WSCorrs', SData.Value);
          writeCSVFile('WCovars', CVData.Value);

    end;}
  end;

  procedure BivariateCorrelations;
  var
    C1, C2: Integer;
    r, adjr, sr, rho, adjrho, srho: Double;
    Col1, Col2: TVector;
  begin
    C1 := ReadPositiveInteger('Input number of first variable: ', nCol);
    C2 := ReadPositiveInteger('Input number of second variable: ', nCol);
    SetLength(Col1, nRow);
    SetLength(Col2, nRow);
    for j := 0 to nRow - 1 do begin
      Col1[j] := WData.Value[j, C1];
      Col2[j] := WData.Value[j, C2];
    end;
    r := PearsonRCorrelation(Col1, Col2);
    sr := PearsonRCorrelationSignificance(r, nCol);
    if (nRow > 2) then
      adjr := 1 - ((1 - r * r) * (nRow - 1) / (nRow - 2))
    else
      adjr := 0.0;
    writeln('Variables ', C1, ' and ', C2);
    writeln('Pearsons R = ', r: 8: 10);
    writeln('Pearsons R Significance (Two-tailed) = ', sr: 8: 10);
    writeln('Pearsons R-Squared = ', r * r: 8: 10);
    writeln('Pearsons Adjusted R = ', adjr: 8: 10);
    writeln('Pearsons Adjusted R-Squared = ', adjr * adjr: 8: 10);
    rho := SpearmanCorrelation(Col1, Col2);
    srho := SpearmanCorrelationSignificance(rho, nCol);
    readln;
  end;

  procedure DescriptiveStatsOne;
  var
    i, V, V1: integer;
    q1, q2, q3: double;
    WMean, WSD, WVariance, WSkewness, WKurtosis, WMin, WMax, WRange,
    WMidrange, WMode, WMedian, WTrimedian, WHighMedian, WLowMedian,
    WHDQ, WSRS, WSEoM, WIQR, WHLE, WSn, WQn, WCV, L1, L2, L3, L4, WTukey,
    WHuberEst, HMean, GMean, CoV, MADM, MAD, MDM, MDD, GiniCoef, HerfIndex: double;
    Col: TVector;
  begin
    ClrScr;
    writeln('Descriptive Statistics on One Variable');
    writeln;
    Write('Input number of first variable: ');
    readln(V);      // V1 is the Data.Value in the range 1..NCol
    V1 := V - 1;    // V1 is the Data.Value in the range 0..NCol-1
    SetLength(Col, nRow);
    for i := 0 to nRow - 1 do Col[i] := WData.Value[i, V1];
{   write(' ', col[0],' ', col[1], ' ',col[2]); write(' ', col[3],' ', col[4]);readln; }

    {Measures of Central Tendency and Dispersion}
    writeln;
    Writeln('Parametric Statistics');
    WMean := Mean(Col);
    writeln('Arithmetic Mean: ', SmartFloat(WMean));
    Write('Geometric Mean: ');
    GMean := GeometricMean(Col);
    Writeln(GMean: 0: 12);
    Write('Harmonic Mean: ');
    HMean := HarmonicMean(Col);
    Writeln(HMean: 0: 12);
    Write('Median Deviation from Mean: ');
    MADM := MedianAbsoluteDeviationFromMean(Col);
    Writeln(MADM);
    Writeln(CoV);
    WVariance := Variance(Col);
    WSD := StandardDeviation(Col);
    Write('Coefficient of Variation: ');
    CoV := CoefficientOfVariation(WMean, WSD);
    writeln('Standard Deviation: ', SmartFloat( WSD));
    writeln('Variance: ', SmartFloat( WVariance));
    WSkewness := Skewness(Col);
    writeln('Skewness: ', SmartFloat( WSkewness));
    WKurtosis := Kurtosis(Col);
    writeln('Kurtosis: ', SmartFloat( WKurtosis));
    WMin := MinValue(Col);
    WMax := MaxValue(Col);
    WRange := WMax - WMin;
    WMidrange := WRange / 2;
    writeln('Minimum: ', SmartFloat( WMin));
    writeln('Maximum: ', SmartFloat( WMax));
    writeln('Range: ', SmartFloat( WRange));
    writeln('Midrange: ', SmartFloat( WMidrange));
    MAD := MeanAbsoluteDeviation(Col);
    Writeln('Mean Absolute Deviation: ', SmartFloat( MAD));
    MDD := MedianAbsoluteDeviationFromMean(Col);
    Writeln('Median Deviation from Mean: ', SmartFloat( MDD));

    {Non-Parametric}

    Quartiles(Col, Q1, Q2, Q3);
    writeln;
    Writeln('Non-Parametric Statistics');
    Writeln('Median: ', SmartFloat( WMedian));
    WMode := Mode(Col);
    writeln('Mode: ', SmartFloat( WMode));
    WTrimedian := Trimedian(Col);
    Writeln('Trimedian: ', SmartFloat( WTrimedian));
    WSn := SnEstimator(Col);
    Writeln('Rousseeuw–Croux Sn Estimator: ', SmartFloat( WSn));
    WQn := QnEstimator(Col);
    Writeln('Rousseeuw–Croux Qn Estimator: ', SmartFloat( WQn));
    WSRS := SpearmanRankSkewness(Col);
    writeln('Spearman''s Rank Skewness: ', SmartFloat( WSRS));
    WIQR := InterquartileRange(Col);
    writeln('Interquartile Range: ', SmartFloat( WIQR));
    WHDQ := HarrisDavisQuantile(Col, 4);  //params
    writeln('Harris Davis Quantile: ', SmartFloat( WHDQ));
    WHLE := HodgesLehmann(Col);
    writeln('Hodges Lehmann Estimator: ', SmartFloat( WHLE));
    WHighMedian := HighMedian(Col);
    Writeln('High Median: ', SmartFloat( WHighMedian));
    WHighMedian := LowMedian(Col);
    Writeln('Low Median: ', SmartFloat( WLowMedian));
    WSEoM := StandardErrorOfMedian(Col);
    writeln('Standard Error of Median: ', SmartFloat( WSEoM));
    WHuberEst := HuberEstimator(Col);
    writeln('Huber Estimator: ', SmartFloat( WHuberEst));
    NonParametricLMoments(Col, L1, L2, L3, L4);
    Writeln('L-scale (2nd moment): ', SmartFloat( L2));
    Writeln('L-skewness (3rd moment): ', SmartFloat( L3));
    Writeln('L-kurtosis (4th moment): ', SmartFloat( L4));

    {Measures of Equality}
    writeln;
    Writeln('Measures of Equality');
    Write('Gini Coefficient: ');
    GiniCoef := GiniCoefficient(Col);
    Writeln(GiniCoef));
    readln;
    //Other
    {WTBW := TukeyBiweight(Col);
    writeln('Tukey Biweight: ', SmartFloat( TBW));  }
  end;

  procedure DescriptiveStatsMultiple;
  var
    i, j, k, V1, O1: word;
    WSD, WMean, WVariance, WSkewness, WKurtosis, WMedian, MinVal, MaxVal: double;
    Col: TVector;
    DData.Value: TMatrix;
  begin
    V1 := nCol - 1;
    O1 := nRow - 1;
    SetLength(DData.Value, nCol, 7);       // DData.Value is V1 rows by 6 cols, 0..nCol, 0..6
    for i := 0 to nCol - 1 do begin      // loop thru vars
      SetLength(Col, nRow);
      for j := 0 to nRow - 1 do    // each var, loop thru obs
        Col[j] := WData.Value[j,i];            // and get a col for that var
      {writeln('Col  '); for k := 0 to O1 do write(' ', Col[k]); readln;}
      // compute stats
      DData.Value[i, 0] := Mean(Col);
      DData.Value[i, 1] := StandardDeviation(Col);
      DData.Value[i, 2] := Skewness(Col);
      DData.Value[i, 3] := Kurtosis(Col);
      DData.Value[i, 4] := Median(Col);
      DData.Value[i, 5] := MinValue(Col);
      DData.Value[i, 6] := MaxValue(Col);
    end;
    ClrScr;
    writeln('Descriptive Statistics for Multiple Variables');
    writeln;
    writeln('               Mean         Std Dev      Skewness    Kurtosis      Median       Minimum      Maximum');
    for i := 0 to nCol - 1 do begin
      Write('Col ', i + 1, '     ');
      for k := 0 to 6 do
        Write(DData.Value[i, k]: 12: 5, ' ');
      writeln;
    end;
    readln;
  end;

  procedure ObtainData.Value;
  var
    WData.ValueInput: string;
  begin
    writeln(Title);
    writeln('Obtain Data.Value');
    writeln;
    repeat
      Write('Enter number of columns (variables): ');
      readln(nCol);
      if (nCol <= 0) or (nCol > 100) then writeln('Number is too big or too small.');
    until (nCol > 0) and (nCol < 100);
    repeat
      Write('Enter number of rows (observations): ');
      readln(nRow);
      if (nCol <= 0) or (nCol > 1000) then writeln('Number is too big or too small.');
    until (nRow > 0) and (nRow < 1001);
    writeln;
    SetLength(WData.Value, nRow, nCol);
    writeln('Options: ');
    writeln('  C Create Random Data.Value');
    writeln('  T Test Data.Value, 12 rows x 5 columns');
    writeln('  R Read Data.Value');
    writeln('  M Manually Enter Data.Value');
    write('>');
    readln(WData.ValueInput);

    case WData.ValueInput of
      'C', 'c': begin
        for i := 0 to nRow - 1 do for j := 0 to nCol - 1 do begin
          WData.Value[i][j] := RandomNormal(0, 10);
          writeln(i, ' ', j, ' ', WData.Value[i, j]);
        end;
      end;
      'T', 't': begin
        nRow := 12;
        nCol := 5;
      end;
      'R', 'r': begin
      end;
      'M', 'm': begin
        {ManuallyInputData.Value(nRow, nCol, WData.Value);}
        Write('Enter each piece of Data.Value, follow by <Enter>');
      end;
    end;

  end;

 procedure LinearAlgebra;
  var

    IMatrix: TMatrix;

    EV: TVector;
    EVl: double;
    EVV: byte;
  begin
    writeln('Options: ');
    writeln('  1 Basic Matrix Analysis');
    writeln('  2 Print Data.Value');
    writeln('  H Help');
    writeln('  X Exit');
    readln(Winput);
    {case WInput of
      '1': begin
        PrintMatrix(nRow, nCol, WData.Value);
         CreateVector(EV, nCol, 0.0);
            CopyMatrix(WData.Value, IMatrix);
            InverseMatrix(IMatrix);
            writeln('Inverse Matrix');
            PrintMatrix(nRow, nCol, IMatrix);
            Write('Determinant: ');
            writeln(Determinante(WData.Value));
            Write('Trace: ');
            writeln(MatrixTrace(WData.Value));
          end;
      '2': PrintMatrix(nRow, nCol, WData.Value);
      'H', 'h': writeln('No help available.');
      'X', 'x': Exit;
    end;
  end;}
end;

begin
  while True do begin
    ClrScr;
    writeln(Title);
    writeln('Main Page');
    writeln;
    repeat
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
      write('>');
      readln(WInput);
      WInput := UpCase(WInput);
      case WInput of
        'O', 'o': ObtainData.Value;    //if (nRow<1) or (nCol<1) then writeln('Too little Data.Value to proceed.);
        'D', 'd': DescriptiveStatsOne;
        'M', 'm': DescriptiveStatsMultiple;
        'E', 'e': CorrelationMatrix;
        'C', 'c': MultCorrelations;
        'B', 'b': BivariateCorrelations;
        'L', 'l': LinearAlgebra;
        'P', 'p': PrintMatrix(nRow, nCol, WData.Value);
        'H', 'h': writeln('No help available.');
        'X', 'x': halt;
      end;
    until (WInput = 'X') or (WInput = 'x');
  end;
end.
