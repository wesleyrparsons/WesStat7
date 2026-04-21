unit Univariate;

{$mode ObjFPC}{$H+}{$I proprietary.txt}

interface

{ WesStat7 September 24, 2025, by Wesley R. Parsons, wespar@bellouth.net, www.wespar.com.}

uses
  BasicFuncs,
  Crt,
  DataDisplay,
  DataManager,
  Globals,
  SysUtils,
  WesUnit;

procedure DisplayUnivariateStatistics(const Col: WVectorType);

implementation

procedure DisplayUnivariateStatistics(const Col: WVectorType);
var
  WAMean, WMPD, WGMean, WHMean,
  WRMS, WCubic,
  WVariance, WPVariance, WSD, WPSD,
  WCoV, WPCoV, WSEoM, WPSEoM,
  WSkewness, WMedcouple, WPSkewness,
  WKurtosis, WPKurtosis, WEKurtosis, WEPKurtosis, WBMC, WSI, WHE,
  WMAD, WGiniCoef, WLMCoef,
  WShannonn, WShannonb, WTheil,
  WMedian, WMedPD, WHighMedian, WLowMedian, WTrimedian, WMode,
  WHuberEst, WHLE, WTukey,
  WMinm, WMaxm, WRange, WMidrange,
  WQ1, WQ2, WQ3, WMidhinge, WQ3Q1R, WIQR, WSIQR,
  WP9010R, WP9505R, WP7525S, WP9505S, WRR,
  WQCoefDisp, WYule,
  WSn, WQn,
  WMADM, WMedADM, WMADMed, WMedADMed,
  WSEoMed, WSRS, WBowley, WKelly, WL1, WL2, WL3, WL4,
  P10, P90,  P5, P95,  P1, P99,
  PSD1Bot, PSD1Top, PSD2Bot, PSD2Top, PSD3Bot, PSD3Top: FloatType;
  ZCol: WVectorType;
begin

  // Take Col from parameter.
  if Length(Col.Value) = 0 then begin
    writeln('Error in Univariate Statistics: selected variable is empty.');
    Exit;
  end;

{Measures of Central Tendency and Dispersion}
  ClrScr;
  if not Ext then writeln('Basic Univariate Statistics');
  writeln('Variable ', Col.Name, ' with ', Length(Col.Value), ' observations.');
  writeln('' : Tab div 2, 'Parametric Statistics');
  writeln('': (Tab div 2) - 2, 'Measures of Central Tendency');

  WAMean := ArithmeticMean(Col);
  writeTab('Arithmetic Mean: ', WAMean);
  WMPD := MeanPairwiseDistance(Col);
  writeTabln('Mean Pairwise Distance: ', WMPD);

  if Ext then begin
    WGMean := GeometricMean(Col);
    writeTab('Geometric Mean: ', WGMean);
    WHMean := HarmonicMean(Col);
    writeTabln('Harmonic Mean: ', WHMean);
   end;

  if Ext then begin
    WRMS := RootMeanSquare(Col);
    writeTab('Root Mean Square: ', WRMS);
    WCubic := CubicMean(Col);
    writeTabln('Cubic Mean: ', WCubic);
  end;

  writeln;
  writeln('': Tab div 2, 'Measures of Dispersion');
  if not Ext then begin
    WVariance := Variance(Col, Sample);
    writeTab('Variance: ', WVariance);
    WSD := StandardDeviation(Col, Sample);
    writeTabln('Standard Deviation: ', WSD);
  end
  else begin
    WVariance := Variance(Col, Sample);
    writeTab('Sample Variance: ', WVariance);
    WSD := StandardDeviation(Col, True);
    writeTabln('Sample Standard Deviation: ', WSD);
    WPVariance := Variance(Col, False);
    writeTab('Population Variance: ', WPVariance);
    WPSD := StandardDeviation(Col, False);
    writeTabln('Population Standard Deviation: ', WPSD);
  end;

  if not Ext then begin
    WCoV := CoefficientOfVariation(Col, Sample);
    writeTabln('Coefficient of Variation: ', WCoV);
  end
  else begin
    WCoV := CoefficientOfVariation(Col, True);
    writeTab('Sample Coefficient of Variation: ', WCoV);
    WPCoV := CoefficientOfVariation(Col, False);
    writeTabln('Population Coefficient of Variation: ', WPCoV);
  end;

  WMAD := MeanAbsoluteDeviation(Col);
  writeTabln('Mean Absolute Deviation: ', WMAD);

  if not Ext then begin
    WSEoM := StandardErrorOfMean(Col, Sample);
    writeTabln('Standard Error of Mean: ', WSEoM);
  end
  else begin
    WSEoM := StandardErrorOfMean(Col, True);
    writeTabln('Sample Standard Error of Mean: ', WSEoM);
    WPSEoM := StandardErrorOfMean(Col, False);
    writeTabln('Population Standard Error of Mean: ', WPSEoM);
  end;

  if not Ext then begin
    WSkewness := Skewness(Col, Sample);
    writeTab('Skewness: ', WSkewness);
  end
  else begin
    WSkewness := Skewness(Col, True);
    writeTab('Sample Skewness: ', WSkewness);
    WPSkewness := Skewness(Col, False);
    writeTabln('Population Skewness: ', WPSkewness);
  end;

  WMedcouple := Medcouple(Col);
  if Ext then writeTabln('Medcouple: ', WMedcouple);

  if not Ext then begin
    WKurtosis := Kurtosis(Col, True);
    writeTabln('Kurtosis: ', WKurtosis);
    WEKurtosis := ExcessKurtosis(Col, False);
    writeTabln('Excess Kurtosis: ', WEKurtosis);
  end
  else begin
    WKurtosis := Kurtosis(Col, True);
    writeTab('Sample Kurtosis: ', WKurtosis);
    WPKurtosis := Kurtosis(Col, False);
    writeTabln('Population Kurtosis: ', WPKurtosis);
    WEKurtosis := ExcessKurtosis(Col, True);
    writeTab('Sample Excess Kurtosis: ', WEKurtosis);
    WEPKurtosis := ExcessKurtosis(Col, False);
    writeTabln('Population Excess Kurtosis: ', WEPKurtosis);
    WBMC := BimodalityCoefficient(Col);
    writeTabln('Bimodality Coefficient: ', WBMC);
    WSI := SymmetryIndex(Col);
    writeTabln('Symmetry Index: ', WSI);
    WHE := HillEstimator(Col, Trunc(Sqrt(nRow)));  // Also do as an advanced statistic
    writeTabln('Hill Estimator: ', WHE);
  end;
  if Ext then begin
    PercentOutliers(Col, 0.1, P10, P90);
    PercentOutliers(Col, 0.05, P5, P95);
    PercentOutliers(Col, 0.01, P1, P99);
    writeln('Percent outliers at 10% and 90%:      ', P10 : 3: 5, '  ', P90: 3 : 5);
    writeln('Percent outliers at  5% and 95%:      ', P5 : 3: 5, '  ', P95: 3 : 5);
    writeln('Percent outliers at  1% and 99%:      ', P1 : 3: 5, '  ', P99: 3 : 5);
    SetLength(ZCol.Value, Length(Col.Value));
    DeepCopyVector(Col, ZCol);
    StandardizeColumn(ZCol);
    PercentOutliers(ZCol, 1.0 - SD1, PSD1Bot, PSD1Top);
    PercentOutliers(ZCol, 1.0 - SD2, PSD2Bot, PSD2Top);
    PercentOutliers(ZCol, 1.0 - SD3, PSD3Bot, PSD3Top);
    writeln('Percent outliers at -1 SD and 1 SD:   ', PSD1Bot : 3: 5, '  ', PSD1Top : 3: 5);
    writeln('Percent outliers at -2 SD and 2 SD:   ', PSD2Bot : 3: 5, '  ', PSD2Top : 3: 5);
    writeln('Percent outliers at -3 SD and 3 SD:   ', PSD3Bot : 3: 5, '  ', PSD3Top : 3: 5);
  end;

  writeln;
  writeln('' : (Tab div 2) + 2, 'Measures of Equality');
  WGiniCoef := GiniCoefficient(Col);
  writeTab('Gini Coefficient: ', WGiniCoef);
  WTheil := TheilIndex(Col);
  writeTabln('Theil Index: ', WTheil);
  // Separate Entropy Routine using NumBinsFromStruges(X).
  if Ext then begin
    WLMCoef := LorenzMuenzner(Col);
    writeTabln('LorenzMuenzner Coefficient: ', WLMCoef);
    WShannonn := ShannonEntropy(Col, Bins);
    WShannonb := WShannonn / Ln2;
    writeTab('Shannon Entropy (nats): ', WShannonn);
    writeTabln('Shannon Entropy (bits): ', WShannonb);
  end;
  Pause;

  {Non-Parametric}
  writeln;
  writeln('' : Tab div 2, 'Non-Parametric Statistics');
  writeln('' : (Tab div 2) - 2, 'Measures of Central Tendency');
  WMedian := Median(Col);
  writeTab('Median: ', WMedian);
  WMedPD := MedianPairwiseDistance(Col);
  writeTabln('Median Pair Distance: ', WMedPD);

  WHighMedian := HighMedian(Col);
  writeTab('High Median: ', WHighMedian);
  WLowMedian := LowMedian(Col);
  writeTabln('Low Median: ', WLowMedian);
  WTrimedian := Trimedian(Col);
  writeTabln('Trimedian: ', WTrimedian);

  WMode := Mode(Col);
  writeTabln('Mode: ', WMode);

  WHuberEst := HuberEstimator(Col);
  if Ext then writeTabln('Huber Estimator: ', WHuberEst);

  WHLE := HodgesLehmann(Col);
  if Ext then writeTab('Hodges Lehmann Estimator: ', WHLE);

  WTukey := TukeyBiweight(Col);
  if Ext then writeTabln('Tukey Biweight: ', WTukey);

  {Measures of Range}
  writeln;
  writeln('' : (Tab div 2) + 3, 'Measures of Range');
  WMinm := MinmValue(Col);
  WMaxm := MaxmValue(Col);
  WRange := WMaxm - WMinm;
  WMidrange := Midrange(Col);
  writeTab('Minimum: ', WMinm);
  writeTabln('Maximum: ', WMaxm);
  writeTab('Range: ', WRange);
  writeTabln('Midrange: ', WMidrange);

  Quartiles(Col, WQ1, WQ2, WQ3);
  if Ext then writeTab('First Quartile: ', WQ1);
  if Ext then writeTabln('Third Quartile: ', WQ3);

  WMidhinge := Midhinge(WQ1, WQ3);
  writeTab('Midhinge: ', WMidhinge);
  WQ3Q1R := Q3Q1Ratio(WQ1, WQ3);
  writeTabln('Q3/Q1 Ratio: ', WQ3Q1R);

  WIQR := InterquartileRange(Col);
  writeTabln('Interquartile Range: ', WIQR);
  WSIQR := SemiInterquartileRange(WQ1, WQ3);
  writeTabln('Semi-Interquartile Range: ', WSIQR);

  if Ext then begin
    WP9010R := P90P10Ratio(Col);
    writeTab('P90-P10 Ratio: ', WP9010R);
    WP9505R := P95P5Ratio(Col);
    writeTabln('P95-P5 Ratio: ', WP9505R);
    WP7525S := P75P25Spread(Col);
    writeTab('P75-P25 Spread: ', WP7525S);
    WP9505S := P95P5Spread(Col);
    writeTabln('P95-P5 Spread: ', WP9505S);
    WRR := RangeRatio(Col);
    writeTabln('Range Ratio: ', WRR);
  end;

  WQCoefDisp := QuartileCoefficientOfDispersion(WQ1, WQ3);
  writeTabln('Quartile Coefficient of Dispersion: ', WQCoefDisp);

  WYule := YuleSkewness(Col, WQ1, WQ3);
  if Ext then writeTabln('Yule''s Coefficient of Skewness: ', WYule);

  WSn := SnEstimator(Col);
  if Ext then writeTab('Rousseeuw Croux Sn Estimator: ', WSn);
  WQn := QnEstimator(Col);
  if Ext then writeTabln('Rousseeuw Croux Qn Estimator: ', WQn);

  {Measures of Scale}
  writeln;
  writeln('' : (Tab div 2) + 3, 'Measures of Scale');
  WMADM := MeanAbsoluteDeviationFromMean(Col);
  if Ext then writeTabln('Mean Absolute Deviation from Mean: ', WMADM);

  WMedADM := MedianAbsoluteDeviationFromMean(Col);
  if Ext then writeTab('Median Absolute Deviation from Mean: ', WMedADM);

  WMADMed := MeanAbsoluteDeviationFromMedian(Col);
  if Ext then writeTabln('Mean Absolute Deviation from Median: ', WMADMed);

  WMedADMed := MedianAbsoluteDeviationFromMedian(Col);
  if Ext then writeTab('Median Absolute Deviation from Median: ', WMedADMed);

  WSEoMed := StandardErrorOfMedian(Col);
  if Ext then writeTabln('Standard Error of Mean: ', WSEoMed);

  WSRS := SpearmanRankSkewness(Col);
  if Ext then writeTabln('Spearman''s Rank Skewness: ', WSRS);

  WBowley := BowleySkewness(WQ1, WQ2, WQ3);
  if Ext then writeTab('Bowley''s Skewness: ', WBowley);
  WKelly := KellySkewness(Col);
  if Ext then writeTabln('Kelly''s Skewness: ', WKelly);

  if Ext then begin
    PercentOutliersMAD(Col, 1.0, PSD1Bot, PSD1Top);
    PercentOutliersMAD(Col, 2.0, PSD2Bot, PSD2Top);
    PercentOutliersMAD(Col, 3.0, PSD3Bot, PSD3Top);
    writeln('MAD percent outliers at -1 SD and 1 SD:   ', PSD1Bot : 3: 5, '  ', PSD1Top : 3: 5);
    writeln('MAD percent outliers at -2 SD and 2 SD:   ', PSD2Bot : 3: 5, '  ', PSD2Top : 3: 5);
    writeln('MAD percent outliers at -3 SD and 3 SD:   ', PSD3Bot : 3: 5, '  ', PSD3Top : 3: 5);
  end;

  LMoments(Col, WL1, WL2, WL3, WL4);
  writeTabln('L-Scale: ', WL2);
  writeTabln('L-Skewness: ', WL3);
  writeTabln('L-Kurtosis: ', WL4);

  DisplayErrorStack;
  if FirstPass then
    SavePartialData('Save univariate statistics to a file? (y/n) ', Col, @DisplayUnivariateStatistics);
  Pause;
end;

end.
