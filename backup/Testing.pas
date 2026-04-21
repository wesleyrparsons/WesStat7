unit Testing;

{$mode ObjFPC}{$H+}{$I proprietary.txt}

interface

{ WesStat7 September 24, 2025, by Wesley R. Parsons, wespar@bellouth.net, www.wespar.com.}

uses
  BasicFuncs,
  Correlations,
  Globals;

procedure TestRoutines;

implementation

procedure TestRoutines;
begin
  writeln('Testing Code');
  writeln;

  begin
    writeln('Gamma Function, Results from  ChatGPT/Wolfram');
    writeln('Arg            Correct Value                    Function Value');
    writeln('0.0                               NaN',          Gamma(0.0): 25: 20);
    writeln('0.1            9.51350769866873128580',          Gamma(0.1): 25: 20);
    writeln('0.2            4.59084371199880305320',          Gamma(0.2): 25: 20);
    writeln('0.3            2.99156898768759062870',          Gamma(0.3): 25: 20);
    writeln('0.5            1.77245385090551602730',          Gamma(0.5): 25: 20);
    writeln('0.7            1.29805533264755778570',          Gamma(0.7): 25: 20);
    writeln('0.9            1.06862870211931935490',          Gamma(0.9): 25: 20);
    writeln('0.99           1.00587171236088241200',          Gamma(0.99): 25: 20);
    writeln('1.00           1.00000000000000000000',          Gamma(1.0): 25: 20);
    writeln('1.5            0.88622692545275801365',          Gamma(1.5): 25: 20);
    writeln('2.0            1.00000000000000000000',          Gamma(2.0): 25: 20);
    writeln('2.5            1.32934038817913702050',          Gamma(2.5): 25: 20);
    writeln('3.0            2.00000000000000000000',          Gamma(3.0): 25: 20);
    writeln('5.0           24.00000000000000000000',          Gamma(5.0): 25: 20);
    writeln('10.0      362880.00000000000000000000',          Gamma(10.0): 25: 20);
    Pause;
  end;

  begin
    writeln('Standard Normal CDF, Results from ChatGPT/Wolfram');
    writeln('Arg            Correct Value                    Function Value');
    writeln('-3.0           0.00134989803163010350',          CDFNormal(-3.0): 25: 20);
    writeln('-2.5           0.00620966532577613290',          CDFNormal(-2.5): 25: 20);
    writeln('-2.0           0.02275013194817919520',          CDFNormal(-2.0): 25: 20);
    writeln('-1.5           0.06680720126885808520',          CDFNormal(-1.5): 25: 20);
    writeln('-1.0           0.15865525393145705140',          CDFNormal(-1.0): 25: 20);
    writeln('-0.5           0.30853753872598689610',          CDFNormal(-0.5): 25: 20);
    writeln('0.0            0.50000000000000000000',          CDFNormal(0.0): 25: 20);
    writeln('0.5            0.69146246127401310390',          CDFNormal(0.5): 25: 20);
    writeln('1.0            0.84134474606854294860',          CDFNormal(1.0): 25: 20);
    writeln('1.5            0.93319279873114191480',          CDFNormal(1.5): 25: 20);
    writeln('2.0            0.97724986805182080480',          CDFNormal(2.0): 25: 20);
    writeln('2.5            0.99379033467422386710',          CDFNormal(2.5): 25: 20);
    writeln('3.0            0.99865010196836989650',          CDFNormal(3.0): 25: 20);
    Pause;
  end;

  begin
    writeln('Log Gamma Function, Results from ChatGPT/Wolfram');
    writeln('Arg            Correct Value                    Function Value');
    writeln('0.0                NaN',                         nLnGamma(0.0): 25: 20);
    writeln('0.1            2.25271265173420595987',          nLnGamma(0.1): 25: 20);
    writeln('0.2            1.52406382243078456480',          nLnGamma(0.2): 25: 20);
    writeln('0.3            1.09579799481807556540',          nLnGamma(0.3): 25: 20);
    writeln('0.5            0.57236494292470008707',          nLnGamma(0.5): 25: 20);
    writeln('0.7            0.26086724653166651439',          nLnGamma(0.7): 25: 20);
    writeln('0.9            0.03898427592308337150',          nLnGamma(0.9): 25: 20);
    writeln('1.0            0.00000000000000000000',          nLnGamma(1.0): 25: 20);
    writeln('1.5           -0.12078223763524522235',          nLnGamma(1.5): 25: 20);
    writeln('2.0            0.00000000000000000000',          nLnGamma(2.0): 25: 20);
    writeln('2.5            0.28468287047291915963',          nLnGamma(2.5): 25: 20);
    writeln('3.0            0.69314718055994530942',          nLnGamma(3.0): 25: 20);
    writeln('5.0            3.17805383034794561965',          nLnGamma(5.0): 25: 20);
    writeln('10.0          12.80182748008146961121',          nLnGamma(10.0): 25: 20);
    Pause;
  end;


  begin
    writeln('Incomplete Beta Function, Results from ChatGPT/Wolfram, IncompleteBetaCF');
    writeln('x   a   b   Correct Value            Function Value');
    writeln('0.1 0.5 0.5 0.20483276469913350160', IncompleteBetaCF(0.1, 0.5, 0.5):25:20);
    writeln('0.5 0.5 0.5 0.50000000000000000000', IncompleteBetaCF(0.5, 0.5, 0.5):25:20);
    writeln('0.9 0.5 0.5 0.79516723530086619309', IncompleteBetaCF(0.9, 0.5, 0.5):25:20);
    writeln('0.3 1.0 1.0 0.30000000000000000000', IncompleteBetaCF(0.3, 1.0, 1.0):25:20);
    writeln('0.5 1.0 2.0 0.75000000000000000000', IncompleteBetaCF(0.5, 1.0, 2.0):25:20);
    writeln('0.5 2.0 1.0 0.25000000000000000000', IncompleteBetaCF(0.5, 2.0, 1.0):25:20);
    writeln('0.5 2.0 2.0 0.50000000000000000000', IncompleteBetaCF(0.5, 2.0, 2.0):25:20);
    writeln('0.4 2.5 3.0 0.41236100688595678232', IncompleteBetaCF(0.4, 2.5, 3.0):25:20);
    writeln('0.7 5.0 2.0 0.42017500000000000000', IncompleteBetaCF(0.7, 5.0, 2.0):25:20);
    writeln('0.5 10.0 10.0 0.50000000000000000000', IncompleteBetaCF(0.5, 10.0, 10.0):25:20);
    Pause;
  end;

  begin
    writeln('Incomplete Beta Function, Results from ChatGPT/Wolfram, IncompleteBetaSimpson');
    writeln('x   a   b   Correct Value            Function Value');
    writeln('0.1 0.5 0.5 0.20483276469913350160', IncompleteBetaSimpson(0.1, 0.5, 0.5):25:20);
    writeln('0.5 0.5 0.5 0.50000000000000000000', IncompleteBetaSimpson(0.5, 0.5, 0.5):25:20);
    writeln('0.9 0.5 0.5 0.79516723530086619309', IncompleteBetaSimpson(0.9, 0.5, 0.5):25:20);
    writeln('0.3 1.0 1.0 0.30000000000000000000', IncompleteBetaSimpson(0.3, 1.0, 1.0):25:20);
    writeln('0.5 1.0 2.0 0.75000000000000000000', IncompleteBetaSimpson(0.5, 1.0, 2.0):25:20);
    writeln('0.5 2.0 1.0 0.25000000000000000000', IncompleteBetaSimpson(0.5, 2.0, 1.0):25:20);
    writeln('0.5 2.0 2.0 0.50000000000000000000', IncompleteBetaSimpson(0.5, 2.0, 2.0):25:20);
    writeln('0.4 2.5 3.0 0.41236100688595678232', IncompleteBetaSimpson(0.4, 2.5, 3.0):25:20);
    writeln('0.7 5.0 2.0 0.42017500000000000000', IncompleteBetaSimpson(0.7, 5.0, 2.0):25:20);
    writeln('0.5 10.0 10.0 0.50000000000000000000', IncompleteBetaSimpson(0.5, 10.0, 10.0):25:20);
    Pause;
  end;

  begin
    writeln('Student''s T Two-Tailed P-Values, Results from Wolfram');
    writeln('t        df        Correct Value                     Function Value');
    writeln('0.0      1.0       1.00000000000000000000',          TDistPValue(0.0, 1): 25: 20);
    writeln('0.0      5.0       1.00000000000000000000',          TDistPValue(0.0, 5): 25: 20);
    writeln('0.0     10.0       1.00000000000000000000',          TDistPValue(0.0, 10): 25: 20);
    writeln('1.0      1.0       0.50000000000000000000',          TDistPValue(1.0, 1): 25: 20);
    writeln('1.0      5.0       0.36321741467523708634',          TDistPValue(1.0, 5): 25: 20);
    writeln('1.0     10.0       0.34188768150831530742',          TDistPValue(1.0, 10): 25: 20);
    writeln('2.0      1.0       0.33333333333333333333',          TDistPValue(2.0, 1): 25: 20);
    writeln('2.0      5.0       0.10193947917070512917',          TDistPValue(2.0, 5): 25: 20);
    writeln('2.0     10.0       0.07366730167714486524',          TDistPValue(2.0, 10): 25: 20);
    writeln('3.0      1.0       0.20000000000000000000',          TDistPValue(3.0, 1): 25: 20);
    writeln('3.0      5.0       0.02807350705994634381',          TDistPValue(3.0, 5): 25: 20);
    writeln('3.0     10.0       0.01332817248871519571',          TDistPValue(3.0, 10): 25: 20);
    writeln('5.0      5.0       0.00303162638123516738',          TDistPValue(5.0, 5): 25: 20);
    writeln('5.0     10.0       0.00053734005360532506',          TDistPValue(5.0, 10): 25: 20);
    writeln('10.0     10.0      0.00000000000000001714',          TDistPValue(10.0, 10): 25: 20);
    writeln('10.0     100.0     0.00000000000000000000',          TDistPValue(10.0, 100): 25: 20);
    Pause;
  end;

end;
end.

