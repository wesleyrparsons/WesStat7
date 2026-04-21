unit Multivariate;

{$mode ObjFPC}{$H+}{$I proprietary.txt}

interface

{ WesStat7 September 24, 2025, by Wesley R. Parsons, wespar@bellouth.net, www.wespar.com.}

uses
  Crt,
  BasicFuncs,
  DataManager,
  Globals,
  SysUtils;

procedure MultivariateStatistics(const WVar: IVectorType);

implementation

procedure MultivariateStatistics(const WVar: IVectorType);
var
  i, k: Integer;
begin

  // Print out the headers.
  k := Width + 2;
  writeln('Descriptive Statistics for Multiple Variables');
  writeln('Observations: ', nRow);
  writeln;
  writeln('Variable', ' ' : k - 8, 'Mean' : k, 'Std Dev' : k,
          'Median' : k, 'Minimum' : k, 'Maximum' : k);

  // Print out the data.
  for i := 0 to Length(WVar) - 1 do begin
    Write(WData[WVar[i]].Name, ' ' : k - Length(WData[WVar[i]].Name));
    Write(WData[WVar[i]].Mean : k : Precision);
    Write(WData[WVar[i]].StdDev : k : Precision);
    Write(WData[WVar[i]].Median : k : Precision);
    Write(WData[WVar[i]].Min : k : Precision);
    Write(WData[WVar[i]].Max : k : Precision);
    Writeln;
  end;

  // Final touches
  if FirstPass then
    SavePartialData('Save multivariate statistics to a file? (y/n) ', WVar, @MultivariateStatistics);
  Pause;
end;

end.
