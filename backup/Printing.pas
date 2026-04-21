unit Printing;
//Wesley R. Parsons,
//wespar@bellouth.net, wespar.com, Miami, Florida, 2025.

interface

uses
  Matrix;

procedure PrintCorrelations(m: word; CData.Value, SData.Value: MatrixTyp; WInput: string);
procedure PrintCovariances(m: word; CData.Value: MatrixTyp);
procedure PrintMatrix(n, m: word; MData.Value: MatrixTyp);

implementation

procedure PrintCorrelations(m: word; CData.Value, SData.Value: MatrixTyp; WInput: string);
var
  i, j: word;
begin
  writeln('Correlation Matrix   Columns: ', m);
  Write('      ');
  for j := 1 to m do Write('Var', j: 3, '    ');
  writeln;
  for i := 1 to m do
  begin
    Write('Var', i: 3, '    ');
    for j := 1 to m do
    begin
      Write(GetMatrixElement(CData.Value, i, j): 10: 6);
    end;
    writeln;
    if WInput = '2' then
      for j := 1 to m do
      begin
        Write('(', GetMatrixElement(SData.Value, i, j): 8: 4, ')');
      end;
    writeln;
  end;
end;

procedure PrintCovariances(m: word; CData.Value: MatrixTyp);
var
  i, j: word;
begin
  writeln('Covariance Matrix   Columns: ', m);
  Write('      ');
  for j := 1 to m do Write('Var', j: 3, '    ');
  writeln;
  for i := 1 to m do
  begin
    Write('Var', i: 3, '    ');
    for j := 1 to m do
    begin
      Write(GetMatrixElement(CData.Value, i, j): 10: 3);
    end;
    writeln;
  end;
end;

procedure PrintMatrix(n, m: word; MData.Value: MatrixTyp);
var
  i, j: word;
begin
  writeln('Matrix   Columns: ', m, '  Rows: ', n);
  for j := 1 to m do Write('   Var ', j: 4, '  ');
  writeln;
  for i := 1 to n do
  begin
    for j := 1 to m do
    begin
      Write(GetMatrixElement(MData.Value, i, j): 10: 5);
    end;
    writeln;
  end;
  readln;
end;

end.
