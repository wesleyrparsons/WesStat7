unit Check;

{$mode ObjFPC}{$H+}{$I proprietary.txt}

interface

{ WesStat7 September 24, 2025, by Wesley R. Parsons, wespar@bellouth.net, www.wespar.com.}

uses
  DataDisplay,
  DataManager,
  Globals;

implementation

// This unit includes procedures to check routines.

// Below is not being used and is out of date.
procedure CheckDeleteRow(var WData: WMatrixType);
begin
  writeln('Delete Row ', 1);
  DisplayMatrixM('', WData, T1);
  DeleteRow(1);
  DisplayMatrixM('', WData, T1);

  writeln('Delete Row ', 5);
  DisplayMatrixM('', WData, T1);
  DeleteRow(5);
  DisplayMatrixM('', WData, T1);

  writeln('Delete Row ', 10);
  DisplayMatrixM('', WData, T1);
  DeleteRow(10);
  DisplayMatrixM('', WData, T1);
end;

end.

