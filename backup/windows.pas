unit Windows;

{$mode ObjFPC}{$H+}

interface

uses
  SysUtils;

implementation

uses
  Windows;

procedure ResizeConsole(cols, rows: Integer);
var
  hConsole: THandle;
  bufferSize: TCoord;
  windowSize: TSmallRect;
begin
  hConsole := GetStdHandle(STD_OUTPUT_HANDLE);

  // Set buffer size
  bufferSize.X := cols;
  bufferSize.Y := rows;
  SetConsoleScreenBufferSize(hConsole, bufferSize);

  // Set window size
  windowSize.Left := 0;
  windowSize.Top := 0;
  windowSize.Right := cols - 1;
  windowSize.Bottom := rows - 1;
  SetConsoleWindowInfo(hConsole, True, windowSize);
end;

end.

