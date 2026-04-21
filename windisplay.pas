unit WinDisplay;

{$mode ObjFPC}{$H+}{$I proprietary.txt}

interface

{ WesStat7 September 24, 2025, by Wesley R. Parsons, wespar@bellouth.net, www.wespar.com.}

procedure ResizeConsoleBuffer;

implementation

uses
  Windows,
  SysUtils;

procedure ResizeConsoleBuffer;
var
  ConsoleHandle: THandle;
  BufferSize: TCoord;
  WindowRect: TSmallRect;
  ConsoleScreenBufferInfo: TConsoleScreenBufferInfo;
begin
  ConsoleHandle := GetStdHandle(STD_OUTPUT_HANDLE);
  if ConsoleHandle = INVALID_HANDLE_VALUE then
    Exit;
  // Set screen buffer size (width: 120 chars, height: 999 lines for scrollback)
  BufferSize.X := 200;  // Width in characters (adjust to your preference)
  BufferSize.Y := 100;  // Height in lines (increases scrollback)
  SetConsoleScreenBufferSize(ConsoleHandle, BufferSize);
  // Optional: Resize the visible window to show more lines (e.g., 50 visible lines)
  // Get current window info first
  GetConsoleScreenBufferInfo(ConsoleHandle, ConsoleScreenBufferInfo);
  WindowRect.Left := 0;
  WindowRect.Top := 0;
  WindowRect.Right := 199;  // Width - 1
  WindowRect.Bottom := 99;  // Height - 1
  SetConsoleWindowInfo(ConsoleHandle, True, WindowRect);
end;

end.



