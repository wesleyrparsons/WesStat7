unit ErrorHandling;

{$mode ObjFPC}{$H+}

interface

type
  String80Type = String[80];
  MessageType = array[0..50] of String80Type;


var
  iES: Integer = 0;
  ErrorStack: MessageType;

procedure DisplayErrorStack(const ErrorStack: MessageType; var iES: Integer);
procedure AddToErrorStack(const Mess: String80Type);

implementation

procedure DisplayErrorStack(const ErrorStack: MessageType; var iES: Integer);
var
  k: Integer;
begin
  writeln('Errors in procedures:');
  for k := 0 to iES - 1 do
    writeln(ErrorStack[k]);
  iES := 0;
end;

procedure AddToErrorStack(const Mess: String80Type);
var k: Integer;
begin
  if iES > 0 then
    for k := 0 to iES - 1 do
      if Mess = ErrorStack[k] then
        Exit;
  ErrorStack[iES] := Mess;
  Inc(iES);
end;

end.

