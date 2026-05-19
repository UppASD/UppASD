!include "EnvVarUpdate.nsh"
!include "LogicLib.nsh"
!include "WinMessages.nsh"
!define StrStr "!insertmacro StrStr"
 
!macro StrStr ResultVar String SubString
  Push `${String}`
  Push `${SubString}`
  Call StrStr
  Pop `${ResultVar}`
!macroend
 
Function StrStr
/*After this point:
  ------------------------------------------
  $R0 = SubString (input)
  $R1 = String (input)
  $R2 = SubStringLen (temp)
  $R3 = StrLen (temp)
  $R4 = StartCharPos (temp)
  $R5 = TempStr (temp)*/
 
  ;Get input from user
  Exch $R0
  Exch
  Exch $R1
  Push $R2
  Push $R3
  Push $R4
  Push $R5
 
  ;Get "String" and "SubString" length
  StrLen $R2 $R0
  StrLen $R3 $R1
  ;Start "StartCharPos" counter
  StrCpy $R4 0
 
  ;Loop until "SubString" is found or "String" reaches its end
  loop:
    ;Remove everything before and after the searched part ("TempStr")
    StrCpy $R5 $R1 $R2 $R4
 
    ;Compare "TempStr" with "SubString"
    StrCmp $R5 $R0 done
    ;If not "SubString", this could be "String"'s end
    IntCmp $R4 $R3 done 0 done
    ;If not, continue the loop
    IntOp $R4 $R4 + 1
    Goto loop
  done:
 
/*After this point:
  ------------------------------------------
  $R0 = ResultVar (output)*/
 
  ;Remove part before "SubString" on "String" (if there has one)
  StrCpy $R0 $R1 `` $R4
 
  ;Return output to user
  Pop $R5
  Pop $R4
  Pop $R3
  Pop $R2
  Pop $R1
  Exch $R0
FunctionEnd

Section -Post
  ; Diagnostic: write current PATH and the candidate addition to a temp file
  ; so we can debug PATH length issues on target machines.
  ; Read current user PATH (HKCU) and the string we plan to add.
  ReadRegStr $R0 HKCU "Environment" "PATH"
  StrCpy $R1 "$INSTDIR\\bin"

  ; Compute lengths
  StrLen $R2 $R0
  StrLen $R3 $R1
  IntOp $R4 $R2 + $R3

  ; Open debug file in TEMP
  FileOpen $R5 "$TEMP\\uppasd_path_debug.txt" w
  FileWrite $R5 "=== UppASD PATH debug ===\r\n"
  FileWrite $R5 "Current PATH (HKCU):\r\n"
  FileWrite $R5 "$R0\r\n"
  FileWrite $R5 "Current PATH length: $R2\r\n"
  FileWrite $R5 "Candidate to add: $R1\r\n"
  FileWrite $R5 "Candidate length: $R3\r\n"
  FileWrite $R5 "Resulting length (sum): $R4\r\n"
  FileWrite $R5 "(Note: system limits ~= 32766)\r\n"
  FileClose $R5

  ; Now perform the actual PATH update (original behavior)
  Push "$INSTDIR"
  Call AddToPath
SectionEnd

