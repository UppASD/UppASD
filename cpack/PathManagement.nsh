!include "EnvVarUpdate.nsh"
!include "LogicLib.nsh"
!include "WinMessages.nsh"

; Function to add installation directory to PATH
Function AddToPath
  Exch $0
  Push $1
  Push $2
  Push $3
  
  ; Get the current PATH
  ReadRegStr $1 HKLM "SYSTEM\CurrentControlSet\Control\Session Manager\Environment" "PATH"
  
  ; Check if our path is already in PATH
  Push "$1"
  Push "$0"
  Call StrStr
  Pop $2
  StrCmp $2 "" 0 AddToPath_done
  
  ; Add our path to PATH
  StrCmp $1 "" AddToPath_no_semicolon
    StrCpy $1 "$1;$0"
    Goto AddToPath_write
  AddToPath_no_semicolon:
    StrCpy $1 "$0"
  
  AddToPath_write:
  WriteRegStr HKLM "SYSTEM\CurrentControlSet\Control\Session Manager\Environment" "PATH" "$1"
  
  ; Broadcast the change
  SendMessage ${HWND_BROADCAST} ${WM_WININICHANGE} 0 "STR:Environment" /TIMEOUT=5000
  
  AddToPath_done:
  Pop $3
  Pop $2
  Pop $1
  Pop $0
FunctionEnd

; Function to remove installation directory from PATH
Function un.RemoveFromPath
  Exch $0
  Push $1
  Push $2
  Push $3
  Push $4
  Push $5
  Push $6
  
  ; Get the current PATH
  ReadRegStr $1 HKLM "SYSTEM\CurrentControlSet\Control\Session Manager\Environment" "PATH"
  
  ; Convert to proper format for StrRep
  StrCpy $5 $1 1 -1 ; Get the last character
  StrCmp $5 ";" 0 +2
    StrCpy $1 $1 -1 ; Remove trailing semicolon
  
  Push $1
  Push "$0;"
  Call un.StrRep
  Pop $2
  
  Push $2
  Push ";$0"
  Call un.StrRep
  Pop $3
  
  Push $3
  Push $0
  Call un.StrRep
  Pop $4
  
  ; Write the modified PATH back
  WriteRegStr HKLM "SYSTEM\CurrentControlSet\Control\Session Manager\Environment" "PATH" "$4"
  
  ; Broadcast the change
  SendMessage ${HWND_BROADCAST} ${WM_WININICHANGE} 0 "STR:Environment" /TIMEOUT=5000
  
  Pop $6
  Pop $5
  Pop $4
  Pop $3
  Pop $2
  Pop $1
  Pop $0
FunctionEnd

; String replacement function
Function un.StrRep
  Exch $R4 ; $R4 = Replacement String
  Exch
  Exch $R3 ; $R3 = String to replace (needle)
  Exch 2
  Exch $R1 ; $R1 = String to do replacement in (haystack)
  Push $R2 ; Replaced haystack
  Push $R5 ; Len (needle)
  Push $R6 ; len (haystack)
  Push $R7 ; Scratch reg
  StrCpy $R2 ""
  StrLen $R5 $R3
  StrLen $R6 $R1
  loop:
    StrCpy $R7 $R1 $R5
    StrCmp $R7 $R3 found
    StrCpy $R7 $R1 1 ; - optimization can be removed if U know len needle=1
    StrCpy $R2 "$R2$R7"
    StrCpy $R1 $R1 $R6 1
    StrLen $R6 $R1
    IntCmp $R6 $R5 loop loop done
  found:
    StrCpy $R2 "$R2$R4"
    StrCpy $R1 $R1 $R6 $R5
    StrLen $R6 $R1
    IntCmp $R6 $R5 loop loop done
  done:
    StrCpy $R3 $R2$R1
    Pop $R7
    Pop $R6
    Pop $R5
    Pop $R2
    Pop $R1
    Pop $R4
    Exch $R3
FunctionEnd

; String search function (StrStr)
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
