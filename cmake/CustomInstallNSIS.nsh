!include "EnvVarUpdate.nsh"

Function AddToPath
    Push "$INSTDIR"
    Call AddEnvVar "PATH"
FunctionEnd

Section -Post
    Call AddToPath
SectionEnd

