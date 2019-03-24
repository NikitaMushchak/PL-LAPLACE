;Include Modern UI
  !include "MUI2.nsh"
  
;General
  Name "Planar3D"
  OutFile "install.exe"
  
  Unicode true
  
  InstallDir "$LOCALAPPDATA\Planar3D"
  
  InstallDirRegKey HKCU "Software\Planar3D" ""

  ;Request application privileges for Windows Vista
  RequestExecutionLevel user

;Variables
  Var StartMenuFolder

;Interface Settings
  !define MUI_ABORTWARNING

;Pages
  !insertmacro MUI_PAGE_WELCOME
  !insertmacro MUI_PAGE_COMPONENTS
  !insertmacro MUI_PAGE_DIRECTORY
  
  ;Start Menu Folder Page Configuration
  !define MUI_STARTMENUPAGE_REGISTRY_ROOT "HKCU" 
  !define MUI_STARTMENUPAGE_REGISTRY_KEY "Planar3D" 
  !define MUI_STARTMENUPAGE_REGISTRY_VALUENAME "Start Menu Folder"
  
  !insertmacro MUI_PAGE_STARTMENU Application $StartMenuFolder
  
  !insertmacro MUI_PAGE_INSTFILES
  !insertmacro MUI_PAGE_FINISH
  
  !insertmacro MUI_UNPAGE_WELCOME
  !insertmacro MUI_UNPAGE_CONFIRM
  !insertmacro MUI_UNPAGE_INSTFILES
  !insertmacro MUI_UNPAGE_FINISH
  
;--------------------------------
;Languages
 
  !insertmacro MUI_LANGUAGE "Russian"

Function .onInit
  StrCpy $language ${LANG_RUSSIAN} ;Always force Russian
  System::Call setupapi::IsUserAdmin()i.r0 
  StrCmp $0 "0" 0 +2
  MessageBox MB_OK|MB_ICONINFORMATION "Для корректной установки библиотек могут требоваться права администратора!"
FunctionEnd

;Installer Sections
Section "Planar3D" SecPlanar3D
  SetDetailsPrint none

  CreateDirectory "$INSTDIR\Results"
  SetOutPath "$INSTDIR"  
  File "Builds\Planar3D\*"
  SetOutPath "$INSTDIR\Include"
  File "Builds\Planar3D\Include\*"
  SetOutPath "$INSTDIR\InitialConditions"
  File "Builds\Planar3D\InitialConditions\*"
  SetOutPath "$INSTDIR\mpl-data"
  File "Builds\Planar3D\mpl-data\*"
  SetOutPath "$INSTDIR\mpl-data\fonts\afm"
  File "Builds\Planar3D\mpl-data\fonts\afm\*"
  SetOutPath "$INSTDIR\mpl-data\fonts\pdfcorefonts"
  File "Builds\Planar3D\mpl-data\fonts\pdfcorefonts\*"
  SetOutPath "$INSTDIR\mpl-data\fonts\ttf"
  File "Builds\Planar3D\mpl-data\fonts\ttf\*"
  SetOutPath "$INSTDIR\mpl-data\images"
  File "Builds\Planar3D\mpl-data\images\*"
  SetOutPath "$INSTDIR\mpl-data\sample_data"
  File "Builds\Planar3D\mpl-data\sample_data\*"
  SetOutPath "$INSTDIR\mpl-data\sample_data\axes_grid"
  File "Builds\Planar3D\mpl-data\sample_data\axes_grid\*"
  SetOutPath "$INSTDIR\mpl-data\stylelib"
  File "Builds\Planar3D\mpl-data\stylelib\*"  
  SetOutPath "$INSTDIR\pytz\zoneinfo"
  File "Builds\Planar3D\pytz\zoneinfo\*"
  SetOutPath "$INSTDIR\pytz\zoneinfo\Africa"
  File "Builds\Planar3D\pytz\zoneinfo\Africa\*"
  SetOutPath "$INSTDIR\pytz\zoneinfo\America"
  File "Builds\Planar3D\pytz\zoneinfo\America\*"
  SetOutPath "$INSTDIR\pytz\zoneinfo\America\Argentina"
  File "Builds\Planar3D\pytz\zoneinfo\America\Argentina\*"
  SetOutPath "$INSTDIR\pytz\zoneinfo\America\Indiana"
  File "Builds\Planar3D\pytz\zoneinfo\America\Indiana\*"
  SetOutPath "$INSTDIR\pytz\zoneinfo\America\Kentucky"
  File "Builds\Planar3D\pytz\zoneinfo\America\Kentucky\*"
  SetOutPath "$INSTDIR\pytz\zoneinfo\America\North_Dakota"
  File "Builds\Planar3D\pytz\zoneinfo\America\North_Dakota\*"
  SetOutPath "$INSTDIR\pytz\zoneinfo\Antarctica"
  File "Builds\Planar3D\pytz\zoneinfo\Antarctica\*"
  SetOutPath "$INSTDIR\pytz\zoneinfo\Arctic"
  File "Builds\Planar3D\pytz\zoneinfo\Arctic\*"
  SetOutPath "$INSTDIR\pytz\zoneinfo\Asia"
  File "Builds\Planar3D\pytz\zoneinfo\Asia\*"
  SetOutPath "$INSTDIR\pytz\zoneinfo\Atlantic"
  File "Builds\Planar3D\pytz\zoneinfo\Atlantic\*"
  SetOutPath "$INSTDIR\pytz\zoneinfo\Australia"
  File "Builds\Planar3D\pytz\zoneinfo\Australia\*"
  SetOutPath "$INSTDIR\pytz\zoneinfo\Brazil"
  File "Builds\Planar3D\pytz\zoneinfo\Brazil\*"
  SetOutPath "$INSTDIR\pytz\zoneinfo\Canada"
  File "Builds\Planar3D\pytz\zoneinfo\Canada\*"
  SetOutPath "$INSTDIR\pytz\zoneinfo\Chile"
  File "Builds\Planar3D\pytz\zoneinfo\Chile\*"
  SetOutPath "$INSTDIR\pytz\zoneinfo\Etc"
  File "Builds\Planar3D\pytz\zoneinfo\Etc\*"
  SetOutPath "$INSTDIR\pytz\zoneinfo\Europe"
  File "Builds\Planar3D\pytz\zoneinfo\Europe\*"
  SetOutPath "$INSTDIR\pytz\zoneinfo\Indian"
  File "Builds\Planar3D\pytz\zoneinfo\Indian\*"
  SetOutPath "$INSTDIR\pytz\zoneinfo\Mexico"
  File "Builds\Planar3D\pytz\zoneinfo\Mexico\*"
  SetOutPath "$INSTDIR\pytz\zoneinfo\Pacific"
  File "Builds\Planar3D\pytz\zoneinfo\Pacific\*"
  SetOutPath "$INSTDIR\pytz\zoneinfo\US"
  File "Builds\Planar3D\pytz\zoneinfo\US\*"
  SetOutPath "$INSTDIR\qt4_plugins\accessible"
  File "Builds\Planar3D\qt4_plugins\accessible\*"
  SetOutPath "$INSTDIR\qt4_plugins\codecs"
  File "Builds\Planar3D\qt4_plugins\codecs\*"
  SetOutPath "$INSTDIR\qt4_plugins\graphicssystems"
  File "Builds\Planar3D\qt4_plugins\graphicssystems\*"
  SetOutPath "$INSTDIR\qt4_plugins\iconengines"
  File "Builds\Planar3D\qt4_plugins\iconengines\*"
  SetOutPath "$INSTDIR\qt4_plugins\imageformats"
  File "Builds\Planar3D\qt4_plugins\imageformats\*"
  SetOutPath "$INSTDIR\tcl"
  File "Builds\Planar3D\tcl\*"
  SetOutPath "$INSTDIR\tcl\encoding"
  File "Builds\Planar3D\tcl\encoding\*"
  SetOutPath "$INSTDIR\tcl\http1.0"
  File "Builds\Planar3D\tcl\http1.0\*"
  SetOutPath "$INSTDIR\tcl\msgs"
  File "Builds\Planar3D\tcl\msgs\*"
  SetOutPath "$INSTDIR\tcl\opt0.4"
  File "Builds\Planar3D\tcl\opt0.4\*"
  SetOutPath "$INSTDIR\tcl\tzdata"
  File "Builds\Planar3D\tcl\tzdata\*"
  SetOutPath "$INSTDIR\tcl\tzdata\Africa"
  File "Builds\Planar3D\tcl\tzdata\Africa\*"
  SetOutPath "$INSTDIR\tcl\tzdata\America"
  File "Builds\Planar3D\tcl\tzdata\America\*"
  SetOutPath "$INSTDIR\tcl\tzdata\America\Argentina"
  File "Builds\Planar3D\tcl\tzdata\America\Argentina\*"
  SetOutPath "$INSTDIR\tcl\tzdata\America\Indiana"
  File "Builds\Planar3D\tcl\tzdata\America\Indiana\*"
  SetOutPath "$INSTDIR\tcl\tzdata\America\Kentucky"
  File "Builds\Planar3D\tcl\tzdata\America\Kentucky\*"
  SetOutPath "$INSTDIR\tcl\tzdata\America\North_Dakota"
  File "Builds\Planar3D\tcl\tzdata\America\North_Dakota\*"
  SetOutPath "$INSTDIR\tcl\tzdata\Antarctica"
  File "Builds\Planar3D\tcl\tzdata\Antarctica\*"
  SetOutPath "$INSTDIR\tcl\tzdata\Arctic"
  File "Builds\Planar3D\tcl\tzdata\Arctic\*"
  SetOutPath "$INSTDIR\tcl\tzdata\Asia"
  File "Builds\Planar3D\tcl\tzdata\Asia\*"
  SetOutPath "$INSTDIR\tcl\tzdata\Atlantic"
  File "Builds\Planar3D\tcl\tzdata\Atlantic\*"
  SetOutPath "$INSTDIR\tcl\tzdata\Australia"
  File "Builds\Planar3D\tcl\tzdata\Australia\*"
  SetOutPath "$INSTDIR\tcl\tzdata\Brazil"
  File "Builds\Planar3D\tcl\tzdata\Brazil\*"
  SetOutPath "$INSTDIR\tcl\tzdata\Canada"
  File "Builds\Planar3D\tcl\tzdata\Canada\*"
  SetOutPath "$INSTDIR\tcl\tzdata\Chile"
  File "Builds\Planar3D\tcl\tzdata\Chile\*"
  SetOutPath "$INSTDIR\tcl\tzdata\Etc"
  File "Builds\Planar3D\tcl\tzdata\Etc\*"
  SetOutPath "$INSTDIR\tcl\tzdata\Europe"
  File "Builds\Planar3D\tcl\tzdata\Europe\*"
  SetOutPath "$INSTDIR\tcl\tzdata\Indian"
  File "Builds\Planar3D\tcl\tzdata\Indian\*"
  SetOutPath "$INSTDIR\tcl\tzdata\Mexico"
  File "Builds\Planar3D\tcl\tzdata\Mexico\*"
  SetOutPath "$INSTDIR\tcl\tzdata\Pacific"
  File "Builds\Planar3D\tcl\tzdata\Pacific\*"
  SetOutPath "$INSTDIR\tcl\tzdata\SystemV"
  File "Builds\Planar3D\tcl\tzdata\SystemV\*"
  SetOutPath "$INSTDIR\tcl\tzdata\US"
  File "Builds\Planar3D\tcl\tzdata\US\*"
  SetOutPath "$INSTDIR\tk"
  File "Builds\Planar3D\tk\*"
  SetOutPath "$INSTDIR\tk\images"
  File "Builds\Planar3D\tk\images\*"
  SetOutPath "$INSTDIR\tk\msgs"
  File "Builds\Planar3D\tk\msgs\*"
  SetOutPath "$INSTDIR\tk\ttk"
  File "Builds\Planar3D\tk\ttk\*"
  
  SetOutPath "$INSTDIR"
  
  ;Store installation folder
  WriteRegStr HKCU "Software\Planar3D" "" $INSTDIR
  
  !insertmacro MUI_STARTMENU_WRITE_BEGIN Application
    
    ;Create shortcuts
    CreateDirectory "$SMPROGRAMS\$StartMenuFolder"
    CreateShortcut "$SMPROGRAMS\$StartMenuFolder\Planar3DGUI.lnk" "$INSTDIR\planar3DGUI.exe"
    CreateShortcut "$SMPROGRAMS\$StartMenuFolder\Planar3D.lnk" "$INSTDIR\planar3D.exe"
    CreateShortcut "$SMPROGRAMS\$StartMenuFolder\InitialConditions.lnk" "$INSTDIR\InitialConditions\"
    CreateShortcut "$SMPROGRAMS\$StartMenuFolder\Results.lnk" "$INSTDIR\Results\"
  
  !insertmacro MUI_STARTMENU_WRITE_END
SectionEnd

Section "ParallelPlanar3D" SecParallelPlanar3D
  SetDetailsPrint none

  CreateDirectory "$INSTDIR\Results"
  SetOutPath "$INSTDIR"  
  File "Builds\ParallelPlanar3D\*"
  SetOutPath "$INSTDIR\InitialConditions"
  File "Builds\ParallelPlanar3D\InitialConditions\*"
  
  SetOutPath "$INSTDIR"
  
  ;Store installation folder
  WriteRegStr HKCU "Software\ParallelPlanar3D" "" $INSTDIR
   
  !insertmacro MUI_STARTMENU_WRITE_BEGIN Application
    
    ;Create shortcuts
    CreateDirectory "$SMPROGRAMS\$StartMenuFolder"
    CreateShortcut "$SMPROGRAMS\$StartMenuFolder\ParallelPlanar3D.lnk" "$INSTDIR\parallelPlanar3D.exe"
    CreateShortcut "$SMPROGRAMS\$StartMenuFolder\InitialConditions.lnk" "$INSTDIR\InitialConditions\"
    CreateShortcut "$SMPROGRAMS\$StartMenuFolder\Results.lnk" "$INSTDIR\Results\"
  
  !insertmacro MUI_STARTMENU_WRITE_END
SectionEnd

Section "Библиотеки" SecLibs
  SetDetailsPrint none 
  SetOutPath "$SYSDIR"
  File "Builds\Libs\*"
SectionEnd


Section "Деинсталлятор" SecUninstaller
  SetDetailsPrint none
  
  ;Create uninstaller
  SetOutPath "$INSTDIR"
  WriteUninstaller "$INSTDIR\uninstall.exe"
   
  SetOutPath "$INSTDIR"
  
  !insertmacro MUI_STARTMENU_WRITE_BEGIN Application
    
    ;Create shortcuts
    CreateDirectory "$SMPROGRAMS\$StartMenuFolder"
    CreateShortcut "$SMPROGRAMS\$StartMenuFolder\Uninstall.lnk" "$INSTDIR\uninstall.exe"
  
  !insertmacro MUI_STARTMENU_WRITE_END
SectionEnd

;Descriptions
  ;Language strings
  LangString DESC_SecPlanar3D ${LANG_RUSSIAN} "Программа расчета геометрии трещины ГРП в однородной среде с учетом контраста напряжений, последовательный алгоритм"
  LangString DESC_SecParallelPlanar3D ${LANG_RUSSIAN} "Программа расчета геометрии трещины ГРП в однородной среде с учетом контраста напряжений, параллельный алгоритм"
  LangString DESC_SecLibs ${LANG_RUSSIAN} "Библиотеки, обеспечивающие работу параллельного алгоритма ParallelPlanar3D"
  LangString DESC_SecUninstaller ${LANG_RUSSIAN} "Данная программа производит удаление ранее установленных компонент Planar3D и ParallelPlanar3D"
  
  ;Assign language strings to sections
  !insertmacro MUI_FUNCTION_DESCRIPTION_BEGIN
    !insertmacro MUI_DESCRIPTION_TEXT ${SecPlanar3D} $(DESC_SecPlanar3D)
    !insertmacro MUI_DESCRIPTION_TEXT ${SecParallelPlanar3D} $(DESC_SecParallelPlanar3D)
    !insertmacro MUI_DESCRIPTION_TEXT ${SecLibs} $(DESC_SecLibs)
    !insertmacro MUI_DESCRIPTION_TEXT ${SecUninstaller} $(DESC_SecUninstaller)
  !insertmacro MUI_FUNCTION_DESCRIPTION_END
  
  ShowInstDetails nevershow
  ShowUninstDetails nevershow

;Uninstaller Section
Section "Uninstall"
  SetDetailsPrint none

  Delete "$INSTDIR\*"
  Delete "$INSTDIR\Include\*"
  Delete "$INSTDIR\InitialConditions\*"
  Delete "$INSTDIR\Results\*"
  Delete "$INSTDIR\mpl-data\*"
  Delete "$INSTDIR\mpl-data\fonts\*"
  Delete "$INSTDIR\mpl-data\fonts\afm\*"
  Delete "$INSTDIR\mpl-data\fonts\pdfcorefonts\*"
  Delete "$INSTDIR\mpl-data\fonts\ttf\*"
  Delete "$INSTDIR\mpl-data\images\*"
  Delete "$INSTDIR\mpl-data\sample_data\*"
  Delete "$INSTDIR\mpl-data\sample_data\axes_grid\*"
  Delete "$INSTDIR\mpl-data\stylelib\*"
  Delete "$INSTDIR\pytz\zoneinfo\*"
  Delete "$INSTDIR\pytz\zoneinfo\Africa\*"
  Delete "$INSTDIR\pytz\zoneinfo\America\*"
  Delete "$INSTDIR\pytz\zoneinfo\America\Argentina\*"
  Delete "$INSTDIR\pytz\zoneinfo\America\Indiana\*"
  Delete "$INSTDIR\pytz\zoneinfo\America\Kentucky\*"
  Delete "$INSTDIR\pytz\zoneinfo\America\North_Dakota\*"
  Delete "$INSTDIR\pytz\zoneinfo\Antarctica\*"
  Delete "$INSTDIR\pytz\zoneinfo\Arctic\*"
  Delete "$INSTDIR\pytz\zoneinfo\Asia\*"
  Delete "$INSTDIR\pytz\zoneinfo\Atlantic\*"
  Delete "$INSTDIR\pytz\zoneinfo\Australia\*"
  Delete "$INSTDIR\pytz\zoneinfo\Brazil\*"
  Delete "$INSTDIR\pytz\zoneinfo\Canada\*"
  Delete "$INSTDIR\pytz\zoneinfo\Chile\*"
  Delete "$INSTDIR\pytz\zoneinfo\Etc\*"
  Delete "$INSTDIR\pytz\zoneinfo\Europe\*"
  Delete "$INSTDIR\pytz\zoneinfo\Indian\*"
  Delete "$INSTDIR\pytz\zoneinfo\Mexico\*"
  Delete "$INSTDIR\pytz\zoneinfo\Pacific\*"
  Delete "$INSTDIR\pytz\zoneinfo\US\*"
  Delete "$INSTDIR\qt4_plugins\accessible\*"
  Delete "$INSTDIR\qt4_plugins\codecs\*"
  Delete "$INSTDIR\qt4_plugins\graphicssystems\*"
  Delete "$INSTDIR\qt4_plugins\iconengines\*"
  Delete "$INSTDIR\qt4_plugins\imageformats\*"
  Delete "$INSTDIR\qt4_plugins\*"
  Delete "$INSTDIR\tcl\encoding\*"
  Delete "$INSTDIR\tcl\http1.0\*"
  Delete "$INSTDIR\tcl\msgs\*"
  Delete "$INSTDIR\tcl\opt0.4\*"
  Delete "$INSTDIR\tcl\tzdata\*"
  Delete "$INSTDIR\tcl\tzdata\Africa\*"
  Delete "$INSTDIR\tcl\tzdata\America\*"
  Delete "$INSTDIR\tcl\tzdata\America\Argentina\*"
  Delete "$INSTDIR\tcl\tzdata\America\Indiana\*"
  Delete "$INSTDIR\tcl\tzdata\America\Kentucky\*"
  Delete "$INSTDIR\tcl\tzdata\America\North_Dakota\*"
  Delete "$INSTDIR\tcl\tzdata\Antarctica\*"
  Delete "$INSTDIR\tcl\tzdata\Arctic\*"
  Delete "$INSTDIR\tcl\tzdata\Asia\*"
  Delete "$INSTDIR\tcl\tzdata\Atlantic\*"
  Delete "$INSTDIR\tcl\tzdata\Australia\*"
  Delete "$INSTDIR\tcl\tzdata\Brazil\*"
  Delete "$INSTDIR\tcl\tzdata\Canada\*"
  Delete "$INSTDIR\tcl\tzdata\Chile\*"
  Delete "$INSTDIR\tcl\tzdata\Etc\*"
  Delete "$INSTDIR\tcl\tzdata\Europe\*"
  Delete "$INSTDIR\tcl\tzdata\Indian\*"
  Delete "$INSTDIR\tcl\tzdata\Mexico\*"
  Delete "$INSTDIR\tcl\tzdata\Pacific\*"
  Delete "$INSTDIR\tcl\tzdata\SystemV\*"
  Delete "$INSTDIR\tcl\tzdata\US\*"
  Delete "$INSTDIR\tcl\*"
  Delete "$INSTDIR\tk\images\*"
  Delete "$INSTDIR\tk\msgs\*"
  Delete "$INSTDIR\tk\ttk\*"
  Delete "$INSTDIR\tk\*"
  RMDir /r "$INSTDIR\Include"
  RMDir /r "$INSTDIR\InitialConditions"
  RMDir /r "$INSTDIR\Results"
  RMDir /r "$INSTDIR\mpl-data\fonts\afm"
  RMDir /r "$INSTDIR\mpl-data\fonts\pdfcorefonts"
  RMDir /r "$INSTDIR\mpl-data\fonts\ttf"
  RMDir /r "$INSTDIR\mpl-data\fonts"
  RMDir /r "$INSTDIR\mpl-data\images"
  RMDir /r "$INSTDIR\mpl-data\sample_data\axes_grid"
  RMDir /r "$INSTDIR\mpl-data\sample_data"
  RMDir /r "$INSTDIR\mpl-data\stylelib"
  RMDir /r "$INSTDIR\mpl-data"
  RMDir /r "$INSTDIR\pytz\zoneinfo"
  RMDir /r "$INSTDIR\pytz\zoneinfo\Africa"
  RMDir /r "$INSTDIR\pytz\zoneinfo\America\Argentina"
  RMDir /r "$INSTDIR\pytz\zoneinfo\America\Indiana"
  RMDir /r "$INSTDIR\pytz\zoneinfo\America\Kentucky"
  RMDir /r "$INSTDIR\pytz\zoneinfo\America\North_Dakota"
  RMDir /r "$INSTDIR\pytz\zoneinfo\America"
  RMDir /r "$INSTDIR\pytz\zoneinfo\Antarctica"
  RMDir /r "$INSTDIR\pytz\zoneinfo\Arctic"
  RMDir /r "$INSTDIR\pytz\zoneinfo\Asia"
  RMDir /r "$INSTDIR\pytz\zoneinfo\Atlantic"
  RMDir /r "$INSTDIR\pytz\zoneinfo\Australia"
  RMDir /r "$INSTDIR\pytz\zoneinfo\Brazil"
  RMDir /r "$INSTDIR\pytz\zoneinfo\Canada"
  RMDir /r "$INSTDIR\pytz\zoneinfo\Chile"
  RMDir /r "$INSTDIR\pytz\zoneinfo\Etc"
  RMDir /r "$INSTDIR\pytz\zoneinfo\Europe"
  RMDir /r "$INSTDIR\pytz\zoneinfo\Indian"
  RMDir /r "$INSTDIR\pytz\zoneinfo\Mexico"
  RMDir /r "$INSTDIR\pytz\zoneinfo\Pacific"
  RMDir /r "$INSTDIR\pytz\zoneinfo\US"
  RMDir /r "$INSTDIR\pytz"
  RMDir /r "$INSTDIR\qt4_plugins\accessible"
  RMDir /r "$INSTDIR\qt4_plugins\codecs"
  RMDir /r "$INSTDIR\qt4_plugins\graphicssystems"
  RMDir /r "$INSTDIR\qt4_plugins\iconengines"
  RMDir /r "$INSTDIR\qt4_plugins\imageformats"
  RMDir /r "$INSTDIR\qt4_plugins"
  RMDir /r "$INSTDIR\tcl\encoding"
  RMDir /r "$INSTDIR\tcl\http1.0"
  RMDir /r "$INSTDIR\tcl\msgs"
  RMDir /r "$INSTDIR\tcl\opt0.4"
  RMDir /r "$INSTDIR\tcl\tzdata\Africa"
  RMDir /r "$INSTDIR\tcl\tzdata\America\Argentina"
  RMDir /r "$INSTDIR\tcl\tzdata\America\Indiana"
  RMDir /r "$INSTDIR\tcl\tzdata\America\Kentucky"
  RMDir /r "$INSTDIR\tcl\tzdata\America\North_Dakota"
  RMDir /r "$INSTDIR\tcl\tzdata\America"
  RMDir /r "$INSTDIR\tcl\tzdata\Antarctica"
  RMDir /r "$INSTDIR\tcl\tzdata\Arctic"
  RMDir /r "$INSTDIR\tcl\tzdata\Asia"
  RMDir /r "$INSTDIR\tcl\tzdata\Atlantic"
  RMDir /r "$INSTDIR\tcl\tzdata\Australia"
  RMDir /r "$INSTDIR\tcl\tzdata\Brazil"
  RMDir /r "$INSTDIR\tcl\tzdata\Canada"
  RMDir /r "$INSTDIR\tcl\tzdata\Chile"
  RMDir /r "$INSTDIR\tcl\tzdata\Etc"
  RMDir /r "$INSTDIR\tcl\tzdata\Europe"
  RMDir /r "$INSTDIR\tcl\tzdata\Indian"
  RMDir /r "$INSTDIR\tcl\tzdata\Mexico"
  RMDir /r "$INSTDIR\tcl\tzdata\Pacific"
  RMDir /r "$INSTDIR\tcl\tzdata\SystemV"
  RMDir /r "$INSTDIR\tcl\tzdata\US"
  RMDir /r "$INSTDIR\tcl\tzdata"
  RMDir /r "$INSTDIR\tcl"
  RMDir /r "$INSTDIR\tk\images"
  RMDir /r "$INSTDIR\tk\msgs"
  RMDir /r "$INSTDIR\tk\ttk"
  RMDir /r "$INSTDIR\tk"  
  Delete "$INSTDIR\uninstall.exe"  
  RMDir "$INSTDIR"
  
  !insertmacro MUI_STARTMENU_GETFOLDER Application $StartMenuFolder
  
  Delete "$SMPROGRAMS\$StartMenuFolder\Planar3DGUI.lnk"
  Delete "$SMPROGRAMS\$StartMenuFolder\Planar3D.lnk"
  Delete "$SMPROGRAMS\$StartMenuFolder\ParallelPlanar3D.lnk"
  Delete "$SMPROGRAMS\$StartMenuFolder\InitialConditions.lnk"
  Delete "$SMPROGRAMS\$StartMenuFolder\Results.lnk"
  Delete "$SMPROGRAMS\$StartMenuFolder\Uninstall.lnk"
  RMDir "$SMPROGRAMS\$StartMenuFolder"

  DeleteRegKey /ifempty HKCU "Software\Planar3D"
  DeleteRegKey /ifempty HKCU "Software\ParallelPlanar3D"
SectionEnd