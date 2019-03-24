REM Набор команд для эмуляции синтаксиса Unix-консоли
REM
REM Пример использования: 
REM 	разместить bat-файл (к примеру, по адресу C:\goUnix.bat),
REM 	в файле запуска консоли Visual Studio vcvarsall.bat 
REM 	(\Community\VC\Auxiliary\Build),
REM 	добавить строчку запуска скрипта @call  "C:\goUnix.bat"
REM
REM Старобинский Егор, 2018 год
@echo off

REM Создание aliases для Unix-команд
echo Setting enviroment...
DOSKEY ls=dir /B $*
DOSKEY cp=copy $*
DOSKEY mv=move $*
DOSKEY rm=del $*
DOSKEY rm-r=rmdir /S /Q $*
DOSKEY bc=calc
DOSKEY cat=powershell get-content $*
DOSKEY grep=findstr $*
DOSKEY pwd=cd
DOSKEY clear=cls
DOSKEY man=help $*
DOSKEY diff=powershell compare-object (get-content $1) (get-content $2)

REM Команда быстрой компиляции
DOSKEY vv=cl /EHsc /Ox /favor:%PROCESSOR_ARCHITECTURE% /Feapp.exe $*

REM Переход в рабочую папку
cd %HOMEPATH%\Documents