set DIR=".\temp\makeahmap-%VER%"
xcopy %DIR%\* .\ /A /EXCLUDE:config.ini /E /Y
copy /Y %DIR%\config.ini .\config_reference.ini

del update_script.bat cleanup_script.bat
del /Q makeahmap_cpu.exe

copy /Y %DIR%\cleanup_script.bat .\temp\
start cmd /c %DIR%\cleanup_script.bat
