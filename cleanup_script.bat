rem Hack make script wait for 1 second
ping -n 2 127.0.0.1 > nul
del /Q .\temp\makeahmap-*.zip
for /d %%i in (".\temp\*") do rmdir /Q /S %%i