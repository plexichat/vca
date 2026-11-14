@echo off
setlocal
cd /d "%~dp0"
python "gui\vector_router_gui.py" %*
endlocal