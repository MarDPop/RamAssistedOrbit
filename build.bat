@echo off

set arg1=%1

setlocal ENABLEDELAYEDEXPANSION

set compiler=g++ -std=c++17 
set debugFlags=-g3 -Wall -Wextra -DDEBUG
set optFlags=-O2 -march=native -DRELEASE
set flags=
set dir=

if "%arg1%" == "release" (
    set flags=%optFlags%
    set dir=Release
) else (
    set flags=%debugFlags%
    set dir=Debug
)

for %%f in (src/*.cpp) do (
    set stmt=%compiler% %flags% -I ./lib -o obj/%%f.o -c src/%%f
    echo !stmt!
    !stmt!
)

setlocal ENABLEDELAYEDEXPANSION

set c= 
for %%f in (obj/*.o) do (
    set c=!c! obj/%%f
)

set stmt=%compiler% %flags% -I ./lib !c! -o bin/%dir%/sim.exe
echo !stmt!
!stmt!

echo built executable

del obj\*.o

REM start "" "python scripts/plots.py"