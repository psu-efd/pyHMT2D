@ECHO OFF

pushd %~dp0

REM Command file for Sphinx documentation

if "%SPHINXBUILD%" == "" (
	set SPHINXBUILD=sphinx-build
)
set SOURCEDIR=source
set BUILDDIR=build

REM Convert README.md to README.rst (no need; changed to rst)
REM pandoc ../README.md --from markdown --to rst -o source/README.rst

REM echo Conversion of ../README.md to source/README.rst Done!

REM Make a copy of the project README.rst to docs/source/
copy ..\README.rst .\source\README_temp.rst

REM modify the location of the logo image
powershell -Command "(gc .\source\README_temp.rst) -replace 'logo/', '_static/images/' | Out-File -encoding UTF8 .\source\README.rst"

REM remove the temp file
del .\source\README_temp.rst

echo Copied project README.rst to source directory.

if "%1" == "" goto help

%SPHINXBUILD% >NUL 2>NUL
if errorlevel 9009 (
	echo.
	echo.The 'sphinx-build' command was not found. Make sure you have Sphinx
	echo.installed, then set the SPHINXBUILD environment variable to point
	echo.to the full path of the 'sphinx-build' executable. Alternatively you
	echo.may add the Sphinx directory to PATH.
	echo.
	echo.If you don't have Sphinx installed, grab it from
	echo.http://sphinx-doc.org/
	exit /b 1
)

%SPHINXBUILD% -M %1 %SOURCEDIR% %BUILDDIR% %SPHINXOPTS% %O%
goto end

:help
%SPHINXBUILD% -M help %SOURCEDIR% %BUILDDIR% %SPHINXOPTS% %O%

:end
popd
