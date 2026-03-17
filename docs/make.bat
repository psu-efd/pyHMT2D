@ECHO OFF

pushd %~dp0

REM Command file for building pyHMT2D Sphinx documentation (including API pages)

if "%SPHINXBUILD%" == "" (
	set SPHINXBUILD=sphinx-build
)
set SOURCEDIR=source
set BUILDDIR=build

REM The top-level README.rst is already maintained in reStructuredText format.
REM Copy it into docs/source/ and adjust the logo path so it works in the docs build.
copy ..\README.rst .\source\README_temp.rst

REM modify the location of the logo image
powershell -Command "(gc .\source\README_temp.rst) -replace 'logo/', '_static/images/' | Out-File -encoding UTF8 .\source\README.rst"

REM remove the temp file
del .\source\README_temp.rst

echo Updated README.rst in source directory for Sphinx build.

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
