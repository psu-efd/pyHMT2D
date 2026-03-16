@ECHO OFF

REM Command for convert pyHMT2D's README.md to README.rst (for Sphinx documentation)

REM mdToRst does not work well
REM python mdToRst.py ../README.md > source/README.rst

pandoc ../README.md --from markdown --to rst -o source/README.rst

echo Conversion of ../README.md to source/README.rst Done!