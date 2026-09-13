@echo off
rem `make examine` for a repository with no Makefile: regenerate the module map and check the claim ledger.
python "%~dp0s26\examine.py" %*
