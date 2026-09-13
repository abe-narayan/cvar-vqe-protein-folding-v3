#!/bin/sh
# `make examine` for a repository with no Makefile: regenerate the module map and check the claim ledger.
exec python "$(dirname "$0")/s26/examine.py" "$@"
